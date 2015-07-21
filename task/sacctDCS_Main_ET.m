function [data,timeStamps] = sacctDCS_Main_ET(xp,startAtLeg,startAtBlock)

try
    %% Setup
    
    for i = 1:xp.nLegs
        data(i).experiment = xp.experiment;
        data(i).codename = xp.codename;
        data(i).task = xp.task;
        data(i).subject = xp.subject;
        data(i).date = xp.date;
        data(i).targetSide = zeros(xp.nBlocks(i),xp.nTrials);
        data(i).leg = xp.legNames{i};
    end
    
    AssertOpenGL; %normally in PsychDefaultSetup
    KbName('UnifyKeyNames'); %normally in PsychDefaultSetup
    Screen('Preference', 'VisualDebuglevel', 3); % to prevent bright splash screens
    
    %define integer luminance values (0-255)
    whiteInt = WhiteIndex(xp.screenNum);
    blackInt = BlackIndex(xp.screenNum);
    grayInt = round(GrayIndex(xp.screenNum));
    [windowPtr, windowRect]=Screen('OpenWindow',xp.screenNum,grayInt); % open a window, color it grey
    [centerX, centerY] = RectCenter(windowRect); % get center of window
    Screen('BlendFunction', windowPtr, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); % set alpha blending to enable smooth, anti-aliased dots
    
    ifi = Screen('GetFlipInterval', windowPtr); %query the inter-frame interval, used for timing
    slack = 0.5*ifi; %headroom of half a frame, to give PTB some time to process Flip commands
    
    Priority(MaxPriority(windowPtr)); % raise the priority of this process to OS to max, for best timing
    GetSecs; % warm up mex functions that can be slow on first call
    HideCursor; % hide the mouse cursor from view
    
    %% Stimuli
    
    %%%TEXT%%%
    textStart = 'Press any key to begin.';
    textBreak = 'Half-way through the current block.\nPlease do not move your head; the task will resume shortly!';
    textBlock = 'You may take a short break now. Press any key to resume';
    textLeg = 'Please wait for the experimenter.';
    textEnd = 'Experiment complete!';
    
    %%%TARGETS%%%
    targetSize = dva2pix(xp.targetSize,[],xp.screenRes,xp.screenDim,xp.screenDist); %convert target size to pixels
    targetEcc = dva2pix(xp.targetEcc,[],xp.screenRes,xp.screenDim,xp.screenDist); %convert target eccentricity to pixels
    
    %%%TIMING%%%
    
    %initialize a structure for keeping stimulus onset timeStamps
    
    for i = 1:xp.nLegs
        timeStamps(i).transDur = round(xp.transTime / ifi) * ifi;
        timeStamps(i).fixDur = zeros(xp.nBlocks(i),xp.nTrials);
        timeStamps(i).targetDur = zeros(xp.nBlocks(i),xp.nTrials);
        timeStamps(i).fix =  zeros(xp.nBlocks(i),xp.nTrials);
        timeStamps(i).target = zeros(xp.nBlocks(i),xp.nTrials);
        timeStamps(i).leg  = xp.legNames{i};
    end
    
    
    %%%EyeLink%%%
    %Specify coordinates to draw to EyeLink PC screen
    stimCoords(1,:) = [centerX-targetEcc centerY];
    stimCoords(2,:) = [centerX, centerY];
    stimCoords(3,:) = [centerX+targetEcc centerY];
    
    %% Experiment loop
    
    DrawFormattedText(windowPtr, textStart, 'center', 'center', blackInt); % draw start text
    Screen('Flip', windowPtr); % flip to screen
    KbStrokeWait; % wait for a key press
    
    for iLeg = startAtLeg:xp.nLegs
        for iBlock = startAtBlock:xp.nBlocks(iLeg)
            
            targetSide = repmat([-1;1],xp.nTrials/2,1);
            targetSide = Shuffle(targetSide);
            
            ISI = xp.fixTime(1) + (xp.fixTime(2) - xp.fixTime(1)).*rand(xp.nTrials,2);
            ISI = round(ISI / ifi) * ifi;
            timeStamps(iLeg).fixDur(iBlock,:) = ISI(:,1);
            timeStamps(iLeg).targetDur(iBlock,:) = ISI(:,2);
            
            %Setup the eye tracker
            EDFname = sacctDCS_ELconfig(windowPtr, xp, iLeg, iBlock, stimCoords);
            
            %Run timer
            countDownTimer(windowPtr,'center','center',5);
            
            %Fixation
            Screen('DrawDots', windowPtr, [centerX centerY], targetSize, xp.targetColor,[],2); % draw the middle dot
            tFixOnset = Screen('Flip', windowPtr);
            Eyelink('Message', sprintf('leg %i block %i started at %f', iLeg, iBlock, tFixOnset));
            
            for iTrial = 1:xp.nTrials
                
                Eyelink('Message', sprintf('trial %i started at %f', iTrial, GetSecs));
                
                %Target
                Screen('DrawDots', windowPtr, [centerX+targetSide(iTrial)*targetEcc centerY], targetSize, xp.targetColor,[],2); % draw a lateral dot
                tTargetOnset = Screen('Flip', windowPtr, tFixOnset + ISI(iTrial,1) - slack);
                Eyelink('Message', sprintf('trial %i phase %i started at %f', iTrial, 1, tTargetOnset));
                
                %Fixation
                Screen('DrawDots', windowPtr, [centerX centerY], targetSize, xp.targetColor,[],2);
                tFixOnset = Screen('Flip', windowPtr, tTargetOnset + ISI(iTrial,2) - slack);
                Eyelink('Message', sprintf('trial %i phase %i started at %f', iTrial, 2, tFixOnset));
                
                %Check for keypresses
                [~,~,keyCode] = KbCheck;
                if keyCode(KbName(xp.abortKey)) % abort the experiment
                    error('Manual quit')
                end
                
                %Store data
                data(iLeg).targetSide(iBlock,iTrial) = targetSide(iTrial);
                timeStamps(iLeg).fix(iBlock,iTrial) = tFixOnset;
                timeStamps(iLeg).target(iBlock,iTrial) = tTargetOnset;
                 
                WaitSecs(xp.saccadeTime); % wait for saccade before sending messages
                
                Eyelink('Message', sprintf('trial %i stopped at %f', iTrial, GetSecs));
                WaitSecs(0.001); % wait a bit, as the eyelink is not always able to write many messages in a short interval
                Eyelink('Message', sprintf('trial %i parameter fixation-target interval : %s', iTrial, ISI(iTrial,1)));
                Eyelink('Message', sprintf('trial %i parameter target-fixation interval : %s', iTrial, ISI(iTrial,2)));
                Eyelink('Message', sprintf('trial %i parameter target saccade direction : %i', iTrial, targetSide(iTrial)));
                Eyelink('Message', sprintf('trial %i parameter fixation saccade direction : %i', iTrial, -1*targetSide(iTrial)));
                
                if ~mod(iTrial,floor(xp.nTrials/(xp.breaksPerBlock+1))+1) % if it's time for a break
                    DrawFormattedText(windowPtr, textBreak, 'center', 'center', blackInt,[],[],[],2); % draw break text
                    Screen('Flip', windowPtr, tFixOnset + timeStamps(iLeg).transDur - slack); % flip one second after final fixation
                    Eyelink('Message', sprintf('leg %i block %i paused at %f', iLeg, iBlock, GetSecs));
                    WaitSecs(5);
                    countDownTimer(windowPtr,'center','center',5);
                    
                    Screen('DrawDots', windowPtr, [centerX centerY], targetSize, xp.targetColor,[],2); % draw the middle dot
                    tFixOnset = Screen('Flip', windowPtr);
                    Eyelink('Message', sprintf('leg %i block %i resumed at %f', iLeg, iBlock, GetSecs));
                end
                 
            end %trial loop
            
            if iBlock < xp.nBlocks(iLeg)
                
                DrawFormattedText(windowPtr, textBlock, 'center', 'center', blackInt); % draw pause text
                Screen('Flip', windowPtr, tFixOnset + timeStamps(iLeg).transDur - slack); % flip one second after final fixation
            elseif iBlock == xp.nBlocks(iLeg) && iLeg < xp.nLegs
                DrawFormattedText(windowPtr, textLeg, 'center', 'center', blackInt); % draw pause text
                Screen('Flip', windowPtr, tFixOnset + timeStamps(iLeg).transDur - slack);
            else
                DrawFormattedText(windowPtr, textEnd, 'center', 'center', blackInt); % draw pause text
                Screen('Flip', windowPtr, tFixOnset + timeStamps(iLeg).transDur - slack);
            end
             Eyelink('Message', sprintf('leg %i block %i stopped at %f', iLeg, iBlock, GetSecs));
             
                %save back-up of data so far
                filename = [xp.backupFolder xp.codename '_' xp.subject '_' xp.tDCS '_' xp.legNames{iLeg} '_' ...
                    datestr(now, 'yyyy-mm-dd_HH-MM-SS') '_block_' int2str(iBlock) '_trial_' int2str(iTrial) '.mat'];
                save(filename,'data', 'xp', 'timeStamps')
                
                % save ET-file
                fullEDFname = sprintf([xp.dataFolder '%s_%s_%s_%s_block%i_%s.edf'], xp.codename, xp.subject, xp.tDCS, data(iLeg).leg, iBlock, datestr(now, 'yyyy-mm-dd_HH-MM-SS'));
                Eyelink('CloseFile');
                Eyelink('WaitForModeReady', 500); % helps to avoid errors in retrieving files
                try
                    status = Eyelink('ReceiveFile', EDFname, fullEDFname); %this collects the file from the eyelink
                    disp(status);
                    disp(['File ' fullEDFname ' saved to disk']);
                catch
                    warning(['File ' fullEDFname ' not saved to disk']);
                end
                KbStrokeWait;
                 
        end % block loop
   
    end % leg loop

    
    %% Save data and close
    filename = [xp.dataFolder xp.codename '_' xp.subject '_' xp.tDCS '_' xp.task '.mat'];
    save(filename,'data', 'xp', 'timeStamps');
    
    sca; %Screen('CloseAll'), to regain control of monitor
    ShowCursor; % show mouse again
    Priority(0); % restore priority
    Eyelink('command', 'clear_screen 0'); % clear tracker display
    
catch err
    
    if exist('iTrial' , 'var')
        %save back-up of data so far
        filename = [xp.backupFolder xp.codename '_' xp.subject '_' xp.tDCS '_' xp.task '_' ...
            datestr(now, 'yyyy-mm-dd_HH-MM-SS') '_leg_' int2str(iLeg) '_block_' int2str(iBlock) '_trial_' int2str(iTrial) '.mat'];
        save(filename,'data', 'xp', 'timeStamps')
        
        % save ET-file
        fullEDFname = sprintf([xp.backupFolder '%s_%s_%s_block%i_%s.edf'], xp.codename, xp.subject, data(iLeg).leg, iBlock, datestr(now, 'yyyy-mm-dd_HH-MM-SS'));
        Eyelink('CloseFile');
        Eyelink('command', 'clear_screen 0'); % clear tracker display
        Eyelink('WaitForModeReady', 500); % helps to avoid errors in retrieving files
        try
            status = Eyelink('ReceiveFile', EDFname, fullEDFname); %this collects the file from the eyelink
            disp(status);
            disp(['File ' fullEDFname ' saved to disk']);
        catch
            warning(['File ' fullEDFname ' not saved to disk']);
        end
    end
    
    sca;
    ShowCursor;
    Priority(0);
    rethrow(err); % issue the caught error
end

