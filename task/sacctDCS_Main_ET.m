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
    textBreak = 'Blink to your hearts content, but please do not move your head!\n Look at the dot for 2 seconds, then press any key to resume.';
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
    
    breakTrials = zeros(1,xp.breaksPerBlock);
    
    for i = 1:xp.breaksPerBlock
        breakTrials(i) = floor(xp.nTrials/(xp.breaksPerBlock+1)*i);
    end
    
    %%%EyeLink%%%
    %Specify coordinates to draw to EyeLink PC screen
    stimCoords(1,:) = [centerX-targetEcc centerY];
    stimCoords(2,:) = [centerX, centerY];
    stimCoords(3,:) = [centerX+targetEcc centerY];
    
    %Setup the eye tracker
    [ELdefaults, ELconfig] = sacctDCS_ELconfig(stimCoords);
    
    %% Experiment loop
    
    DrawFormattedText(windowPtr, textStart, 'center', 'center', blackInt); % draw start text
    Screen('Flip', windowPtr); % flip to screen
    KbStrokeWait; % wait for a key press
    
    for iLeg = startAtLeg:xp.nLegs
        for iBlock = startAtBlock:xp.nBlocks(iLeg)
            
            %%%RANDOMIZE%%%
            
            %%%Randomize hemifield%%%
            targetSide = repmat([-1;1],xp.nTrials/2,1);
            targetSide = Shuffle(targetSide);
            
            %%%Randomize ISI%%%
            lowBound = xp.fixTime(1);
            upBound = xp.fixTime(2);
            Mu = xp.mu - lowBound;
            % get cumulative probalities from exponential cdf
            lowBoundP = expcdf(lowBound, Mu);
            upBoundP = expcdf(upBound , Mu);
            % draw uniformly from this range
            probDraw = lowBoundP + (upBoundP - lowBoundP) .* rand(xp.nTrials, 2);
            % evaluate the exponential funxtion at the drawn probabilities
            ISI = expinv(probDraw, Mu);
            % round to ifi and store
            ISI = round(ISI / ifi) * ifi;
            timeStamps(iLeg).fixDur(iBlock,:) = ISI(:,1);
            timeStamps(iLeg).targetDur(iBlock,:) = ISI(:,2);
            
            %Start EyeLink calibration and recording
            EDFname = sacctDCS_ELrecord(windowPtr, xp, iLeg, iBlock);
            
            %Draw start screen
            DrawFormattedText(windowPtr, textBreak, 'center', centerY*0.8, blackInt,[],[],[],2);
            Screen('DrawDots', windowPtr, [centerX centerY], targetSize, xp.driftColor,[],2);
            Screen('Flip', windowPtr); 
            KbStrokeWait;
            
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
                
                if ismember(iTrial, breakTrials) % if it's time for a break
                    DrawFormattedText(windowPtr, textBreak, 'center', centerY*0.8, blackInt,[],[],[],2); % draw pause text
                    Screen('DrawDots', windowPtr, [centerX centerY], targetSize, xp.driftColor,[],2);
                    Screen('Flip', windowPtr, tISIonset(1) + timeStamps(iLeg).transDur - slack); % flip one second after final fixation
                    Eyelink('Message', sprintf('leg %i block %i paused at %f', iLeg, iBlock, GetSecs));
                    KbStrokeWait;
                    
                    Screen('DrawDots', windowPtr, [centerX centerY], targetSize, xp.targetColor,[],2); % draw the middle dot
                    tFixOnset = Screen('Flip', windowPtr);
                    Eyelink('Message', sprintf('leg %i block %i resumed at %f', iLeg, iBlock, GetSecs));
                end
                 
            end %trial loop
            
            if iBlock < xp.nBlocks(iLeg)
                breakText = textBlock;
            elseif iBlock == xp.nBlocks(iLeg) && iLeg < xp.nLegs
                breakText = textLeg;
            else
                breakText = textEnd;
            end
            DrawFormattedText(windowPtr, breakText, 'center', centerY*0.8, blackInt); % draw pause text
            Eyelink('Message', sprintf('leg %i block %i stopped at %f', iLeg, iBlock, GetSecs));
             
                %save back-up of data so far
                filename = [xp.backupFolder xp.codename '_' xp.subject '_' xp.tDCS '_' xp.legNames{iLeg} '_' ...
                    datestr(now, 'yyyy-mm-dd_HH-MM-SS') '_leg_' int2str(iLeg) '_block_' int2str(iBlock) '_trial_' int2str(iTrial) '.mat'];
                save(filename,'data', 'xp', 'timeStamps')
                
                % save ET-file
                sacctDCS_ELsave(xp,EDFname,iBlock,iLeg)
                
                KbStrokeWait;
                 
        end % block loop
   
    end % leg loop

    
    %% Save data and close
    filename = [xp.dataFolder xp.codename '_' xp.subject '_' xp.tDCS '_' xp.task '.mat'];
    save(filename,'data', 'xp', 'timeStamps');
    
    sca; %Screen('CloseAll'), to regain control of monitor
    ShowCursor; % show mouse again
    Priority(0); % restore priority
    sacctDCS_ELwrapup(ELdefaults) %restore eyetracker defaults
    
catch err
    
    if exist('iTrial' , 'var')
        %save back-up of data so far
        filename = [xp.backupFolder xp.codename '_' xp.subject '_' xp.tDCS '_' xp.task '_' ...
            datestr(now, 'yyyy-mm-dd_HH-MM-SS') '_leg_' int2str(iLeg) '_block_' int2str(iBlock) '_trial_' int2str(iTrial) '.mat'];
        save(filename,'data', 'xp', 'timeStamps')
        
        % save ET-file
        sacctDCS_ELsave(xp,EDFname,iBlock,iLeg)
        
    end
    
    sca;
    ShowCursor;
    Priority(0);
    sacctDCS_ELwrapup(ELdefaults) %restore eyetracker defaults
    rethrow(err); % issue the caught error
end

