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
        data(i).saccLatency = nan(xp.nBlocks(i),xp.nTrials,2);
        data(i).saccAccuracy = nan(xp.nBlocks(i),xp.nTrials,2);
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
    textFB = ['Your average reaction time was %.f milliseconds.\n' ...
        'The endpoint of your eye movements were off by %.2f millimeters on average.\n ' ...
        'Great job; try and see if you can get these scores down even further!'];
    
    %%%TARGETS%%%
    targetSize = dva2pix(xp.targetSize,[],xp.screenRes,xp.screenDim,xp.screenDist); %convert target size to pixels
    targetEcc = dva2pix(xp.targetEcc,[],xp.screenRes,xp.screenDim,xp.screenDist); %convert target eccentricity to pixels
    
    %%%TIMING%%%
    
    %initialize a structure for keeping stimulus onset timeStamps
    
    for i = 1:xp.nLegs
        timeStamps(i).transDur = round(xp.transTime / ifi) * ifi;
        timeStamps(i).maxSaccadeTime = round(xp.maxSaccadeTime / ifi) * ifi;
        timeStamps(i).fixDur = zeros(xp.nBlocks(i),xp.nTrials);
        timeStamps(i).targetDur = zeros(xp.nBlocks(i),xp.nTrials);
        timeStamps(i).fix =  zeros(xp.nBlocks(i),xp.nTrials);
        timeStamps(i).target = zeros(xp.nBlocks(i),xp.nTrials);
        timeStamps(i).leg  = xp.legNames{i};
        timeStamps(i).saccTarget = zeros(xp.nBlocks(i),xp.nTrials);
        timeStamps(i).saccFix = zeros(xp.nBlocks(i),xp.nTrials);
        
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
    maxSaccDev = dva2pix(xp.maxSaccDev, [] ,xp.screenRes,xp.screenDim,xp.screenDist);
    
    %Setup the eye tracker
    [ELdefaults, ELconfig] = sacctDCS_ELconfig(stimCoords);
    
    %% Experiment loop
    
    DrawFormattedText(windowPtr, textStart, 'center', 'center', blackInt); % draw start text
    Screen('Flip', windowPtr); % flip to screen
    fprintf('Waiting for subject to respond to: "%s"\n', textStart); % print what subject sees to command window
    KbStrokeWait; % wait for a key press
    fprintf('Subject indicated he/she''s ready for calibration\n');
    
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
            fprintf('Subject is taking a short break\n');
            KbStrokeWait;
            fprintf('Subject started leg %i block %i\n', iLeg, iBlock);
            
            %Fixation
            Screen('DrawDots', windowPtr, [centerX centerY], targetSize, xp.targetColor,[],2); % draw the middle dot
            tISIonset(1) = Screen('Flip', windowPtr);
            Eyelink('Message', sprintf('leg %i block %i started at %f', iLeg, iBlock, tISIonset(1)));
            
            
            for iTrial = 1:xp.nTrials
                
                Eyelink('Message', sprintf('trial %i started at %f', iTrial, GetSecs));
                
                %Target
                Screen('DrawDots', windowPtr, [centerX+targetSide(iTrial)*targetEcc centerY], targetSize, xp.targetColor,[],2); % draw a lateral dot
                tISIonset(2) = Screen('Flip', windowPtr, tISIonset(1) + ISI(iTrial,1) - slack);
                Eyelink('Message', sprintf('trial %i phase %i started at %f', iTrial, 1, tISIonset(2)));
                stimTime = Eyelink('TrackerTime');
                timeStamps(iLeg).target(iBlock,iTrial) = tISIonset(2);

                % Check for saccades
                saccMade = false;
                while GetSecs < tISIonset(2) + timeStamps(iLeg).maxSaccadeTime % check for saccades % check untill the minimum ISI, minus 2 frames to have some slack for the next flip
                    evtype = Eyelink('GetNextDataType'); % check identity of current event
                    if evtype==ELconfig.ENDSACC % if it is an end saccade event
                         evt = Eyelink('GetFloatData', evtype); % fetch data for the event
                         if sqrt(((centerX+targetSide(iTrial)*targetEcc) - evt.genx)^2 + (centerY - evt.geny)^2) < maxSaccDev
                             tISIonset(2) = GetSecs;
                             saccMade = true;
                             break
                         end
                    end
                end
                
                if saccMade % if there was a saccade/we were able to get it
                    data(iLeg).saccLatency(iBlock,iTrial,1) = evt.sttime - stimTime*1000; %calculate latency
                    data(iLeg).saccAccuracy(iBlock,iTrial,1) = sqrt(((centerX+targetSide(iTrial)*targetEcc) - evt.genx)^2 + (centerY - evt.geny)^2); %calculate error in pixels
                else
                    tISIonset(2) = tISIonset(2) + timeStamps(iLeg).maxSaccadeTime;
                end
                timeStamps(iLeg).saccTarget(iBlock,iTrial) = tISIonset(2);
                
                %Fixation
                Screen('DrawDots', windowPtr, [centerX centerY], targetSize, xp.targetColor,[],2);
                tISIonset(1) = Screen('Flip', windowPtr, tISIonset(2) + ISI(iTrial,2) - slack);
                Eyelink('Message', sprintf('trial %i phase %i started at %f', iTrial, 2, tISIonset(1)));
                stimTime = Eyelink('TrackerTime');
                timeStamps(iLeg).fix(iBlock,iTrial) = tISIonset(1);

                % Check for saccades
                saccMade = false;
                while GetSecs < tISIonset(1) + timeStamps(iLeg).maxSaccadeTime % check for saccades
                    evtype = Eyelink('GetNextDataType'); % check identity of current event
                    if evtype==ELconfig.ENDSACC % if it is an end saccade event
                        evt = Eyelink('GetFloatData', evtype); % fetch data for the event
                        if sqrt((centerX - evt.genx)^2 + (centerY - evt.geny)^2) < maxSaccDev
                            tISIonset(1) = GetSecs;
                            saccMade = true;
                            break
                        end
                    end
                end
                
                if saccMade % if there was a saccade/we were able to get it
                    data(iLeg).saccLatency(iBlock,iTrial,2) = evt.sttime - stimTime*1000; %calculate latency
                    data(iLeg).saccAccuracy(iBlock,iTrial,2) = sqrt((centerX - evt.genx)^2 + (centerY - evt.geny)^2); %calculate error in pixels
                else
                    tISIonset(1) = tISIonset(1) + timeStamps(iLeg).maxSaccadeTime;
                end
                timeStamps(iLeg).saccFix(iBlock,iTrial) = tISIonset(1);
                
                %Check for keypresses
                [~,~,keyCode] = KbCheck;
                if keyCode(KbName(xp.abortKey)) % abort the experiment
                    error('Manual quit')
                end
                
                data(iLeg).targetSide(iBlock,iTrial) = targetSide(iTrial);
                
                Eyelink('Message', sprintf('trial %i stopped at %f', iTrial, GetSecs));
                WaitSecs(0.001); % wait a bit, as the eyelink is not always able to write many messages in a short interval
                Eyelink('Message', sprintf('trial %i parameter fixation-target interval : %f', iTrial, ISI(iTrial,1)));
                Eyelink('Message', sprintf('trial %i parameter target-fixation interval : %f', iTrial, ISI(iTrial,2)));
                Eyelink('Message', sprintf('trial %i parameter target saccade direction : %i', iTrial, targetSide(iTrial)));
                Eyelink('Message', sprintf('trial %i parameter fixation saccade direction : %i', iTrial, -1*targetSide(iTrial)));
                WaitSecs(0.001);
                Eyelink('Message', sprintf('trial %i parameter target saccade latency : %f', iTrial, data(iLeg).saccLatency(iBlock,iTrial,1)));
                Eyelink('Message', sprintf('trial %i parameter fixation saccade latency : %f', iTrial, data(iLeg).saccLatency(iBlock,iTrial,2)));
                Eyelink('Message', sprintf('trial %i parameter target saccade accuracy : %f', iTrial, data(iLeg).saccAccuracy(iBlock,iTrial,1)));
                Eyelink('Message', sprintf('trial %i parameter fixation saccade accuracy : %f', iTrial, data(iLeg).saccAccuracy(iBlock,iTrial,2)));
                
                if ismember(iTrial, breakTrials) % if it's time for a break
                    DrawFormattedText(windowPtr, textBreak, 'center', centerY*0.8, blackInt,[],[],[],2); % draw pause text
                    Screen('DrawDots', windowPtr, [centerX centerY], targetSize, xp.driftColor,[],2);
                    Screen('Flip', windowPtr, tISIonset(1) + timeStamps(iLeg).transDur - slack); % flip one second after final fixation
                    Eyelink('Message', sprintf('leg %i block %i paused at %f', iLeg, iBlock, GetSecs));
                    fprintf('Subject is taking a short break\n');
                    KbStrokeWait;
                    fprintf('Subject resumed leg %i block %i\n', iLeg, iBlock);
                    
                    Screen('DrawDots', windowPtr, [centerX centerY], targetSize, xp.targetColor,[],2); % draw the middle dot
                    tISIonset(1) = Screen('Flip', windowPtr);
                    Eyelink('Message', sprintf('leg %i block %i resumed at %f', iLeg, iBlock, GetSecs));
                end
                 
            end %trial loop
            
            if iBlock < xp.nBlocks(iLeg)
                breakText = textBlock;
                fprintf('Block %i of leg %i finished. Subject can take a longer break now.\n', iBlock, iLeg);
            elseif iBlock == xp.nBlocks(iLeg) && iLeg < xp.nLegs
                breakText = textLeg;
                fprintf('Leg %i finished. Please go in to the subject room!\n', iLeg);
            else
               breakText = textEnd;
               fprintf('Experiment completed.\n');
            end
            
             DrawFormattedText(windowPtr, breakText, 'center', centerY*0.8, blackInt); % draw pause text
             % select data from saccades to lateral targets only for feedback
             latFBdata = data(iLeg).saccLatency(iBlock,:,1); 
             accFBdata = data(iLeg).saccAccuracy(iBlock,:,1);
             FB = sprintf(textFB, mean(latFBdata(~isnan(latFBdata))), mean(accFBdata(~isnan(accFBdata))./(xp.screenRes(1)/xp.screenDim(1))*10) ); %format feedback text
             DrawFormattedText(windowPtr, FB, 'center', centerY*0.6, blackInt,[],[],[],2); % draw feedback text
             Screen('Flip', windowPtr, tISIonset(1) + timeStamps(iLeg).transDur - slack);
             Eyelink('Message', sprintf('leg %i block %i stopped at %f', iLeg, iBlock, GetSecs));
             
                %save back-up of data so far
                filename = [xp.backupFolder xp.codename '_' xp.subject '_' xp.tDCS '_' xp.legNames{iLeg} '_' ...
                    datestr(now, 'yyyy-mm-dd_HH-MM-SS') '_leg_' int2str(iLeg) '_block_' int2str(iBlock) '_trial_' int2str(iTrial) '.mat'];
                save(filename,'data', 'xp', 'timeStamps')
                
                % save ET-file
                sacctDCS_ELsave(xp,EDFname,iBlock,iLeg)
                
                KbStrokeWait;
                fprintf('Subject indicated he/she''s ready for calibration\n');
                 
        end % block loop
    startAtBlock = 1; % even if the current leg was started at a differen block, for the next leg, should always start at block 1 again
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

