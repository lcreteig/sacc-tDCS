function [data,timeStamps] = sacctDCS_Main_noET(xp,placeHolderFlag,overlap)

try
    %% Setup
    
    if nargin < 2
        placeHolderFlag = true;
    elseif nargin < 3
        overlap = 0; % in seconds
    end
    
    for i = 1:length(xp.nLegs)
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
    textBreak = 'You may take a short break now. Press any key to resume';
    textLeg = 'Please wait for the experimenter.';
    textEnd = 'Experiment complete!';
    
    %%%TARGETS%%%
    targetSize = dva2pix(xp.targetSize,[],xp.screenRes,xp.screenDim,xp.screenDist); %convert target size to pixels
    targetEcc = dva2pix(xp.targetEcc,[],xp.screenRes,xp.screenDim,xp.screenDist); %convert target eccentricity to pixels
    
    %%%PLACEHOLDERS%%%
    if placeHolderFlag
        placeSize = dva2pix(xp.placeSize,[],xp.screenRes,xp.screenDim,xp.screenDist);
        placeHolder(:,1) = CenterRectOnPoint([0 0 placeSize placeSize], centerX-targetEcc, centerY); % left placeholder
        placeHolder(:,2) = CenterRectOnPoint([0 0 placeSize placeSize], centerX, centerY); % center placeholder
        placeHolder(:,3) = CenterRectOnPoint([0 0 placeSize placeSize], centerX+targetEcc, centerY); % right placeholder
    end
    
    %%%TIMING%%%
    
    %initialize a structure for keeping stimulus onset timeStamps
    timeStamps.fixDur = round(xp.fixTime / ifi) * ifi;
    
    for i = 1:xp.nLegs
        timeStamps(i).fix =  zeros(xp.nBlocks(i),xp.nTrials);
        timeStamps(i).target = zeros(xp.nBlocks(i),xp.nTrials);
        timeStamps(i).leg  = xp.legNames{i};
    end
    
    if overlap
        overlap = round(overlap / ifi) * ifi;
    end
    
    %% Experiment loop
    
    DrawFormattedText(windowPtr, textStart, 'center', 'center', blackInt); % draw start text
    Screen('Flip', windowPtr); % flip to screen
    KbStrokeWait; % wait for a key press
    
    for iLeg = 1:xp.nLegs
        for iBlock = 1:xp.nBlocks(iLeg)
            
            targetSide = repmat([-1;1],xp.nTrials/2,1);
            targetSide = Shuffle(targetSide);
            
            ISI = xp.fixTime(1) + (xp.fixTime(2) - xp.fixTime(1)).*rand(xp.nTrials,2);
            ISI = round(ISI / ifi) * ifi;
            
            %Fixation
            if placeHolderFlag
                Screen('FrameRect', windowPtr, xp.placeColor, placeHolder); % draw the place holders
            end
            Screen('DrawDots', windowPtr, [centerX centerY], targetSize, xp.targetColor,[],2); % draw the middle dot
            tFixOnset = Screen('Flip', windowPtr);
            
            for iTrial = 1:xp.nTrials
                
                %Target
                if overlap
                    if placeHolderFlag
                        Screen('FrameRect', windowPtr, xp.placeColor, placeHolder);
                    end
                    Screen('DrawDots', windowPtr, [centerX centerY], targetSize, xp.targetColor,[],2); % draw a lateral dot
                    Screen('DrawDots', windowPtr, [centerX+targetSide(iTrial)*targetEcc centerY], targetSize, xp.targetColor,[],2); % draw a lateral dot
                    tTargetOnset = Screen('Flip', windowPtr, tFixOnset + ISI(iTrial,2) - slack);
                end
                
                if placeHolderFlag
                    Screen('FrameRect', windowPtr, xp.placeColor, placeHolder);
                end
                Screen('DrawDots', windowPtr, [centerX+targetSide(iTrial)*targetEcc centerY], targetSize, xp.targetColor,[],2); % draw a lateral dot
                tTargetOnset = Screen('Flip', windowPtr, tFixOnset + ISI(iTrial,2) + overlap - slack);
                
                %Fixation
                if overlap
                if placeHolderFlag
                    Screen('FrameRect', windowPtr, xp.placeColor, placeHolder);
                end
                Screen('DrawDots', windowPtr, [centerX+targetSide(iTrial)*targetEcc centerY], targetSize, xp.targetColor,[],2); % draw a lateral dot
                Screen('DrawDots', windowPtr, [centerX centerY], targetSize, xp.targetColor,[],2);
                tFixOnset = Screen('Flip', windowPtr, tTargetOnset +ISI(iTrial,1));
                end
                
                if placeHolderFlag
                    Screen('FrameRect', windowPtr, xp.placeColor, placeHolder);
                end
                Screen('DrawDots', windowPtr, [centerX centerY], targetSize, xp.targetColor,[],2);
                tFixOnset = Screen('Flip', windowPtr, tTargetOnset + ISI(iTrial,1) + overlap - slack);
                
                %Check for keypresses
                [~,~,keyCode] = KbCheck;
                if keyCode(KbName(xp.abortKey)) % abort the experiment
                    error('Manual quit')
                end
                
                %Store data
                data(iLeg).targetSide(iBlock,iTrial) = targetSide(iTrial);
                timeStamps(iBlock,iTrial).fix = tFixOnset;
                timeStamps(iBlock,iTrial).target = tTargetOnset;
                
            end %trial loop
            
            if iBlock < xp.nBlocks(iLeg)
            
            DrawFormattedText(windowPtr, textBreak, 'center', 'center', blackInt); % draw pause text
            Screen('Flip', windowPtr);
            %save back-up of data so far
            filename = [xp.backupFolder xp.codename '_' xp.subject '_' xp.tDCS '_' xp.task '_' ...
                datestr(now, 'yyyy-mm-dd_HH-MM-SS') '_leg_' int2str(iLeg) '_block_' int2str(iBlock) '_trial_' int2str(iTrial) '.mat'];
            save(filename,'data', 'xp', 'timeStamps')
            KbStrokeWait;
            
            end
            
        end % block loop
        
        DrawFormattedText(windowPtr, textLeg, 'center', 'center', blackInt); % draw pause text
        Screen('Flip', windowPtr);
        KbStrokeWait;
        
    end % leg loop
    
    DrawFormattedText(windowPtr, textEnd, 'center', 'center', blackInt); % draw pause text
    Screen('Flip', windowPtr);
    pause(5);
    
    %% Save data and close
    filename = [xp.dataFolder xp.codename '_' xp.subject '_' xp.tDCS '_' xp.task '.mat'];
    save(filename,'data', 'xp', 'timeStamps');
    
    sca; %Screen('CloseAll'), to regain control of monitor
    ShowCursor; % show mouse again
    Priority(0); % restore priority
    
catch err
    
    %save back-up of data so far
    filename = [xp.backupFolder xp.codename '_' xp.subject '_' xp.tDCS '_' xp.task '_' ...
                datestr(now, 'yyyy-mm-dd_HH-MM-SS') '_leg_' int2str(iLeg) '_block_' int2str(iBlock) '_trial_' int2str(iTrial) '.mat'];
    save(filename,'data', 'xp', 'timeStamps')
    
    sca;
    ShowCursor; 
    Priority(0);
    rethrow(err); % issue the caught error
end

