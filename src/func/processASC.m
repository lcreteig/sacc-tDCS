function saccData = processASC(ASCfile,MATfile,leg,block,microSaccCrit)
% PROCESSEDF: Extract saccade measures from Eyelink Data Files that were 
% converted to ascii text
%
% Usage: saccData = processEDF(ASCfile,MATfile,leg,block,microSaccCrit)
%
% Inputs:
%
% -ASCfile: string with full path to .asc file containing eye data
% -MATfile: string with full path to .mat file with experiment metadata
% -leg: string specifying task section with respect to stimulation 
%       ('pre', 'tDCS', 'post') 
% -block: number of task block (each recalibration of the eyelink)
% -microSaccCrit: number in degrees of visual angle; saccades with an 
%                 amplitude smaller than this criterion will be marked as 
%                 microsaccades, and not extracted for further analysis.
%
% Outputs: 
%
% -saccData: structure with trial metadata and saccade measures

%% Setup variables

load(MATfile); % load the .mat file containing experiment metadata
fid = fopen(ASCfile, 'r'); %open ascii for read access

%Initialize answer variables
saccData.trialNum = zeros(2*xp.nTrials,1);
saccData.trialNum(1:2:end) = 1:xp.nTrials; saccData.trialNum(2:2:end) = 1:xp.nTrials;
saccData.lateral = repmat([1 0]', xp.nTrials,1); % code saccades away from center as 1 (always first saccade), saccades toward center as 0 (always second saccade)
saccData.direction = zeros(2*xp.nTrials,1);
saccData.direction(1:2:end) = data(strcmpi(leg, xp.legNames)).targetSide(block,:); % fill in direction for saccades away from center (1 = to the right, -1 = to the left)
saccData.direction(2:2:end) = -1*data(strcmpi(leg, xp.legNames)).targetSide(block,:); % direction for the following saccade (towards center) is always opposite
saccData.amplitude = nan(2*xp.nTrials,1);
saccData.startDev = nan(2*xp.nTrials,1);
saccData.endDev.x = nan(2*xp.nTrials,1);
saccData.endDev.y = nan(2*xp.nTrials,1);
saccData.latency = nan(2*xp.nTrials,1);
saccData.drift.x = nan(2*xp.nTrials,1);
saccData.drift.y = nan(2*xp.nTrials,1);


%Initialize helper variables
phaseStart = false;
blinkStart = false;
breakStart = false;
rowCounter = 1;

%Stimuli info
centerDot = xp.screenRes/2; % [x y] coordinates of middle target/fixation point (= center of monitor);
targetEcc = dva2pix(xp.targetEcc,[], xp.screenRes, xp.screenDim, xp.screenDist); %how many x-pixels the targets were away from center

%% Parse the text file; extract saccade measures

while ~feof(fid) %loop untill we reach the end of the file
    tline = fgetl(fid); % read line, move pointer to next line
    
    if regexp(tline, 'leg \d block \d started|resumed') % if current line marks end of a break
        breakStart = true;
    end
    
    if breakStart & regexp(tline, '^EFIX') % if this is the first "end fixation" event
       EFIXdata = sscanf(tline, ['%*s %*s' repmat('%f', 1,6)]); %skip the first two strings, then read the next 6 numbers
       breakStart = false; % reset to be able to write out data for next break
    end
    
    if regexp(tline, 'trial \d+ phase \d') % if current line marks the onset of a stimulus
        if phaseStart % if we already passed a line like this, apparently there was no saccade for this phase
            blinkStart = false; % reset the blink memory variable
            rowCounter = rowCounter+1; %increment counter, to write data for new trial/phase
        else
            phaseStart = true; % set variable to remember this later
        end
         stimOnset = sscanf(tline, '%*s%i',1); %read only the first integer you encounter (which is the timestamp)
    end
    
    if phaseStart & regexp(tline, '^SBLINK') % if current line marks a blink
        blinkStart = true;  % set variable to remember this later
    end
    
    if phaseStart & regexp(tline, '^ESACC') & ~blinkStart  % if current line marks the end of a saccade, which was not preceded by a blink
        ESACCdata = sscanf(tline, ['%*s %*s' repmat('%f', 1,9)]); %skip the first two strings, then read the next nine numbers
        
        %Calculate saccadic amplitude.
        %a = subtract the x-coordinates of the starting and ending points of the saccade;
        %b = subtract the y-coordinates of the starting and ending points of the saccade;
        %c^2 = a^2 + b^2; saccadic amplitude = square root of c.
        amplPix = sqrt((ESACCdata(4) - ESACCdata(6))^2 + (ESACCdata(5) - ESACCdata(7))^2);
        amplDVA = pix2dva(amplPix, [], xp.screenRes, xp.screenDim, xp.screenDist); % convert from pixels to degrees of visual angle
        
        if amplDVA < microSaccCrit % if the current saccade is a microsaccade
            continue % move on with the loop (read next line of data file)
        else % if this is the large saccade we're looking for, analyze it
            
            % store the amplitude
            saccData.amplitude(rowCounter) = amplDVA;
            
            % Compute x- and y- offset between average fixation position and
            % actual stimulus location, to optionally do drift correction later.
            saccData.drift.x(rowCounter) = pix2dva(EFIXdata(4)-centerDot(1),[],xp.screenRes, xp.screenDim, xp.screenDist);
            [~,saccData.drift.y(rowCounter)] = pix2dva([],EFIXdata(5)-centerDot(2),xp.screenRes, xp.screenDim, xp.screenDist);
            
            %Determine the [x y] coordinates where the present saccade should end (i.e., where the dot was presented)
            %and where the present saccade should have started (i.e., where participants should have been fixating)
            diffFromCenter = saccData.direction(rowCounter)*targetEcc; % distance from target to fixation for the current trial
            startPoint = [centerDot(1) - ~saccData.lateral(rowCounter)*diffFromCenter, centerDot(2)];
            endPoint = [centerDot(1) + saccData.lateral(rowCounter)*diffFromCenter, centerDot(2)];
            
            % Calculate vector / Euclidian distance between actual saccade start point and stimulus coordinates 
            % (i.e. how much participants were off when fixating).
            saccData.startDev(rowCounter) = pix2dva(sqrt((startPoint(1) - ESACCdata(4))^2 + (startPoint(2) - ESACCdata(5))^2),[], ...
                xp.screenRes, xp.screenDim, xp.screenDist);
            
            % Calculate x- and y- coordinates of deviation between saccade end point and stimulus coordinates. 
            % Note that here we don't calculate the vector yet, so we can still correct it / do different kinds of analyses
            saccData.endDev.x(rowCounter) = pix2dva(ESACCdata(6)-endPoint(1),[],xp.screenRes, xp.screenDim, xp.screenDist);
            [~,saccData.endDev.y(rowCounter)] = pix2dva([],ESACCdata(7)-endPoint(2), xp.screenRes, xp.screenDim, xp.screenDist);
            
            saccData.latency(rowCounter) = ESACCdata(1)-stimOnset; %calculate saccade latency: difference between starting times of saccade and stimulus
            
            phaseStart = false;
            blinkStart = false;
            rowCounter = rowCounter+1;
        end
        
    end

end
fclose(fid); %close the asci file
