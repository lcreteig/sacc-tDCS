%% Parameters

summary_func = @nanmedian; % calculate mean or median of saccade measures

microSaccCrit = 1.5; % saccades with smaller amplitude are considered microsaccades and not extracted for analysis
crit.saccStart = 1.8; % exclude saccades whose start point is more than this degrees away from fixation
crit.saccEnd = 8; % exclude saccades whose end point is more than this degrees away from the target
crit.saccLatencyFast = 50; % exclude saccades with latency faster than 50 ms
crit.saccLatencySlow = 400; % exclude saccades with latency slower than 50 ms

%% I/O

if ismac % if we're working on a mac
    driveLetter = '/Volumes/research$/'; % set appropriate file prefix
    converterName = 'edf2asc-mac'; % name of the utility that converts files from .edf (European Data Format) to .asc (ASCII plain text)
elseif ispc % if we're working on a pc
    driveLetter = 'Y:\'; % assume BigBrother is mapped to the Z drive
    converterName = 'edf2asc.exe';
end

dataDir = fullfile(driveLetter, 'reteig', 'sacc-tDCS', 'data', filesep); % directory where data is stored
srcDir = fullfile(driveLetter, 'reteig', 'sacc-tDCS', 'src', filesep); % directory where code is located
addpath(genpath(srcDir)); % add the scripts in this directory to Matlab's path
converterPath = fullfile(srcDir, 'bin', converterName);
converterInput = '-ns'; % optional flag to .edf file converter; '-ns' means Not to include Samples but only events to keep the file size down

subFolders = dir(fullfile(dataDir, 'S*'));
rawFolder = 'raw'; % sub-directory containing raw Eyelink Data Files (.edf) and .mat files for each subject
subjects = {subFolders.name};
stimulation = {'anodal', 'cathodal'}; % type of stimulation over the FEF
tDCScodes = {'C', 'I'; 'H', 'A'}; % in subject 1-10, setting 'C' on the device was mapped to anodal; for cathodal or 'I' the wires were reversed 
% In subject 11-15, 'A' was mapped to cathodal; for anhodal or 'H' the wires were reversed 
legs = {'pre', 'tDCS', 'post'};

fid = fopen(fullfile(dataDir, ['saccLatency' '_' func2str(summary_func) '.txt']),'w'); % open text file for writing
printHeaders = true; % headers still need to be printed
%% Process files

for iSub = [subjects(1) subjects] % loop over first subject twice, just to print headers to data table
    
    if printHeaders
        fprintf(fid,'subject\t');  % Write out subject column header
    else
        fprintf(fid,'%s\t',iSub{:}); % Write row header
    end
    
    % get correspondance between stimulation type and tDCS setting
    if str2double(iSub{:}(2:end)) <= 10 % for the first 10 subjects
        mapping = 1;
    elseif str2double(iSub{:}(2:end)) > 10 && str2double(iSub{:}(2:end)) <= 15 % for subjects 11-15
        mapping = 2;
    end
    
    for iStim = tDCScodes(mapping,:)
        
        MATfile = fullfile(dataDir, iSub{:}, rawFolder, ['sacc-tDCS_' iSub{:} '_' iStim{:} '_main.mat']); % filename of the accompanying .mat file, containing info about the experiment
        stimType = stimulation(strcmp(iStim, tDCScodes(mapping,:))); % get name of stimulation type used in current session, for printing data table
        
        for iLeg = legs
            
            files2load = dir(fullfile(dataDir, iSub{:}, rawFolder, ['sacc-tDCS_' iSub{:} '_' iStim{:} '_' iLeg{:} '*.edf'])); % names of files from all blocks
            
            for iBlock = 1:length(files2load)
                
                if printHeaders % if headers are still to be printed
                    % print headers
                    fprintf(fid, '%s\t', [stimType{:} '_' iLeg{:} '_' num2str(iBlock) '_' 'lateral' '_' 'left']);
                    fprintf(fid, '%s\t', [stimType{:} '_' iLeg{:} '_' num2str(iBlock) '_' 'lateral' '_' 'right']);
                    fprintf(fid, '%s\t', [stimType{:} '_' iLeg{:} '_' num2str(iBlock) '_' 'center' '_' 'left']);
                    fprintf(fid, '%s\t', [stimType{:} '_' iLeg{:} '_' num2str(iBlock) '_' 'center' '_' 'right']);
                    continue % skip to next column
                end
                
                fileName = files2load(iBlock).name; %file name of current block
                EDFfile = fullfile(dataDir, iSub{:}, rawFolder, fileName); % full path to .edf file
                
                %% Convert EDF (if necessary)
                
                ASCfile = fullfile(dataDir, iSub{:}, [fileName(1:end-4) '_' converterInput '.asc']); % full path to .asc file
                if ~exist(ASCfile, 'file') % if there does not exist an ASCII version of this file yet
                    %call the converter utility to make one (executes in terminal / DOS); remove samples ('ns') and only keep events
                    system(['"' converterPath '"' ' ' converterInput ' ' '"' EDFfile '"']);
                    movefile([EDFfile(1:end-4), '.asc'], ASCfile); % move the newly created file out of the raw directory; change name (not sure if this can be done with converter directly)
                end
                
                %% Parse the text file and extract saccade measures
                
                saccData = processEDF(ASCfile,MATfile,iLeg{:},iBlock,microSaccCrit);
                
                %% Clean data
                
                saccStartIdx = saccData.startDev < crit.saccStart; % trials with saccade start point smaller than criterion
                saccEndIdx = saccData.endDev < crit.saccEnd; % trials with saccade end point smaller than criterion
                saccLatencyIdx = saccData.latency > crit.saccLatencyFast & saccData.latency < crit.saccLatencySlow; % trials with saccade latency neither too fast nor too slow
                
                corrSaccs = saccStartIdx & saccEndIdx & saccLatencyIdx; % combine all to get indices of trials for analysis
                fprintf('%i saccades of %i (%i %%) included for analysis in file %s\n', sum(corrSaccs), length(corrSaccs), round(sum(corrSaccs)/length(corrSaccs)*100), fileName);
                                
                %% Print data to table
                
                fprintf(fid,'%g\t', summary_func(saccData.latency(corrSaccs & saccData.lateral == 1 & saccData.direction == -1))); % saccades to a LATERAL target, LEFT of current fixation position
                fprintf(fid,'%g\t', summary_func(saccData.latency(corrSaccs & saccData.lateral == 1 & saccData.direction == 1)));  % saccades to a LATERAL target, RIGHT of current fixation position
                fprintf(fid,'%g\t', summary_func(saccData.latency(corrSaccs & saccData.lateral == 0 & saccData.direction == -1))); % saccades to a CENTRAL target, LEFT of current fixation position
                fprintf(fid,'%g\t', summary_func(saccData.latency(corrSaccs & saccData.lateral == 0 & saccData.direction == 1))); % saccades to a CENTRAL target, RIGHT of current fixation position
                
            end
        end
    end
    printHeaders = false; % stop printing headers after first outer loop
    fprintf(fid,'\n'); % go to next line (subject)
end
fclose(fid);