%% Parameters
microSaccCrit = 1.5; % saccades with smaller amplitude are considered microsaccades and not extracted for analysis

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

fid = fopen(fullfile(dataDir, ['saccLatency' '.csv']),'w'); % open text file for writing
% print column headers
fprintf(fid, '%s,', 'subject');
fprintf(fid, '%s,', 'stimulation');
fprintf(fid, '%s,', 'leg');
fprintf(fid, '%s,', 'block');
fprintf(fid, '%s,', 'trial');
fprintf(fid, '%s,', 'type');
fprintf(fid, '%s,', 'direction');
fprintf(fid, '%s,', 'deviation.start');
fprintf(fid, '%s,', 'deviation.end');
fprintf(fid, '%s,', 'amplitude');
fprintf(fid, '%s,', 'latency');
fprintf(fid,'\n'); % go to next line

%% Process files

for iSub = [subjects(1) subjects] % loop over first subject twice, just to print headers to data table

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
                
                % initialize cells for saccade info
                saccType = cell(length(saccData.lateral),1);
                saccDirection = cell(length(saccData.direction),1);
                
                % convert numerical condition code to strings
                [saccType{saccData.lateral == 1}] = deal('lateral');
                [saccType{saccData.lateral == 0}] = deal('center');
                [saccDirection{saccData.direction == -1}] = deal('left');
                [saccDirection{saccData.direction == 1}] = deal('right');
                
                %% Print data to table

                for iTrial = 1:length(saccData.trialNum)

                    fprintf(fid,'%s,%s,%s,%i,%i,%s,%s,%g,%g,%g,%g\n', ...
                        iSub{:}, ... % print subject ID
                        stimType{:}, ... % print stimulation type
                        iLeg{:}, ... % print pre/during/post tDCS
                        iBlock, ... % print block number
                        saccData.trialNum(iTrial), ... % print current trial number
                        saccType{iTrial}, ... % print whether saccade was to side or center
                        saccDirection{iTrial}, ... % print whether saccade direction was right or left
                        saccData.startDev(iTrial), ... % print deviation of saccade start point to fixation
                        saccData.endDev(iTrial), ... % print deviation from saccade end point to target location
                        saccData.amplitude(iTrial), ... % print saccade amplitude (length of straight line between start and end point)
                        saccData.latency(iTrial) ... % print saccade latency (time from target onset to start of saccade)
                    );

                end
   
            end
        end
    end
 
end
fclose(fid);
