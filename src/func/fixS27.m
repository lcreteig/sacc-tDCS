%% Merge all data from S27 into one file again

%% Prepare master file for session B
clear all
rootFolder = 'FILL IN PATH HERE';
sessionFolder = fullfile(rootFolder,'data','S27','session_B_orig');
load(fullfile(sessionFolder,'sacc-tDCS_S27_B_main.mat'));
timeStampsAll = timeStamps;
dataAll = data;

clear timeStamps data

%% Leg 1, block 1
leg = 1;
origBlock = 1;
newBlock = 1;

folderName = 'pre_block1';
fileName = 'sacc-tDCS_S27_B_practice.mat';

[timeStampsAll,dataAll] = replaceData(leg,origBlock,newBlock,sessionFolder,folderName,fileName,timeStampsAll,dataAll);

%% Save file for session B
timeStamps = timeStampsAll;
data = dataAll;
save(fullfile('/Volumes','students$','reteig students','FEF-tDCS','data','S27','sacc-tDCS_S27_B_main.mat'),'timeStamps','data','xp');

%% Prepare master file for session K
clear all
sessionFolder = fullfile(rootFolder,'data','S27','session_K_orig');
load(fullfile(sessionFolder,'sacc-tDCS_S27_K_main.mat'));
timeStampsAll = timeStamps;
dataAll = data;

clear timeStamps data

%% Leg 1, blocks 1,2,3 and leg 2, blocks 1,2,3
for iLeg = 1:2
    for iBlock = 1:3
leg = iLeg;
origBlock = iBlock;
newBlock = iBlock;

folderName = 'pre-tDCS';
fileName = 'sacc-tDCS_S27_K_main_2017-06-16_12-53-16_leg_3_block_1_trial_1.mat';

[timeStampsAll,dataAll] = replaceData(leg,origBlock,newBlock,sessionFolder,folderName,fileName,timeStampsAll,dataAll);
    end
end

%% Save file for session K
timeStamps = timeStampsAll;
data = dataAll;
save(fullfile(rootFolder, 'data','S27','sacc-tDCS_S27_K_main.mat'),'timeStamps','data','xp');

%%
function [timeStampsAll,dataAll] = replaceData(leg,origBlock,newBlock,rootFolder,folderName,fileName,timeStampsAll,dataAll)
load(fullfile(rootFolder,folderName,fileName));
timeFields = fieldnames(timeStampsAll);
for k = 1:length(timeFields)
    val = timeStamps(leg).(timeFields{k});
            if isnumeric(val) && ~isscalar(val) 
            timeStampsAll(leg).(timeFields{k})(newBlock,:) = val(origBlock,:);
            end
end

dataFields = fieldnames(dataAll);
for k = 1:length(dataFields)
    val = data(leg).(dataFields{k});
            if isnumeric(val) && ~isscalar(val) 
            dataAll(leg).(dataFields{k})(newBlock,:,:) = val(origBlock,:,:);
            end
end

end