%% Merge all data from S18 session K into one file again

%% Prepare master file
rootFolder = fullfile('/Volumes','students$','reteig students','FEF-tDCS','data','S18','session_K_orig');
load(fullfile(rootFolder,'post','sacc-tDCS_S18_K_main.mat'));
timeStampsAll = timeStamps;
dataAll = data;

clear timeStamps data

%% Leg 1, block 1
leg = 1;
origBlock = 1;
newBlock = 1;

folderName = 'pre_block1';
fileName = 'sacc-tDCS_S18_K_practice.mat';

[timeStampsAll,dataAll] = replaceData(leg,origBlock,newBlock,rootFolder,folderName,fileName,timeStampsAll,dataAll);

%% Leg 1, block 2
leg = 1;
origBlock = 2;
newBlock = 2;

folderName = 'pre_block2_3-tDCS_block1_2';
fileName = 'sacc-tDCS_S18_K_pre_2017-04-04_16-50-28_leg_1_block_2_trial_120.mat';

[timeStampsAll,dataAll] = replaceData(leg,origBlock,newBlock,rootFolder,folderName,fileName,timeStampsAll,dataAll);

%% Leg 1, block 3
leg = 1;
origBlock = 3;
newBlock = 3;

folderName = 'pre_block2_3-tDCS_block1_2';
fileName = 'sacc-tDCS_S18_K_pre_2017-04-04_16-54-41_leg_1_block_3_trial_120.mat';

[timeStampsAll,dataAll] = replaceData(leg,origBlock,newBlock,rootFolder,folderName,fileName,timeStampsAll,dataAll);

%% Leg 2, block 1
leg = 2;
origBlock = 2;
newBlock = 1;

folderName = 'pre_block2_3-tDCS_block1_2';
fileName = 'sacc-tDCS_S18_K_tDCS_2017-04-04_17-02-54_leg_2_block_2_trial_120.mat';

[timeStampsAll,dataAll] = replaceData(leg,origBlock,newBlock,rootFolder,folderName,fileName,timeStampsAll,dataAll);

%% Leg 2, block 2
leg = 2;
origBlock = 3;
newBlock = 2;

folderName = 'pre_block2_3-tDCS_block1_2';
fileName = 'sacc-tDCS_S18_K_tDCS_2017-04-04_17-06-53_leg_2_block_3_trial_120.mat';

[timeStampsAll,dataAll] = replaceData(leg,origBlock,newBlock,rootFolder,folderName,fileName,timeStampsAll,dataAll);

%% Leg 2, block 3
leg = 2;
origBlock = 1;
newBlock = 3;

folderName = 'tDCS_block3';
fileName = 'sacc-tDCS_S18_K_tDCS_2017-04-04_17-13-26_leg_2_block_1_trial_120.mat';

[timeStampsAll,dataAll] = replaceData(leg,origBlock,newBlock,rootFolder,folderName,fileName,timeStampsAll,dataAll);

%% Save
timeStamps = timeStampsAll;
data = dataAll;
save(fullfile('/Volumes','students$','reteig students','FEF-tDCS','data','S18','sacc-tDCS_S18_K_main.mat'),'timeStamps','data','xp');

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