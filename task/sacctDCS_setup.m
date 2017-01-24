function sacctDCS_setup

prompt={'Subject ID:', 'tDCS code:', 'Task?(1=practice, 2=main)', 'Screen distance (cm):', ...
    'Environment (1=L1.09, 2=L1.01, 3=mac, 4=pc)', 'Eyetracker? (1=yes, 0=no)', 'Start at leg:', 'Start at block:'};
name='sacctDCS';
numlines=1;
defaultanswer={'sID','A','1','73','1','1','1','1'};
answer=inputdlg(prompt,name,numlines,defaultanswer);


%%
environment = {'L1.09', 'L1.01', 'mac', 'pc'};
task = {'practice', 'main'};

xp = sacctDCS_getParams(answer{1},answer{2},environment{str2double(answer{5})},task{str2double(answer{3})},str2double(answer{4}));

startAtLeg = str2double(answer{7});
startAtBlock = str2double(answer{8});

placeHolderFlag = false;
overlap = 0; 
expISI = 1;

if str2double(answer{6})
    [data,timeStamps] = sacctDCS_Main_ET(xp,startAtLeg,startAtBlock);
else
    [data,timeStamps] = sacctDCS_Main_noET(xp,placeHolderFlag,overlap,expISI,startAtLeg,startAtBlock);
end

assignin('base', 'xp', xp);
assignin('base', 'data', data);
assignin('base', 'timeStamps', timeStamps);