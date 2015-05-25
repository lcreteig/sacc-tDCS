function sacctDCS_setup

prompt={'Subject ID:', 'Task?(1=practice, 2=main)', 'Screen distance (cm):', ...
    'Environment (1=lab, 2=mac, 3=pc)', 'Eyetracker? (1=yes, 0=no)'};
name='sacctDCS';
numlines=1;
defaultanswer={'sID','1','62','1','1'};
answer=inputdlg(prompt,name,numlines,defaultanswer);


%%
environment = {'lab', 'mac', 'pc'};
task = {'practice', 'main'};

xp = sacctDCS_getParams(answer{1},environment{str2double(answer{4})},task{str2double(answer{2})},str2double(answer{3}));

if str2double(answer{5})
else
    placeHolderFlag = false;
    overlap = 0.100;
    [data,timeStamps] = sacctDCS_Main_noET(xp,placeHolderFlag,overlap);
end

assignin('base', 'xp', xp);
assignin('base', 'data', data);
assignin('base', 'timeStamps', timeStamps);