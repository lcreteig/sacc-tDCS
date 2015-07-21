function [xp] = sacctDCS_getParams(subjectID,tDCScode,environment,task,screenDist)

% GETPARAMS_SACCTDCS Set the eXperimental Parameters for sacc-tDCS.

% leonreteig@gmail.com

%ADD elconfig
%ADD eyetracking
%ADD gaze contingent

%% I/O

switch environment
    
      case 'L1.09'
        xp.dataFolder = fullfile('D:', 'USERS', 'Leon', 'sacc-tDCS', 'data', filesep);
        xp.backupFolder = fullfile('D:', 'USERS', 'Leon', 'sacc-tDCS', 'data', 'backup', filesep);
        xp.taskFolder = fullfile('D:', 'USERS', 'Leon', 'sacc-tDCS', 'task', filesep);
        addpath(genpath(xp.taskFolder));
    
    case 'L1.01'
        xp.dataFolder = fullfile('/Users', 'test', 'Documents', 'Leon', 'sacc-tDCS', 'data', filesep);
        xp.backupFolder = fullfile('/Users', 'test', 'Documents', 'Leon', 'sacc-tDCS', 'data', 'backup', filesep);
        xp.taskFolder = fullfile('/Users', 'test', 'Documents', 'Leon', 'sacc-tDCS', 'task', filesep);
        addpath(genpath(xp.taskFolder));
        
    case 'pc'
        driveletter = 'Y:'; % default for PCs
        xp.dataFolder = fullfile(driveletter, 'reteig', 'sacc-tDCS', 'data', filesep);
        xp.backupFolder = fullfile(driveletter, 'reteig', 'sacc-tDCS', 'data', 'backup', filesep);
        xp.taskFolder = fullfile(driveletter, 'reteig',  'sacc-tDCS', 'task', filesep);
        
    case 'mac'
        driveletter = '/Volumes/research$/'; % default for macs
        xp.dataFolder = fullfile(driveletter, 'reteig', 'sacc-tDCS', 'data', filesep);
        xp.backupFolder = fullfile(driveletter, 'reteig', 'sacc-tDCS', 'data', 'backup', filesep);
        xp.taskFolder = fullfile(driveletter, 'reteig',  'sacc-tDCS', 'task', filesep);
end

xp.environment = environment;
xp.experiment = 'pro-saccade task for frontal eye field tDCS';
xp.codename = 'sacc-tDCS'; %prepended to all saved files
xp.subject = subjectID;
xp.tDCS = tDCScode;
xp.task = task;
xp.date = datestr(now); %current date and time

%% Experiment

switch task
    
    case 'practice' 
        xp.nLegs = 1;
        xp.legNames = {'practice'};
        xp.nBlocks = 1;
        xp.nTrials = 140;
        xp.breaksPerBlock = 1;
        
    case 'main' 
        xp.nLegs = 3;
        xp.legNames = {'pre'; 'tDCS'; 'post'};
        xp.nBlocks = [3 3 6]; % blocks per leg
        xp.nTrials = 140; % trials per block
        xp.breaksPerBlock = 1; %number of short, timed breaks in each block
end

%% Screen

switch environment
    
       case 'L1.09'
         xp.screenNum = 2;
         xp.screenDim = [51 28.5]; % screen dimensions [width height] in cm.
         xp.screenRes = [1920 1080]; % screen dimensions [width height] in pixels.
         xp.screenRefresh = 120;
    
    case 'L1.01'
        xp.screenNum = 0;
        xp.screenDim = [47.4 29.6]; % screen dimensions [width height] in cm.
        xp.screenRes = [1600 1200]; % screen dimensions [width height] in pixels.
        xp.screenRefresh = 60;
        
    case {'pc', 'mac'}
        xp.screenNum = max(Screen('Screens')); % select most outer screen
        [screenDim(1), screenDim(2)] = Screen('DisplaySize',xp.screenNum); % query screen size
        xp.screenDim = screenDim/10; % convert to cm (original is mm)
        defSettings = Screen('Resolution', xp.screenNum); %query screen res, store result
        xp.screenRes = [defSettings.width defSettings.height];
        xp.screenRefresh = defSettings.hz;
               
end

xp.screenDist = screenDist; % subject-screen distance in cm.
Screen('Resolution',xp.screenNum,xp.screenRes(1),xp.screenRes(2),xp.screenRefresh); % set monitor to desired settings

%% Stimuli

% N.B. all of the size parameters should be specified in degrees of visual
% angle (except where noted). They are converted to pixels in the main code
% using the accompanying 'dva2pix' function.

%%%TEXT%%%
xp.FontSize = 18; % in points
xp.Font = 'Arial';
xp.TextStyle = 0; %0 = normal, 1 = bold, 2 = italic.

%%%TARGETS%%%
xp.targetEcc = 8; % displacement along horizontal meridian
xp.targetSize = 0.5;
xp.targetColor = 0; % single contrast value (0-255)

%%%PLACEHOLDERS%%%
xp.placeSize = xp.targetSize*1.25;
xp.placeColor = xp.targetColor;

%% Timing

% Presentation times in seconds
xp.saccadeTime = 0.300;
xp.fixTime = [0.300 0.700]+xp.saccadeTime; % [lowerbound upperbound], each trial is random number in between
xp.transTime = 1.000;% fixation at first and last trial of every block

%% Keys
xp.abortKey = 'ESCAPE';

