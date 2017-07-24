function [pixX,pixY] = dva2pix(dvaX, dvaY, screenRes, screenDim, screenDist, screenNum)
%
% DVA2PIX Convert degrees of visual angle (dva) to pixels. If the screen
% resolution (screenRes) or screen dimensions (screenDim) are not specified,
% DVA2PIX will try to gather this information using Psychtoolbox if it is
% installed.
%
%
% Usage:
% [pixX, pixY] = dva2pix(dvaX, dvaY)
% [pixX pixY] = dva2pix(dvaX, dvaY, [screenRes], [screenDim], [screenDist], [screenNum])
%
% Inputs:
%
% -DVAX:            requested stimulus width in degrees of visual angle.
%                   If only desiring an x-value, dvaY can be left empty and
%                   it's corresponding output pixY can be supressed.
%
% -DVAY:            requested stimulus height in degrees of visual angle.
%                   If only desiring a y-value, dvaX can be left empty and
%                   it's corresponding output pixX can be supressed.
%
% -SCREENRES:       screen resolution in pixels [width, height].
% (optional)        Default (when input is []): uses the Psychtoolbox
%                   command "Screen('Resolution, screenNum)" to obtain the
%                   resolution.
%
% -SCREENDIM:       screen dimension in centimeters [width, height].
%                   Default (when input is []): uses the Psychtoolbox
%                   command "Screen('DisplaySize', screenNum)" to obtain
%                   the dimensions. Note that this output is not guaranteed
%                   to be correct, so always physically measure yourself
%                   and/or look at manufacturer's the website!
%
% -SCREENDIST:      distance from subject to screen in centimeters.
%  (optional)       Default (when input is []) is 90 cm.
%
% -SCREENNUM:       pointer to the monitor as defined in Psychtoolbox.
% (optional)        Default (when input is []): uses the Psychtoolbox
%                   command "max(Screen('Screens'))" to obtain the screen
%                   number.
%
% Outputs:
%
% -PIXX:            stimulus width in pixels.
% -PIXY:            stimulus height in pixels.

%% Parse input

narginchk(1,6)
nargoutchk(0,2)

if      ( ~exist('screenDim','var') || isempty(screenDim) ) && ...
        ( ~exist('screenRes','var') || isempty(screenRes) ) ...
        
        try
        Psychtoolboxroot
        catch err
            if strcmp(err.identifier, 'MATLAB:UndefinedFunction')
                error('MATLAB:dva2pix', ['Looks like Psychtoolbox is not installed' ...
                    'on this machine\n. In that case, specify the screenRes and screenDim' ...
                    'arguments!'])
            end
        end
end

if ~exist('dvaX','var') || isempty(dvaX)
    dvaX = NaN; % if not specified, assign something so the code doesn't jam;
end

if ~exist('dvaY','var') || isempty(dvaY)
    dvaY = NaN; % if not specified, assign something so the code doesn't jam;
end

if      ( ~exist('screenNum','var') || isempty(screenNum) ) && ...
        ( ~exist('screenDim','var') || isempty(screenDim) ) && ...
        ( ~exist('screenRes','var') || isempty(screenRes) ) ...
        
    screenNum = max(Screen('Screens'));
end

if ~exist('screenDist','var') || isempty(screenDist)
    screenDist = 90;
end

if ~exist('screenRes','var') || isempty(screenRes)
    resolution = Screen('Resolution',screenNum);
    screenRes = [resolution.width resolution.height]; % convert to input format of this function
end

if ~exist('screenDim','var') || isempty(screenDim)
    [screenDim(1), screenDim(2)] = Screen('DisplaySize',screenNum);
    screenDim = screenDim/10; % convert to cm, since original output is in mm
end

%% Process

pixPerCm = screenRes ./ screenDim;

cmX = 2 * tan(dvaX/2 * pi/180) * screenDist;
cmY = 2 * tan(dvaY/2 * pi/180) * screenDist;

pixX = round(cmX .* pixPerCm(1));
pixY = round(cmY .* pixPerCm(2));

