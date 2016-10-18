function [dvaX,dvaY] = pix2dva(pixX, pixY, screenRes, screenDim, screenDist, screenNum)
%
% PIX2DVA Convert pixels to degrees of visual angle (dva). If the screen
% resolution (screenRes) or screen dimensions (screenDim) are not specified,
% PIX2DVA will try to gather this information using Psychtoolbox if it is
% installed.
%
%
% Usage:
% [dvaX,dvaY] = pix2dva(pixX, pixY)
% [dvaX,dvaY] = pix2dva(pixX, pixY, screenRes, screenDim, screenDist, screenNum)
%
% Inputs:
%
% -pixX:            requested stimulus width in pixels.
%                   If only desiring an x-value, pixY can be left empty and
%                   it's corresponding output dvaY can be supressed.
%
% -pixY:            requested stimulus height in pixels.
%                   If only desiring a y-value, pixX can be left empty and
%                   it's corresponding output dvaX can be supressed.
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
% -DVAX:            stimulus width in degrees of visual angle.
% -DVAY:            stimulus height in degrees of visual angle.

%% Parse input

error(nargchk(1,6,nargin))
error(nargoutchk(0,2,nargout))

if      ( ~exist('screenDim','var') || isempty(screenDim) ) && ...
        ( ~exist('screenRes','var') || isempty(screenRes) ) ...
        
        try
        Psychtoolboxroot
        catch err
            if strcmp(err.identifier, 'MATLAB:UndefinedFunction')
                error('MATLAB:pix2dva', ['Looks like Psychtoolbox is not installed' ...
                    'on this machine\n. In that case, specify the screenRes and screenDim' ...
                    'arguments!'])
            end
        end
end

if ~exist('pixX','var') || isempty(pixX)
    pixX = NaN; % if not specified, assign something so the code doesn't jam;
end

if ~exist('pixY','var') || isempty(pixY)
    pixY = NaN; % if not specified, assign something so the code doesn't jam;
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

cmX = pixX ./ pixPerCm(1); 
cmY = pixY ./ pixPerCm(2);

dvaX = ( 2 * atan((cmX/2) /screenDist) * (180/pi) );
dvaY = ( 2 * atan((cmY/2) /screenDist) * (180/pi) );

