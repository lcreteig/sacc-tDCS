function [ELdefaults, ELconfig] = sacctDCS_ELconfig(stimCoords)

ELconfig = EyelinkInitDefaults; % set standard parameters

% Initialize connection with Eyelink
if ~EyelinkInit(0)
    error('Eyelink initialization aborted.\n');
end

% First query some EL parameters as they are now, and store 
[~, ELdefaults.config] = Eyelink('ReadFromTracker', 'elcl_select_configuration');
[~, ELdefaults.sRate] = Eyelink('ReadFromTracker', 'sample_rate');
[~, ELdefaults.autoCalib] = Eyelink('ReadFromTracker', 'enable_automatic_calibration');
[~, ELdefaults.calibType] = Eyelink('ReadFromTracker', 'calibration_type');
[~, ELdefaults.eye] = Eyelink('ReadFromTracker', 'active_eye');
[~, ELdefaults.parseType] = Eyelink('ReadFromTracker', 'recording_parse_type');
[~, ELdefaults.parseConfig] = Eyelink('ReadFromTracker', 'select_parser_configuration');
[~, ELdefaults.filter] = Eyelink('ReadFromTracker', 'heuristic_filter');
[~, ELdefaults.pupilTrack] = Eyelink('ReadFromTracker', 'use_ellipse_fitter');
[~, ELdefaults.pupilType] = Eyelink('ReadFromTracker', 'pupil_size_diameter');
[~, ELdefaults.fileSampleData] = Eyelink('ReadFromTracker', 'file_sample_data');
[~, ELdefaults.fileEventFilter] = Eyelink('ReadFromTracker', 'file_event_filter');

% Now set these parameters to the desired values through the link

%Eyelink('command', 'simulation_screen_distance = %d', xp.screenDist*10); or properly set PHYSICAL.INI instead
Eyelink('command', 'elcl_select_configuration = MTABLER'); % use the desktop level configuration (with illuminator on the right)
Eyelink('command', 'sample_rate = 1000'); % sample at 1k Hz
Eyelink('command', 'enable_automatic_calibration = YES'); %automatically move to next calibration point (this can become unset for some reason sometimes)
Eyelink('command', 'calibration_type = HV9'); %set 9-point calibration
Eyelink('command', 'active_eye = RIGHT'); % or track the dominant eye instead
Eyelink('command', 'recording_parse_type = GAZE'); % use gaze data to process events
Eyelink('command', 'select_parser_configuration = 0'); %determine thresholds for on-line parsing of eye data into events. %0 matches "cognitive" configuration in manual (more conservative); 1 is psychophysical (more sensitive)
Eyelink('command', 'heuristic_filter = 1 1'); % set filtering for file and link data to STD level 
Eyelink('command', 'use_ellipse_fitter = YES'); %to convert pupil area to diameter
Eyelink('command', 'pupil_size_diameter = YES'); %to convert pupil area to diameter

%"file_sample_data" specifies what types of samples will be written to the EDF file
Eyelink('command', 'file_sample_data = LEFT,RIGHT,GAZE,GAZERES,AREA,STATUS'); %DEFAULT: LEFT,RIGHT,GAZE,GAZERES,AREA,STATUS
%"file_event_data" sets data in events written to the EDF file
%"file_event_filter" specifies what types of events will be written to the EDF file
Eyelink('command', 'file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON'); %DEFAULT: LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON
%"link_sample_data" specifies which types of samples will be available online through the link
%"link_event_data" sets data in events sent through the link
%"link_event_filter" specifies which types of samples will be available online through the link

% Draw shapes to host pc at saccade target locations, for online viewing of gaze cursor
Eyelink('command', 'set_idle_mode'); % must be offline to draw to EyeLink screen
Eyelink('command', 'clear_screen 7'); % clear tracker display, fill with color
Eyelink('command', 'draw_cross %d %d 0', stimCoords(1,1), stimCoords(1,2)); %draw left target location
Eyelink('command', 'draw_cross %d %d 0', stimCoords(2,1), stimCoords(2,2)); %draw center target location
Eyelink('command', 'draw_cross %d %d 0', stimCoords(3,1), stimCoords(3,2)); %draw right target location