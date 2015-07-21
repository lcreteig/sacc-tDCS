function [EDFname] = sacctDCS_ELconfig(windowPtr, xp, leg, block, stimCoords)

el = EyelinkInitDefaults(windowPtr); % set standard parameters

%modify some to match present experiment
el.msgfontcolour    = BlackIndex(xp.screenNum);
el.imgtitlecolour   = BlackIndex(xp.screenNum);
el.calibrationtargetcolour = WhiteIndex(xp.screenNum);
el.backgroundcolour = round(GrayIndex(xp.screenNum));

EyelinkUpdateDefaults(el);
if ~EyelinkInit(0)
    error('Eyelink initialization aborted.\n');
end

[~, vs]  = Eyelink('GetTrackerVersion');
fprintf('Running experiment on a ''%s'' tracker.\n', vs);

% Set some parameters of the eyetracker through the link
Eyelink('command', ['add_file_preamble_text ''' xp.codename '_'  xp.subject '_' num2str(leg) '_' num2str(block) '''']); %write experiment details to start of file
Eyelink('Command', 'simulation_screen_distance = %d', 100)
Eyelink('command', 'enable_automatic_calibration = YES'); %automatically move to next calibration point (this can become unset for some reason sometimes)
Eyelink('command', 'calibration_type = HV9'); %set 9-point calibration
Eyelink('command', 'active_eye = LEFT'); %record left eye
Eyelink('command', 'recording_parse_type = GAZE'); % use gaze data to process events
%determine thresholds for on-line parsing of eye data into events
Eyelink('command', 'set_parser_configuration = 0'); %0 matches "cognitive" configuration in manual (more conservative); 1 is psychophysical (more sensitive)
Eyelink('command', 'pupil_size_diameter = YES'); %to convert pupil area to diameter
%"file_sample_data" specifies what types of samples will be written to the EDF file
Eyelink('command', 'file_sample_data = LEFT,RIGHT,GAZE,GAZERES,AREA,STATUS'); %DEFAULT: LEFT,RIGHT,GAZE,GAZERES,AREA,STATUS
%"file_event_data" sets data in events written to the EDF file
%"file_event_filter" specifies what types of events will be written to the EDF file
Eyelink('command', 'file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON'); %DEFAULT: LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON

%"link_sample_data" specifies which types of samples will be available online through the link
%"link_event_data" sets data in events sent through the link
%"link_event_filter" specifies which types of samples will be available online through the link

% Open EDF file for recording data from Eyelink
% Filename can be 8 characters (letters or numbers) max!
EDFname = sprintf('%sL%iB%i.edf', xp.subject, leg, block);
stat = Eyelink('Openfile', EDFname);
if stat~=0
    error('Cannot create EDF file');
end

% Calibrate the eye tracker
EyelinkDoTrackerSetup(el);

% Draw shapes to host pc at saccade target locations, for online viewing of gaze cursor
Eyelink('command', 'set_idle_mode'); % must be offline to draw to EyeLink screen
Eyelink('command', 'clear_screen 7'); % clear tracker display
Eyelink('command', 'draw_cross %d %d 0', stimCoords(1,1), stimCoords(1,2)); %draw left target location
Eyelink('command', 'draw_cross %d %d 0', stimCoords(2,1), stimCoords(2,2)); %draw center target location
Eyelink('command', 'draw_cross %d %d 0', stimCoords(3,1), stimCoords(3,2)); %draw right target location

% Send screen info
Eyelink('message', 'DISPLAY_COORDS %ld %ld %ld %ld', 0, 0, xp.screenRes(1)-1, xp.screenRes(2)-1);
Eyelink('message', sprintf('degrees per pixel %f',  2 * tan(1/2 * pi/180) * xp.screenDist * (xp.screenRes(1) / xp.screenDim(1))));

% Start recording eye position 
if Eyelink('IsConnected')~=1 % Make sure we're still connected.
    error('Eyelink connection lost')
end;
Eyelink('command', 'set_idle_mode'); %short pause so that the tracker can finish the mode transition; might crash otherwise
WaitSecs(0.05);
Eyelink('StartRecording');
WaitSecs(0.1); % record a few samples before start displaying, otherwise a few msec of data might be lost
Eyelink('message', 'start recording Eyelink'); % mark start of recording in data file

end


%BELOW IS FROM change.M AND EyelinkEventExample.m

% we are disabling samples transfered to the display but we're still
% getting events. This speeds things up a little since in this
% example we do not care about samples but about events from the
% tracker
%Eyelink('StartRecording',1,1,0,1);

%   % check for endsaccade events
%         if Eyelink('isconnected') == el.dummyconnected % in dummy mode use mousecoordinates
%             [x,y,button] = GetMouse(window);
%             evt.type=el.ENDSACC;
%             evt.genx=x;
%             evt.geny=y;
%             evtype=el.ENDSACC;
%         else % check for events
%             evtype=Eyelink('getnextdatatype');
%         end
%         if evtype==el.ENDSACC		% if the subject finished a saccade check if it fell on an object
%             if Eyelink('isconnected') == el.connected % if we're really measuring eye-movements
%                 evt = Eyelink('getfloatdata', evtype); % get data
%             end
%             % check if saccade landed on an object
%             choice=-1;
%             noobject=0;
%             i=1;
%             while 1
%                 if 1==IsInRect(evt.genx,evt.geny, object(i).rect )
%                     choice=i;
%                     break;
%                 end
%                 i=i+1;
%                 if i>length(object)
%                     noobject=1;
%                     break;
%                 end
%             end
%             if lastchoice>0 && (choice~=lastchoice || noobject==1) % toggle object color
%                 if object(lastchoice).on==1 % restore screen
%                     Screen('CopyWindow', buffer, window, object(lastchoice).rect, object(lastchoice).rect);
%                     object(lastchoice).on=0;
%                     lastchoice=-1;
%                     doflip=1;
%                 end
%             end
%             if choice>0 && choice~=lastchoice % toggle object color
%                 if object(choice).on==0 % toggle object on screen
%                     Screen('CopyWindow', altbuffer, window, object(choice).rect, object(choice).rect);
%                     object(choice).on=1;
%                     doflip=1;
%                 end
%                 lastchoice=choice;
%             end
%             if doflip==1
%                 Screen('Flip',  window, [], 1);
%                 doflip=0;
%             end
%         end % saccade?


% if dummymode==0
%     error=Eyelink('CheckRecording');
%     if(error~=0)
%         disp('Error in Recording');
%         break;
%     end
%     % we need to loop over this a few times ( 30 is
%     % randomly chosen) so that we do not miss any events
%     % and to prevent any buffer overflow
%     for j=1:30
%         evtype = Eyelink('GetNextDataType');
%         if evtype == el.FIXUPDATE
%             if Eyelink('isconnected') == el.connected % if we're really measuring eye-movements
%                 evt = Eyelink('getfloatdata', evtype);% get data
%                 
%                 % only process if its the desired eye
%                 if evt.eye == eye_used
%                     % send msg with details of fixation
%                     % update event
%                     Eyelink('message', 'Fixupdate: avg_x %d, y %d, dur %d',evt.gavx, evt.gavy, evt.entime-evt.sttime);
%                     
%                     % determine if gaze values are within
%                     % interest region and if gaze has been
%                     % maintained over 300 ms. This method
%                     % allows for saccades as long as they
%                     % are withing interest area
%                     if infixationWindow(evt.gavx,evt.gavy)
%                         totalFixTime = totalFixTime + 50;
%                         if totalFixTime >= 300
%                             break;
%                         end
%                     else % broke fixation reset time
%                         totalFixTime = 0;
%                     end
%                 end
%             else
%                 disp('Eyelink disconnected!');
%             end
%         end
%     end %end for
% else