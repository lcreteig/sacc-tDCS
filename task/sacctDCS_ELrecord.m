function EDFname = sacctDCS_ELrecord(windowPtr, xp, leg, block)

el = EyelinkInitDefaults(windowPtr); % set standard parameters

%modify some to match present experiment
el.msgfontcolour    = BlackIndex(xp.screenNum);
el.imgtitlecolour   = BlackIndex(xp.screenNum);
el.calibrationtargetcolour = WhiteIndex(xp.screenNum);
el.backgroundcolour = round(GrayIndex(xp.screenNum));
el.calibrationtargetsize = 2;
el.calibrationtargetwidth = 0.5;

EyelinkUpdateDefaults(el);

% Open EDF file for recording data from Eyelink
% Filename can be 8 characters (letters or numbers) max!
EDFname = sprintf('%sL%iB%i.edf', xp.subject, leg, block);
stat = Eyelink('Openfile', EDFname);
if stat~=0
    error('Cannot create EDF file');
end
Eyelink('command', ['add_file_preamble_text ''' xp.codename '_'  xp.subject '_' num2str(leg) '_' num2str(block) '''']); %write experiment details to start of file

% Calibrate the eye tracker
EyelinkDoTrackerSetup(el);

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