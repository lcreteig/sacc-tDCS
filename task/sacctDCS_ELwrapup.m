function sacctDCS_ELwrapup(ELdefaults)

Eyelink('command', 'clear_screen 0'); % clear tracker display

% Reset parameters that were changed to their initial value (see sacctDCS_ELconfig)

Eyelink('command', ['elcl_select_configuration = ' ELdefaults.config]);
Eyelink('command', ['sample_rate = ' ELdefaults.sRate]);
Eyelink('command', ['enable_automatic_calibration = ' ELdefaults.autoCalib]);
Eyelink('command', ['calibration_type = ' ELdefaults.calibType]); 
Eyelink('command', ['active_eye = ' ELdefaults.eye]); 
Eyelink('command', ['recording_parse_type = ' ELdefaults.parseType]); 
Eyelink('command', ['select_parser_configuration = ' ELdefaults.parseConfig]);
Eyelink('command', ['heuristic_filter = ' ELdefaults.filter]); 
Eyelink('command', ['use_ellipse_fitter = ' ELdefaults.pupilTrack]);
Eyelink('command', ['pupil_size_diameter = ' ELdefaults.pupilType]);
Eyelink('command', ['file_sample_data = ' ELdefaults.fileSampleData]); 
Eyelink('command', ['file_event_filter = ' ELdefaults.fileEventFilter]);