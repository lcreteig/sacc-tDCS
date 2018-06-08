These functions run the experimental task using [Psychtoolbox-3](http://psychtoolbox.org/) from within MATLAB.

# Main task files

* `sacctDCS_setup.m`: Top-level file that allows for experimenter inputs (participant ID, whether eyetracking is used, etc.). Run only this file to start the task.
* `sacctDCS_getParams.m`: Definitions of experimental parameters (stimulus size, trial counts, etc). Called in `sacctDCS_setup.m`
* `sacctDCS_Main_ET.m`: Main code for running the task, with Eye Tracking. Called in `sacctDCS_setup.m`
* `sacctDCS_Main_noET.m`: Main code for running the task, without Eye Tracking (for testing purposes). Called in `sacctDCS_setup.m`

`dva2pix.m` converts Degrees of Visual Angle (used for specifying task parameters in `sacctDCS_getParams.m`) to Pixels.

# Eye tracker control files

All of these are called from within `sacctDCS_Main_ET.m`

* `sacctDCS_ELconfig.m`: Commands for setting up the Eyelink with the desired parameters from within MATLAb
* `sacctDCS_ELrecord.m`: Commands for opening data files and starting the recording
* `sacctDCS_ELsave.m`: Commands for saving the recorded data
* `sacctDCS_ELwrapup.m`: Commands for resetting the Eyelink parameters to how they were before `sacctDCS_ELconfig.m` was run
