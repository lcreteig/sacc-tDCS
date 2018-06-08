# `/func`

Functions or scripts particular to this project:

* `fixS18.m`: Quick MATLAB script to fix the raw data files for S18.
* `fixS27.m`: Quick MATLAB script to fix the raw data files for S27.
* `mni_coords.sh`: BASH shell script that computes MNI coordinates for the frontal eye field coordinates in native space. Called in the `notebooks/frontal_eye_field.md` R notebook.
* `fef_rois.sh`: BASH shell script that creates regions of interest for the frontal eye field of each participant in MNI space. Was originally used for visualizing frontal eye field location, but now deprecated in favor of `plot_coords.m`. Called in the `notebooks/frontal_eye_field.md` R notebook.
* `plot_coords.m` MATLAB script and function (latter retrieved from <https://github.com/rordenlab/spmScripts/blob/master/nii_coord2spheres.m>) to visualize MNI frontal eye field coordinates of each participant with [Surf Ice](https://www.nitrc.org/projects/surfice/) software.
* `processASC.m`: MATLAB function that extracts saccade measures from the `.asc` files. Called in the main `analysis.m` script.

# `/bin`

Binary code (only executable, not source). These files are proprietary and should be obtained from [SR Research](http://download.sr-support.com/dispdoc/page25.html) directly.

* `edf2asc.exe`: Windows Utility for converting Eyelink `.edf` files as they come from the eye tracker to plain-text `.asc` files.
* `edf2asc-mac`: same as `edf2asc.exe`, but for macs.

# `/lib`

Source code _not_ particular to this project:

* `dva2pix.m`: MATLAB function for converting Degrees of Visual Angle to Pixels.
* `pix2dva.m`: MATLAB function for converting Pixels to Degrees of Visual Angle.
* `inclusionBF.R`: R function for calculating inclusion Bayes factors, as in the [JASP software](https://jasp-stats.org/)
