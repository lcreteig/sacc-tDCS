# `/func`

Functions or scripts particular to this project:

* `fixS18.m`: Quick MATLAB script to fix the raw data files for S18. `See data/S18/raw/notes/S18_notes.md` for details.
* `fixS27.m`: Quick MATLAB script to fix the raw data files for S27. `See data/S27/raw/notes/S27_notes.md` for details.
* `mni_coords.sh`: BASH shell script that computes MNI coordinates for the frontal eye field coordinates in native space. Loads `data/FEF_coords_native.csv` and creates `data/FEF_coords_MNI.csv`. Called in the `notebooks/frontal_eye_field.Rmd` R notebook.
* `fef_rois.sh`: BASH shell script that creates regions of interest for the frontal eye field of each participant in MNI space. Was originally used for visualizing frontal eye field location, but now deprecated in favor of `plot_coords.m` Called in the `notebooks/frontal_eye_field.Rmd` R notebook.
* `plot_coords.m` MATLAB script and function (latter retrieved from <https://github.com/rordenlab/spmScripts/blob/master/nii_coord2spheres.m>) to visualize MNI frontal eye field coordinates of each participant with [Surf Ice](https://www.nitrc.org/projects/surfice/) software. Loads `data/FEF_coords_MNI_excl.csv` and creates `data/FEF.node`.
* `processASC.m`: MATLAB function that extracts saccade measures from the `.asc` files (Eyelink event data; see `data/README.md`). Called in the main `analysis.m` script.

# `/bin`

Binary code (only executable, not source). These files are proprietary and should be obtained from [SR Research](http://download.sr-support.com/dispdoc/page25.html) directly.

* `edf2asc.exe`: Windows Utility for converting Eyelink `.edf` files as they come from the eye tracker to plain-text `.asc` files.
* `edf2asc-mac`: same as `edf2asc.exe`, but for macs.

# `/lib`

Source code _not_ particular to this project:

* `dva2pix.m`: MATLAB function for converting Degrees of Visual Angle to Pixels.
* `pix2dva.m`: MATLAB function for converting Pixels to Degrees of Visual Angle.
* `inclusionBF.R`: R function for calculating inclusion Bayes factors, as in the [JASP software](https://jasp-stats.org/)
