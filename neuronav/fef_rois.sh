#!/bin/sh
# First argument: input csv file with MNI coordinates (as created by mni_coords.sh)
# Second argument: output NIfTI file with FEF ROIs for each subject
#
# e.g.
# bash fef_rois.sh FEF_coords_MNI.csv FEF_ROIs.nii.gz

inFile=$1
outFile=$2

# Create mask of MNI template filled with zeroes
fslmaths ${FSLDIR}/data/standard/MNI152_T1_1mm_brain -mul 0 ${outFile}

COUNTER=1 # subject number

# read all columns of file, line by line
# "|| [ -n "$MNI_Z" ]" makes sure last line is also read even though there is no newline character after last column
while IFS=";", read -r subject folder scan MNI_X MNI_Y MNI_Z || [ -n "$MNI_Z" ]
do
  subjectFolder=neuronav/${folder}
  writeFolder=${subjectFolder}/registration

  # Get MNI coordinates in voxel space (instead of mm), for defining ROIs
  read -r vox_X vox_Y vox_Z <<<$(echo "$MNI_X $MNI_Y $MNI_Z" | img2stdcoord -img ${FSLDIR}/data/standard/MNI152_T1_1mm_brain -mm)

  # Define mask for each subject
  fslmaths ${FSLDIR}/data/standard/MNI152_T1_1mm_brain -mul 0 -add 1 -roi $vox_X 1 $vox_Y 1 $vox_Z 1 0 1 ${writeFolder}/${folder}_FEF_vox.nii.gz -odt float # Create mask with FEF voxel set to "1"
  fslmaths ${writeFolder}/${folder}_FEF_vox.nii.gz -kernel sphere 1 -fmean ${writeFolder}/${folder}_FEF_sphere.nii.gz -odt float # inflate to sphere with radius 1 mm
  fslmaths ${writeFolder}/${folder}_FEF_sphere.nii.gz -bin ${writeFolder}/${folder}_FEF_sphere.nii.gz #binarize (set all values to 1)
  fslmaths ${writeFolder}/${folder}_FEF_sphere.nii.gz -mul $COUNTER ${writeFolder}/${folder}_FEF_sphere.nii.gz -odt float # multiply values in sphere with subject number

  # Add mask for this subject to group file
  fslmaths ${outFile} -add ${writeFolder}/${folder}_FEF_sphere.nii.gz ${outFile}

COUNTER=$[$COUNTER +1]
done < <(tail -n +2 $inFile) # skip reading the header line
