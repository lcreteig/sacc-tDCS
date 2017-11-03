#!/bin/sh
# First argument: input csv file with native coordinates
# Second argument: output csv file with MNI coordinates
# Third argument (optional): "all" to also (redo) strip skull (BET) and register to standard (FLIRT)
# e.g.
# bash mni_coords.sh FEF_coords_native.csv FEF_coords_MNI.csv all

inFile=$1
outFile=$2
bet_flirt=$3

echo "subject;folder;scan;MNI_X;MNI_Y;MNI_Z" > $outFile # header for new file

# read all columns of file, line by line
# "|| [ -n "$native_Z" ]" makes sure last line is also read even though there is no newline character after last column
while IFS=";", read -r subject folder scan native_X native_Y native_Z || [ -n "$native_Z" ]
do
    subjectFolder=neuronav/${folder}
    readFolder=${subjectFolder}/from_neuronav
    writeFolder=${subjectFolder}/registration

    mkdir -p "$writeFolder" # if it doesn't exist yet, create directory

    if [ $bet_flirt = "all" ]; then

      # Reorient to match orientation of MNI template (if not already)
      fslreorient2std ${readFolder}/${scan}.nii ${writeFolder}/${scan}_reorient.nii.gz

      # BET (strip the skull)
      bet ${writeFolder}/${scan}_reorient.nii.gz ${writeFolder}/${scan}_reorient_brain.nii.gz

      # FLIRT (register to MNI)
      flirt -in ${writeFolder}/${scan}_reorient_brain.nii.gz -ref ${FSLDIR}/data/standard/MNI152_T1_1mm_brain -out ${writeFolder}/${scan}_reorient_brain_toMNI1mm.nii.gz -omat ${writeFolder}/${folder}_native2mNI.mat
    fi

    # Convert coordinates to MNI
    MNI_mm=$(echo "$native_X $native_Y $native_Z" | img2stdcoord -img ${writeFolder}/${scan}_reorient.nii.gz -std ${FSLDIR}/data/standard/MNI152_T1_1mm -xfm ${writeFolder}/${folder}_native2mNI.mat)

    # Append MNI coordinates to file
    echo "$subject;$folder;$scan;${MNI_mm//  /;}" >> $outFile # replace (double) spaces in MNI coords with ";"

done < <(tail -n +2 $inFile) # skip reading the header line
