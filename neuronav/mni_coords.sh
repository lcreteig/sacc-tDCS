#!/bin/sh

#SUBJECT FOLDER   SCAN

# S01: FFM_03     FFM_03-20150504-0005-WIPsT13DTFE_P25_S2_3mSENSE.nii.gz
# S02: P03_sc1    P03_sc1503_K1-20120514-0003-WIPsT13DTFE_P25_S2_3mSENSE.nii.gz
# S03: Unclear markers
# S04: DL         DL_140514_10.nii.gz
# S05: MB_012     MB_012-20141118-0011-WIPsT13DTFE3mSENSE.nii.gz
# S06: BM         BM_T1.nii.gz
# S07: NA         NA_T1.nii.gz
# S08: DE         DE_T1.nii.gz
# S09: JW         T1.nii.gz
# S10: sc2658     sc2658-20140203-0003-WIPsT13DTFE_P25_S2_3mSENSE.nii
# S11: AU         o20140312_160938WIPsT13DTFE_P25_S2_3mSENSEs003a001.nii.gz
# S12: JS         T1.nii.gz
# S13: LN         LN.nii
# S14: P06_sc1507 P06_sc1507_K1-20120515-0003-WIPsT13DTFE_P25_S2_3mSENSE.nii.gz
# S15: MB         t1_1mm.nii.gz
# S16: SOAT_21    SOAT_21-20161214-0003-WIPT1w_1mm_4m40_S1SENSE.nii
# S17: SOAT_33    SOAT_33-20161230-0003-WIPT1w_1mm_4m40_S1SENSE.nii
# S18: SOAT_11    SOAT_11-20161208-0003-WIPT1w_1mm_4m40_S1SENSE.nii
# S19: SOAT_26    SOAT_26-20161215-0003-WIPT1w_1mm_4m40_S1SENSE.nii
# S20: pp20_PH_1  PP20_PH_1correct-20170117-0006-WIPsT1W_3D_TFE_32channelSENSE.nii

# S01 and S16 are to be excluded!

homeDir=/Volumes/research$/reteig/sacc-tDCS/neuronav # working directory

# list of folders in directory
declare -a folders=("P03_sc1"
"DL"
"MB_012"
"BM"
"NA"
"DE"
"JW"
"sc2658"
"AU"
"JS"
"LN"
"P06_sc1507"
"MB"
"SOAT_33"
"SOAT_11"
"SOAT_26"
"pp20_PH_1"
)

# get length the array
numsubs=${#folders[@]}

# list of scans (1 per folder)
declare -a scans=("P03_sc1503_K1-20120514-0003-WIPsT13DTFE_P25_S2_3mSENSE"
"DL_140514_10"
"MB_012-20141118-0011-WIPsT13DTFE3mSENSE"
"BM_T1_reorient"
"NA_T1_reorient"
"T1_reorient"
"T1_reorient"
"sc2658-20140203-0003-WIPsT13DTFE_P25_S2_3mSENSE"
"o20140312_160938WIPsT13DTFE_P25_S2_3mSENSEs003a001"
"T1_reorient"
"LN"
"P06_sc1507_K1-20120515-0003-WIPsT13DTFE_P25_S2_3mSENSE"
"t1_1mm"
"SOAT_33-20161230-0003-WIPT1w_1mm_4m40_S1SENSE"
"SOAT_11-20161208-0003-WIPT1w_1mm_4m40_S1SENSE"
"SOAT_26-20161215-0003-WIPT1w_1mm_4m40_S1SENSE"
"PP20_PH_1correct-20170117-0006-WIPsT1W_3D_TFE_32channelSENSE"
)

for (( i=0; i<${numsubs}; i++ )); # loop over subjects
do
  cd ${homeDir}/${folders[$i]} # change directory to current folder
  T1=${scans[$i]} # name of current scan

  # BET (strip the skull)
  bet $T1 ${T1}_brain.nii.gz

  # FLIRT (register to MNI)
  flirt -in ${T1}_brain.nii.gz -ref ${FSLDIR}/data/standard/MNI152_T1_1mm_brain -out ${T1}_brain_toMNI1mm.nii.gz -omat native2mNI.mat

  # Convert coord to MNI
  cat FEF_nativeCoords.txt | img2stdcoord -img $T1 -std ${FSLDIR}/data/standard/MNI152_T1_1mm -xfm native2mNI.mat >> ${homeDir}/FEF_MNIcoords.txt

done
