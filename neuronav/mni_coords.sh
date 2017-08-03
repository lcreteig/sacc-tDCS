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
# S16: SOAT_21    SOAT_21-20161214-0003-WIPT1w_1mm_4m40_S1SENSE.nii.gz
# S17: SOAT_33    SOAT_33-20161230-0003-WIPT1w_1mm_4m40_S1SENSE.nii
# S18: SOAT_11    SOAT_11-20161208-0003-WIPT1w_1mm_4m40_S1SENSE.nii
# S19: Unclear markers
# S20: pp20_PH_1  PP20_PH_1correct-20170117-0006-WIPsT1W_3D_TFE_32channelSENSE.nii
# S21: sub-026    sub-026Active-CurioVal-20170421-0002-WIPT1anat_acq-3m_T1wSENSE.nii.gz
# S22: sub-024    sub-024Active-CurioVal-20170421-0002-WIPT1anat_acq-3m_T1wSENSE.nii.gz
# S23: experiment aborted
# S24: sub-018    sub-018Active-CurioVal-20170411-0002-WIPT1anat_acq-3m_T1wSENSE.nii.gz
# S25: sub-052    sub-052Passive-CurioVal-20170430-0002-WIPT1anat_acq-3m_T1wSENSE.nii.gz
# S26: sub-003    sub-003Active-Curioval-20170407-0002-WIPT1anat_acq-3m_T1wSENSE.nii.gz
# S27: sub-040    sub-040Passive-CurioVal-20170424-0002-WIPT1anat_acq-3m_T1wSENSE.nii.gz
# S28: sub-051    sub-051Passive-CurioVal-20170425-0002-WIPT1anat_acq-3m_T1wSENSE.nii.gz
# S29: BCS034     BCS034_UvA_T1.nii
# S30: BCS038     BCS038_UvA_T1.nii
# S31: did not show up
# S32: sub-010    sub-010Active-Curioval-20170409-0002-WIPT1anat_acq-3m_T1wSENSE.nii.gz
# S33: sub-016    sub-016Active-CurioVal-20170411-0002-WIPT1anat_acq-3m_T1wSENSE.nii.gz

homeDir=/Volumes/research$/reteig/sacc-tDCS/neuronav # working directory

# list of folders in directory
declare -a folders=(
"FFM_03"
"P03_sc1"
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
"SOAT_21"
"SOAT_33"
"SOAT_11"
"pp20_PH_1"
"sub-026"
"sub-024"
"sub-018"
"sub-052"
"sub-003"
"sub-040"
"sub-051"
"BCS034"
"BCS038"
"sub-010"
"sub-016"
)

# get length the array
numsubs=${#folders[@]}

# list of scans (1 per folder)
declare -a scans=(
"FFM_03-20150504-0005-WIPsT13DTFE_P25_S2_3mSENSE.nii.gz"
"P03_sc1503_K1-20120514-0003-WIPsT13DTFE_P25_S2_3mSENSE"
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
"SOAT_21-20161214-0003-WIPT1w_1mm_4m40_S1SENSE"
"SOAT_33-20161230-0003-WIPT1w_1mm_4m40_S1SENSE"
"SOAT_11-20161208-0003-WIPT1w_1mm_4m40_S1SENSE"
"PP20_PH_1correct-20170117-0006-WIPsT1W_3D_TFE_32channelSENSE"
"sub-026Active-CurioVal-20170421-0002-WIPT1anat_acq-3m_T1wSENSE"
"sub-024Active-CurioVal-20170421-0002-WIPT1anat_acq-3m_T1wSENSE"
"sub-018Active-CurioVal-20170411-0002-WIPT1anat_acq-3m_T1wSENSE"
"sub-052Passive-CurioVal-20170430-0002-WIPT1anat_acq-3m_T1wSENSE"
"sub-003Active-Curioval-20170407-0002-WIPT1anat_acq-3m_T1wSENSE"
"sub-040Passive-CurioVal-20170424-0002-WIPT1anat_acq-3m_T1wSENSE"
"sub-051Passive-CurioVal-20170425-0002-WIPT1anat_acq-3m_T1wSENSE"
"BCS034_UvA_T1"
"BCS038_UvA_T1"
"sub-010Active-Curioval-20170409-0002-WIPT1anat_acq-3m_T1wSENSE"
"sub-016Active-CurioVal-20170411-0002-WIPT1anat_acq-3m_T1wSENSE"
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
