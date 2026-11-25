#!/bin/bash
filepath="../Outputs/L3FillHistOutputs/1123_PtRef50to600_PVfilter/"
#outfilepath="/afs/cern.ch/user/x/xirong/CMSSW_14_1_9_patch2/src/residualanalysis/Outputs/L3FillHistOutputs/0914MergedHists_allPD_rougherbins/"
outfilepath="$filepath/Merged/"

hadd -f "$outfilepath/AK4_photonjet_Data.root" "$filepath"/AK4_photonjet_HP*.root
hadd -f "$outfilepath/Ak4_photonjet_MC.root" "$filepath"/AK4_photonjet_MC*.root
