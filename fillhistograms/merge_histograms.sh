#!/bin/bash
filepath="/afs/cern.ch/user/x/xirong/CMSSW_14_1_9_patch2/src/residualanalysis/Outputs/L3FillHistOutputs/0914Histograms_allPD_rougherbins/"
outfilepath="/afs/cern.ch/user/x/xirong/CMSSW_14_1_9_patch2/src/residualanalysis/Outputs/L3FillHistOutputs/0914MergedHists_allPD_rougherbins/"

#hadd -f "$outfilepath/AK4_photonjet_Data.root" "$filepath"/AK4_photonjet_HP*.root
hadd -f "$outfilepath/Ak4_photonjet_MC.root" "$filepath"/AK4_photonjet_MC*.root
