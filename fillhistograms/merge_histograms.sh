#!/bin/bash
filepath="/eos/cms/store/group/phys_heavyions/xirong/L3Residuals/Histograms/"
outfilepath="/eos/cms/store/group/phys_heavyions/xirong/L3Residuals/MergedHistograms/"

hadd -f "$outfilepath/AK4_photojet_HP0_0827Merged.root" "$filepath"/AK4_photonjet_HP0_batch*_of_46.root
hadd -f "$outfilepath/Ak4_photojet_MC_0827Merged.root" "$filepath"/AK4_photonjet_MC_batch*_of_7.root
