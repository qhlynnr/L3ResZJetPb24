###
#!/bin/sh

#MCFile=/eos/cms/store/group/phys_heavyions/xirong/L3Residuals/MergedHistograms/Ak4_photojet_MC_0827Merged.root
#DataFile=/eos/cms/store/group/phys_heavyions/xirong/L3Residuals/MergedHistograms/AK4_photojet_HP0_0827Merged.root


##MCFile=/eos/cms/store/group/phys_heavyions/bharikri/JetMinPOG/L3ResPhotonJet/2026_05_26_QCDPhoton_V1L2Relative.root
##DataFile=/eos/cms/store/group/phys_heavyions/bharikri/JetMinPOG/L3ResPhotonJet/2026_05_26_2024ppRefHP_full_V1L2Relative.root

MCFile=/afs/cern.ch/user/x/xirong/CMSSW_14_1_9_patch2/src/residualanalysis/Outputs/L3FillHistOutputs/0907MergedHists/AK4_photonjet_MC.root
DataFile=/afs/cern.ch/user/x/xirong/CMSSW_14_1_9_patch2/src/residualanalysis/Outputs/L3FillHistOutputs/0907MergedHists/AK4_photonjet_HP0.root
OutFolderName=/afs/cern.ch/user/x/xirong/CMSSW_14_1_9_patch2/src/residualanalysis/Outputs/DeriveL3Outpus/
mkdir -p $OutFolderName
OutFileName=$OutFolderName/0908DeriveL3.root
root -l -b -q "deriveL3_from_photonjet.C(\"$MCFile\",\"$DataFile\",\"$OutFileName\",true,5,true,false)"
