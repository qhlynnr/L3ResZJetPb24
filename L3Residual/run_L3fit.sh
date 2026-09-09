###
#!/bin/sh

#InFile=/eos/cms/store/group/phys_heavyions/xirong/L3Residuals/0831L3DerivedHistogram.root
#OutFolder=/eos/cms/store/group/phys_heavyions/xirong/L3Residuals/0831Fits/

InFile=/afs/cern.ch/user/x/xirong/CMSSW_14_1_9_patch2/src/residualanalysis/Outputs/DeriveL3Outpus/0909DeriveL3.root
OutFolder=/afs/cern.ch/user/x/xirong/CMSSW_14_1_9_patch2/src/residualanalysis/Outputs/DoFitOutput/0909DoFit_newrange/
mkdir -p OutFolder
OutFile=L3_0901_Fits_
root -l -b -q "dofits_L3.C(\"$InFile\",
               30,
               700,
               \"$OutFile\",
	       \"$OutFolder\",
               true,
               \"2024ppRef\",
               \"pp Reference\",
               true,
               true,
               5,
               0.15,
	       0.35)"

