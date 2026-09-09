#!/bin/sh

# Directory this script lives in, so dofits_L3.C resolves regardless of cwd.
MacroDir="$(cd "$(dirname "$0")" && pwd)"

InFile=/afs/cern.ch/user/x/xirong/CMSSW_14_1_9_patch2/src/residualanalysis/Outputs/DeriveL3Outpus/0908DeriveL3.root
OutFolder=/afs/cern.ch/user/x/xirong/CMSSW_14_1_9_patch2/src/residualanalysis/Outputs/DoFitOutput/0908DoFit_newcode
OutFile=L3_0901_Fits_

SampleTypes=photonjet
PtRange=
RunLabel=2024ppRef
LumiLabel="pp Reference"
RefAlphaBin=5
FitAlphaMin=0.15
FitAlphaMax=0.35

# Real L2Residual file used to lay out the L2L3Residual export (the macro's
# built-in default path is stale/missing the 2024PPRefUpToDateCorrections
# subfolder, so this must be passed explicitly).
L2ResidualFile=/afs/cern.ch/user/x/xirong/CMSSW_14_1_9_patch2/src/residualanalysis/fillhistograms/jecfiles/2024PPRefUpToDateCorrections/Prompt24HIpp_V1_DATA_L2Residual_AK4PF.txt

mkdir -p "$OutFolder"

root -l -b -q "${MacroDir}/dofits_L3.C(\"$InFile\", \
               \"$SampleTypes\", \
               \"$PtRange\", \
               \"$OutFile\", \
               \"$RunLabel\", \
               \"$LumiLabel\", \
               $RefAlphaBin, \
               $FitAlphaMin, \
               $FitAlphaMax, \
               \"$L2ResidualFile\", \
               \"$OutFolder\")"
