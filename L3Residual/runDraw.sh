#!/bin/bash

INPUT="/afs/cern.ch/user/x/xirong/CMSSW_14_1_9_patch2/src/residualanalysis/Outputs/DeriveL3Outpus//0908DeriveL3.root"
OUTDIR="/afs/cern.ch/user/x/xirong/CMSSW_14_1_9_patch2/src/residualanalysis/Outputs/PlotsFromDeriveL3/0908Derivel3/"

mkdir -p $OUTDIR

root -l -b -q "DrawL3ResidualBalance.C(\"$INPUT\",\"$OUTDIR\")"
