#!/bin/bash
# Shared output locations for the PhotonJet L3 residual analysis.
#
# Sourced by both run.sh (which writes the test output and does the hadd) and
# run_photonjet_job.sh (which writes the Condor batch outputs). Keep it here and
# only here - when these two disagree, `./run.sh merge` silently finds nothing
# because the jobs wrote somewhere else.

#OUTPUT_DIR="/eos/cms/store/group/phys_heavyions/xirong/L3Residuals/L3FillHistogramOutputs/0907Histograms/"
#MERGED_DIR="/eos/cms/store/group/phys_heavyions/xirong/L3Residuals/L3FillHistogramOutputs/0907MergedHistograms/"

OUTPUT_DIR="/afs/cern.ch/user/x/xirong/CMSSW_14_1_9_patch2/src/residualanalysis/Outputs/L3FillHistOutputs/0907Histograms/"
MERGED_DIR="/afs/cern.ch/user/x/xirong/CMSSW_14_1_9_patch2/src/residualanalysis/Outputs/L3FillHistOutputs/0907MergedHists/"
LOG_DIR="/eos/cms/store/group/phys_heavyions/xirong/L3Residuals/L3FillHistogram/logs/0907/"

mkdir -p $OUTPUT_DIR
mkdir -p $MERGED_DIR
