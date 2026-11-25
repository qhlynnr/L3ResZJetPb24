#!/bin/sh
# L3ResidualBharad pipeline driver.
#
#   ./run.sh derive   run only deriveL3_from_photonjet.C
#   ./run.sh fit      run only dofits_L3.C (needs an existing DeriveOutFile)
#   ./run.sh all      run both, in order (default if no argument given)

# Directory this script lives in, so the macros resolve regardless of cwd.
MacroDir="$(cd "$(dirname "$0")" && pwd)"
Mode="${1:-all}"

# ------------------------------------------------------------------
# Step 1: deriveL3_from_photonjet.C - build the pTref (and, if the
# *_jetpt profiles happen to be present upstream in fillhistograms,
# jetpt) response inputs from the merged MC/data histogram files.
# ------------------------------------------------------------------
InputDir=/afs/cern.ch/user/x/xirong/CMSSW_14_1_9_patch2/src/residualanalysis_old/Outputs/L3FillHistOutputs/1123_PtRef50to600_PVfilter/Merged/
McFile=$InputDir/Ak4_photonjet_MC.root
DataFile=$InputDir/AK4_photonjet_Data.root
DeriveOutFile=/afs/cern.ch/user/x/xirong/CMSSW_14_1_9_patch2/src/residualanalysis_old/Outputs/DeriveL3Outputs/L3Residual_1123_PtRef50to600_PVfilter_eta1p3_DeriveL3.root
UseAbs=true
UseWideAbs=true

# ------------------------------------------------------------------
# Step 2: dofits_L3.C - fit the pTref series and export the L3/L2L3
# text payloads (calls L3Res.C then createL2L3ResTextFile.C).
# ------------------------------------------------------------------
InFile="$DeriveOutFile"
OutFolder=/afs/cern.ch/user/x/xirong/CMSSW_14_1_9_patch2/src/residualanalysis_old/Outputs/DoFitOutput/L3Residual1123_PtRef50to600_PVfilter_eta1p3/
OutFile=L3Res_Fits

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
L2ResidualFile=/afs/cern.ch/user/x/xirong/CMSSW_14_1_9_patch2/src/residualanalysis_old/fillhistograms/jecfiles/2024PPRefUpToDateCorrections/Prompt24HIpp_V1_DATA_L2Residual_AK4PF.txt

# '+' (ACLiC) is required, not optional, on both steps: running these plain
# interpreted crashes Cling with a segfault (llvm::orc JIT) once the
# vsjetpt-optional code path in deriveL3_from_photonjet.C/L3Res.C is
# exercised (which it is by default now, since fillhistograms doesn't
# currently produce the *_jetpt profiles). Compiled mode does not hit this.

run_derive() {
    mkdir -p "$(dirname "$DeriveOutFile")"
    root -l -b -q "${MacroDir}/deriveL3_from_photonjet.C+(\"$McFile\", \
                   \"$DataFile\", \
                   \"$DeriveOutFile\", \
                   $RefAlphaBin, \
                   $UseAbs, \
                   $UseWideAbs)"
}

run_fit() {
    if [ ! -s "$InFile" ]; then
        echo "ERROR: ${InFile} does not exist or is empty - run './run.sh derive' first."
        exit 1
    fi
    mkdir -p "$OutFolder"
    root -l -b -q "${MacroDir}/dofits_L3.C+(\"$InFile\", \
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
}

case "$Mode" in
  derive) run_derive ;;
  fit)    run_fit ;;
  all)    run_derive && run_fit ;;
  *)
    echo "ERROR: unknown mode \"$Mode\" - expected \"derive\", \"fit\", or \"all\""
    exit 1
    ;;
esac
