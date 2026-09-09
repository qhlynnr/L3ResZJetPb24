#!/bin/bash
# Condor job wrapper for the PhotonJet L3 residual analysis.
# Arguments: INPUT OUTPUT_TAG BATCH_INDEX TOTAL_BATCHES FILES_PER_JOB JET_TREE IS_MC INPUT_TYPE
#
# INPUT is a filelist file or a local/EOS directory, depending on INPUT_TYPE
# ("filelist" or "directory" - see run.sh). Each job reads a contiguous slice:
#   skipFiles = BATCH_INDEX * FILES_PER_JOB   (computed inside analyse_PhotonJet)
# and writes  <OUTPUT_DIR>/<OUTPUT_TAG>_batch<N>_of_<M>.root
#
# IS_MC must be the literal "true" or "false" - it is pasted straight into the
# ROOT call as the analyse_PhotonJet isMC argument. run.sh sets it from MODE.

set -u

FILELIST=$1
OUTPUT_TAG=$2
BATCH_INDEX=$3
TOTAL_BATCHES=$4
FILES_PER_JOB=$5
JET_TREE=${6:-ak4PFJetAnalyzerSDZcut1/t}
IS_MC=${7:-false}
INPUT_TYPE=${8:-filelist}

# ---- paths ----
CMSSW_DIR="/afs/cern.ch/user/x/xirong/CMSSW_14_1_9_patch2/src"
ANALYSIS_DIR="${CMSSW_DIR}/residualanalysis/fillhistograms"
# OUTPUT_DIR / LOG_DIR come from paths.sh, shared with run.sh so that
# `./run.sh merge` looks where the jobs actually wrote.
source "${ANALYSIS_DIR}/paths.sh"
# ---------------

if [ "$IS_MC" != "true" ] && [ "$IS_MC" != "false" ]; then
    echo "ERROR: IS_MC must be \"true\" or \"false\", got \"$IS_MC\""
    exit 1
fi

mkdir -p "$OUTPUT_DIR" "$LOG_DIR"
exec >> "${LOG_DIR}/photonjet_${OUTPUT_TAG}_batch${BATCH_INDEX}.log" 2>&1

echo "=========================================="
echo "Start        : $(date)"
echo "Input        : $FILELIST"
echo "Input type   : $INPUT_TYPE"
echo "Output tag   : $OUTPUT_TAG"
echo "Batch        : $BATCH_INDEX / $TOTAL_BATCHES  ($FILES_PER_JOB files)"
echo "Jet tree     : $JET_TREE"
echo "isMC         : $IS_MC"
echo "Host         : $(hostname)"

if [ "$INPUT_TYPE" = "filelist" ]; then
    # The forests are read over xrootd, so the job needs the grid proxy that
    # Condor transferred in (x509userproxy in the submit file).
    if [ -n "${X509_USER_PROXY:-}" ] && [ -f "${X509_USER_PROXY}" ]; then
        echo "Proxy        : ${X509_USER_PROXY} ($(voms-proxy-info -timeleft 2>/dev/null)s left)"
    else
        echo "ERROR: no X509_USER_PROXY in the job environment - xrootd reads will fail."
        exit 1
    fi
fi

# Condor's execute nodes are bare AlmaLinux9, but this area is
# el8_amd64_gcc12 - on AlmaLinux9's newer glibc, ACLiC's clang/cling can't
# parse /usr/include/dirent.h ("use of undeclared identifier 'opendir'")
# and segfaults. Run the CMSSW setup + ROOT call inside the matching EL8
# container instead. (cmsenv/scramv1 alone does not fix this - it only sets
# environment variables, it doesn't change which OS headers get compiled
# against.)
INNER_SCRIPT="${ANALYSIS_DIR}/.inner_${OUTPUT_TAG}_batch${BATCH_INDEX}.sh"
cat > "$INNER_SCRIPT" <<EOF
#!/bin/bash
set -u
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd "${CMSSW_DIR}" || { echo "ERROR: cannot cd to ${CMSSW_DIR}"; exit 1; }
eval \$(scramv1 runtime -sh)
cd "${ANALYSIS_DIR}" || { echo "ERROR: cannot cd to ${ANALYSIS_DIR}"; exit 1; }

# Two container-specific fixes, both found by reproducing a Condor job's
# exact failure interactively (see git log for the diagnosis):
#
# 1. Inside this container, R__LOAD_LIBRARY's bare-name dlopen of
#    histograms_C.so/eventhistograms_C.so unreliably fails to find them via
#    cwd even though cwd is right here - unlike plain interactive
#    'root -l -b -q' outside the container, which never had this problem.
#    Putting cwd on LD_LIBRARY_PATH makes the resolution explicit instead of
#    relying on cwd being in ROOT's dynamic search path.
export LD_LIBRARY_PATH="${ANALYSIS_DIR}:\${LD_LIBRARY_PATH}"
#
# 2. This container image's /etc/nsswitch.conf lists "sss" first for passwd/
#    group, but doesn't actually ship libnss_sss.so.2 (sssd-client isn't
#    installed here). Something in CMSSW's library init code (triggered when
#    histograms_C.so/eventhistograms_C.so are dlopen'd) does a passwd/group
#    lookup, glibc tries "sss" first per nsswitch.conf, dlopen fails, and
#    Cling's DynamicLibraryManager treats that dlerror() as fatal to the
#    whole load - aborting root with exit 1 before analyse_PhotonJet ever
#    runs. A harmless empty stub satisfies the dlopen() (its symbols are
#    never actually called for anything analyse_PhotonJet does) and clears
#    the error. nssstub/libnss_sss.so.2 in this directory is that stub.
export LD_LIBRARY_PATH="${ANALYSIS_DIR}/nssstub:\${LD_LIBRARY_PATH}"

# Signature: (input, outputfiletag, isMC, checkjetid, inputType, maxFiles,
#             maxEvents, outputDir, batchIndex, totalBatches, jetPath)
root -l -b -q "analyse_PhotonJet.cc+(\"${FILELIST}\", \"${OUTPUT_TAG}\", ${IS_MC}, false, \"${INPUT_TYPE}\", ${FILES_PER_JOB}, -1, \"${OUTPUT_DIR}/\", ${BATCH_INDEX}, ${TOTAL_BATCHES}, \"${JET_TREE}\")"
EOF

/cvmfs/cms.cern.ch/common/cmssw-el8 --command-to-run "bash ${INNER_SCRIPT}"
STATUS=$?
rm -f "$INNER_SCRIPT"

echo "root exit status: $STATUS"
echo "Finished     : $(date)"
exit $STATUS
