#!/bin/bash
# PhotonJet L3 residual analysis - driver.
#
# Every mode takes an optional sample argument, "data" (default) or "mc":
#
#   ./run.sh test [mc]      one batch locally (~40 files) to sanity-check the setup
#   ./run.sh jobs [mc]      write the job list only
#   ./run.sh submit [mc]    write the job list and condor_submit it
#   ./run.sh merge [mc]     hadd the per-batch outputs into one file per dataset
#   ./run.sh status [mc]    condor_q summary
#
# data = the five 2024 PP Ref HardProbes streams, isMC=false.
# mc   = the 2024 PP Ref QCD pThat-15to1200 Pythia8 sample (EOS-mounted,
#        nbarnett's production), isMC=true.
#
# Input: each sample picks whichever of input_config.h's two file-discovery
# modes matches where its forest actually lives, via INPUT_TYPE below:
#   filelist  - a text file of xrootd URLs, for forests reachable only over
#               xrootd (no local/EOS mount there). Lists are made by
#               ./make_filelists.sh; ./make_filelists_mc.sh is the old
#               Vanderbilt QCDPhoton xrootd list generator, unused by default
#               now that mc uses directory mode - keep it around only if you
#               need to fall back to that sample.
#   directory - a local/EOS path, recursed with POSIX opendir()/stat(). No
#               filelist needed; override DATA_BASE_DIR/MC_BASE_DIR to repoint.
# Both data and mc use directory mode (EOS) by default.
#
# mc was switched from the Vanderbilt QCDPhoton pThat-30to9999 xrootd sample
# to this one because that sample's "weight" branch was uniformly 1.0 (a
# single un-stitched, unweighted pThat>30 inclusive bin) - low-pthat events
# migrating into high-photon-pT bins made photonjet_balance collapse from
# ~0.9 to ~0.2 over 60-300 GeV, an unusable input for the L3 fit. This QCD
# pThat-15to1200 sample has a properly pthat-dependent weight and gives a
# flat, physical balance-vs-pT shape instead. It's a general QCD sample
# rather than a dedicated photon+jet one, so photon+jet selection efficiency
# is much lower (~0.13% vs ~0.6% of events) - budget more files accordingly.
#
# Turn off `debug` in analyse_PhotonJet.cc before a real run: debug caps the
# event loop at 1000 events and renames the output to *_test.root.

set -u

MODE=${1:-help}
SAMPLE=${2:-data}

FILES_PER_JOB=${FILES_PER_JOB:-200}
JET_TREE=${JET_TREE:-ak4PFJetAnalyzer/t}
OUTPUT_TAG_BASE=${OUTPUT_TAG_BASE:-AK4_photonjet}

FILELIST_DIR="filelists"
# Base of the local/EOS-mounted data forest, for INPUT_TYPE=directory.
# Temporary location - the permanent one will be a similar local path,
# so this only needs a one-line update (or override) when it moves.
DATA_BASE_DIR=${DATA_BASE_DIR:-"/eos/cms/store/group/phys_heavyions/nbarnett/HiForest/DATA/2024ppref"}
# Base of the local/EOS-mounted MC forest (QCD pThat-15to1200), for
# INPUT_TYPE=directory. Recursed automatically (the 0000/0001/0002/...
# subdirectories under it are walked by GetFilesFromDirectory).
MC_BASE_DIR=${MC_BASE_DIR:-"/eos/cms/store/group/phys_heavyions/nbarnett/HiForest/MC/2024ppref/QCD_pThat-15to1200_TuneCP5_5p36TeV_pythia8/crab_foresting_2024ppref_MC_v2/260526_133929"}
# OUTPUT_DIR / MERGED_DIR / LOG_DIR are shared with run_photonjet_job.sh.
source ./paths.sh

# ------------------------------------------------------------ sample config
# Each sample defines:
#   IS_MC      passed through to analyse_PhotonJet's isMC argument
#   INPUT_TYPE "filelist" or "directory" - see the header comment above
#   JOBFILE    the Condor job list this sample writes
#   SUBFILE    the submit file that reads that job list
#   DATASETS   "<path>:<output tag suffix>" pairs - path is a filelist
#              basename (filelist mode) or an absolute directory (directory
#              mode)
#   MAKELIST   the script that generates the filelists (filelist mode only)

# Resolve the crab/timestamp subdir under DATA_BASE_DIR/PPRefHardProbes<N>,
# same idea as make_filelists.sh's xrdfs walk but over a local/EOS mount.
# Echoes the resolved directory, or nothing if it can't be found.
resolve_stream_dir() {
    local base="${DATA_BASE_DIR}/PPRefHardProbes$1"
    local crabdir tsdir
    crabdir=$(find "$base" -mindepth 1 -maxdepth 1 -type d -name 'crab_*' 2>/dev/null | sort | head -1)
    [ -z "$crabdir" ] && return
    tsdir=$(find "$crabdir" -mindepth 1 -maxdepth 1 -type d 2>/dev/null | sort | head -1)
    echo "$tsdir"
}

case "$SAMPLE" in
  data)
    IS_MC="false"
    INPUT_TYPE="directory"
    JOBFILE="photonjet_jobs.txt"
    SUBFILE="condor_photonjet.sub"
    MAKELIST="./make_filelists.sh"
    # One entry per HardProbes stream (all five by default). Override to
    # pick a subset, e.g.
    #   STREAMS="0" ./run.sh submit
    STREAMS=${STREAMS:-"0 1 2 3 4"}
    DATASETS=()
    for s in $STREAMS; do
        dir=$(resolve_stream_dir "$s")
        if [ -z "$dir" ]; then
            echo "WARNING: no directory found for stream ${s} under ${DATA_BASE_DIR}/PPRefHardProbes${s}"
        fi
        DATASETS+=("${dir}:HP${s}")
    done
    if [ ${#DATASETS[@]} -eq 0 ]; then
        echo "ERROR: STREAMS is empty - nothing to do"
        exit 1
    fi
    ;;
  mc)
    IS_MC="true"
    INPUT_TYPE="directory"
    JOBFILE="photonjet_mc_jobs.txt"
    SUBFILE="condor_photonjet_mc.sub"
    MAKELIST="./make_filelists_mc.sh"
    # One sample, not five streams. See MC_BASE_DIR above to repoint it.
    DATASETS=("${MC_BASE_DIR}:MC")
    ;;
  help|-h|--help)
    sed -n '2,25p' "$0"
    exit 0
    ;;
  *)
    echo "ERROR: unknown sample \"$SAMPLE\" - expected \"data\" or \"mc\""
    exit 1
    ;;
esac

mkdir -p "$OUTPUT_DIR" "$MERGED_DIR"

# ---------------------------------------------------------------- helpers

# "path:suffix" -> input path: the path itself (directory mode) or
# FILELIST_DIR/path.txt (filelist mode)
input_path() {
    local key="${1%%:*}"
    if [ "$INPUT_TYPE" = "directory" ]; then
        echo "$key"
    else
        echo "${FILELIST_DIR}/${key}.txt"
    fi
}
# "path:suffix" -> output tag
out_tag()   { echo "${OUTPUT_TAG_BASE}_${1##*:}"; }
# "path:suffix" -> label for messages (the suffix tag in directory mode,
# since the path itself is a long absolute directory there)
label() {
    if [ "$INPUT_TYPE" = "directory" ]; then
        echo "${1##*:}"
    else
        echo "${1%%:*}"
    fi
}

check_proxy() {
    if ! voms-proxy-info -exists 2>/dev/null; then
        echo "ERROR: no valid grid proxy. The forests are xrootd-only, so jobs"
        echo "       cannot read anything without one. Run:"
        echo "         export X509_USER_PROXY=~/private/x509up_u\$(id -u)"
        echo "         voms-proxy-init -voms cms -valid 192:00 -out \"\$X509_USER_PROXY\""
        exit 1
    fi
    case "$(voms-proxy-info -path 2>/dev/null)" in
        /tmp/*)
            echo "WARNING: proxy lives in /tmp, which is local to this lxplus node."
            echo "         The Condor schedd (bigbird102) is a different host and"
            echo "         can't see it there - submitted jobs will hold with"
            echo "         'Transfer input files failure ... No such file or directory'."
            echo "         Regenerate it on AFS instead:"
            echo "           export X509_USER_PROXY=~/private/x509up_u\$(id -u)"
            echo "           voms-proxy-init -voms cms -valid 192:00 -out \"\$X509_USER_PROXY\""
            echo
            ;;
    esac
    echo "Proxy OK: $(voms-proxy-info -timeleft)s remaining"
}

check_inputs() {
    local missing=0
    for d in "${DATASETS[@]}"; do
        local p; p=$(input_path "$d")
        if [ "$INPUT_TYPE" = "directory" ]; then
            if [ -z "$p" ] || [ ! -d "$p" ]; then
                echo "ERROR: missing directory for $(label "$d") (${p:-<unresolved>})"
                missing=1
            fi
        else
            if [ ! -s "$p" ]; then
                echo "ERROR: missing ${p}"
                missing=1
            fi
        fi
    done
    if [ $missing -eq 1 ]; then
        if [ "$INPUT_TYPE" = "directory" ] && [ "$SAMPLE" = "data" ]; then
            echo "Check DATA_BASE_DIR (${DATA_BASE_DIR}) and its crab/timestamp subdirectories."
        elif [ "$INPUT_TYPE" = "directory" ]; then
            echo "Check MC_BASE_DIR (${MC_BASE_DIR})."
        else
            echo "Generate the filelists first:  ${MAKELIST}"
        fi
        exit 1
    fi
}

check_debug_flag() {
    if grep -qE '^\s*bool\s+debug\s*=\s*true' analyse_PhotonJet.cc; then
        echo "WARNING: debug = true in analyse_PhotonJet.cc"
        echo "         The event loop is capped at 1000 events and output is"
        echo "         renamed *_test.root. Set it to false for a real run."
        echo
    fi
}

# Every Condor job execs "analyse_PhotonJet.cc+" (see run_photonjet_job.sh) in
# this same shared AFS directory. ACLiC only rebuilds a .so when its source is
# newer, so if histograms.C/eventhistograms.C/analyse_PhotonJet.cc were edited
# since the last build, the .so's are stale when "submit" fires - and every
# job launched at once independently notices that and races to rebuild the
# same analyse_PhotonJet_cc.so/histograms_C.so, stepping on each other's
# rootcling dictionary output (leftover analyse_PhotonJet_cc_ACLiC_dict.cxx_
# tmp_* files are the fingerprint of that race). A stale histograms_C.so is
# worse than just slow: if histograms.h changed but that .so wasn't rebuilt,
# it disagrees with analyse_PhotonJet_cc.so about the histograms class's
# member layout, so h->someNewMember reads garbage instead of the real
# pointer - a histogram silently isn't filled/written instead of crashing.
# Building all three serially, once, before any job starts avoids both.
compile_lib() {
    echo "Compiling histograms.C, eventhistograms.C, analyse_PhotonJet.cc once"
    echo "(serially, before any Condor job can race on the shared .so files)..."
    local src log
    for src in histograms.C eventhistograms.C analyse_PhotonJet.cc; do
        log=$(root -l -b -q -e ".L ${src}+" 2>&1)
        if echo "$log" | grep -qE ' error:| Error in '; then
            echo "$log"
            echo "ERROR: compiling ${src} failed - see above. Not submitting."
            exit 1
        fi
    done
    echo "Compile OK."
    echo
}

# ------------------------------------------------------------ job list gen

write_jobs() {
    check_inputs
    # NOTE: no "# input output_tag ..." header line in $JOBFILE. HTCondor's
    # "queue ...  from <file>" does NOT skip '#'-prefixed lines the way a
    # shell script would - it reads the header as a real (garbage) data row
    # and queues an extra bogus job for it. Columns, for reference:
    #   input output_tag batch_index total_batches files_per_job jet_tree is_mc input_type
    : > "$JOBFILE"

    local grand_files=0 grand_jobs=0
    for d in "${DATASETS[@]}"; do
        local input tag nfiles nbatches
        input=$(input_path "$d")
        tag=$(out_tag "$d")
        if [ "$INPUT_TYPE" = "directory" ]; then
            nfiles=$(find "$input" -name '*.root' | wc -l)
        else
            nfiles=$(grep -c . "$input")
        fi
        # ceiling division
        nbatches=$(( (nfiles + FILES_PER_JOB - 1) / FILES_PER_JOB ))

        for (( b=0; b<nbatches; b++ )); do
            echo "${input} ${tag} ${b} ${nbatches} ${FILES_PER_JOB} ${JET_TREE} ${IS_MC} ${INPUT_TYPE}" >> "$JOBFILE"
        done

        echo "  $(label "$d"): ${nfiles} files -> ${nbatches} jobs"
        grand_files=$(( grand_files + nfiles ))
        grand_jobs=$(( grand_jobs + nbatches ))
    done

    echo
    echo "Sample: ${SAMPLE} (isMC=${IS_MC})"
    echo "Total: ${grand_files} files -> ${grand_jobs} jobs, ${FILES_PER_JOB} files/job"
    echo "Wrote ${JOBFILE}"
}

# ------------------------------------------------------------------ modes

case "$MODE" in

  test)
    [ "$INPUT_TYPE" = "filelist" ] && check_proxy
    check_debug_flag
    check_inputs
    d="${DATASETS[0]}"
    input=$(input_path "$d")
    NTEST=${NTEST:-100}

    echo "Local test: first ${NTEST} files of $(label "$d")"
    echo "Input type: ${INPUT_TYPE}"
    echo "Jet tree  : ${JET_TREE}"
    echo "isMC      : ${IS_MC}"
    echo "Output    : ${OUTPUT_DIR}/"
    echo

    # batchIndex 0 with maxFiles=NTEST reads exactly the first NTEST files.
    # The _TEST tag keeps this out of the way of the real batch outputs, and out
    # of what `merge` picks up.
    root -l -b -q "analyse_PhotonJet.cc+(\"${input}\", \"${OUTPUT_TAG_BASE}_TEST_${SAMPLE}\", ${IS_MC}, true, \"${INPUT_TYPE}\", ${NTEST}, -1, \"${OUTPUT_DIR}/\", 0, 1, \"${JET_TREE}\")"
    ;;

  jobs)
    write_jobs
    ;;

  submit)
    check_debug_flag
    compile_lib
    write_jobs
    echo
    if [ "$INPUT_TYPE" = "filelist" ]; then
        check_proxy
        # Condor copies the proxy from its current path; make sure the submit
        # file points at the proxy this shell is actually using.
        proxy_path=$(voms-proxy-info -path)
        if ! grep -q "x509userproxy *= *${proxy_path}" "$SUBFILE"; then
            echo "NOTE: ${SUBFILE} x509userproxy does not match ${proxy_path}"
            echo "      Update it before submitting if the jobs fail on permissions."
            echo
        fi
    fi
    read -p "Submit $(grep -vc '^#' "$JOBFILE") ${SAMPLE} jobs? [y/N] " ans
    if [ "$ans" = "y" ] || [ "$ans" = "Y" ]; then
        condor_submit "$SUBFILE"
    else
        echo "Not submitted. Job list is in ${JOBFILE}."
    fi
    ;;

  merge)
    mkdir -p "$MERGED_DIR"
    for d in "${DATASETS[@]}"; do
        tag=$(out_tag "$d")
        parts=( ${OUTPUT_DIR}/${tag}_batch*.root )
        if [ ! -e "${parts[0]}" ]; then
            echo "$(label "$d"): no outputs yet, skipping"
            continue
        fi
        echo "$(label "$d"): merging ${#parts[@]} files"
        hadd -f "${MERGED_DIR}/${tag}.root" "${parts[@]}"
    done
    echo
    echo "Merged outputs in ${MERGED_DIR}"
    ;;

  status)
    condor_q
    echo
    echo "Completed batch outputs (${SAMPLE}):"
    for d in "${DATASETS[@]}"; do
        tag=$(out_tag "$d")
        n=$(ls -1 ${OUTPUT_DIR}/${tag}_batch*.root 2>/dev/null | wc -l)
        exp=$(grep -c " ${tag} " "$JOBFILE" 2>/dev/null || echo "?")
        echo "  $(label "$d"): ${n} / ${exp}"
    done
    ;;

  *)
    sed -n '2,25p' "$0"
    ;;
esac
