#!/bin/bash
# Enumerate the 2024 PP Ref HardProbes forests over xrootd and write one
# filelist per stream into filelists/.
#
# These files live ONLY on the Vanderbilt SE - there is no local/EOS mount -
# so the analysis has to read them as xrootd URLs. input_config.h's
# DIRECTORY mode uses POSIX opendir()/stat(), which cannot see them; FILELIST
# mode is the only workable input type. This script produces those lists.
#
# Needs a valid grid proxy: voms-proxy-init -voms cms -valid 192:00

set -u

REDIR="root://xrootd-se31-vanderbilt.sites.opensciencegrid.org/"
BASE="/store/user/bharikri/Forest/GammaJetSubstructure/2024_PP_Data/2026_01_31_PPRefHardProbes0_EGMcorr_PhoPFiso_ChargedConst_main"

OUTDIR="filelists"
mkdir -p "$OUTDIR"

if ! voms-proxy-info -exists 2>/dev/null; then
    echo "ERROR: no valid grid proxy. Run: voms-proxy-init -voms cms -valid 192:00"
    exit 1
fi

STREAMS=${STREAMS:-"0 1 2 3 4"}

for s in $STREAMS; do
    stream="PPRefHardProbes${s}"
    out="${OUTDIR}/${stream}.txt"
    : > "$out"

    echo "=== $stream ==="

    # crab_* dir, then the timestamp dir underneath it (differs per stream)
    crabdir=$(xrdfs "$REDIR" ls "${BASE}/${stream}" 2>/dev/null | head -1)
    if [ -z "$crabdir" ]; then
        echo "  WARNING: no crab dir found, skipping"
        continue
    fi
    tsdir=$(xrdfs "$REDIR" ls "$crabdir" 2>/dev/null | head -1)
    if [ -z "$tsdir" ]; then
        echo "  WARNING: no timestamp dir found, skipping"
        continue
    fi

    # numbered subdirs 0000, 0001, ...
    for sub in $(xrdfs "$REDIR" ls "$tsdir" 2>/dev/null); do
        n=$(xrdfs "$REDIR" ls "$sub" 2>/dev/null \
            | grep '\.root$' \
            | sed "s|^|${REDIR}|" \
            | tee -a "$out" | wc -l)
        echo "  $(basename "$sub"): $n files"
    done

    echo "  -> $(wc -l < "$out") files in $out"
done

echo
echo "Done. Totals:"
wc -l "${OUTDIR}"/*.txt
