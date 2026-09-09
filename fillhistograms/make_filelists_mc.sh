#!/bin/bash
# Enumerate the 2024 PP Ref QCDPhoton MC forest over xrootd and write
# filelists/PPRefQCDPhotonMC.txt.
#
# Same story as make_filelists.sh for data: these files live only on the
# Vanderbilt SE, so the analysis reads them as xrootd URLs and has to run in
# FILELIST mode (input_config.h's DIRECTORY mode uses POSIX opendir()/stat()
# and cannot see them).
#
# Unlike the data, the MC is one sample rather than five HardProbes streams, so
# this writes a single list covering all the numbered subdirs (0000, 0001, ...).
#
# Needs a valid grid proxy: voms-proxy-init -voms cms -valid 192:00

set -u

REDIR="root://xrootd-se31-vanderbilt.sites.opensciencegrid.org/"
BASE="/store/user/bharikri/Forest/GammaJetSubstructure/2024_PP_MC/2026_04_05_Pythia8_2024_ppRef_QCDPhoton/QCD-Photon_pThat-30to9999_TuneCP5_5p36TeV_pythia8/crab_2026_04_05_Pythia8_2024_ppRef_QCDPhoton30to9999_EGMcorr_PhoPFiso_ChargedConst/260405_205020"

OUTDIR="filelists"
OUT="${OUTDIR}/PPRefQCDPhotonMC.txt"

mkdir -p "$OUTDIR"

if ! voms-proxy-info -exists 2>/dev/null; then
    echo "ERROR: no valid grid proxy. Run: voms-proxy-init -voms cms -valid 192:00"
    exit 1
fi

: > "$OUT"

echo "=== QCDPhoton pThat-30to9999 MC ==="

# The numbered subdirs are listed explicitly rather than globbed: only 0000 and
# 0001 hold this sample. Add to the list if more show up later.
for sub in 0000 0001; do
    n=$(xrdfs "$REDIR" ls "${BASE}/${sub}" 2>/dev/null \
        | grep '\.root$' \
        | sed "s|^|${REDIR}|" \
        | tee -a "$OUT" | wc -l)
    if [ "$n" -eq 0 ]; then
        echo "  WARNING: ${sub}: no files found"
    else
        echo "  ${sub}: ${n} files"
    fi
done

total=$(wc -l < "$OUT")
if [ "$total" -eq 0 ]; then
    echo "ERROR: no files found at all - check the proxy and the BASE path."
    exit 1
fi

echo "  -> ${total} files in ${OUT}"
