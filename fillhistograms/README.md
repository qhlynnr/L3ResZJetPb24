# PhotonJet L3 Residual Analysis - `run.sh`

`run.sh` is the driver for filling PhotonJet histograms over the 2024 PP Ref
HardProbes data streams and the matching QCDPhoton MC sample, both locally
(quick test) and via HTCondor (full production).

## Setup

```bash
cd /afs/cern.ch/user/x/xirong/CMSSW_14_1_9_patch2/src
cmsenv
cd residualanalysis/fillhistograms
```

Compile the analysis code before running anything. `compile.C` builds
`histograms.C` and `eventhistograms.C` first, then `analyse_PhotonJet.cc`
last, since its dictionary is built against those two - building it out of
order (or editing one of the earlier files without recompiling all three)
leaves a stale `.so` that crashes every Condor job with a cling dictionary
redefinition segfault (root exit status 129):

```bash
root -l -b -q 'compile.C(true)'   # true = clean rebuild (recommended after any edit)
```

Re-run `compile.C` any time you change `histograms.C`, `eventhistograms.C`,
or `analyse_PhotonJet.cc`.

## Two samples

Every `run.sh` mode takes an optional second argument, the sample:

| Sample | `run.sh <mode>` | Description |
|---|---|---|
| `data` (default) | `./run.sh <mode>` | The five 2024 PP Ref HardProbes streams, `isMC=false`, read from an EOS-mounted directory. |
| `mc` | `./run.sh <mode> mc` | The 2024 PP Ref QCD `pThat-15to1200` Pythia8 sample, `isMC=true`, read from an EOS-mounted directory. |

Both samples use `INPUT_TYPE=directory` - no grid proxy or filelist needed.
Data walks `DATA_BASE_DIR/PPRefHardProbes<N>/crab_*/<timestamp>/`; MC walks
`MC_BASE_DIR` directly (override either to repoint). `make_filelists.sh` /
`make_filelists_mc.sh` and the xrootd/proxy path only matter as a fallback,
if a forest isn't reachable via an EOS mount.

MC was switched from the Vanderbilt `QCDPhoton pThat-30to9999` xrootd sample
to this one: that sample's `weight` branch was uniformly `1.0` (a single
un-stitched, unweighted `pThat>30` inclusive bin), so low-pthat events
migrated into high-photon-pT bins and `photonjet_balance` collapsed from
~0.9 to ~0.2 over 60-300 GeV - unusable as L3 fit input. This sample has a
properly pthat-dependent `weight` and gives a flat, physical balance-vs-pT
shape instead. It's a general QCD sample rather than a dedicated photon+jet
one, so photon+jet selection efficiency is much lower (~0.13% vs ~0.6% of
events) - budget more files/statistics accordingly.

### Grid proxy / filelists (fallback only, xrootd input)

```bash
export X509_USER_PROXY=~/private/x509up_u$(id -u)
voms-proxy-init -voms cms -valid 192:00 -out "$X509_USER_PROXY"
```

Generate it on AFS (`~/private/...`), not the default `/tmp` location - the
Condor schedd runs on a different host and can't see a proxy left in
`/tmp` on lxplus, which makes submitted jobs hold with a transfer-input
failure. `run.sh` checks for this and warns. `./run.sh submit`/`test`/`jobs`
will tell you if `make_filelists.sh` (data) or `make_filelists_mc.sh` (the
old Vanderbilt MC sample) is needed - only relevant if a forest isn't
reachable via an EOS mount and has to fall back to filelist mode.

## Modes

```bash
./run.sh test   [mc]     # one batch locally (~40-100 files) to sanity-check the setup
./run.sh jobs   [mc]     # write the Condor job list only, don't submit
./run.sh submit [mc]     # write the job list and condor_submit it (asks for confirmation)
./run.sh merge  [mc]     # hadd the per-batch outputs into one file per dataset
./run.sh status [mc]     # condor_q summary + per-dataset completed/expected batch counts
```

### `test`

Runs `analyse_PhotonJet.cc+` locally and directly (no Condor) over the first
`NTEST` files (default 100, override with `NTEST=<n>`) of the first dataset
in the sample. Output is tagged `*_TEST_<sample>.root` so it never collides
with real batch outputs and is excluded from `merge`. Use this first to
confirm the setup (proxy, filelists/directory, jet tree, `isMC`) works
before submitting a full job.

```bash
./run.sh test          # data, 100 files
NTEST=20 ./run.sh test mc
```

### `jobs`

Writes the Condor job list (`photonjet_jobs.txt` for data,
`photonjet_mc_jobs.txt` for MC) without submitting - useful for inspecting
the batch split first.

### `submit`

Writes the job list, then prompts `Submit N <sample> jobs? [y/N]` before
calling `condor_submit`. For data, jobs default to `STREAMS="0"` (only
HardProbes stream 0); override to include more streams, e.g.:

```bash
STREAMS="0 1 2 3 4" ./run.sh submit
./run.sh submit mc
```

Each batch is `FILES_PER_JOB` files (default 200, override with
`FILES_PER_JOB=<n>`). Jobs run via `run_photonjet_job.sh`, which wraps the
ROOT call in the `cmssw-el8` container (Condor's execute nodes are bare
AlmaLinux9, whose newer glibc breaks ACLiC's directory-listing headers under
`el8_amd64_gcc12`), and logs to `LOG_DIR` (see `paths.sh`).

### `merge`

`hadd`s each dataset's `*_batch*.root` outputs (from `OUTPUT_DIR`) into one
file per dataset in `MERGED_DIR`. Skips datasets with no outputs yet.

### `status`

Shows `condor_q` plus, per dataset, how many batch output files exist versus
how many jobs were queued (from the job list file).

## Output locations

Defined once in `paths.sh` and shared between `run.sh` and
`run_photonjet_job.sh` - do not let the two disagree, or `./run.sh merge`
will silently find nothing:

- `OUTPUT_DIR` - per-batch `.root` outputs and the local `test` output
- `MERGED_DIR` - `merge`'s hadd'd output, one file per dataset
- `LOG_DIR` - per-batch Condor job logs

## Before a real (non-test) run

Make sure `debug = false` in `analyse_PhotonJet.cc`. `run.sh` checks this
and warns on `test`/`submit`, but won't stop you - with `debug = true` the
event loop is capped at 1000 events and the output is renamed
`*_test.root`.

## Useful overrides

| Variable | Default | Effect |
|---|---|---|
| `STREAMS` | `"0"` | Which HardProbes streams to include (data only) |
| `FILES_PER_JOB` | `200` | Files per Condor batch |
| `NTEST` | `100` | Files read by `test` mode |
| `JET_TREE` | `ak4PFJetAnalyzer/t` | Jet tree path inside the forest |
| `DATA_BASE_DIR` | `/eos/cms/.../2024ppref` | Base of the EOS-mounted data forest |
| `MC_BASE_DIR` | `/eos/cms/.../QCD_pThat-15to1200_.../260526_133929` | Base of the EOS-mounted MC forest |
| `OUTPUT_TAG_BASE` | `AK4_photonjet` | Prefix for output file tags |
