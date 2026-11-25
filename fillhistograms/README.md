# PhotonJet L3 Residual Analysis - `run.sh`

`run.sh` is the driver for filling PhotonJet histograms over the 2024 PP Ref
HardProbes data streams and the 2024 PP Ref QCD MC sample, both locally
(quick test) and via HTCondor (full production).

## Setup

```bash
cd /afs/cern.ch/user/x/xirong/CMSSW_14_1_9_patch2/src
cmsenv
cd residualanalysis_old/fillhistograms
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
calling `condor_submit`. For data, jobs default to all five HardProbes
streams (`STREAMS="0 1 2 3 4"`); override to pick a subset, e.g.:

```bash
STREAMS="0" ./run.sh submit
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

## Event selection (`analyse_PhotonJet.cc`)

Applied in this order. Line numbers refer to `analyse_PhotonJet.cc` unless
another file is named.

### Inputs

| Item | Value |
|---|---|
| Jet tree | `ak4PFJetAnalyzer/t` (AK4 PF, R = 0.4) |
| Photon tree | `ggHiNtuplizer/EventTree` (also provides `rho` for L1) |
| Event tree | `hiEvtAnalyzer/HiTree` (`vz`, `hiBin`, MC `weight`/`pthat`) |
| Event weight | data: 1; MC: forest `weight` branch (pthat-dependent) |

### Trigger and event filters

- **Data:** intended trigger is `HLT_PPRefGEDPhoton30_v6`, but the branch read
  is currently commented out and `HLT_Photon30` is hard-coded to 1, so **no
  trigger requirement is applied** right now. The 50 GeV photon threshold
  sits well above the 30 GeV HLT threshold.
- **MC:** no trigger requirement.
- **Event filters (data only):** `pprimaryVertexFilter == 1` from
  `skimanalysis/HltTree`. The MC forest has no skim tree, so it isn't applied
  to MC. Part of run 387445 (the first files of each HardProbes stream) has
  no reconstructed tracks (CHF = 0 for every jet) and contains noise events
  with 100+ jets; the filter removes all of those, and the run's remaining
  events have normal tracking. About 64-76% of events pass overall. No |vz|
  cut is applied.
- **Jet multiplicity guard:** events with `nref > MAXJETS` (250, `settings.h`)
  are skipped before the jet arrays are read, and counted in the job log.
- **Basic event requirements:** at least one photon (`nPho >= 1`), at least
  one jet (`nref >= 1`), leading jet *uncorrected* pT >= 15 GeV (`jtptmin`
  in `settings.h`).

### Jet energy corrections (`settings.h`, `REDOJES = 1`)

All jets are re-corrected from raw pT before any jet selection, using the
`Prompt24HIpp_V1` files in `jecfiles/2024PPRefUpToDateCorrections/`:

- L1FastJet (MC file), uses jet area and `rho`
- L2Relative (MC file)
- L2Residual (data only)

No L3 residual is applied (that is what this analysis measures). No JER
smearing is applied to MC.

### Photon selection

Kinematic preselection, applied to every photon; the highest-pT survivor is
the tag photon:

| Cut | Value |
|---|---|
| E<sub>T</sub> | >= 50 GeV |
| \|η\| | <= 1.3 (barrel) |
| σ<sub>iηiη</sub> | >= 0.002 (anomalous-signal / spike rejection) |

MC only: the photon must be gen-matched (`pho_genMatchedIndex != -1`) to a
generator photon (`mcPID == 22`) whose mother has \|PID\| <= 22 or is
unknown (-999), with generator calorimeter isolation `mcCalIsoDR04 < 3` GeV.
This selects prompt, isolated photons and rejects photons from hadron decays.

ID and isolation, applied to the leading photon only (event rejected if it
fails, the next photon is not tried):

| Cut | Value |
|---|---|
| H/E | <= 0.2 |
| σ<sub>iηiη</sub> (`phoSigmaIEtaIEta_2012`) | <= 0.021 |
| Charged PF isolation (`pfcIso3subUEec`) | <= 2.0 GeV |
| Neutral PF isolation (`pfnIso3subUEec`) | <= 2.0 GeV |
| Photon PF isolation (`pfpIso3subUEec`) | <= 2.0 GeV |

The isolations are ΔR < 0.3 cones with UE subtraction.

### Away-side jets

A jet is an away-side jet if, after JEC:

| Cut | Value |
|---|---|
| Corrected p<sub>T</sub> | >= 15 GeV |
| Δφ(γ, jet) | >= 2π/3 (2.094) |
| ΔR(γ, jet) | >= 0.4 |

There is no jet \|η\| cut at this stage. Away-side jets are sorted by
corrected pT:

- **Probe jet:** the leading away-side jet. At least one is required.
- **Jet ID on the probe:** the probe must pass `JetSelection_pp`
  (`JetSelection.h`), otherwise the event is rejected. Jet ID is checked
  after sorting, so a failing jet is never replaced by the next one. Jobs run
  with `checkjetid = true`. The second away-side jet is *not* required to
  pass jet ID.

`JetSelection_pp` (2024 JetID13p6TeV recommendation), with
cm = charged-hadron + electron + muon multiplicity:

| \|η\| range | Requirement |
|---|---|
| <= 2.6 | NEF < 0.9, MUF < 0.8, CEF < 0.8, NHF < 0.99, CHF > 0.01, cm > 0 |
| 2.6 – 2.7 | NHF < 0.90, NEF < 0.99, MUF < 0.8, CEF < 0.8, cm > 0 |
| 2.7 – 3.0 | NHF < 0.99 |
| 3.0 – 5.0 | NEF < 0.9 |

`JetSelection_pp` also rejects jets in a hot cell of the veto map.

### Jet veto map

Map: `jetvetomap_all` from `jecfiles/Winter24Prompt24_2024BCDEFGHI.root`
(`applyjetvetomap = true`). The event is rejected if the photon, the probe
jet, or the second away-side jet (if present) falls in a vetoed (η, φ) cell.

### Analysis variables

| Variable | Definition |
|---|---|
| p<sub>T</sub><sup>ref</sup> (`ptavgtp`, the histogram x-axis) | photon E<sub>T</sub> |
| Balance (response) | p<sub>T</sub><sup>probe</sup> / p<sub>T</sub><sup>γ</sup> |
| α | p<sub>T</sub> of the 2nd away-side jet / p<sub>T</sub><sup>γ</sup> (0 if there is only one away-side jet) |

Only away-side jets enter α. Jets on the photon side (Δφ < 2π/3) are not
counted.

### Binning (`histograms.h`)

- **p<sub>T</sub><sup>ref</sup> (`ptforjec`):** 50, 60, 70, 80, 90, 100, 125,
  150, 175, 200, 300, 400, 500, 600 GeV
- **α thresholds (`alphavalues`):** 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40,
  0.45, 0.50. The 3D profiles are filled cumulatively: an event with
  α < threshold goes into every bin with that threshold or a larger one.
- **Output η folders (`etaedges` in `settings.h`), selected on probe-jet η:**
  `eta_-2.5_2.5` and `eta_-5.2_5.2`. The 3D profiles are only in
  `eta_-5.2_5.2`. Inside them the probe jet η is binned with `wetarange`
  (signed), `wabsetarange` (\|η\|), `etarange` (narrow) or `dwabsetarange`
  (single 0–1.3 bin).
- **Centrality folder:** `hibin_-1.0_0.0` only (pp).

### Disabled options

- `checkvalidjet = false`: the old JEC-validity η/pT window is not applied.
- `debug = false`: required for production. See [Before a real run](#before-a-real-non-test-run).

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
| `STREAMS` | `"0 1 2 3 4"` | Which HardProbes streams to include (data only) |
| `FILES_PER_JOB` | `200` | Files per Condor batch |
| `NTEST` | `100` | Files read by `test` mode |
| `JET_TREE` | `ak4PFJetAnalyzer/t` | Jet tree path inside the forest |
| `DATA_BASE_DIR` | `/eos/cms/.../2024ppref` | Base of the EOS-mounted data forest |
| `MC_BASE_DIR` | `/eos/cms/.../QCD_pThat-15to1200_.../260526_133929` | Base of the EOS-mounted MC forest |
| `OUTPUT_TAG_BASE` | `AK4_photonjet` | Prefix for output file tags |
