# Residual Analysis Framework

Framework for JEC residual corrections (L2, L3) and JER scale factors for CMS Heavy Ion physics.

## Prerequisites

- CMSSW environment (tested with CMSSW_15_1_0_patch3)
- ROOT (provided by CMSSW)
- HTCondor (for batch processing)

## Quick Start

### Setup Environment
```bash
cd /eos/home-b/bharikri/lxplus_private/EGamma/CMSSW_15_1_0_patch3/src
cmsenv
cd /eos/home-b/bharikri/lxplus_private/EGamma/residualanalysis/fillhistograms
```

### Compile Support Classes
```bash
root -l -b -q compile.C
```

### Run Analysis

#### Photon+jet L3 residual workflow (2024 ppRef, updated 2026-11-23)

Full details of `run.sh`, the event selection and the binning are in
[`fillhistograms/README.md`](fillhistograms/README.md).

**1. Fill histograms** (`fillhistograms/`)

- Edit `paths.sh`: set `OUTPUT_DIR` for this production. `MERGED_DIR`
  (`$OUTPUT_DIR/Merged/`) and `LOG_DIR` (`$OUTPUT_DIR/Log/`) follow from it.
- Edit `run.sh` so the MC and data input paths match your samples, and set
  `ANALYSIS_DIR` in `run_photonjet_job.sh` to your checkout.
- Test, then submit:
  ```bash
  ./run.sh test mc
  ./run.sh test data
  ./run.sh submit mc
  ./run.sh submit data        # all five HardProbes streams by default
  STREAMS="0" ./run.sh submit # or a subset
  ```
  `submit` compiles `histograms.C`, `eventhistograms.C` and
  `analyse_PhotonJet.cc` once, serially, before queuing jobs, so the Condor
  jobs don't race to rebuild the shared `.so` files.
- Jobs now run with jet ID on (`checkjetid = true`), the data-only
  `pprimaryVertexFilter`, and `MAXJETS = 250` (the old 50 segfaulted on
  noisy run-387445 events).
- pT<sub>ref</sub> binning: 50, 60, 70, 80, 90, 100, 125, 150, 175, 200,
  300, 400, 500, 600 GeV. The histograms also carry `*_jetpt` profiles
  (balance binned in leading-jet pT).
- To resubmit only failed batches, list them in `photonjet_jobs_rerun.txt`
  and `condor_submit condor_photonjet_rerun.sub`.

**2. Merge** with `./merge_histograms.sh` (set `filepath` to `OUTPUT_DIR`;
writes `AK4_photonjet_Data.root` and `Ak4_photonjet_MC.root` to
`$filepath/Merged/`).

**3. Derive and fit L3** (`L3ResidualNew/`)

Edit the paths at the top of `run.sh` (`InputDir`, `DeriveOutFile`,
`OutFolder`, `L2ResidualFile`), then:
```bash
./run.sh all      # or: ./run.sh derive / ./run.sh fit
```
`derive` now also writes `balance_vsabseta_*` histograms (MC/Data balance
and Data/MC vs fine |η| per pT bin), and `L3Res.C` plots them, plus the
response vs leading-jet pT next to the pT<sub>ref</sub> series. Fit text
files land in `Outputs/DoFitOutput/<tag>/L3Res_Fits/textfiles/`.

**4. Sanity plots** (optional): data vs MC photon pT and Δφ(γ, jet),
shape-normalised, with ratio panels:
```bash
root -l -b -q 'plotSanity.C("Data.root", "MC.root", "outdir")'
```


#### Single File (Legacy Era Mode)
```bash
root -l -b -q 'analyse_PhotonJet.cc("PHOTONHP", "output_tag", false, false)'
```

#### Directory Mode (Multiple Files)
```bash
root -l -b -q 'analyse_PhotonJet.cc("/path/to/directory", "output_tag", false, false, "directory", 100, 10000, "/output/dir")'
```

#### Filelist Mode
```bash
root -l -b -q 'analyse_PhotonJet.cc("/path/to/filelist.txt", "output_tag", false, false, "filelist", -1, -1, "/output/dir")'
```

### Batch Processing

See batch/README.md for HTCondor batch submission.

```bash
python3 batch/submit_condor.py \
    --era PHOTONHP_FULL \
    --input /path/to/input/directory \
    --output-dir /path/to/output \
    --files-per-job 50 \
    --analysis PhotonJet
```

## Analysis Scripts

| Script | Purpose | Output |
|--------|---------|--------|
| analyse.cc | Dijet analysis | L2 residual histograms |
| analyse_PhotonJet.cc | Photon+Jet analysis | L3 residual histograms |
| analyse_JER.cc | JER scale factor | JER SF histograms |
| dofits_L3.C | Fit L3 residuals, optional raw MC/Data and alpha-extrapolation plots | L3 text/plots |

### dofits_L3 helper flags
- `plotRawResponses` (default false): also draw MC and Data balance vs pT (no ratio/fits) in `L3fits_<tag>/raw/`.
- `saveAlphaExtrap` (default false): build MC/Data/ratio vs alpha for each (pT,|eta|), fit a linear alpha→0 extrapolation, and save PNGs in `L3fits_<tag>/alpha_extrap/`.

Example:
```bash
root -l -b -q 'dofits_L3.C("L3_derived.root", 30, 150, "L3Res_photonjet", false, "2024ppRef", "pp Reference", true, true)'
```

## Function Parameters

```cpp
void analyse_PhotonJet(
    string input,           // Era name, file path, directory, or filelist
    string outputfiletag,   // Tag for output filename
    bool isMC,              // Is Monte Carlo
    bool checkjetid,        // Apply jet ID cuts
    string inputType,       // "era", "file", "directory", "filelist"
    int maxFiles,           // Max files to process (-1 = all)
    int maxEvents,          // Max events to process (-1 = all)
    string outputDir,       // Output directory
    int batchIndex,         // Batch job index (-1 = non-batch)
    int totalBatches        // Total number of batches
)
```

## Directory Structure

```
residualanalysis/
├── fillhistograms/     # Main analysis code
│   ├── analyse*.cc     # Analysis scripts
│   ├── histograms.h/C  # Histogram definitions
│   ├── input_config.h  # Multi-file input utilities
│   ├── chain_builder.h # TChain construction
│   └── jecfiles/       # JEC correction files
├── batch/              # Batch submission
├── L2Residual/         # L2 derivation
├── L3Residual/         # L3 derivation
└── JER/                # JER analysis
```

## Configuration

- configurations.h - Era-to-file mapping (legacy, still supported)
- settings.h - Global settings (JEC files, binning, cuts)

At the moment direct balance method with tag-and-probing a dijet system is used.

## Required steps:
```
cd fillhistograms/jecfiles/
wget https://github.com/cms-jet/JECDatabase/blob/master/jet_veto_maps/Summer23BPixPrompt23/Summer23BPixPrompt23_RunD_v1.root
wget https://indico.cern.ch/event/1592138/contributions/6711482/attachments/3141457/5575822/Spring23_HI_V1_MC_L2Relative_AK4PF.txt
```

## Quick start

-> First: root -l compile.C 
-> 

Then, to run the analysis:

1. analyse.cc fills histograms you need.
Flags:

At the moment this is run with hadded ntuples per dat set (MC, HP0, HP1, HP2, ZB0, etc.)

You can do a test run with simply root -l analyse.cc

-> to run locally: run.sh

1. To fill histograms for L2residuals:

- Code is in directory L2Residual.

-> call this code separarely for different datasets. you need cmsenv to be able to call functions to apply jec (and jer sf). shouldn't really matter which cmssw you use.
-> in settings.h you can manually set which JEC files to use. The JEC (and JER SF) files are always read from the path set up in here.


When you have the output from the former
2. run deriveL2_from3D.C

3. Fit the response ratios vs. alpha: dofits.C

4. plotresponses.C can be used to plot the responses

5. Produce txt files:


JET pT/eta/phi RESOLUTION AND JER SF:

To fill histograms for JER SF (the tag-and-probe conditions are slightly different): turn on flag "" in analyse.cc
-> you need to apply the L2 residual JEC before deriving the SF, remember check that too


When you have the outputs from that, you can run MC-checks:
--- JER/MCJER.C -> pt resolution
--- JER/MCJPR.C -> eta/phi resolition (wip?)
--- JER/MCRESP.C -> this plots the MC response <pT(reco)/pT(gen)>
--- doTxtMCJER.C -> print txt files of resolution fit parameters

For SF there are a couple of scripts:
1. These are alternatives:
--- JERSF_fits.C -> Extracts resoluiton by fitting gaussian to the dijet asymmetry distributions
--- JERSF_RMS.C -> Extracts resolution from trunct RMS of dijet asymmetry distributions
2. Do fits against alpha:
--- JERSF_fits_vsalpha.C -> this needs as an input the output from the previous step
3. Produce txt files
... JERSF_printtxt.C

To look at trigger turn-ons:
(this is has been used to merge results from different datasets/triggers, for 2023 it has been zero bias and hard probes datasets)

Plot things:
