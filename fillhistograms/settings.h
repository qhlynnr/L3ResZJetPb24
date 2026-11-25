// Some general parameters

#define REDOJES 1

// 2024 PP Reference Corrections (up-to-date)
string l1file = "jecfiles/2024PPRefUpToDateCorrections/Prompt24HIpp_V1_MC_L1FastJet_AK4PF.txt";
string l2relative = "jecfiles/2024PPRefUpToDateCorrections/Prompt24HIpp_V1_MC_L2Relative_AK4PF.txt";
string l2residual = "jecfiles/2024PPRefUpToDateCorrections/Prompt24HIpp_V1_DATA_L2Residual_AK4PF.txt";
string l2file = "";  // Don't apply L2Residual when measuring L3 residuals!

// JER - not needed for L3 PhotonJet residual analysis
string resolutionFile = "";
string scaleFactorFile = "";

// Jet veto map - how to implement?
TString vetomapFile = "jecfiles/";


// Size of the per-event jet arrays. Some 2024 ppRef HardProbes files (the
// first ~200 of each stream) have events with nref up to ~195, which
// overflowed the old value of 50 and segfaulted the batch-0 jobs.
int MAXJETS = 250;

float akradius = 0.4;
Float_t jtptmin = 15.;


float hibins[] = {-1, 0.}; //, 0., 5000.};
const int nhibins = sizeof(hibins)/sizeof(hibins[0])-1;

//float ptedges[] = {80., 100., 120., 140., 80., 140., 100., 5000., 0., 5000};
float ptedges[] = {0., 5000};
const int nptbins = sizeof(ptedges)/sizeof(ptedges[0])-1;

//float etaedges[] = {-5.2, -3.9, -2.6, -1.3, 0, 1.3, 2.6, 3.9, 5.2, -5.2, 5.2};
float etaedges[] = {-2.5, 2.5, -5.2, 5.2};
const int netabins = sizeof(etaedges)/sizeof(etaedges[0])-1;

int hltedges[] = {0, 40, 60, 80, 100, 120};
const int nhltbins = sizeof(hltedges)/sizeof(hltedges[0])-1;

