// Histograms that are filled per-event.


#include "eventhistograms.h"
eventhistograms::eventhistograms(TDirectory *dir, bool ismc) {
  TDirectory *curdir = gDirectory;
  bool enter = dir->cd();
  assert(enter);
  this->dir = dir;
  
  this->isMC = ismc;
  
  // Event quality histograms
  event_vz = new TH1D("vz", "vz;;", 100, -20, 20);
  event_hiBin = new TH1D("hiBin", "hiBin;;", 120, 0, 120);
  
  // y-axis is log10(weight), not weight - see histograms.C for why (a plain
  // 0-0.1 linear axis missed the entire weight distribution).
  if (isMC)  event_pthatwsgenweight = new TH2D("pthatvsgenwt","",200,0,1500,140,-10,4);
  
}
void eventhistograms::Write() {
  dir->cd();
  dir->Write();
  }
