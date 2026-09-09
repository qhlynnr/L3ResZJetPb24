// Plot the histograms produced by deriveL3_from_photonjet.C
//
// Usage:
//   root -l -b -q 'plotL3_photonjet.C("L3_derived_TEST0825.root", "plots_L3residual")'

#include "TFile.h"
#include "TKey.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLine.h"
#include "TString.h"
#include "TRegexp.h"
#include "TSystem.h"
#include "TStyle.h"
#include <iostream>
#include <vector>
#include <set>

using namespace std;

TFile *_f = nullptr;
TString _outdir;

TH1D* get1D(TString name) {
  TH1D* h = (TH1D*)_f->Get(name);
  if (!h) cout << "  [missing] " << name << endl;
  return h;
}

// Overlay MC (blue) vs Data (red) on top, MC/Data ratio on bottom.
void drawOverlayWithRatio(TH1D* hmc, TH1D* hdata, TH1D* hratio,
                           TString title, TString xtitle, TString outname) {
  if (!hmc || !hdata) return;

  TH1D* ratio = hratio;
  bool ownRatio = false;
  if (!ratio) {
    ratio = (TH1D*)hmc->Clone(outname + "_ratio");
    ratio->Divide(hdata);
    ownRatio = true;
  }

  TCanvas c("c", "c", 700, 700);
  TPad *p1 = new TPad("p1", "p1", 0, 0.32, 1, 1.0);
  TPad *p2 = new TPad("p2", "p2", 0, 0.0, 1, 0.32);
  p1->SetBottomMargin(0.02);
  p2->SetTopMargin(0.03);
  p2->SetBottomMargin(0.32);
  p1->Draw();
  p2->Draw();

  p1->cd();
  hmc->SetLineColor(kBlue + 1);
  hmc->SetMarkerColor(kBlue + 1);
  hmc->SetMarkerStyle(20);
  hdata->SetLineColor(kRed + 1);
  hdata->SetMarkerColor(kRed + 1);
  hdata->SetMarkerStyle(21);

  hmc->SetTitle(title + ";;Balance");
  hmc->GetXaxis()->SetLabelSize(0);
  hmc->Draw("E1");
  hdata->Draw("E1 SAME");

  TLegend *leg = new TLegend(0.65, 0.72, 0.88, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(hmc, "MC", "lep");
  leg->AddEntry(hdata, "Data", "lep");
  leg->Draw();

  p2->cd();
  ratio->SetTitle(";" + xtitle + ";MC / Data");
  ratio->SetLineColor(kBlack);
  ratio->SetMarkerColor(kBlack);
  ratio->SetMarkerStyle(20);
  ratio->GetYaxis()->SetRangeUser(0.5, 1.5);
  ratio->GetYaxis()->SetNdivisions(505);
  ratio->GetXaxis()->SetLabelSize(0.09);
  ratio->GetYaxis()->SetLabelSize(0.09);
  ratio->GetXaxis()->SetTitleSize(0.11);
  ratio->GetYaxis()->SetTitleSize(0.11);
  ratio->GetYaxis()->SetTitleOffset(0.45);
  ratio->Draw("E1");

  TLine *line = new TLine(ratio->GetXaxis()->GetXmin(), 1.0,
                           ratio->GetXaxis()->GetXmax(), 1.0);
  line->SetLineStyle(2);
  line->Draw();

  c.SaveAs(_outdir + "/" + outname + ".pdf");
  c.SaveAs(_outdir + "/" + outname + ".png");

  if (ownRatio) delete ratio;
}

// Overlay a family of ratio histograms (e.g. one per alpha cut) with a color gradient.
void drawFamily(vector<TH1D*> hists, vector<TString> labels, TString title,
                 TString xtitle, TString ytitle, TString outname) {
  if (hists.empty()) return;
  TCanvas c("c", "c", 700, 600);
  TLegend *leg = new TLegend(0.65, 0.55, 0.88, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);

  int colors[] = {kBlue+2, kAzure+1, kCyan+2, kGreen+2, kOrange+1,
                   kOrange+7, kRed+1, kMagenta+1, kViolet+7, kBlack};
  int ncolors = sizeof(colors) / sizeof(colors[0]);

  for (size_t i = 0; i < hists.size(); ++i) {
    TH1D* h = hists[i];
    if (!h) continue;
    h->SetLineColor(colors[i % ncolors]);
    h->SetMarkerColor(colors[i % ncolors]);
    h->SetMarkerStyle(20);
    h->SetTitle(title + ";" + xtitle + ";" + ytitle);
    h->GetYaxis()->SetRangeUser(0.5, 1.5);
    h->Draw(i == 0 ? "E1" : "E1 SAME");
    leg->AddEntry(h, labels[i], "lep");
  }
  TLine *line = new TLine(hists[0]->GetXaxis()->GetXmin(), 1.0,
                           hists[0]->GetXaxis()->GetXmax(), 1.0);
  line->SetLineStyle(2);
  line->Draw();
  leg->Draw();

  c.SaveAs(_outdir + "/" + outname + ".pdf");
  c.SaveAs(_outdir + "/" + outname + ".png");
}

void draw2D(TString name, TString outname) {
  TH2D* h = (TH2D*)_f->Get(name);
  if (!h) { cout << "  [missing] " << name << endl; return; }
  TCanvas c("c", "c", 700, 600);
  c.SetRightMargin(0.15);
  h->Draw("COLZ");
  c.SaveAs(_outdir + "/" + outname + ".pdf");
  c.SaveAs(_outdir + "/" + outname + ".png");
}

void plotL3_photonjet(TString infile = "L3_derived_TEST0825.root",
                       TString outdir = "plots_L3residual") {
  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);

  _f = TFile::Open(infile);
  if (!_f || _f->IsZombie()) {
    cout << "ERROR: cannot open " << infile << endl;
    return;
  }
  _outdir = outdir;
  gSystem->mkdir(_outdir, true);

  // ---- 1) Balance vs |eta_jet| per photon-pT bin, with L3 residual ratio ----
  cout << "\n=== Balance vs eta, per pT bin ===" << endl;
  set<TString> ptalphaTags;
  TIter next(_f->GetListOfKeys());
  TKey *key;
  while ((key = (TKey*)next())) {
    TString name = key->GetName();
    if (name.BeginsWith("mc_balance_pt")) {
      TString tag = name;
      tag.Remove(0, TString("mc_balance_pt").Length());
      ptalphaTags.insert(tag);
    }
  }
  for (auto tag : ptalphaTags) {
    TH1D* hmc   = get1D("mc_balance_pt" + tag);
    TH1D* hdata = get1D("dt_balance_pt" + tag);
    TH1D* hrat  = get1D("ratio_pt" + tag);
    drawOverlayWithRatio(hmc, hdata, hrat,
                          "Balance vs |#eta_{jet}|, p_{T}^{#gamma} bin " + tag,
                          "|#eta_{jet}|", "balance_vs_eta_pt" + tag);
  }

  // ---- 2) Balance & L3 residual vs photon pT, one plot per alpha cut ----
  cout << "\n=== Balance vs photon pT, per alpha cut ===" << endl;
  vector<TH1D*> ratioFamily_pt;
  vector<TString> ratioLabels_pt;
  for (int a = 1; a <= 20; ++a) {
    TH1D* hmc   = (TH1D*)_f->Get(Form("balance_vsphotonpt_mc_alpha%d", a));
    TH1D* hdata = (TH1D*)_f->Get(Form("balance_vsphotonpt_data_alpha%d", a));
    TH1D* hrat  = (TH1D*)_f->Get(Form("ratio_vsphotonpt_alpha%d", a));
    if (!hmc || !hdata) break;
    drawOverlayWithRatio(hmc, hdata, hrat,
                          Form("Balance vs Photon p_{T} (alpha cut %d)", a),
                          "Photon p_{T} (GeV)", Form("balance_vs_photonpt_alpha%d", a));
    if (hrat) {
      ratioFamily_pt.push_back(hrat);
      ratioLabels_pt.push_back(Form("alpha bin %d", a));
    }
  }
  drawFamily(ratioFamily_pt, ratioLabels_pt,
             "L3 Residual (MC/Data) vs Photon p_{T}, all alpha cuts",
             "Photon p_{T} (GeV)", "MC / Data", "ratio_vs_photonpt_allalpha");

  // ---- 3) Balance & ratio vs "derived jet pT" family (NOTE: see caveat below) ----
  cout << "\n=== Balance vs derived jet pT, per alpha cut ===" << endl;
  vector<TH1D*> ratioFamily_jet;
  vector<TString> ratioLabels_jet;
  for (int a = 1; a <= 20; ++a) {
    TH1D* hmc   = (TH1D*)_f->Get(Form("balance_vsjetpt_mc_alpha%d", a));
    TH1D* hdata = (TH1D*)_f->Get(Form("balance_vsjetpt_data_alpha%d", a));
    TH1D* hrat  = (TH1D*)_f->Get(Form("ratio_vsjetpt_alpha%d", a));
    if (!hmc || !hdata) break;
    drawOverlayWithRatio(hmc, hdata, hrat,
                          Form("Balance vs Derived Jet p_{T} (alpha cut %d)", a),
                          "Derived Jet p_{T} (GeV)",
                          Form("balance_vs_jetpt_alpha%d", a));
    if (hrat) {
      ratioFamily_jet.push_back(hrat);
      ratioLabels_jet.push_back(Form("alpha bin %d", a));
    }
  }
  drawFamily(ratioFamily_jet, ratioLabels_jet,
             "L3 Residual (MC/Data) vs Derived Jet p_{T}",
             "Derived Jet p_{T} (GeV)", "MC / Data",
             "ratio_vs_jetpt_allalpha");

  // ---- 4) 2D summary maps ----
  cout << "\n=== 2D summary maps ===" << endl;
  draw2D("balanceMap_mc", "balanceMap_mc");
  draw2D("balanceMap_data", "balanceMap_data");
  draw2D("l3resMap", "l3resMap");

  cout << "\nDone. Plots written to: " << _outdir << "/" << endl;

  _f->Close();
}
