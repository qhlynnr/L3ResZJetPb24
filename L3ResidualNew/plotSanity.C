// Purpose: data vs MC sanity plots for the photon+jet selection, from the
// merged fillhistograms outputs:
//   1. photon pT        (photon_pt), drawn 40-300 GeV
//   2. delta phi(gamma, probe jet)  (photonjet_dphi), drawn 0-pi; the
//      away-side jet definition requires delta phi > 2pi/3, so below that
//      it is empty
//
// Data and MC are each normalised per jet, 1/N_jets dN/dx (every selected
// event has exactly one probe jet, so N_jets = number of selected events).
// This is a shape comparison only - the MC weight is a generator weight, not
// a luminosity weight. Each plot has a data/MC ratio panel, linear y axis,
// and is saved as .png and .pdf.
//
// Usage:
//   root -l -b -q 'plotSanity.C("Data.root", "MC.root", "outdir")'
//   root -l -b -q 'plotSanity.C("Data.root", "MC.root", "outdir", "eta_-2.5_2.5")'

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"

#include <iostream>

using namespace std;

// Draw one data/MC comparison with a ratio panel and save it.
static void drawComparison(TH1D *hData, TH1D *hMC, const TString &xTitle,
                           const TString &label, const TString &yTitle,
                           double xMin, double xMax, const TString &outBase) {
  // 1/N_jets dN/dx: unit integral, divided by bin width
  hData->Scale(1.0 / hData->Integral(), "width");
  hMC->Scale(1.0 / hMC->Integral(), "width");
  // Keep the range inside the histogram, or SetRangeUser pulls in the
  // overflow bin.
  xMax = min(xMax, hData->GetXaxis()->GetXmax());
  const double xMaxRange = xMax - 1e-6 * (xMax - xMin);
  hData->GetXaxis()->SetRangeUser(xMin, xMaxRange);
  hMC->GetXaxis()->SetRangeUser(xMin, xMaxRange);

  TCanvas c("c_" + outBase, "", 800, 800);
  TPad top("top", "", 0, 0.30, 1, 1);
  TPad bot("bot", "", 0, 0, 1, 0.30);
  top.SetBottomMargin(0.02);
  top.SetLeftMargin(0.14);
  bot.SetTopMargin(0.03);
  bot.SetBottomMargin(0.32);
  bot.SetLeftMargin(0.14);
  top.Draw();
  bot.Draw();

  // ---- top: distributions
  top.cd();
  const int dataColor = kOrange + 7;
  const int mcColor = kAzure + 2;
  hData->SetMarkerStyle(kFullCircle);
  hData->SetMarkerSize(0.9);
  hData->SetMarkerColor(dataColor);
  hData->SetLineColor(dataColor);
  hMC->SetMarkerStyle(kFullSquare);
  hMC->SetMarkerSize(0.8);
  hMC->SetMarkerColor(mcColor);
  hMC->SetLineColor(mcColor);

  const double ymax = max(hData->GetMaximum(), hMC->GetMaximum());
  hMC->SetMaximum(ymax * 1.35);
  hMC->SetMinimum(0);
  hMC->SetTitle("");
  hMC->GetYaxis()->SetTitle(yTitle);
  hMC->GetYaxis()->SetTitleSize(0.05);
  hMC->GetYaxis()->SetTitleOffset(1.3);
  hMC->GetYaxis()->SetLabelSize(0.045);
  hMC->GetXaxis()->SetLabelSize(0);
  hMC->Draw("PE");
  hData->Draw("PE SAME");

  TLegend leg(0.60, 0.72, 0.88, 0.88);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextSize(0.04);
  leg.AddEntry(hData, "Data", "PE");
  leg.AddEntry(hMC, "MC", "PE");
  leg.Draw();

  TLatex tex;
  tex.SetNDC();
  tex.SetTextSize(0.045);
  tex.DrawLatex(0.16, 0.92, "#bf{CMS} #it{Preliminary}");
  tex.SetTextSize(0.035);
  tex.DrawLatex(0.18, 0.84, "pp Reference 2024, #gamma+jet");
  tex.DrawLatex(0.18, 0.79, label);

  // ---- bottom: data/MC
  bot.cd();
  TH1D *hRatio = (TH1D *)hData->Clone("ratio_" + outBase);
  hRatio->Divide(hMC);
  hRatio->SetTitle("");
  hRatio->SetMinimum(0.5);
  hRatio->SetMaximum(1.5);
  hRatio->GetYaxis()->SetTitle("Data / MC");
  hRatio->GetYaxis()->SetNdivisions(505);
  hRatio->GetYaxis()->SetTitleSize(0.11);
  hRatio->GetYaxis()->SetTitleOffset(0.55);
  hRatio->GetYaxis()->SetLabelSize(0.10);
  hRatio->GetXaxis()->SetTitle(xTitle);
  hRatio->GetXaxis()->SetTitleSize(0.12);
  hRatio->GetXaxis()->SetLabelSize(0.10);
  hRatio->SetMarkerColor(kBlack);
  hRatio->SetLineColor(kBlack);
  hRatio->Draw("PE");
  TLine one(xMin, 1, xMax, 1);
  one.SetLineStyle(kDashed);
  one.SetLineColor(kGray + 2);
  one.Draw();

  c.SaveAs(outBase + ".png");
  c.SaveAs(outBase + ".pdf");
}

void plotSanity(TString dataFile, TString mcFile, TString outDir = ".",
                TString etaFolder = "eta_-5.2_5.2",
                TString hibinFolder = "hibin_-1.0_0.0") {
  gStyle->SetOptStat(0);
  gSystem->mkdir(outDir, kTRUE);

  TFile *fData = TFile::Open(dataFile);
  TFile *fMC = TFile::Open(mcFile);
  if (!fData || fData->IsZombie() || !fMC || fMC->IsZombie()) {
    cout << "ERROR: cannot open " << dataFile << " or " << mcFile << endl;
    return;
  }

  const TString dir = hibinFolder + "/" + etaFolder + "/";
  // Probe-jet eta range of this folder, for the plot label.
  float etaLo = 0, etaHi = 0;
  sscanf(etaFolder.Data(), "eta_%f_%f", &etaLo, &etaHi);
  const TString label =
      Form("p_{T}^{#gamma} > 50 GeV, |#eta^{#gamma}| < 1.3, %.1f < #eta^{jet} < %.1f",
           etaLo, etaHi);

  struct Plot {
    const char *hist, *xTitle, *tag;
    const char *yTitle;
    int rebin;
    double xMin, xMax; // drawn x range
  };
  const Plot plots[] = {
      {"photon_pt", "p_{T}^{#gamma} (GeV)", "photon_pt",
       "1/N_{jets} dN/dp_{T}^{#gamma} (GeV^{-1})", 1, 40, 300},
      {"photonjet_dphi", "#Delta#phi(#gamma, jet)", "photonjet_dphi",
       "1/N_{jets} dN/d#Delta#phi", 1, 0, 3.1416},
  };

  for (const Plot &p : plots) {
    TH1D *hData = (TH1D *)fData->Get(dir + p.hist);
    TH1D *hMC = (TH1D *)fMC->Get(dir + p.hist);
    if (!hData || !hMC) {
      cout << "ERROR: " << dir + p.hist << " missing in data or MC" << endl;
      continue;
    }
    hData = (TH1D *)hData->Clone(Form("%s_data", p.tag));
    hMC = (TH1D *)hMC->Clone(Form("%s_mc", p.tag));
    if (p.rebin > 1) {
      hData->Rebin(p.rebin);
      hMC->Rebin(p.rebin);
    }
    drawComparison(hData, hMC, p.xTitle, label, p.yTitle, p.xMin, p.xMax,
                   outDir + "/sanity_" + p.tag + "_" + etaFolder);
  }
}
