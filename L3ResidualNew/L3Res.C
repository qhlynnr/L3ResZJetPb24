// Purpose: fit the active L3 residual response workflow in pTref.
//
// Step 1. Read one or more derived ROOT files from deriveL3_from_photonjet.C.
// Step 2. Fit the alpha dependence in every pTref bin to extract kFSR and the
//         alpha->0 corrected response for each input independently.
// Step 3. Combine the selected pTref points from all inputs and run the
//         final full direct pTref fit used by the export step.
// Step 4. Save the fit ROOT file, metadata, and validation plots needed by
//         createL2L3ResTextFile.C.

#include "TAxis.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TNamed.h"
#include "TObjArray.h"
#include "TPad.h"
#include "TParameter.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TStyle.h"
#include "TSystem.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "CMS_lumi.C"
#include "L3ResCommon.h"
#include "tdrStyle.C"

using namespace std;

static void useCmsPlotStyle(const string &lumiLabel) {
  setTDRStyle();
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  writeExtraText = true;
  extraText = "Preliminary";
  lumi_sqrtS = Form("%s, #sqrt{s} = 5.36 TeV", lumiLabel.c_str());
}

static void drawCmsLabel(TCanvas &canvas) {
  canvas.cd();
  CMS_lumi(&canvas, 0, 0);
  canvas.Modified();
  canvas.Update();
}

static void styleFrameAxes(TH1D *frame, double yTitleOffset = 1.35) {
  if (!frame)
    return;
  frame->GetXaxis()->SetTitleSize(0.046);
  frame->GetYaxis()->SetTitleSize(0.046);
  frame->GetXaxis()->SetLabelSize(0.034);
  frame->GetYaxis()->SetLabelSize(0.034);
  frame->GetXaxis()->SetTitleOffset(1.05);
  frame->GetYaxis()->SetTitleOffset(yTitleOffset);
}

static TString balanceAxisLabel(const PreparedInput &prepared) {
  return Form("B_{%s}", prepared.label.Data());
}

static TString responseAxisLabel(const PreparedInput &prepared) {
  return Form("R_{%s}", prepared.label.Data());
}

static TString alphaRatioAxisLabel(const PreparedInput &prepared) {
  return Form("r_{%s}(#alpha)", prepared.label.Data());
}

// Purpose: describe which pT axis a per-input plot is drawn against, so the
// same drawing code serves both the pTref series and the vsjetpt series.
struct PlotAxis {
  TString symbol; // axis symbol, e.g. "p_{T}^{ref}"
  TString tag;    // file-name tag, e.g. "ptref" (-> *_vsptref.png)
};

static PlotAxis ptRefAxis() { return {"p_{T}^{ref}", "ptref"}; }
static PlotAxis jetPtAxis() { return {"p_{T}^{jet}", "jetpt"}; }

struct SequentialFitResult {
  TF1 *f0 = nullptr;
  TF1 *f1 = nullptr;
  TF1 *f2 = nullptr;
  TF1 *f3 = nullptr;
  TF1 *f4 = nullptr;
  TF1 *finalFit = nullptr;
  int finalStatus = -1;
  bool valid = false;
};

static double fitChi2Ndf(const TF1 *fit) {
  if (!fit)
    return -1.0;
  const int ndf = fit->GetNDF();
  if (ndf <= 0)
    return -1.0;
  return fit->GetChisquare() / ndf;
}

static bool fitNeedsReview(const TF1 *fit, int fitStatus,
                           double chi2Limit = 5.0) {
  if (fitStatus != 0)
    return true;
  const double chi2PerNdf = fitChi2Ndf(fit);
  return (chi2PerNdf > chi2Limit);
}

static const char *finalPtRefFitExpr() {
  return "1./([0]+[1]/x+[2]*log(x)/x+[3]*(pow(x/[4],[5])-1)/(pow(x/"
         "[4],[5])+1)+[6]*pow(x,-0.3051)+[7]*x)";
}

static void initFinalPtRefFit(TF1 *fit, double ySeed = 1.0) {
  if (!fit)
    return;

  const double safeSeed = (ySeed > 0.0 ? ySeed : 1.0);
  fit->SetParameters(0.9722 / safeSeed, 0.7944, 2.14069, 0.10229, 10.52, 1.5550,
                     -0.72222, -2.25e-06);
  fit->SetParLimits(0, 0.5, 2.0);
  fit->SetParLimits(4, 1.0, 100.0);
  fit->SetParLimits(5, 0.1, 5.0);
}

// Purpose: keep the historical sequential pTref fits for diagnostics, but use
// the full direct pTref model as the final shared fit.
static SequentialFitResult runSequentialTemplateFit(TGraphErrors *graph,
                                                    double fitMin,
                                                    double fitMax,
                                                    const TString &fitTag) {
  SequentialFitResult result;
  if (!graph || graph->GetN() < 3 || fitMax <= fitMin)
    return result;

  const double ySeed = graphMeanY(graph);
  const double oneOverXLimit = std::max(0.5, 0.40 * fitMin / 10.0);

  result.f0 =
      new TF1(Form("f_ptref_f0_%s", fitTag.Data()), "[0]", fitMin, fitMax);
  result.f1 = new TF1(Form("f_ptref_f1_%s", fitTag.Data()),
                      "[0]+[1]*log10(0.01*x)", fitMin, fitMax);
  result.f2 = new TF1(Form("f_ptref_f2_%s", fitTag.Data()),
                      "[0]+[1]*log10(0.01*x)+[2]/(x/10.)", fitMin, fitMax);
  result.f3 =
      new TF1(Form("f_ptref_f3_%s", fitTag.Data()),
              "[0]+[1]*log10(0.01*x)+[2]*pow(log10(0.01*x),2)", fitMin, fitMax);
  result.f4 =
      new TF1(Form("f_ptref_f4_%s", fitTag.Data()),
              "[0]+[1]*log10(0.01*x)+[2]*pow(log10(0.01*x),2)+[3]/(x/10.)",
              fitMin, fitMax);
  result.finalFit = new TF1(Form("f_ptref_final_%s", fitTag.Data()),
                            finalPtRefFitExpr(), fitMin, fitMax);

  result.f0->SetParameter(0, ySeed);
  result.f0->SetParLimits(0, 0.5, 1.5);
  graph->Fit(result.f0, "QRN");

  result.f1->SetParameters(result.f0->GetParameter(0), -0.01);
  result.f1->SetParLimits(0, 0.5, 1.5);
  result.f1->SetParLimits(1, -0.20, 0.20);
  graph->Fit(result.f1, "QRN");

  result.f2->SetParameters(result.f1->GetParameter(0),
                           result.f1->GetParameter(1), 0.02);
  result.f2->SetParLimits(0, 0.5, 1.5);
  result.f2->SetParLimits(1, -0.20, 0.20);
  result.f2->SetParLimits(2, -oneOverXLimit, oneOverXLimit);
  graph->Fit(result.f2, "QRN");

  result.f3->SetParameters(result.f1->GetParameter(0),
                           result.f1->GetParameter(1), 0.005);
  result.f3->SetParLimits(0, 0.5, 1.5);
  result.f3->SetParLimits(1, -0.20, 0.20);
  result.f3->SetParLimits(2, -0.15, 0.15);
  graph->Fit(result.f3, "QRN");

  result.f4->SetParameters(
      result.f3->GetParameter(0), result.f3->GetParameter(1),
      result.f3->GetParameter(2), result.f2->GetParameter(2));
  result.f4->SetParLimits(0, 0.5, 1.5);
  result.f4->SetParLimits(1, -0.20, 0.20);
  result.f4->SetParLimits(2, -0.15, 0.15);
  result.f4->SetParLimits(3, -oneOverXLimit, oneOverXLimit);
  graph->Fit(result.f4, "QRN");

  initFinalPtRefFit(result.finalFit, ySeed);
  result.finalStatus = (int)graph->Fit(result.finalFit, "QRN");
  result.finalStatus = (int)graph->Fit(result.finalFit, "QRN");

  result.valid = std::isfinite(result.finalFit->GetParameter(0)) &&
                 std::isfinite(result.finalFit->GetParameter(7));
  return result;
}

void L3Res(TString inFileL3Derived = "L3_derived.root",
           TString sampleTypesCSV = "", TString inputPtRangesCSV = "",
           string outfilename = "L3Res_photonjet",
           string runLabel = "2024ppRef", string lumiLabel = "pp 480.4 pb^{-1}",
           int refAlphaBin = 5, double fitAlphaMin = 0.0,
           double fitAlphaMax = 0.4, string outBaseDir = "L3Residual") {

  useCmsPlotStyle(lumiLabel);

  const string outfolder = outBaseDir + "/" + outfilename;
  const string pdfFolder = outfolder + "/pdf";
  gSystem->mkdir(outfolder.c_str(), kTRUE);
  gSystem->mkdir(pdfFolder.c_str(), kTRUE);

  vector<TString> inputFiles = splitCsvTokens(inFileL3Derived);
  if (inputFiles.empty())
    inputFiles.push_back(inFileL3Derived);

  vector<SampleType> sampleTypes;
  string sampleError;
  if (!parseSampleTypesCSV(sampleTypesCSV, inputFiles.size(), sampleTypes,
                           sampleError)) {
    cout << "ERROR: " << sampleError << endl;
    return;
  }

  const vector<TString> fitWindowTokens = splitCsvTokens(inputPtRangesCSV);
  const int colors[] = {kBlue + 1,    kRed + 1,    kGreen + 2,
                        kMagenta + 2, kOrange + 7, kCyan + 2};
  const int markers[] = {kFullCircle,  kFullSquare,       kFullTriangleUp,
                         kFullDiamond, kFullTriangleDown, kFullCross};
  const int colorCount = sizeof(colors) / sizeof(colors[0]);
  const int markerCount = sizeof(markers) / sizeof(markers[0]);

  TFile *outputFile =
      new TFile(Form("%s/%s_fit.root", outfolder.c_str(), outfilename.c_str()),
                "RECREATE");

  // Purpose: clone a required histogram and stop the macro with a clear error
  // if the contract is broken.
  auto cloneRequiredHist = [](TFile *inputFile, const TString &name,
                              const TString &context) -> TH1D * {
    TH1D *hist = inputFile ? (TH1D *)inputFile->Get(name) : nullptr;
    if (!hist) {
      cout << "ERROR: Missing required histogram '" << name << "' in "
           << context << endl;
      return nullptr;
    }
    TH1D *clone = (TH1D *)hist->Clone(Form("%s_clone", name.Data()));
    clone->SetDirectory(nullptr);
    return clone;
  };

  // Purpose: load the full alpha series for a required histogram family.
  auto loadSeries = [&](TFile *inputFile, const TString &prefix,
                        const TString &context,
                        bool required = true) -> map<int, TH1D *> {
    map<int, TH1D *> series;
    for (int alphaBin = 1; alphaBin <= 50; ++alphaBin) {
      TH1D *sourceHist =
          inputFile
              ? (TH1D *)inputFile->Get(Form("%s%d", prefix.Data(), alphaBin))
              : nullptr;
      TH1D *hist = nullptr;
      if (sourceHist) {
        hist = (TH1D *)sourceHist->Clone(
            Form("%s%d_clone", prefix.Data(), alphaBin));
        hist->SetDirectory(nullptr);
      }
      if (!hist) {
        if (alphaBin == 1 && required) {
          cout << "ERROR: Missing required histogram '"
               << Form("%s%d", prefix.Data(), alphaBin) << "' in " << context
               << endl;
          series.clear();
        }
        break;
      }
      series[alphaBin] = hist;
    }
    return series;
  };

  auto normalizeSeriesByReference =
      [&](const map<int, TH1D *> &ratioSeries, int referenceAlphaBin,
          const TString &prefix, const TString &context) -> map<int, TH1D *> {
    map<int, TH1D *> normalizedSeries;
    auto refIt = ratioSeries.find(referenceAlphaBin);
    if (refIt == ratioSeries.end() || !refIt->second) {
      cout << "ERROR: Missing required reference alpha histogram for "
              "normalized series '"
           << prefix << "' in " << context << endl;
      return normalizedSeries;
    }

    TH1D *referenceHist = refIt->second;
    for (const auto &entry : ratioSeries) {
      TH1D *hist = (TH1D *)entry.second->Clone(
          Form("%s%d_clone", prefix.Data(), entry.first));
      hist->SetDirectory(nullptr);

      for (int bin = 1; bin <= hist->GetNbinsX(); ++bin) {
        const double value = hist->GetBinContent(bin);
        const double error = hist->GetBinError(bin);
        const double referenceValue = referenceHist->GetBinContent(bin);
        const double referenceError = referenceHist->GetBinError(bin);

        if (!std::isfinite(value) || !std::isfinite(referenceValue) ||
            !(referenceValue > 0.0) || !(value > 0.0)) {
          hist->SetBinContent(bin, 0.0);
          hist->SetBinError(bin, 0.0);
          continue;
        }

        if (entry.first == referenceAlphaBin) {
          hist->SetBinContent(bin, 1.0);
          hist->SetBinError(bin, 0.0);
          continue;
        }

        // ROOT tracks the per-histogram bin errors, but not the covariance
        // between cumulative-alpha ratios from different histograms. When an
        // older derived file does not already store ratio_norm_*, use a
        // conservative independent numerator/denominator propagation here.
        const double normalizedValue = value / referenceValue;
        const double relativeError = sqrt(
            pow(error / value, 2) + pow(referenceError / referenceValue, 2));
        hist->SetBinContent(bin, normalizedValue);
        hist->SetBinError(bin, normalizedValue * relativeError);
      }

      normalizedSeries[entry.first] = hist;
    }

    return normalizedSeries;
  };

  // Purpose: draw the per-alpha ratio overlay used to inspect the stored pTref
  // contract.
  auto saveRatioOverlay = [&](const PreparedInput &prepared,
                              const map<int, TH1D *> &ratioSeries,
                              TH1D *refHist, const vector<double> &alphaEdges,
                              const string &inputPdfFolder,
                              const PlotAxis &axis) {
    TCanvas canvas(
        Form("cRatioOverlay_%s_%s", prepared.token.Data(), axis.tag.Data()),
        Form("cRatioOverlay_%s_%s", prepared.token.Data(), axis.tag.Data()),
        900, 700);
    canvas.SetLogx();

    TH1D *frame = new TH1D(
        Form("hRatioOverlay_%s_%s", prepared.token.Data(), axis.tag.Data()),
        Form(";%s (GeV);%s", axis.symbol.Data(),
             responseAxisLabel(prepared).Data()),
        100, refHist->GetXaxis()->GetBinLowEdge(1),
        refHist->GetXaxis()->GetBinLowEdge(refHist->GetNbinsX() + 1));
    frame->SetMinimum(0.7);
    frame->SetMaximum(1.5);
    frame->GetXaxis()->SetMoreLogLabels(kTRUE);
    frame->GetXaxis()->SetNoExponent(kTRUE);
    styleFrameAxes(frame);
    frame->Draw();

    TLine unity(frame->GetXaxis()->GetXmin(), 1.0, frame->GetXaxis()->GetXmax(),
                1.0);
    unity.SetLineStyle(kDashed);
    unity.SetLineColor(kGray + 1);
    unity.Draw("SAME");

    TLegend legend(0.56, 0.55, 0.88, 0.88);
    legend.SetBorderSize(0);
    legend.SetFillStyle(0);
    legend.SetTextFont(42);
    legend.SetTextSize(0.028);

    int colorIndex = 0;
    for (const auto &entry : ratioSeries) {
      TH1D *hist = entry.second;
      hist->SetMarkerStyle(kFullCircle);
      hist->SetMarkerSize(0.8);
      hist->SetMarkerColor(colors[colorIndex % colorCount]);
      hist->SetLineColor(colors[colorIndex % colorCount]);
      hist->Draw("PE1 SAME");

      const int alphaBin = entry.first;
      const double alphaCut =
          (alphaBin >= 1 && alphaBin < (int)alphaEdges.size())
              ? alphaEdges[alphaBin]
              : alphaBin;
      legend.AddEntry(hist, Form("#alpha < %.2f", alphaCut), "PE");
      ++colorIndex;
    }

    TLatex label;
    label.SetNDC();
    label.SetTextFont(42);
    label.SetTextSize(0.032);
    label.DrawLatex(0.18, 0.86, prepared.label.Data());
    if (prepared.haveEtaRange) {
      label.DrawLatex(
          0.18, 0.81,
          Form("%.3f < #eta < %.3f", prepared.etaMin, prepared.etaMax));
    }
    label.DrawLatex(0.18, 0.76,
                    Form("%s = B_{%s}^{Data} / B_{%s}^{MC}",
                         responseAxisLabel(prepared).Data(),
                         prepared.label.Data(), prepared.label.Data()));
    legend.Draw();
    drawCmsLabel(canvas);
    canvas.SaveAs(Form("%s/L3Res_%s_ratio_overlay_vs%s.png",
                       inputPdfFolder.c_str(), prepared.token.Data(),
                       axis.tag.Data()));

    delete frame;
  };

  // Purpose: draw the reference-alpha raw MC/Data response used in the alpha
  // fit.
  auto saveRefRawOverlay = [&](const PreparedInput &prepared, TH1D *rawMc,
                               TH1D *rawDt, const string &inputPdfFolder,
                               const PlotAxis &axis) {
    TCanvas canvas(
        Form("cRawOverlay_%s_%s", prepared.token.Data(), axis.tag.Data()),
        Form("cRawOverlay_%s_%s", prepared.token.Data(), axis.tag.Data()), 900,
        700);
    canvas.SetLogx();

    TH1D *frame = new TH1D(
        Form("hRawOverlay_%s_%s", prepared.token.Data(), axis.tag.Data()),
        Form(";%s (GeV);%s", axis.symbol.Data(),
             balanceAxisLabel(prepared).Data()),
        100, rawMc->GetXaxis()->GetBinLowEdge(1),
        rawMc->GetXaxis()->GetBinLowEdge(rawMc->GetNbinsX() + 1));
    double yLow = 0.7;
    double yHigh = 1.02;
    if (axis.tag != "ptref") {
      // Binned in jet pT the balance is not bounded near 1 (it grows well
      // above 1 at high jet pT from the pTref threshold), so fit the range
      // to the filled bins instead.
      yLow = 1e9;
      yHigh = -1e9;
      for (TH1D *hist : {rawMc, rawDt}) {
        for (int bin = 1; bin <= hist->GetNbinsX(); ++bin) {
          const double value = hist->GetBinContent(bin);
          if (!(value > 0.0) || !std::isfinite(value))
            continue;
          yLow = std::min(yLow, value - hist->GetBinError(bin));
          yHigh = std::max(yHigh, value + hist->GetBinError(bin));
        }
      }
      if (yLow > yHigh) {
        yLow = 0.7;
        yHigh = 1.02;
      } else {
        const double pad = 0.1 * (yHigh - yLow);
        yLow -= pad;
        yHigh += 3.0 * pad; // leave room for the labels and legend
      }
    }
    frame->SetMinimum(yLow);
    frame->SetMaximum(yHigh);
    frame->GetXaxis()->SetMoreLogLabels(kTRUE);
    frame->GetXaxis()->SetNoExponent(kTRUE);
    styleFrameAxes(frame);
    frame->Draw();

    rawMc->SetMarkerStyle(kFullCircle);
    rawMc->SetMarkerColor(kBlue + 1);
    rawMc->SetLineColor(kBlue + 1);
    rawMc->Draw("PE1 SAME");

    rawDt->SetMarkerStyle(kFullCircle);
    rawDt->SetMarkerColor(kRed + 1);
    rawDt->SetLineColor(kRed + 1);
    rawDt->Draw("PE1 SAME");

    TLegend legend(0.60, 0.74, 0.88, 0.88);
    legend.SetBorderSize(0);
    legend.SetFillStyle(0);
    legend.SetTextFont(42);
    legend.SetTextSize(0.032);
    legend.AddEntry(rawMc, "MC", "PE");
    legend.AddEntry(rawDt, "Data", "PE");
    legend.Draw();

    TLatex label;
    label.SetNDC();
    label.SetTextFont(42);
    label.SetTextSize(0.032);
    label.DrawLatex(0.18, 0.86, prepared.label.Data());
    label.DrawLatex(0.18, 0.81,
                    Form("#alpha_{ref} < %.2f (bin %d)", prepared.refAlphaValue,
                         refAlphaBin));
    label.DrawLatex(0.18, 0.76,
                    Form("B_{%s}^{Data,MC} = p_{T}^{jet} / p_{T}^{ref}",
                         prepared.label.Data()));
    if (prepared.haveEtaRange) {
      label.DrawLatex(
          0.18, 0.71,
          Form("%.3f < #eta < %.3f", prepared.etaMin, prepared.etaMax));
    }
    drawCmsLabel(canvas);
    canvas.SaveAs(Form("%s/L3Res_%s_refalpha_raw_vs%s.png",
                       inputPdfFolder.c_str(), prepared.token.Data(),
                       axis.tag.Data()));

    delete frame;
  };

  // Purpose: draw the alpha fit in one pTref bin and keep the fitted intercept
  // visible.
  auto saveAlphaFitPlot = [&](const PreparedInput &prepared, int ptBin,
                              double ptLow, double ptHigh, TGraph *graph,
                              TF1 *fitLine, double kfsr, double kfsrError,
                              const string &alphaFolder, const PlotAxis &axis) {
    TCanvas canvas(Form("cAlphaFit_%s_%s_%d", prepared.token.Data(),
                        axis.tag.Data(), ptBin),
                   Form("cAlphaFit_%s_%s_%d", prepared.token.Data(),
                        axis.tag.Data(), ptBin),
                   800, 600);
    TH1D *frame =
        new TH1D(Form("hAlphaFit_%s_%s_%d", prepared.token.Data(),
                      axis.tag.Data(), ptBin),
                 Form(";#alpha;%s", alphaRatioAxisLabel(prepared).Data()), 100,
                 0.0, 0.5);
    frame->SetMinimum(0.97);
    frame->SetMaximum(1.05);
    styleFrameAxes(frame, 1.30);
    frame->Draw();

    TLine unity(0.0, 1.0, 0.5, 1.0);
    unity.SetLineStyle(kDashed);
    unity.SetLineColor(kGray + 1);
    unity.Draw("SAME");

    graph->SetMarkerStyle(kFullCircle);
    graph->SetMarkerColor(kBlue + 1);
    graph->SetLineColor(kBlue + 1);
    graph->Draw("PZ SAME");

    fitLine->SetLineColor(kRed + 1);
    fitLine->SetLineWidth(2);
    fitLine->Draw("SAME");

    TLatex label;
    label.SetNDC();
    label.SetTextFont(42);
    label.SetTextSize(0.035);
    label.DrawLatex(0.18, 0.86, prepared.label.Data());
    label.DrawLatex(0.18, 0.81,
                    Form("%s: %.0f-%.0f GeV", axis.symbol.Data(), ptLow,
                         ptHigh));
    label.DrawLatex(0.18, 0.76,
                    Form("#alpha_{ref} < %.2f, fit %.2f #leq #alpha #leq %.2f",
                         prepared.refAlphaValue, fitAlphaMin, fitAlphaMax));
    label.DrawLatex(0.18, 0.71,
                    Form("%s = %s(#alpha) / %s(#alpha_{ref})",
                         alphaRatioAxisLabel(prepared).Data(),
                         responseAxisLabel(prepared).Data(),
                         responseAxisLabel(prepared).Data()));
    label.DrawLatex(
        0.18, 0.66,
        Form("f(#alpha) = k_{FSR} + m#alpha,  k_{FSR} = %.4f #pm %.4f", kfsr,
             kfsrError));
    label.DrawLatex(0.18, 0.61,
                    Form("m = %.4f #pm %.4f", fitLine->GetParameter(1),
                         fitLine->GetParError(1)));
    drawCmsLabel(canvas);
    // Keep the historical pTref name (..._alphaFit_pt60to65.png); the jet-pT
    // series gets ..._alphaFit_jetpt60to65.png.
    const TString binStem = (axis.tag == "ptref") ? TString("pt") : axis.tag;
    canvas.SaveAs(Form("%s/L3Res_%s_kFSR_alphaFit_%s%.0fto%.0f.png",
                       alphaFolder.c_str(), prepared.token.Data(),
                       binStem.Data(), ptLow, ptHigh));

    delete frame;
  };

  // Purpose: draw the per-input kFSR summary after the alpha fits finish.
  auto saveKfsrSummary = [&](const PreparedInput &prepared, TH1D *kfsrHist,
                             const string &alphaFolder, const PlotAxis &axis) {
    TCanvas canvas(Form("cKfsr_%s_%s", prepared.token.Data(), axis.tag.Data()),
                   Form("cKfsr_%s_%s", prepared.token.Data(), axis.tag.Data()),
                   900, 700);
    canvas.SetLogx();

    TH1D *frame = new TH1D(
        Form("hKfsr_%s_%s", prepared.token.Data(), axis.tag.Data()),
        Form(";%s (GeV);k_{FSR}^{%s}", axis.symbol.Data(),
             prepared.label.Data()),
        100, kfsrHist->GetXaxis()->GetBinLowEdge(1),
        kfsrHist->GetXaxis()->GetBinLowEdge(kfsrHist->GetNbinsX() + 1));
    frame->SetMinimum(0.9);
    frame->SetMaximum(1.1);
    frame->GetXaxis()->SetMoreLogLabels(kTRUE);
    frame->GetXaxis()->SetNoExponent(kTRUE);
    styleFrameAxes(frame);
    frame->Draw();

    TLine unity(frame->GetXaxis()->GetXmin(), 1.0, frame->GetXaxis()->GetXmax(),
                1.0);
    unity.SetLineStyle(kDashed);
    unity.SetLineColor(kGray + 1);
    unity.Draw("SAME");

    kfsrHist->SetMarkerStyle(kFullCircle);
    kfsrHist->SetMarkerColor(kBlue + 1);
    kfsrHist->SetLineColor(kBlue + 1);
    kfsrHist->Draw("PE1 SAME");

    TLatex label;
    label.SetNDC();
    label.SetTextFont(42);
    label.SetTextSize(0.032);
    label.DrawLatex(0.18, 0.86, prepared.label.Data());
    label.DrawLatex(0.18, 0.81,
                    Form("#alpha_{ref} < %.2f, fit %.2f #leq #alpha #leq %.2f",
                         prepared.refAlphaValue, fitAlphaMin, fitAlphaMax));
    drawCmsLabel(canvas);
    canvas.SaveAs(Form("%s/L3Res_%s_kFSR_vs%s.png", alphaFolder.c_str(),
                       prepared.token.Data(), axis.tag.Data()));

    delete frame;
  };

  // Purpose: draw the reference-alpha ratio and the alpha->0 corrected response
  // on the same axis.
  auto saveCorrSummary = [&](const PreparedInput &prepared, TH1D *refHist,
                             TH1D *corrHist, const string &alphaFolder,
                             const PlotAxis &axis) {
    TCanvas canvas(Form("cCorr_%s_%s", prepared.token.Data(), axis.tag.Data()),
                   Form("cCorr_%s_%s", prepared.token.Data(), axis.tag.Data()),
                   900, 700);
    canvas.SetLogx();

    TH1D *frame = new TH1D(
        Form("hCorr_%s_%s", prepared.token.Data(), axis.tag.Data()),
        Form(";%s (GeV);%s", axis.symbol.Data(),
             responseAxisLabel(prepared).Data()),
        100, corrHist->GetXaxis()->GetBinLowEdge(1),
        corrHist->GetXaxis()->GetBinLowEdge(corrHist->GetNbinsX() + 1));
    frame->SetMinimum(0.7);
    frame->SetMaximum(1.5);
    frame->GetXaxis()->SetMoreLogLabels(kTRUE);
    frame->GetXaxis()->SetNoExponent(kTRUE);
    styleFrameAxes(frame);
    frame->Draw();

    TLine unity(frame->GetXaxis()->GetXmin(), 1.0, frame->GetXaxis()->GetXmax(),
                1.0);
    unity.SetLineStyle(kDashed);
    unity.SetLineColor(kGray + 1);
    unity.Draw("SAME");

    refHist->SetMarkerStyle(kOpenCircle);
    refHist->SetMarkerColor(kBlue + 1);
    refHist->SetLineColor(kBlue + 1);
    refHist->Draw("PE1 SAME");

    corrHist->SetMarkerStyle(kFullCircle);
    corrHist->SetMarkerColor(kRed + 1);
    corrHist->SetLineColor(kRed + 1);
    corrHist->Draw("PE1 SAME");

    TLegend legend(0.54, 0.74, 0.88, 0.88);
    legend.SetBorderSize(0);
    legend.SetFillStyle(0);
    legend.SetTextFont(42);
    legend.SetTextSize(0.032);
    legend.AddEntry(refHist, Form("Ref #alpha < %.2f", prepared.refAlphaValue),
                    "PE");
    legend.AddEntry(corrHist, "#alpha #rightarrow 0", "PE");
    legend.Draw();

    TLatex label;
    label.SetNDC();
    label.SetTextFont(42);
    label.SetTextSize(0.032);
    label.DrawLatex(0.18, 0.86, prepared.label.Data());
    label.DrawLatex(0.18, 0.81,
                    Form("%s(0) = k_{FSR} #times %s(#alpha_{ref})",
                         responseAxisLabel(prepared).Data(),
                         responseAxisLabel(prepared).Data()));
    label.DrawLatex(0.18, 0.76,
                    Form("#alpha_{ref} < %.2f, fit %.2f #leq #alpha #leq %.2f",
                         prepared.refAlphaValue, fitAlphaMin, fitAlphaMax));
    drawCmsLabel(canvas);
    canvas.SaveAs(Form("%s/L3Res_%s_corr_alpha0_vs%s.png", alphaFolder.c_str(),
                       prepared.token.Data(), axis.tag.Data()));

    delete frame;
  };

  // Purpose: overlay the alpha->0 corrected response binned in pTref and in
  // leading-jet pT on one axis, to compare the two binnings directly.
  auto saveCorrPtRefVsJetPt = [&](const PreparedInput &prepared,
                                  const string &alphaFolder) {
    TCanvas canvas(Form("cCorrCmp_%s", prepared.token.Data()),
                   Form("cCorrCmp_%s", prepared.token.Data()), 900, 700);
    canvas.SetLogx();

    const double xMin =
        std::min(prepared.corrHist->GetXaxis()->GetBinLowEdge(1),
                 prepared.corrJetHist->GetXaxis()->GetBinLowEdge(1));
    const double xMax = std::max(
        prepared.corrHist->GetXaxis()->GetBinLowEdge(
            prepared.corrHist->GetNbinsX() + 1),
        prepared.corrJetHist->GetXaxis()->GetBinLowEdge(
            prepared.corrJetHist->GetNbinsX() + 1));
    TH1D *frame = new TH1D(
        Form("hCorrCmp_%s", prepared.token.Data()),
        Form(";p_{T}^{ref} or p_{T}^{jet} (GeV);%s(#alpha#rightarrow0)",
             responseAxisLabel(prepared).Data()),
        100, xMin, xMax);
    frame->SetMinimum(0.7);
    frame->SetMaximum(1.5);
    frame->GetXaxis()->SetMoreLogLabels(kTRUE);
    frame->GetXaxis()->SetNoExponent(kTRUE);
    styleFrameAxes(frame);
    frame->Draw();

    TLine unity(xMin, 1.0, xMax, 1.0);
    unity.SetLineStyle(kDashed);
    unity.SetLineColor(kGray + 1);
    unity.Draw("SAME");

    TH1D *ptRefHist = (TH1D *)prepared.corrHist->Clone(
        Form("corrCmp_ptref_%s", prepared.token.Data()));
    ptRefHist->SetDirectory(nullptr);
    ptRefHist->SetMarkerStyle(kFullCircle);
    ptRefHist->SetMarkerColor(kRed + 1);
    ptRefHist->SetLineColor(kRed + 1);
    ptRefHist->Draw("PE1 SAME");

    TH1D *jetHist = (TH1D *)prepared.corrJetHist->Clone(
        Form("corrCmp_jetpt_%s", prepared.token.Data()));
    jetHist->SetDirectory(nullptr);
    jetHist->SetMarkerStyle(kOpenSquare);
    jetHist->SetMarkerColor(kAzure + 2);
    jetHist->SetLineColor(kAzure + 2);
    jetHist->Draw("PE1 SAME");

    TLegend legend(0.54, 0.74, 0.88, 0.88);
    legend.SetBorderSize(0);
    legend.SetFillStyle(0);
    legend.SetTextFont(42);
    legend.SetTextSize(0.032);
    legend.AddEntry(ptRefHist, "vs p_{T}^{ref}", "PE");
    legend.AddEntry(jetHist, "vs p_{T}^{jet}", "PE");
    legend.Draw();

    TLatex label;
    label.SetNDC();
    label.SetTextFont(42);
    label.SetTextSize(0.032);
    label.DrawLatex(0.18, 0.86, prepared.label.Data());
    label.DrawLatex(0.18, 0.81,
                    Form("#alpha_{ref} < %.2f, fit %.2f #leq #alpha #leq %.2f",
                         prepared.refAlphaValue, fitAlphaMin, fitAlphaMax));
    drawCmsLabel(canvas);
    canvas.SaveAs(Form("%s/L3Res_%s_corr_alpha0_ptref_vs_jetpt.png",
                       alphaFolder.c_str(), prepared.token.Data()));

    delete ptRefHist;
    delete jetHist;
    delete frame;
  };

  // Purpose: draw the raw (not alpha->0 extrapolated) reference-alpha
  // response vs |eta|, overlaid for each pTref bin. This reuses the
  // ratio_pt<ptrange>_alpha<cut> histograms deriveL3_from_photonjet.C already
  // writes but L3Res.C otherwise never reads.
  auto saveEtaRawOverlay = [&](const PreparedInput &prepared,
                               const map<int, TH1D *> &etaSeries,
                               TAxis *ptAxisForEta,
                               const string &inputPdfFolder) {
    if (etaSeries.empty())
      return;

    TCanvas canvas(Form("cEtaOverlay_%s", prepared.token.Data()),
                   Form("cEtaOverlay_%s", prepared.token.Data()), 900, 700);

    TH1D *frame = new TH1D(
        Form("hEtaOverlay_%s", prepared.token.Data()),
        Form(";#eta;%s", responseAxisLabel(prepared).Data()), 100,
        prepared.etaMin, prepared.etaMax);
    frame->SetMinimum(0.7);
    frame->SetMaximum(1.5);
    styleFrameAxes(frame);
    frame->Draw();

    TLine unity(prepared.etaMin, 1.0, prepared.etaMax, 1.0);
    unity.SetLineStyle(kDashed);
    unity.SetLineColor(kGray + 1);
    unity.Draw("SAME");

    TLegend legend(0.56, 0.15, 0.88, 0.45);
    legend.SetBorderSize(0);
    legend.SetFillStyle(0);
    legend.SetTextFont(42);
    legend.SetTextSize(0.026);

    int colorIndex = 0;
    for (const auto &entry : etaSeries) {
      TH1D *hist = entry.second;
      hist->SetMarkerStyle(markers[colorIndex % markerCount]);
      hist->SetMarkerSize(0.8);
      hist->SetMarkerColor(colors[colorIndex % colorCount]);
      hist->SetLineColor(colors[colorIndex % colorCount]);
      hist->Draw("PE1 SAME");

      const int ptBin = entry.first;
      legend.AddEntry(hist,
                      Form("p_{T}^{ref}: %.0f-%.0f GeV",
                           ptAxisForEta->GetBinLowEdge(ptBin),
                           ptAxisForEta->GetBinLowEdge(ptBin + 1)),
                      "PE");
      ++colorIndex;
    }

    TLatex label;
    label.SetNDC();
    label.SetTextFont(42);
    label.SetTextSize(0.032);
    label.DrawLatex(0.18, 0.86, prepared.label.Data());
    label.DrawLatex(
        0.18, 0.81,
        Form("#alpha_{ref} < %.2f (bin %d), raw - not #alpha#rightarrow0 "
             "corrected",
             prepared.refAlphaValue, refAlphaBin));
    label.DrawLatex(0.18, 0.76,
                    Form("%s = B_{%s}^{Data} / B_{%s}^{MC}",
                         responseAxisLabel(prepared).Data(),
                         prepared.label.Data(), prepared.label.Data()));
    legend.Draw();
    drawCmsLabel(canvas);
    canvas.SaveAs(Form("%s/L3Res_%s_ratio_raw_vseta.png", inputPdfFolder.c_str(),
                       prepared.token.Data()));

    delete frame;
  };

  // Purpose: balance and Data/MC vs fine |eta| for one pT bin, read from the
  // balance_vsabseta_* histograms deriveL3_from_photonjet.C writes from the
  // fine-|eta| profiles (independent of the |eta| binning used by the fit).
  struct AbsEtaSeries {
    map<int, TH1D *> mc, data, ratio; // keyed by pT bin
    TAxis *ptAxis = nullptr;
  };

  auto loadAbsEtaSeries = [&](TFile *inputFile, const PlotAxis &axis,
                              const char *profileName) -> AbsEtaSeries {
    AbsEtaSeries series;
    TProfile3D *profile = (TProfile3D *)inputFile->Get(profileName);
    if (!profile)
      return series;
    series.ptAxis = (TAxis *)profile->GetXaxis()->Clone(
        Form("absEtaPtAxis_%s", axis.tag.Data()));
    for (int ptBin = 1; ptBin <= series.ptAxis->GetNbins(); ++ptBin) {
      const TString tag = Form("%s%.0fto%.0f", axis.tag.Data(),
                               series.ptAxis->GetBinLowEdge(ptBin),
                               series.ptAxis->GetBinLowEdge(ptBin + 1));
      TH1D *mc = (TH1D *)inputFile->Get("balance_vsabseta_mc_" + tag);
      TH1D *data = (TH1D *)inputFile->Get("balance_vsabseta_data_" + tag);
      TH1D *ratio = (TH1D *)inputFile->Get("ratio_vsabseta_" + tag);
      if (!mc || !data || !ratio)
        continue;
      series.mc[ptBin] = (TH1D *)mc->Clone(Form("%s_clone", mc->GetName()));
      series.data[ptBin] =
          (TH1D *)data->Clone(Form("%s_clone", data->GetName()));
      series.ratio[ptBin] =
          (TH1D *)ratio->Clone(Form("%s_clone", ratio->GetName()));
      series.mc[ptBin]->SetDirectory(nullptr);
      series.data[ptBin]->SetDirectory(nullptr);
      series.ratio[ptBin]->SetDirectory(nullptr);
    }
    return series;
  };

  // Purpose: y-range covering the filled bins (value +- error) of a set of
  // histograms, padded, with a minimum span so flat series stay readable.
  auto filledRange = [](const vector<TH1D *> &hists, double minSpan,
                        double &yLow, double &yHigh) {
    yLow = 1e9;
    yHigh = -1e9;
    for (TH1D *hist : hists) {
      for (int bin = 1; bin <= hist->GetNbinsX(); ++bin) {
        const double value = hist->GetBinContent(bin);
        if (!(value > 0.0) || !std::isfinite(value))
          continue;
        yLow = std::min(yLow, value - hist->GetBinError(bin));
        yHigh = std::max(yHigh, value + hist->GetBinError(bin));
      }
    }
    if (yLow > yHigh) {
      yLow = 0.5;
      yHigh = 1.5;
      return;
    }
    if (yHigh - yLow < minSpan) {
      const double mid = 0.5 * (yLow + yHigh);
      yLow = mid - 0.5 * minSpan;
      yHigh = mid + 0.5 * minSpan;
    }
    const double pad = 0.1 * (yHigh - yLow);
    yLow -= pad;
    yHigh += pad;
  };

  // Purpose: one pT bin - MC/Data balance vs |eta| on top, Data/MC below.
  auto saveAbsEtaBinPlot = [&](const PreparedInput &prepared,
                               const AbsEtaSeries &series, int ptBin,
                               const string &etaFolder, const PlotAxis &axis) {
    TH1D *mc = series.mc.at(ptBin);
    TH1D *data = series.data.at(ptBin);
    TH1D *ratio = series.ratio.at(ptBin);
    const double ptLow = series.ptAxis->GetBinLowEdge(ptBin);
    const double ptHigh = series.ptAxis->GetBinLowEdge(ptBin + 1);
    const double etaMin = mc->GetXaxis()->GetXmin();
    const double etaMax = mc->GetXaxis()->GetXmax();
    const TString id =
        Form("%s_%s%d", prepared.token.Data(), axis.tag.Data(), ptBin);

    TCanvas canvas(Form("cAbsEta_%s", id.Data()), Form("cAbsEta_%s", id.Data()),
                   800, 800);
    TPad top(Form("pTop_%s", id.Data()), "", 0.0, 0.32, 1.0, 1.0);
    TPad bottom(Form("pBot_%s", id.Data()), "", 0.0, 0.0, 1.0, 0.32);
    top.SetBottomMargin(0.02);
    bottom.SetTopMargin(0.03);
    bottom.SetBottomMargin(0.32);
    top.Draw();
    bottom.Draw();

    top.cd();
    double yLow, yHigh;
    filledRange({mc, data}, 0.1, yLow, yHigh);
    yHigh += 0.6 * (yHigh - yLow); // room for labels and legend
    TH1D *frameTop =
        new TH1D(Form("hAbsEtaTop_%s", id.Data()),
                 Form(";;%s", balanceAxisLabel(prepared).Data()), 100, etaMin,
                 etaMax);
    frameTop->SetMinimum(yLow);
    frameTop->SetMaximum(yHigh);
    styleFrameAxes(frameTop);
    frameTop->GetXaxis()->SetLabelSize(0);
    frameTop->GetYaxis()->SetTitleSize(0.060);
    frameTop->GetYaxis()->SetLabelSize(0.045);
    frameTop->GetYaxis()->SetTitleOffset(1.05);
    frameTop->Draw();

    mc->SetMarkerStyle(kOpenCircle);
    mc->SetMarkerColor(kBlue + 1);
    mc->SetLineColor(kBlue + 1);
    mc->Draw("PE1 SAME");
    data->SetMarkerStyle(kFullCircle);
    data->SetMarkerColor(kRed + 1);
    data->SetLineColor(kRed + 1);
    data->Draw("PE1 SAME");

    TLegend legend(0.66, 0.66, 0.90, 0.84);
    legend.SetBorderSize(0);
    legend.SetFillStyle(0);
    legend.SetTextFont(42);
    legend.SetTextSize(0.045);
    legend.AddEntry(mc, "MC", "PE");
    legend.AddEntry(data, "Data", "PE");
    legend.Draw();

    TLatex label;
    label.SetNDC();
    label.SetTextFont(42);
    label.SetTextSize(0.042);
    label.DrawLatex(0.18, 0.84, prepared.label.Data());
    label.DrawLatex(0.18, 0.78, Form("%s: %.0f-%.0f GeV", axis.symbol.Data(),
                                     ptLow, ptHigh));
    label.DrawLatex(0.18, 0.72,
                    Form("#alpha_{ref} < %.2f (bin %d), raw - not "
                         "#alpha#rightarrow0 corrected",
                         prepared.refAlphaValue, refAlphaBin));
    CMS_lumi(&top, 0, 0);

    bottom.cd();
    filledRange({ratio}, 0.1, yLow, yHigh);
    TH1D *frameBot = new TH1D(Form("hAbsEtaBot_%s", id.Data()),
                              ";|#eta_{jet}|;Data/MC", 100, etaMin, etaMax);
    frameBot->SetMinimum(std::min(yLow, 0.98));
    frameBot->SetMaximum(std::max(yHigh, 1.02));
    frameBot->GetXaxis()->SetTitleSize(0.13);
    frameBot->GetXaxis()->SetLabelSize(0.10);
    frameBot->GetXaxis()->SetTitleOffset(1.0);
    frameBot->GetYaxis()->SetTitleSize(0.11);
    frameBot->GetYaxis()->SetLabelSize(0.09);
    frameBot->GetYaxis()->SetTitleOffset(0.55);
    frameBot->GetYaxis()->SetNdivisions(505);
    frameBot->Draw();
    TLine unity(etaMin, 1.0, etaMax, 1.0);
    unity.SetLineStyle(kDashed);
    unity.SetLineColor(kGray + 1);
    unity.Draw("SAME");
    ratio->SetMarkerStyle(kFullCircle);
    ratio->SetMarkerColor(kBlack);
    ratio->SetLineColor(kBlack);
    ratio->Draw("PE1 SAME");

    canvas.cd();
    canvas.SaveAs(Form("%s/L3Res_%s_balance_vsabseta_%s%.0fto%.0f.png",
                       etaFolder.c_str(), prepared.token.Data(),
                       axis.tag.Data(), ptLow, ptHigh));

    delete frameTop;
    delete frameBot;
  };

  // Purpose: Data/MC vs fine |eta| for a few merged pT ranges on one axis.
  // The ranges are merged from the stored fine-|eta| TProfile3Ds (sums of
  // weights), not by averaging the per-bin histograms, so the weighted-MC
  // statistics and errors stay correct. Drawn as step histograms.
  auto saveAbsEtaRatioOverlay = [&](const PreparedInput &prepared,
                                    TFile *inputFile, const TString &profileStem,
                                    int nEtaKeep, const string &inputPdfFolder,
                                    const PlotAxis &axis) {
    TProfile3D *mc3 = (TProfile3D *)inputFile->Get(profileStem + "_mc");
    TProfile3D *dt3 = (TProfile3D *)inputFile->Get(profileStem + "_data");
    if (!mc3 || !dt3) {
      cout << "INFO: " << profileStem << "_{mc,data} missing in "
           << prepared.path << "; skipping the merged-pT |eta| ratio plot vs "
           << axis.tag << "." << endl;
      return;
    }

    // Project to (pT, |eta|) at the reference alpha bin only (cumulative
    // alpha filling: bin refAlphaBin already means alpha < alpha_ref).
    auto projectAtRefAlpha = [&](TProfile3D *profile,
                                 const TString &name) -> TProfile2D * {
      profile->GetZaxis()->SetRange(refAlphaBin, refAlphaBin);
      TProfile2D *projected = profile->Project3DProfile("yx"); // x=pT, y=|eta|
      // MC and Data profiles share one name, so detach and rename right away
      // to keep the next projection from replacing this one.
      projected->SetDirectory(nullptr);
      projected->SetName(name);
      profile->GetZaxis()->SetRange(0, 0);
      return projected;
    };
    TProfile2D *mc2 = projectAtRefAlpha(
        mc3, Form("mcAbsEta_%s_%s", prepared.token.Data(), axis.tag.Data()));
    TProfile2D *dt2 = projectAtRefAlpha(
        dt3, Form("dtAbsEta_%s_%s", prepared.token.Data(), axis.tag.Data()));

    // Requested pT ranges; the last one runs to the end of the pT axis.
    const vector<pair<double, double>> ptRanges = {
        {60., 80.}, {80., 120.}, {120., 300.}, {300., 1e9}};
    TAxis *ptAxis = mc2->GetXaxis();
    nEtaKeep = std::min(nEtaKeep, mc2->GetNbinsY());
    vector<double> etaEdges(nEtaKeep + 1);
    for (int bin = 1; bin <= nEtaKeep + 1; ++bin)
      etaEdges[bin - 1] = mc2->GetYaxis()->GetBinLowEdge(bin);

    TCanvas canvas(
        Form("cAbsEtaRatio_%s_%s", prepared.token.Data(), axis.tag.Data()),
        Form("cAbsEtaRatio_%s_%s", prepared.token.Data(), axis.tag.Data()),
        900, 700);

    // Fixed range (as in ratio_raw_vseta), comparable across inputs.
    TH1D *frame = new TH1D(
        Form("hAbsEtaRatio_%s_%s", prepared.token.Data(), axis.tag.Data()),
        Form(";|#eta_{jet}|;%s", responseAxisLabel(prepared).Data()), 100,
        etaEdges.front(), etaEdges.back());
    frame->SetMinimum(0.6);
    frame->SetMaximum(1.4);
    styleFrameAxes(frame);
    frame->Draw();

    TLine unity(etaEdges.front(), 1.0, etaEdges.back(), 1.0);
    unity.SetLineStyle(kDashed);
    unity.SetLineColor(kGray + 1);
    unity.Draw("SAME");

    TLegend legend(0.58, 0.70, 0.90, 0.90);
    legend.SetBorderSize(0);
    legend.SetFillStyle(0);
    legend.SetTextFont(42);
    legend.SetTextSize(0.030);

    vector<TH1D *> drawn;
    int index = 0;
    for (const auto &range : ptRanges) {
      int firstBin = 0, lastBin = 0;
      for (int bin = 1; bin <= ptAxis->GetNbins(); ++bin) {
        if (ptAxis->GetBinLowEdge(bin) >= range.first - 1e-6 &&
            ptAxis->GetBinUpEdge(bin) <= range.second + 1e-6) {
          if (!firstBin)
            firstBin = bin;
          lastBin = bin;
        }
      }
      if (!firstBin)
        continue;
      const double ptLow = ptAxis->GetBinLowEdge(firstBin);
      const double ptHigh = ptAxis->GetBinUpEdge(lastBin);
      const TString id = Form("%s_%s_%.0fto%.0f", prepared.token.Data(),
                              axis.tag.Data(), ptLow, ptHigh);

      TProfile *mcEta = mc2->ProfileY("mcEta_" + id, firstBin, lastBin);
      TProfile *dtEta = dt2->ProfileY("dtEta_" + id, firstBin, lastBin);
      TH1D *ratio = new TH1D("ratioEta_" + id, "", nEtaKeep, etaEdges.data());
      ratio->SetDirectory(nullptr);
      for (int bin = 1; bin <= nEtaKeep; ++bin) {
        const double mcVal = mcEta->GetBinContent(bin);
        const double dtVal = dtEta->GetBinContent(bin);
        if (!(mcVal > 0.0) || !(dtVal > 0.0) || mcEta->GetBinEntries(bin) <= 0 ||
            dtEta->GetBinEntries(bin) <= 0)
          continue;
        const double value = dtVal / mcVal;
        const double relErr = sqrt(pow(dtEta->GetBinError(bin) / dtVal, 2) +
                                   pow(mcEta->GetBinError(bin) / mcVal, 2));
        ratio->SetBinContent(bin, value);
        ratio->SetBinError(bin, value * relErr);
      }
      delete mcEta;
      delete dtEta;

      const int color = colors[index % colorCount];
      ratio->SetLineColor(color);
      ratio->SetLineWidth(2);
      ratio->SetMarkerSize(0);
      ratio->Draw("HIST SAME"); // bin outline
      ratio->Draw("E1 SAME");   // error bars at bin centres
      legend.AddEntry(ratio,
                      range.second > 1e8
                          ? Form("%s #geq %.0f GeV", axis.symbol.Data(), ptLow)
                          : Form("%.0f #leq %s < %.0f GeV", ptLow,
                                 axis.symbol.Data(), ptHigh),
                      "LE");
      drawn.push_back(ratio);
      ++index;
    }

    TLatex label;
    label.SetNDC();
    label.SetTextFont(42);
    label.SetTextSize(0.032);
    label.DrawLatex(0.18, 0.86, prepared.label.Data());
    label.DrawLatex(0.18, 0.81,
                    Form("%s = B_{%s}^{Data} / B_{%s}^{MC}",
                         responseAxisLabel(prepared).Data(),
                         prepared.label.Data(), prepared.label.Data()));
    label.DrawLatex(0.18, 0.76,
                    Form("#alpha_{ref} < %.2f, raw - not #alpha#rightarrow0",
                         prepared.refAlphaValue));
    legend.Draw();
    frame->Draw("AXIS SAME"); // keep tick marks on top of the outlines
    drawCmsLabel(canvas);
    canvas.SaveAs(Form("%s/L3Res_%s_ratio_vsabseta_%s.png",
                       inputPdfFolder.c_str(), prepared.token.Data(),
                       axis.tag.Data()));

    for (TH1D *hist : drawn)
      delete hist;
    delete mc2;
    delete dt2;
    delete frame;
  };

  // Purpose: all fine-|eta| plots for one pT axis, skipped quietly if the
  // derived file predates the balance_vsabseta_* histograms.
  auto saveAbsEtaPlots = [&](const PreparedInput &prepared, TFile *inputFile,
                             const PlotAxis &axis, const TString &profileStem,
                             const string &inputPdfFolder,
                             const string &etaFolder) {
    AbsEtaSeries series =
        loadAbsEtaSeries(inputFile, axis, (profileStem + "_mc").Data());
    if (series.ratio.empty()) {
      cout << "INFO: balance_vsabseta_*_" << axis.tag << "* missing in "
           << prepared.path << "; skipping the |eta| plots vs "
           << axis.tag << "." << endl;
      return;
    }
    gSystem->mkdir(etaFolder.c_str(), kTRUE);
    for (const auto &entry : series.ratio)
      saveAbsEtaBinPlot(prepared, series, entry.first, etaFolder, axis);
    // Same |eta| range as the per-bin histograms (the fit's range).
    saveAbsEtaRatioOverlay(prepared, inputFile, profileStem,
                           series.ratio.begin()->second->GetNbinsX(),
                           inputPdfFolder, axis);
  };

  vector<PreparedInput> preparedInputs;
  preparedInputs.reserve(inputFiles.size());

  for (size_t inputIndex = 0; inputIndex < inputFiles.size(); ++inputIndex) {
    PreparedInput prepared;
    prepared.sample = sampleTypes[inputIndex];
    prepared.path = inputFiles[inputIndex];
    prepared.label = sampleLabel(prepared.sample);
    prepared.token = sampleToken(prepared.sample, (int)inputIndex + 1);
    prepared.jecTag = sampleTag(prepared.sample);

    const string inputFolder =
        outfolder + "/inputs/" + string(prepared.token.Data());
    const string inputPdfFolder = inputFolder + "/pdf";
    const string alphaFolder = inputFolder + "/alpha_extrap";
    gSystem->mkdir(inputFolder.c_str(), kTRUE);
    gSystem->mkdir(inputPdfFolder.c_str(), kTRUE);
    gSystem->mkdir(alphaFolder.c_str(), kTRUE);

    TFile *inputFile = TFile::Open(prepared.path, "READ");
    if (!inputFile || inputFile->IsZombie()) {
      cout << "ERROR: Cannot open derived input '" << prepared.path << "'."
           << endl;
      return;
    }

    TProfile3D *balance3D = (TProfile3D *)inputFile->Get("balance3D_mc");
    if (!balance3D) {
      cout << "ERROR: Missing required object 'balance3D_mc' in "
           << prepared.path << endl;
      inputFile->Close();
      delete inputFile;
      return;
    }

    vector<double> alphaEdges(balance3D->GetZaxis()->GetNbins() + 1, 0.0);
    for (int bin = 1; bin <= balance3D->GetZaxis()->GetNbins() + 1; ++bin) {
      alphaEdges[bin - 1] = balance3D->GetZaxis()->GetBinLowEdge(bin);
    }
    if (refAlphaBin < 1 || refAlphaBin >= (int)alphaEdges.size()) {
      cout << "ERROR: Reference alpha bin " << refAlphaBin
           << " is outside the input alpha range for " << prepared.path << endl;
      inputFile->Close();
      delete inputFile;
      return;
    }

    prepared.refAlphaValue = alphaEdges[refAlphaBin];
    prepared.etaMin = balance3D->GetYaxis()->GetBinLowEdge(1);
    prepared.etaMax = balance3D->GetYaxis()->GetBinLowEdge(
        balance3D->GetYaxis()->GetNbins() + 1);
    prepared.haveEtaRange = (prepared.etaMax > prepared.etaMin);

    // Purpose: pull the raw (ref-alpha, not alpha->0 extrapolated) response
    // vs |eta| for each pTref bin, if deriveL3_from_photonjet.C wrote them.
    map<int, TH1D *> etaRatioByPt;
    TAxis *ptAxisForEta = balance3D->GetXaxis();
    {
      // Matches the "alpha%.2f" tag deriveL3_from_photonjet.C builds from a
      // float alphaCutValue when naming ratio_pt<range>_alpha<cut>.
      const float alphaCutValueForEta = (float)prepared.refAlphaValue;
      const TString alphaTagForEta = Form("alpha%.2f", alphaCutValueForEta);
      for (int ptBin = 1; ptBin <= ptAxisForEta->GetNbins(); ++ptBin) {
        const TString ptTag = Form("%.0fto%.0f",
                                   ptAxisForEta->GetBinLowEdge(ptBin),
                                   ptAxisForEta->GetBinLowEdge(ptBin + 1));
        const TString histName =
            Form("ratio_pt%s_%s", ptTag.Data(), alphaTagForEta.Data());
        TH1D *sourceHist = (TH1D *)inputFile->Get(histName);
        if (!sourceHist)
          continue;
        TH1D *clone = (TH1D *)sourceHist->Clone(Form("%s_clone", histName.Data()));
        clone->SetDirectory(nullptr);
        etaRatioByPt[ptBin] = clone;
      }
      if (etaRatioByPt.empty()) {
        cout << "INFO: ratio_pt*_alpha* (vs eta) histograms missing in "
             << prepared.path
             << "; skipping the raw JES-vs-eta diagnostic plot." << endl;
      }
    }

    map<int, TH1D *> ratioSeries =
        loadSeries(inputFile, "ratio_vsptref_alpha", prepared.path);
    map<int, TH1D *> jetRatioSeries =
        loadSeries(inputFile, "ratio_vsjetpt_alpha", prepared.path);
    if (ratioSeries.empty()) {
      cout << "ERROR: Missing required ratio_vsptref series in "
           << prepared.path << endl;
      inputFile->Close();
      delete inputFile;
      return;
    }
    prepared.haveJetPt =
        !jetRatioSeries.empty() && jetRatioSeries.count(refAlphaBin) > 0;
    if (!prepared.haveJetPt) {
      cout << "INFO: ratio_vsjetpt_alpha* missing or incomplete in "
           << prepared.path
           << "; skipping the vsjetpt diagnostic series for this input."
           << endl;
    }
    map<int, TH1D *> ratioNormSeries =
        loadSeries(inputFile, "ratio_norm_vsptref_alpha", prepared.path, false);
    if (ratioNormSeries.empty()) {
      cout << "INFO: ratio_norm_vsptref_alpha* missing in " << prepared.path
           << "; building conservative normalized series in memory." << endl;
      ratioNormSeries = normalizeSeriesByReference(
          ratioSeries, refAlphaBin, "ratio_norm_vsptref_alpha", prepared.path);
    }
    map<int, TH1D *> jetRatioNormSeries;
    if (prepared.haveJetPt) {
      jetRatioNormSeries = loadSeries(inputFile, "ratio_norm_vsjetpt_alpha",
                                      prepared.path, false);
      if (jetRatioNormSeries.empty()) {
        cout << "INFO: ratio_norm_vsjetpt_alpha* missing in " << prepared.path
             << "; building conservative normalized series in memory." << endl;
        jetRatioNormSeries =
            normalizeSeriesByReference(jetRatioSeries, refAlphaBin,
                                       "ratio_norm_vsjetpt_alpha", prepared.path);
      }
      if (jetRatioNormSeries.empty()) {
        prepared.haveJetPt = false;
      }
    }
    if (ratioNormSeries.empty()) {
      inputFile->Close();
      delete inputFile;
      return;
    }
    if (!ratioSeries.count(refAlphaBin)) {
      cout << "ERROR: Missing the requested reference alpha bin " << refAlphaBin
           << " in " << prepared.path << endl;
      inputFile->Close();
      delete inputFile;
      return;
    }

    prepared.refRatioHist = cloneRequiredHist(
        inputFile, Form("ratio_vsptref_alpha%d", refAlphaBin), prepared.path);
    if (prepared.haveJetPt) {
      prepared.refJetRatioHist = cloneRequiredHist(
          inputFile, Form("ratio_vsjetpt_alpha%d", refAlphaBin), prepared.path);
      if (!prepared.refJetRatioHist) prepared.haveJetPt = false;
    }
    TH1D *rawMc = cloneRequiredHist(
        inputFile, Form("balance_vsptref_mc_alpha%d", refAlphaBin),
        prepared.path);
    TH1D *rawDt = cloneRequiredHist(
        inputFile, Form("balance_vsptref_data_alpha%d", refAlphaBin),
        prepared.path);
    if (!prepared.refRatioHist || !rawMc || !rawDt) {
      inputFile->Close();
      delete inputFile;
      return;
    }

    // Optional: raw MC/Data balance vs leading-jet pT at the reference alpha.
    TH1D *rawJetMc = nullptr;
    TH1D *rawJetDt = nullptr;
    if (prepared.haveJetPt) {
      TH1D *srcMc = (TH1D *)inputFile->Get(
          Form("balance_vsjetpt_mc_alpha%d", refAlphaBin));
      TH1D *srcDt = (TH1D *)inputFile->Get(
          Form("balance_vsjetpt_data_alpha%d", refAlphaBin));
      if (srcMc && srcDt) {
        rawJetMc = (TH1D *)srcMc->Clone(
            Form("balance_vsjetpt_mc_alpha%d_clone", refAlphaBin));
        rawJetMc->SetDirectory(nullptr);
        rawJetDt = (TH1D *)srcDt->Clone(
            Form("balance_vsjetpt_data_alpha%d_clone", refAlphaBin));
        rawJetDt->SetDirectory(nullptr);
      } else {
        cout << "INFO: balance_vsjetpt_{mc,data}_alpha" << refAlphaBin
             << " missing in " << prepared.path
             << "; skipping the raw balance vs jet pT plot." << endl;
      }
    }

    prepared.fitWindow = parseFitWindow(
        inputIndex < fitWindowTokens.size() ? fitWindowTokens[inputIndex] : "",
        prepared.refRatioHist);
    saveRatioOverlay(prepared, ratioSeries, prepared.refRatioHist, alphaEdges,
                     inputPdfFolder, ptRefAxis());
    saveRefRawOverlay(prepared, rawMc, rawDt, inputPdfFolder, ptRefAxis());
    saveEtaRawOverlay(prepared, etaRatioByPt, ptAxisForEta, inputPdfFolder);
    saveAbsEtaPlots(prepared, inputFile, ptRefAxis(), "balance3D_abseta",
                    inputPdfFolder, inputFolder + "/eta");
    saveAbsEtaPlots(prepared, inputFile, jetPtAxis(),
                    "balance3D_abseta_jetpt", inputPdfFolder,
                    inputFolder + "/eta");
    if (prepared.haveJetPt) {
      saveRatioOverlay(prepared, jetRatioSeries, prepared.refJetRatioHist,
                       alphaEdges, inputPdfFolder, jetPtAxis());
      if (rawJetMc && rawJetDt) {
        saveRefRawOverlay(prepared, rawJetMc, rawJetDt, inputPdfFolder,
                          jetPtAxis());
      }
    }

    prepared.kfsrHist = (TH1D *)prepared.refRatioHist->Clone(
        Form("kFSR_vsptref_%s", prepared.token.Data()));
    prepared.kfsrHist->Reset();
    prepared.kfsrHist->SetDirectory(nullptr);
    prepared.corrHist = (TH1D *)prepared.refRatioHist->Clone(
        Form("corr_vsptref_%s", prepared.token.Data()));
    prepared.corrHist->Reset();
    prepared.corrHist->SetDirectory(nullptr);
    if (prepared.haveJetPt) {
      prepared.kfsrJetHist = (TH1D *)prepared.refJetRatioHist->Clone(
          Form("kFSR_vsjetpt_%s", prepared.token.Data()));
      prepared.kfsrJetHist->Reset();
      prepared.kfsrJetHist->SetDirectory(nullptr);
      prepared.corrJetHist = (TH1D *)prepared.refJetRatioHist->Clone(
          Form("corr_vsjetpt_%s", prepared.token.Data()));
      prepared.corrJetHist->Reset();
      prepared.corrJetHist->SetDirectory(nullptr);
    }

    TDirectory *inputDir =
        outputFile->mkdir(Form("input_%s", prepared.token.Data()));
    if (!inputDir) {
      cout << "ERROR: Could not create output directory for input '"
           << prepared.label << "'." << endl;
      inputFile->Close();
      delete inputFile;
      return;
    }

    inputDir->cd();
    prepared.refRatioHist->Write("ratio_vsptref_refalpha");
    if (prepared.haveJetPt) {
      prepared.refJetRatioHist->Write("ratio_vsjetpt_refalpha");
    }
    rawMc->Write("balance_vsptref_mc_refalpha");
    rawDt->Write("balance_vsptref_data_refalpha");
    if (rawJetMc && rawJetDt) {
      rawJetMc->Write("balance_vsjetpt_mc_refalpha");
      rawJetDt->Write("balance_vsjetpt_data_refalpha");
    }
    for (auto &entry : etaRatioByPt) {
      entry.second->Write(Form("ratio_vseta_ptbin%d_refalpha", entry.first));
    }

    auto fillCorrectedResponse =
        [&](TH1D *refHist, const map<int, TH1D *> &normalizedAxisRatioSeries,
            TH1D *kfsrHist, TH1D *corrHist, const TString &graphStem,
            const PlotAxis &axis) {
          for (int ptBin = 1; ptBin <= refHist->GetNbinsX(); ++ptBin) {
            vector<double> alphaAllX;
            vector<double> alphaAllY;
            vector<double> alphaAllErr;
            vector<double> alphaFitX;
            vector<double> alphaFitY;
            vector<double> alphaFitErr;

            const double refValue = refHist->GetBinContent(ptBin);
            const double refError = refHist->GetBinError(ptBin);
            if (!(refValue > 0.0) || !std::isfinite(refValue)) {
              kfsrHist->SetBinContent(ptBin, 0.0);
              kfsrHist->SetBinError(ptBin, 0.0);
              corrHist->SetBinContent(ptBin, 0.0);
              corrHist->SetBinError(ptBin, 0.0);
              continue;
            }

            for (const auto &entry : normalizedAxisRatioSeries) {
              const int alphaBin = entry.first;
              TH1D *hist = entry.second;
              const double value = hist->GetBinContent(ptBin);
              const double error = hist->GetBinError(ptBin);
              if (!(value > 0.0) || !std::isfinite(value))
                continue;

              const double alphaCenter =
                  0.5 * (alphaEdges[alphaBin - 1] + alphaEdges[alphaBin]);
              alphaAllX.push_back(alphaCenter);
              alphaAllY.push_back(value);
              alphaAllErr.push_back(
                  (error > 0.0 && std::isfinite(error)) ? error : 0.0);

              if (alphaBin == refAlphaBin)
                continue;
              if (alphaCenter < fitAlphaMin || alphaCenter > fitAlphaMax)
                continue;
              if (!(error > 0.0) || !std::isfinite(error))
                continue;
              alphaFitX.push_back(alphaCenter);
              alphaFitY.push_back(value);
              alphaFitErr.push_back(error);
            }

            if (alphaAllX.size() < 2 || alphaFitX.size() < 2) {
              kfsrHist->SetBinContent(ptBin, 0.0);
              kfsrHist->SetBinError(ptBin, 0.0);
              corrHist->SetBinContent(ptBin, 0.0);
              corrHist->SetBinError(ptBin, 0.0);
              continue;
            }

            TGraphErrors *alphaGraph = new TGraphErrors((int)alphaAllX.size());
            for (int point = 0; point < (int)alphaAllX.size(); ++point) {
              alphaGraph->SetPoint(point, alphaAllX[point], alphaAllY[point]);
              alphaGraph->SetPointError(point, 0.0, alphaAllErr[point]);
            }
            alphaGraph->SetName(Form("%s_pt%d", graphStem.Data(), ptBin));

            TGraphErrors fitGraph((int)alphaFitX.size());
            for (int point = 0; point < (int)alphaFitX.size(); ++point) {
              fitGraph.SetPoint(point, alphaFitX[point], alphaFitY[point]);
              fitGraph.SetPointError(point, 0.0, alphaFitErr[point]);
            }
            TF1 alphaLine(Form("fAlpha_%s_pt%d", graphStem.Data(), ptBin),
                          "[0]+[1]*x", fitAlphaMin, fitAlphaMax);
            alphaLine.SetParameters(1.0, 0.0);
            fitGraph.Fit(&alphaLine, "QRN");

            const double kfsr = alphaLine.GetParameter(0);
            const double kfsrError = alphaLine.GetParError(0);
            kfsrHist->SetBinContent(ptBin, kfsr);
            kfsrHist->SetBinError(ptBin, kfsrError);

            if (kfsr > 0.0) {
              corrHist->SetBinContent(ptBin, refValue * kfsr);
              // Keep the reference-ratio uncertainty on the corrected response.
              // The fitted kFSR uncertainty is preserved separately in
              // kfsrHist; multiplying it back in here would double-count the
              // same reference denominator that already appears in the
              // normalized alpha series.
              corrHist->SetBinError(ptBin, refError * kfsr);
            } else {
              corrHist->SetBinContent(ptBin, 0.0);
              corrHist->SetBinError(ptBin, 0.0);
            }

            TF1 *alphaLineOut = (TF1 *)alphaLine.Clone(
                Form("fAlpha_%s_pt%d", graphStem.Data(), ptBin));
            alphaGraph->Write();
            alphaLineOut->Write();

            saveAlphaFitPlot(
                prepared, ptBin, refHist->GetXaxis()->GetBinLowEdge(ptBin),
                refHist->GetXaxis()->GetBinLowEdge(ptBin + 1), alphaGraph,
                alphaLineOut, kfsr, kfsrError, alphaFolder, axis);
          }
        };

    fillCorrectedResponse(prepared.refRatioHist, ratioNormSeries,
                          prepared.kfsrHist, prepared.corrHist,
                          "gAlphaNormPtRef", ptRefAxis());
    if (prepared.haveJetPt) {
      fillCorrectedResponse(prepared.refJetRatioHist, jetRatioNormSeries,
                            prepared.kfsrJetHist, prepared.corrJetHist,
                            "gAlphaNormJetPt", jetPtAxis());
    }

    prepared.kfsrHist->Write("kFSR_vsptref");
    prepared.corrHist->Write("corr_vsptref");
    if (prepared.haveJetPt) {
      prepared.kfsrJetHist->Write("kFSR_vsjetpt");
      prepared.corrJetHist->Write("corr_vsjetpt");
    }
    saveKfsrSummary(prepared, prepared.kfsrHist, alphaFolder, ptRefAxis());
    saveCorrSummary(prepared, prepared.refRatioHist, prepared.corrHist,
                    alphaFolder, ptRefAxis());
    if (prepared.haveJetPt) {
      saveKfsrSummary(prepared, prepared.kfsrJetHist, alphaFolder,
                      jetPtAxis());
      saveCorrSummary(prepared, prepared.refJetRatioHist, prepared.corrJetHist,
                      alphaFolder, jetPtAxis());
      saveCorrPtRefVsJetPt(prepared, alphaFolder);
    }

    prepared.displayGraph = buildWindowedGraph(
        prepared.corrHist, prepared.fitWindow, true,
        markers[inputIndex % markerCount], colors[inputIndex % colorCount]);
    prepared.fitGraph = buildWindowedGraph(
        prepared.corrHist, prepared.fitWindow, false,
        markers[inputIndex % markerCount], colors[inputIndex % colorCount]);
    if (!prepared.fitGraph || prepared.fitGraph->GetN() < 3) {
      cout << "ERROR: Not enough valid pTref points in the selected fit window "
              "for input '"
           << prepared.label << "'." << endl;
      inputFile->Close();
      delete inputFile;
      return;
    }

    prepared.displayGraph->Write("g_ptref_display");
    prepared.fitGraph->Write("g_ptref_fit");
    if (prepared.haveJetPt) {
      FitWindow jetWindow = parseFitWindow("", prepared.refJetRatioHist);
      prepared.displayJetGraph = buildWindowedGraph(
          prepared.corrJetHist, jetWindow, true,
          markers[inputIndex % markerCount], colors[inputIndex % colorCount]);
      prepared.fitJetGraph = buildWindowedGraph(
          prepared.corrJetHist, jetWindow, true,
          markers[inputIndex % markerCount], colors[inputIndex % colorCount]);
      if (!prepared.fitJetGraph || prepared.fitJetGraph->GetN() < 3) {
        cout << "INFO: Not enough valid JetPt points for input '"
             << prepared.label
             << "'; dropping its vsjetpt diagnostic series." << endl;
        prepared.haveJetPt = false;
        prepared.displayJetGraph = nullptr;
        prepared.fitJetGraph = nullptr;
      } else {
        prepared.displayJetGraph->Write("g_jetpt_display");
        prepared.fitJetGraph->Write("g_jetpt_fit");
      }
    }
    outputFile->cd();
    preparedInputs.push_back(prepared);

    for (auto &entry : etaRatioByPt)
      delete entry.second;
    for (auto &entry : ratioSeries)
      delete entry.second;
    for (auto &entry : jetRatioSeries)
      delete entry.second;
    for (auto &entry : ratioNormSeries)
      delete entry.second;
    for (auto &entry : jetRatioNormSeries)
      delete entry.second;
    delete rawMc;
    delete rawDt;
    inputFile->Close();
    delete inputFile;
  }

  if (preparedInputs.empty()) {
    cout << "ERROR: No valid derived inputs were prepared for the combined L3 "
            "fit."
         << endl;
    outputFile->Close();
    delete outputFile;
    return;
  }

  double fitMinOut = preparedInputs.front().fitWindow.min;
  double fitMaxOut = preparedInputs.front().fitWindow.max;
  double xMin =
      preparedInputs.front().refRatioHist->GetXaxis()->GetBinLowEdge(1);
  double xMax = preparedInputs.front().refRatioHist->GetXaxis()->GetBinLowEdge(
      preparedInputs.front().refRatioHist->GetNbinsX() + 1);
  double jetXMin = 1e18;
  double jetXMax = -1e18;
  double etaMinOut = preparedInputs.front().etaMin;
  double etaMaxOut = preparedInputs.front().etaMax;
  bool haveEtaRange = preparedInputs.front().haveEtaRange;

  TGraphErrors *combinedGraph = new TGraphErrors();
  int combinedPoint = 0;
  double yMin = 1e9;
  double yMax = -1e9;
  TGraphErrors *combinedJetGraph = new TGraphErrors();
  int combinedJetPoint = 0;
  bool anyJetPt = false;

  for (size_t inputIndex = 0; inputIndex < preparedInputs.size();
       ++inputIndex) {
    fitMinOut = std::min(fitMinOut, preparedInputs[inputIndex].fitWindow.min);
    fitMaxOut = std::max(fitMaxOut, preparedInputs[inputIndex].fitWindow.max);
    xMin = std::min(
        xMin,
        preparedInputs[inputIndex].refRatioHist->GetXaxis()->GetBinLowEdge(1));
    xMax = std::max(
        xMax,
        preparedInputs[inputIndex].refRatioHist->GetXaxis()->GetBinLowEdge(
            preparedInputs[inputIndex].refRatioHist->GetNbinsX() + 1));
    if (preparedInputs[inputIndex].haveJetPt) {
      anyJetPt = true;
      jetXMin = std::min(
          jetXMin,
          preparedInputs[inputIndex].refJetRatioHist->GetXaxis()->GetBinLowEdge(
              1));
      jetXMax = std::max(
          jetXMax,
          preparedInputs[inputIndex].refJetRatioHist->GetXaxis()->GetBinLowEdge(
              preparedInputs[inputIndex].refJetRatioHist->GetNbinsX() + 1));
    }

    if (preparedInputs[inputIndex].haveEtaRange) {
      etaMinOut = std::min(etaMinOut, preparedInputs[inputIndex].etaMin);
      etaMaxOut = std::max(etaMaxOut, preparedInputs[inputIndex].etaMax);
      haveEtaRange = true;
    }

    for (int point = 0; point < preparedInputs[inputIndex].fitGraph->GetN();
         ++point) {
      combinedGraph->SetPoint(
          combinedPoint, preparedInputs[inputIndex].fitGraph->GetX()[point],
          preparedInputs[inputIndex].fitGraph->GetY()[point]);
      combinedGraph->SetPointError(
          combinedPoint, preparedInputs[inputIndex].fitGraph->GetEX()[point],
          preparedInputs[inputIndex].fitGraph->GetEY()[point]);
      yMin = std::min(yMin, preparedInputs[inputIndex].fitGraph->GetY()[point]);
      yMax = std::max(yMax, preparedInputs[inputIndex].fitGraph->GetY()[point]);
      ++combinedPoint;
    }

    for (int point = 0; point < preparedInputs[inputIndex].displayGraph->GetN();
         ++point) {
      yMin = std::min(yMin,
                      preparedInputs[inputIndex].displayGraph->GetY()[point]);
      yMax = std::max(yMax,
                      preparedInputs[inputIndex].displayGraph->GetY()[point]);
    }

    if (preparedInputs[inputIndex].haveJetPt) {
      for (int point = 0;
           point < preparedInputs[inputIndex].fitJetGraph->GetN(); ++point) {
        combinedJetGraph->SetPoint(
            combinedJetPoint,
            preparedInputs[inputIndex].fitJetGraph->GetX()[point],
            preparedInputs[inputIndex].fitJetGraph->GetY()[point]);
        combinedJetGraph->SetPointError(
            combinedJetPoint,
            preparedInputs[inputIndex].fitJetGraph->GetEX()[point],
            preparedInputs[inputIndex].fitJetGraph->GetEY()[point]);
        ++combinedJetPoint;
      }
    }
  }

  if (combinedGraph->GetN() < 3) {
    cout << "ERROR: The combined fit graph has fewer than three points."
         << endl;
    outputFile->Close();
    delete outputFile;
    delete combinedGraph;
    delete combinedJetGraph;
    return;
  }

  if (anyJetPt && combinedJetGraph->GetN() < 3) {
    cout << "WARNING: The combined JetPt graph has fewer than three points; "
            "dropping the vsjetpt diagnostic output."
         << endl;
    anyJetPt = false;
  }

  const double yLowerMargin = 0.08;
  const double yUpperMargin = 0.30;
  yMin = std::max(0.45, yMin - yLowerMargin);
  yMax = std::min(1.8, yMax + yUpperMargin);

  SequentialFitResult fitResult = runSequentialTemplateFit(
      combinedGraph, fitMinOut, fitMaxOut, outfilename.c_str());
  if (!fitResult.valid || !fitResult.finalFit) {
    cout << "ERROR: The sequential pTref fit failed." << endl;
    outputFile->Close();
    delete outputFile;
    delete combinedGraph;
    return;
  }

  TCanvas finalCanvas("cL3ResFinal", "cL3ResFinal", 900, 700);
  finalCanvas.SetLogx();
  TH1D *frame = new TH1D("hL3ResFinalFrame",
                         ";p_{T}^{ref} (GeV);R_{Data/MC} = B^{Data}/B^{MC}",
                         100, xMin, xMax);
  frame->SetMinimum(yMin);
  frame->SetMaximum(yMax);
  frame->GetXaxis()->SetMoreLogLabels(kTRUE);
  frame->GetXaxis()->SetNoExponent(kTRUE);
  styleFrameAxes(frame);
  frame->Draw();

  TLine unity(xMin, 1.0, xMax, 1.0);
  unity.SetLineStyle(kDashed);
  unity.SetLineColor(kGray + 1);
  unity.Draw("SAME");

  TLegend legend(0.52, 0.68, 0.88, 0.88);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  legend.SetTextFont(42);
  legend.SetTextSize(0.028);

  for (size_t inputIndex = 0; inputIndex < preparedInputs.size();
       ++inputIndex) {
    preparedInputs[inputIndex].displayGraph->Draw("PE SAME");
    legend.AddEntry(preparedInputs[inputIndex].displayGraph,
                    Form("%s (%.0f-%.0f GeV)",
                         preparedInputs[inputIndex].label.Data(),
                         preparedInputs[inputIndex].fitWindow.min,
                         preparedInputs[inputIndex].fitWindow.max),
                    "PE");
  }

  fitResult.finalFit->SetLineColor(kBlack);
  fitResult.finalFit->SetLineWidth(2);
  fitResult.finalFit->Draw("SAME");

  TLatex label;
  label.SetNDC();
  label.SetTextFont(42);
  label.SetTextSize(0.032);
  label.DrawLatex(0.18, 0.86, Form("Run label: %s", runLabel.c_str()));
  if (haveEtaRange) {
    label.DrawLatex(0.18, 0.81,
                    Form("%.3f < #eta < %.3f", etaMinOut, etaMaxOut));
  }
  label.DrawLatex(0.18, 0.76,
                  Form("#alpha_{ref} < %.2f, fit %.2f #leq #alpha #leq %.2f",
                       preparedInputs.front().refAlphaValue, fitAlphaMin,
                       fitAlphaMax));
  label.DrawLatex(0.18, 0.71,
                  Form("Fit window: %.0f-%.0f GeV", fitMinOut, fitMaxOut));
  label.DrawLatex(0.18, 0.66, "Direct global p_{T}^{ref} fit to R_{Data/MC}");
  label.DrawLatex(0.18, 0.61,
                  "C_{L3}(x)=1./(p_{0}+p_{1}/x+p_{2}log(x)/"
                  "x+p_{3}T(x)+p_{6}x^{-0.3051}+p_{7}x)");
  label.DrawLatex(0.18, 0.56,
                  "T(x)=((x/p_{4})^{p_{5}}-1)/((x/p_{4})^{p_{5}}+1)");
  label.DrawLatex(0.18, 0.51,
                  Form("p_{0}=%.4f, p_{1}=%.4f, p_{2}=%.4f, p_{3}=%.4f",
                       fitResult.finalFit->GetParameter(0),
                       fitResult.finalFit->GetParameter(1),
                       fitResult.finalFit->GetParameter(2),
                       fitResult.finalFit->GetParameter(3)));
  label.DrawLatex(0.18, 0.46,
                  Form("p_{4}=%.4f, p_{5}=%.4f, p_{6}=%.4f, p_{7}=%.6f",
                       fitResult.finalFit->GetParameter(4),
                       fitResult.finalFit->GetParameter(5),
                       fitResult.finalFit->GetParameter(6),
                       fitResult.finalFit->GetParameter(7)));
  label.DrawLatex(0.18, 0.41,
                  Form("#chi^{2}/ndf = %.1f/%d",
                       fitResult.finalFit->GetChisquare(),
                       fitResult.finalFit->GetNDF()));
  legend.Draw();
  drawCmsLabel(finalCanvas);
  finalCanvas.SaveAs(
      Form("%s/L3Res_%s_ptref_final.png", pdfFolder.c_str(), runLabel.c_str()));

  const string outputTag = (preparedInputs.size() > 1
                                ? string("combined")
                                : string(preparedInputs.front().jecTag.Data()));
  outputFile->cd();
  combinedGraph->SetName("g_ptref_combined");
  combinedGraph->Write();
  if (anyJetPt) {
    combinedJetGraph->SetName("g_jetpt_combined");
    combinedJetGraph->Write();
  }
  fitResult.f0->Write("f_ptref_f0");
  fitResult.f1->Write("f_ptref_f1");
  fitResult.f2->Write("f_ptref_f2");
  fitResult.f3->Write("f_ptref_f3");
  fitResult.f4->Write("f_ptref_f4");
  fitResult.finalFit->Write("f_ptref_final");

  TNamed runLabelObj("run_label", runLabel.c_str());
  runLabelObj.Write();
  TNamed outputNameObj("output_name", outfilename.c_str());
  outputNameObj.Write();
  TNamed outputTagObj("output_tag", outputTag.c_str());
  outputTagObj.Write();
  TNamed lumiLabelObj("lumi_label", lumiLabel.c_str());
  lumiLabelObj.Write();
  TParameter<int> refAlphaBinObj("ref_alpha_bin", refAlphaBin);
  refAlphaBinObj.Write();
  TParameter<double> refAlphaValueObj("ref_alpha_value",
                                      preparedInputs.front().refAlphaValue);
  refAlphaValueObj.Write();
  TParameter<double> fitMinObj("pt_fit_min", fitMinOut);
  fitMinObj.Write();
  TParameter<double> fitMaxObj("pt_fit_max", fitMaxOut);
  fitMaxObj.Write();
  TParameter<double> etaMinObj("eta_min", haveEtaRange ? etaMinOut : -5.191);
  etaMinObj.Write();
  TParameter<double> etaMaxObj("eta_max", haveEtaRange ? etaMaxOut : 5.191);
  etaMaxObj.Write();
  TParameter<int> inputCountObj("n_inputs", (int)preparedInputs.size());
  inputCountObj.Write();

  outputFile->Close();

  cout << "============================================" << endl;
  cout << "Wrote fit ROOT file: " << outfolder << "/" << outfilename
       << "_fit.root" << endl;
  cout << "Plots saved in: " << pdfFolder << endl;
  cout << "NOTE: L3Res.C stops at the p_{T}^{ref} fit stage." << endl;
  cout << "      Run createL2L3ResTextFile.C on the fit ROOT file, or use "
          "dofits_L3.C, to produce direct p_{T}^{ref} text payloads."
       << endl;
  cout << "============================================" << endl;
}
