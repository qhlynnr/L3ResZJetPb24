// draw_histograms.cpp
//
// Draws every histogram family found in the ROOT file's key listing, with
// axis titles, per-plot titles, and legends. Every histogram (or, where an
// MC/Data pair is naturally overlaid, every pair) gets its own standalone
// canvas and PNG file -- no divided/multi-pad canvases.
//
// The keys fall into 6 naming patterns, each handled by its own function:
//
//   1. mc_balance_pt{X}to{Y}_alpha0.30 / dt_balance_... / l3res_... / ratio_...
//        -> DrawPtBinBalances()
//   2. Balance_vsa_mc_{i}_{j} / Balance_vsa_data_{i}_{j} /
//      L3Res_vsa_{i}_{j} / L3Res_vsa_norm_{i}_{j}   (i = pt bin 1-9, j = eta bin 1-18)
//        -> DrawVsAlpha()
//   3. balance_vsphotonpt_mc_alpha{n} / balance_vsphotonpt_data_alpha{n} /
//      ratio_vsphotonpt_alpha{n} / ratio_norm_vsphotonpt_alpha{n} / ratio_ref_vsphotonpt
//        -> DrawVsPhotonPt()
//   4. balance_vsjetpt_mc_alpha{n} / balance_vsjetpt_data_alpha{n} / ratio_vsjetpt_alpha{n}
//        -> DrawVsJetPt()
//   5. balanceMap_mc / balanceMap_data / l3resMap   (TH2D)
//        -> Draw2DMaps()
//   6. balance3D_mc / balance3D_data (TProfile3D), counts3D_mc / counts3D_data (TH3D)
//        -> Draw3DObjects()
//
// NOTE ON AXIS LABELS: the ROOT file itself carries no axis-title metadata
// for most of these keys, only object titles (mostly generic, e.g. "MC
// balance" for nearly everything). The axis titles set below are the
// physically sensible choice given each family's name (e.g. "vsphotonpt"
// -> x = photon pT, "vsa" -> x = alpha), but adjust kXTitle*/kYTitle*
// below if your analysis defines these differently.
//
// Usage (from a ROOT prompt or command line):
//   root -l -q 'draw_histograms.cpp("input.root","plots")'
// or compiled:
//   root -l
//   .L draw_histograms.cpp+
//   DrawAll("input.root","plots")

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TProfile3D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TString.h>
#include <vector>
#include <iostream>

// ---------------------------------------------------------------------------
// Shared bin/label definitions inferred from the key listing
// ---------------------------------------------------------------------------
TString filetype = ".png";
// The 9 pT bins used for the "pt{X}to{Y}" and vs-alpha (index i) histograms
static const std::vector<TString> kPtBinLabels = {
    "pt60to70", "pt70to80", "pt80to90", "pt90to100", "pt100to120",
    "pt120to140", "pt140to160", "pt160to200", "pt200to300"
};
// Human-readable version of the same bins, used in titles/legends
static const std::vector<TString> kPtBinPretty = {
    "60-70 GeV", "70-80 GeV", "80-90 GeV", "90-100 GeV", "100-120 GeV",
    "120-140 GeV", "140-160 GeV", "160-200 GeV", "200-300 GeV"
};

// The 9 alpha-cut values used for the vsphotonpt / vsjetpt families
static const std::vector<TString> kAlphaLabels = {
    "0.10", "0.15", "0.20", "0.25", "0.30", "0.35", "0.40", "0.45", "0.50"
};

// Number of eta bins used in the vs-alpha (index j) family
static const int kNEtaBins = 18;

// Common axis titles
static const char* kXTitlePhotonPt = "Photon p_{T} [GeV]";
static const char* kXTitleJetPt    = "Jet p_{T} [GeV]";
static const char* kXTitleAlpha    = "#alpha = p_{T,jet2}/p_{T,photon}";
static const char* kXTitleBalance  = "Balance = p_{T,jet}/p_{T,photon}";
static const char* kYTitleBalance  = "Balance";
static const char* kYTitleRatio    = "Ratio (MC/Data)";
static const char* kYTitleRatioNorm = "Normalized ratio";
static const char* kYTitleEntries  = "Events / bin";
static const char* kYTitleL3Res    = "L3 Residual";
static const char* kXTitlePtAvg    = "p_{T,avg} [GeV]";
static const char* kYTitleAbsEta   = "|#eta_{jet}|";

// Helper: fetch a histogram/object by name, warn if missing
template <typename T>
T* GetObj(TFile* f, const TString& name) {
    T* obj = dynamic_cast<T*>(f->Get(name));
    if (!obj) {
        std::cerr << "  [warn] object not found: " << name << std::endl;
    }
    return obj;
}

// Helper: make sure the output directory exists
void EnsureDir(const TString& dir) {
    if (gSystem->AccessPathName(dir)) {
        gSystem->mkdir(dir, kTRUE);
    }
}

// Helper: apply a consistent title/axis-label/color style to a TH1
void StyleHist(TH1* h, const TString& title, const TString& xTitle,
               const TString& yTitle, Color_t color, int lineWidth = 2) {
    if (!h) return;
    h->SetTitle(title);
    h->GetXaxis()->SetTitle(xTitle);
    h->GetYaxis()->SetTitle(yTitle);
    h->GetXaxis()->SetTitleSize(0.045);
    h->GetYaxis()->SetTitleSize(0.045);
    h->SetLineColor(color);
    h->SetLineWidth(lineWidth);
}

// Helper: draw a small legend with a single "name/color" entry.
// Useful even for single-curve pads so a standalone PDF is self-explanatory.
TLegend* MakeLegend(double x1 = 0.62, double y1 = 0.75, double x2 = 0.88, double y2 = 0.88) {
    TLegend* leg = new TLegend(x1, y1, x2, y2);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.04);
    return leg;
}

// Helper: draw one histogram on its own canvas and save it, standalone.
void DrawSingle(TH1* h, const TString& canvasName, const TString& outPath,
                 const TString& title, const TString& xTitle, const TString& yTitle,
                 Color_t color, const TString& legendLabel) {
    if (!h) return;
    TCanvas c(canvasName, canvasName, 800, 650);
    StyleHist(h, title, xTitle, yTitle, color);
    h->Draw("hist");
    TLegend* leg = MakeLegend();
    leg->AddEntry(h, legendLabel, "l");
    leg->Draw();
    c.SaveAs(outPath);
}

// ---------------------------------------------------------------------------
// 1) mc_balance_pt{X}to{Y}_alpha0.30 / dt_balance / l3res / ratio
// ---------------------------------------------------------------------------
void DrawPtBinBalances(TFile* f, const TString& outDir) {
    std::cout << "== DrawPtBinBalances ==" << std::endl;
    TString outDirReal = outDir + "/PtBinBalances/";
    gSystem->mkdir(outDirReal,true);
    for (size_t b = 0; b < kPtBinLabels.size(); ++b) {
        const TString& ptBin    = kPtBinLabels[b];
        const TString& ptPretty = kPtBinPretty[b];
        TString suffix = ptBin + "_alpha0.30";

        TH1D* hMc    = GetObj<TH1D>(f, "mc_balance_" + suffix);
        TH1D* hDt    = GetObj<TH1D>(f, "dt_balance_" + suffix);
        TH1D* hL3res = GetObj<TH1D>(f, "l3res_" + suffix);
        TH1D* hRatio = GetObj<TH1D>(f, "ratio_" + suffix);


        DrawSingle(hMc, "c_mc_" + suffix, outDirReal + "/ptbin_mc_" + suffix + filetype,
                   "MC balance, p_{T} " + ptPretty, kXTitleBalance, kYTitleEntries, kBlue + 1, "MC");

        DrawSingle(hDt, "c_dt_" + suffix, outDirReal + "/ptbin_data_" + suffix + filetype,
                   "Data balance, p_{T} " + ptPretty, kXTitleBalance, kYTitleEntries, kRed + 1, "Data");

        DrawSingle(hL3res, "c_l3res_" + suffix, outDirReal + "/ptbin_l3res_" + suffix + filetype,
                   "L3 residual, p_{T} " + ptPretty, kXTitleBalance, kYTitleL3Res, kGreen + 2,
                   "L3 residual (MC/Data)");

        DrawSingle(hRatio, "c_ratio_" + suffix, outDirReal + "/ptbin_ratio_" + suffix + filetype,
                   "Ratio, p_{T} " + ptPretty, kXTitleBalance, kYTitleRatio, kMagenta + 1, "MC / Data");
    }
}

// ---------------------------------------------------------------------------
// 2) Balance_vsa_mc_{i}_{j} / Balance_vsa_data_{i}_{j} /
//    L3Res_vsa_{i}_{j} / L3Res_vsa_norm_{i}_{j}
// ---------------------------------------------------------------------------
void DrawVsAlpha(TFile* f, const TString& outDir) {
    std::cout << "== DrawVsAlpha ==" << std::endl;
    TString outDirReal = outDir + "/vsAlpha/";
    gSystem->mkdir(outDirReal,true);

    for (size_t i = 1; i <= kPtBinLabels.size(); ++i) {
        const TString& ptPretty = kPtBinPretty[i - 1];
        for (int j = 1; j <= kNEtaBins; ++j) {
            TString tag = TString::Format("%zu_%d", i, j);
            TString binLabel = TString::Format("p_{T} %s, #eta bin %d", ptPretty.Data(), j);
            TString fileTag = "pt" + TString::Format("%zu", i) + "_eta" + TString::Format("%d", j);

            TH1D* hMc      = GetObj<TH1D>(f, "Balance_vsa_mc_" + tag);
            TH1D* hData    = GetObj<TH1D>(f, "Balance_vsa_data_" + tag);
            TH1D* hL3res   = GetObj<TH1D>(f, "L3Res_vsa_" + tag);
            TH1D* hL3norm  = GetObj<TH1D>(f, "L3Res_vsa_norm_" + tag);

            if (!hMc && !hData && !hL3res && !hL3norm) continue; // skip empty combos

            DrawSingle(hMc, "c_vsa_mc_" + tag, outDirReal + "/vsa_mc_" + fileTag + filetype,
                       "MC balance vs #alpha, " + binLabel, kXTitleAlpha, kYTitleBalance, kBlue + 1, "MC");

            DrawSingle(hData, "c_vsa_data_" + tag, outDirReal + "/vsa_data_" + fileTag + filetype,
                       "Data balance vs #alpha, " + binLabel, kXTitleAlpha, kYTitleBalance, kRed + 1, "Data");

            DrawSingle(hL3res, "c_vsa_l3res_" + tag, outDirReal + "/vsa_l3res_" + fileTag + filetype,
                       "L3 residual vs #alpha, " + binLabel, kXTitleAlpha, kYTitleL3Res, kGreen + 2,
                       "L3 residual (MC/Data)");

            DrawSingle(hL3norm, "c_vsa_l3norm_" + tag, outDirReal + "/vsa_l3norm_" + fileTag + filetype,
                       "Normalized L3 residual vs #alpha, " + binLabel, kXTitleAlpha, kYTitleRatioNorm,
                       kOrange + 7, "L3 residual, normalized");
        }
    }
}

// ---------------------------------------------------------------------------
// 3) balance_vsphotonpt_mc_alpha{n} / _data_alpha{n} / ratio_vsphotonpt_alpha{n}
//    ratio_norm_vsphotonpt_alpha{n} / ratio_ref_vsphotonpt
// ---------------------------------------------------------------------------
void DrawVsPhotonPt(TFile* f, const TString& outDir) {
    std::cout << "== DrawVsPhotonPt ==" << std::endl;
    TString outDirReal = outDir + "/vsPhotonPt/";
    gSystem->mkdir(outDirReal,true);

    for (size_t n = 1; n <= kAlphaLabels.size(); ++n) {
        TString idx = TString::Format("alpha%zu", n);
        const TString& aLabel = kAlphaLabels[n - 1];

        TH1D* hMc    = GetObj<TH1D>(f, "balance_vsphotonpt_mc_" + idx);
        TH1D* hData  = GetObj<TH1D>(f, "balance_vsphotonpt_data_" + idx);
        TH1D* hRatio = GetObj<TH1D>(f, "ratio_vsphotonpt_" + idx);
        TH1D* hNorm  = GetObj<TH1D>(f, "ratio_norm_vsphotonpt_" + idx);

        if (hMc && hData) {
            double ymax = std::max(hMc->GetMaximum(), hData->GetMaximum());
            double ymin = std::min(hMc->GetMinimum(), hData->GetMinimum());
            double pad = 0.1 * (ymax - ymin);
            hMc->GetYaxis()->SetRangeUser(ymin - pad, ymax + pad);
        }
        // MC and Data are naturally overlaid on one plot (same quantity, same axes)
        if (hMc || hData) {
            TCanvas c("c_photonpt_" + idx, "photonpt_" + idx, 800, 650);
            TLegend* leg = MakeLegend();
            if (hMc) {
                StyleHist(hMc, "Balance vs photon p_{T} (#alpha < " + aLabel + ")", kXTitlePhotonPt, kYTitleBalance, kBlue + 1);
                hMc->Draw("hist");
                leg->AddEntry(hMc, "MC", "l");
            }
            if (hData) {
                StyleHist(hData, hMc ? TString("") : ("Data balance vs photon p_{T} (#alpha < " + aLabel + ")"),
                           kXTitlePhotonPt, kYTitleBalance, kRed + 1);
                hData->Draw(hMc ? "hist same" : "hist");
                leg->AddEntry(hData, "Data", "l");
            }
            leg->Draw();
            c.SaveAs(outDirReal + "/vsphotonpt_" + idx + filetype);
        }

        DrawSingle(hRatio, "c_photonpt_ratio_" + idx, outDirReal + "/vsphotonpt_ratio_" + idx + filetype,
                   "Balance ratio vs photon p_{T} (#alpha < " + aLabel + ")", kXTitlePhotonPt, kYTitleRatio,
                   kGreen + 2, "MC / Data");

        DrawSingle(hNorm, "c_photonpt_norm_" + idx, outDirReal + "/vsphotonpt_norm_" + idx + filetype,
                   "Normalized ratio vs photon p_{T} (#alpha < " + aLabel + " / ref 0.30)", kXTitlePhotonPt,
                   kYTitleRatioNorm, kMagenta + 1, "Ratio / ref (#alpha<0.30)");
    }

    // Reference ratio (alpha < 0.30), drawn on its own
    TH1D* hRef = GetObj<TH1D>(f, "ratio_ref_vsphotonpt");
    DrawSingle(hRef, "c_photonpt_ref", outDir + "/vsphotonpt_ratio_ref" + filetype,
               "Reference balance ratio vs photon p_{T} (#alpha < 0.30)", kXTitlePhotonPt, kYTitleRatio,
               kBlack, "MC / Data, #alpha<0.30");
}

// ---------------------------------------------------------------------------
// 4) balance_vsjetpt_mc_alpha{n} / _data_alpha{n} / ratio_vsjetpt_alpha{n}
// ---------------------------------------------------------------------------
void DrawVsJetPt(TFile* f, const TString& outDir) {
    TString outDirReal = outDir + "/vsJetPt/";
    gSystem->mkdir(outDirReal,true);

    std::cout << "== DrawVsJetPt ==" << std::endl;
    for (size_t n = 1; n <= kAlphaLabels.size(); ++n) {
        TString idx = TString::Format("alpha%zu", n);
        const TString& aLabel = kAlphaLabels[n - 1];

        TH1D* hMc    = GetObj<TH1D>(f, "balance_vsjetpt_mc_" + idx);
        TH1D* hData  = GetObj<TH1D>(f, "balance_vsjetpt_data_" + idx);
        TH1D* hRatio = GetObj<TH1D>(f, "ratio_vsjetpt_" + idx);

        // MC and Data are naturally overlaid on one plot (same quantity, same axes)
        if (hMc || hData) {
            TCanvas c("c_jetpt_" + idx, "jetpt_" + idx, 800, 650);
            TLegend* leg = MakeLegend();
            if (hMc) {
                StyleHist(hMc, "Balance vs jet p_{T} (#alpha < " + aLabel + ")", kXTitleJetPt, kYTitleBalance, kBlue + 1);
                hMc->Draw("hist");
                leg->AddEntry(hMc, "MC", "l");
            }
            if (hData) {
                StyleHist(hData, hMc ? TString("") : ("Data balance vs jet p_{T} (#alpha < " + aLabel + ")"),
                           kXTitleJetPt, kYTitleBalance, kRed + 1);
                hData->Draw(hMc ? "hist same" : "hist");
                leg->AddEntry(hData, "Data", "l");
            }
            leg->Draw();
            c.SaveAs(outDirReal + "/vsjetpt_" + idx + filetype);
        }

        DrawSingle(hRatio, "c_jetpt_ratio_" + idx, outDirReal + "/vsjetpt_ratio_" + idx + filetype,
                   "Balance ratio vs jet p_{T} (#alpha < " + aLabel + ")", kXTitleJetPt, kYTitleRatio,
                   kGreen + 2, "MC / Data");
    }
}

// ---------------------------------------------------------------------------
// 5) balanceMap_mc / balanceMap_data / l3resMap  (TH2D)
// ---------------------------------------------------------------------------
void Draw2DMaps(TFile* f, const TString& outDir) {
    std::cout << "== Draw2DMaps ==" << std::endl;
    TString outDirReal = outDir + "/2DPlots/";
    gSystem->mkdir(outDirReal,true);

    struct MapDef {
        const char* name;
        const char* title;
        const char* zTitle;
    };
    std::vector<MapDef> maps = {
        {"balanceMap_mc",   "MC balance (#alpha < 0.30)",             "Balance"},
        {"balanceMap_data", "Data balance (#alpha < 0.30)",           "Balance"},
        {"l3resMap",        "L3 residual, MC/Data (#alpha < 0.30)",   "L3 Residual"},
    };

    for (const auto& m : maps) {
        TH2D* h = GetObj<TH2D>(f, m.name);
        if (!h) continue;
        h->SetTitle(m.title);
        h->GetXaxis()->SetTitle(kXTitlePtAvg);
        h->GetYaxis()->SetTitle(kYTitleAbsEta);
        h->GetZaxis()->SetTitle(m.zTitle);
        h->GetXaxis()->SetTitleSize(0.04);
        h->GetYaxis()->SetTitleSize(0.04);
        h->GetZaxis()->SetTitleSize(0.04);

        TCanvas c(TString("c_") + m.name, m.name, 850, 700);
        c.SetRightMargin(0.16);
        h->Draw("colz");
        c.SaveAs(outDirReal + "/" + m.name + filetype);
    }
}

// ---------------------------------------------------------------------------
// 6) balance3D_mc / balance3D_data (TProfile3D), counts3D_mc / counts3D_data (TH3D)
// ---------------------------------------------------------------------------
void Draw3DObjects(TFile* f, const TString& outDir) {
    std::cout << "== Draw3DObjects ==" << std::endl;

    TProfile3D* pMc   = GetObj<TProfile3D>(f, "balance3D_mc");
    TProfile3D* pData = GetObj<TProfile3D>(f, "balance3D_data");
    TH3D*       hMc   = GetObj<TH3D>(f, "counts3D_mc");
    TH3D*       hData = GetObj<TH3D>(f, "counts3D_data");

    // 3D objects don't project well onto a single view by default, so this
    // draws a labeled "box" overview; reslice interactively in a live ROOT
    // session (e.g. Project3D) for detailed 1D/2D looks at fixed alpha/eta.
    auto drawOne = [&](TH1* obj, const char* name, const char* title, Color_t color) {
        if (!obj) return;
        obj->SetTitle(title);
        obj->GetXaxis()->SetTitle(kXTitlePtAvg);
        obj->GetYaxis()->SetTitle(kYTitleAbsEta);
        obj->GetZaxis()->SetTitle(kXTitleAlpha);
        obj->SetLineColor(color);
        obj->SetMarkerColor(color);

        TCanvas c(TString("c_") + name, name, 850, 700);
        obj->Draw("box2z");
        TLegend* leg = MakeLegend(0.60, 0.80, 0.88, 0.90);
        leg->AddEntry(obj, name, "l");
        leg->Draw();
        c.SaveAs(outDir + "/" + TString(name) + filetype);
    };

    drawOne(pMc,   "balance3D_mc",   "MC balance vs p_{T,avg}, |#eta_{jet}|, #alpha", kBlue + 1);
    drawOne(pData, "balance3D_data", "Data balance vs p_{T,avg}, |#eta_{jet}|, #alpha", kRed + 1);
    drawOne(hMc,   "counts3D_mc",    "MC entries vs p_{T,avg}, |#eta_{jet}|, #alpha", kBlue + 2);
    drawOne(hData, "counts3D_data",  "Data entries vs p_{T,avg}, |#eta_{jet}|, #alpha", kRed + 2);
}
void DrawPtBinBalancesOverlay(TFile* f, const TString& outDir) {
    std::cout << "== DrawPtBinBalancesOverlay ==" << std::endl;
    TString outDirReal = outDir + "/PtBinBalances/";
    gSystem->mkdir(outDirReal, true);

    // A distinct color per pT bin (9 bins -> 9 colors)
    static const int kPtColors[] = {
        kBlue + 1, kRed + 1, kGreen + 2, kMagenta + 1, kOrange + 7,
        kCyan + 2, kViolet + 7, kAzure + 7, kSpring + 9
    };

    // One overlay canvas per histogram family: mc, data, l3res, ratio
    struct FamilyDef {
        const char* prefix;      // key prefix, e.g. "mc_balance_"
        const char* canvasTag;   // used in filename
        const char* title;
        const char* yTitle;
    };
    std::vector<FamilyDef> families = {
        {"mc_balance_",  "mc",    "MC balance, all p_{T} bins (#alpha < 0.30)",   kYTitleEntries},
        {"dt_balance_",  "data",  "Data balance, all p_{T} bins (#alpha < 0.30)", kYTitleEntries},
        {"l3res_",       "l3res", "L3 residual, all p_{T} bins (#alpha < 0.30)",  kYTitleL3Res},
        {"ratio_",       "ratio", "Ratio, all p_{T} bins (#alpha < 0.30)",        kYTitleRatio},
    };

    for (const auto& fam : families) {
        TCanvas c(TString("c_overlay_") + fam.canvasTag, fam.canvasTag, 900, 700);
        TLegend* leg = MakeLegend(0.68, 0.55, 0.90, 0.90);

        double ymax = -1e300, ymin = 1e300;
        std::vector<TH1D*> hists;

        // First pass: fetch all bins, track shared y-range
        for (size_t b = 0; b < kPtBinLabels.size(); ++b) {
            TString suffix = kPtBinLabels[b] + "_alpha0.30";
            TH1D* h = GetObj<TH1D>(f, TString(fam.prefix) + suffix);
            hists.push_back(h);
            if (h) {
                ymax = std::max(ymax, h->GetMaximum());
                ymin = std::min(ymin, h->GetMinimum());
            }
        }
        if (ymax < ymin) continue; // nothing found for this family

        double pad = 0.1 * (ymax - ymin);
        bool first = true;
        for (size_t b = 0; b < hists.size(); ++b) {
            TH1D* h = hists[b];
            if (!h) continue;
            Color_t color = kPtColors[b % 9];
            StyleHist(h, fam.title, kXTitleBalance, fam.yTitle, color);
            if (first) {
                h->GetYaxis()->SetRangeUser(ymin - pad, ymax + pad);
                h->Draw("hist");
                first = false;
            } else {
                h->Draw("hist same");
            }
            leg->AddEntry(h, kPtBinPretty[b], "l");
        }
        leg->Draw();
        c.SaveAs(outDirReal + "/overlay_" + fam.canvasTag + "_alpha0.30" + filetype);
    }
}
void DrawVsAlphaOverlay(TFile* f, const TString& outDir) {
    std::cout << "== DrawVsAlphaOverlay ==" << std::endl;
    TString outDirReal = outDir + "/vsAlpha/";
    gSystem->mkdir(outDirReal, true);

    static const int kPtColors[] = {
        kBlue + 1, kRed + 1, kGreen + 2, kMagenta + 1, kOrange + 7,
        kCyan + 2, kViolet + 7, kAzure + 7, kSpring + 9
    };

    struct FamilyDef {
        const char* prefix;    // key prefix, e.g. "Balance_vsa_mc_"
        const char* canvasTag; // used in filename
        const char* title;
        const char* yTitle;
    };
    std::vector<FamilyDef> families = {
        {"Balance_vsa_mc_",   "mc",     "MC balance vs #alpha, all p_{T} bins",             kYTitleBalance},
        {"Balance_vsa_data_", "data",   "Data balance vs #alpha, all p_{T} bins",           kYTitleBalance},
        {"L3Res_vsa_",        "l3res",  "L3 residual vs #alpha, all p_{T} bins",            kYTitleL3Res},
        {"L3Res_vsa_norm_",   "l3norm", "Normalized L3 residual vs #alpha, all p_{T} bins", kYTitleRatioNorm},
    };

    for (int j = 1; j <= kNEtaBins; ++j) {
        for (const auto& fam : families) {
            std::vector<TH1D*> hists;
            double ymax = -1e300, ymin = 1e300;

            for (size_t i = 1; i <= kPtBinLabels.size(); ++i) {
                TString tag = TString::Format("%zu_%d", i, j);
                TH1D* h = GetObj<TH1D>(f, TString(fam.prefix) + tag);
                hists.push_back(h);
                if (h) {
                    ymax = std::max(ymax, h->GetMaximum());
                    ymin = std::min(ymin, h->GetMinimum());
                }
            }
            if (ymax < ymin) continue; // nothing found for this (family, eta bin)

            TString etaTag = TString::Format("eta%d", j);
            TCanvas c(TString("c_vsa_overlay_") + fam.canvasTag + "_" + etaTag,
                      TString(fam.canvasTag) + "_" + etaTag, 900, 700);
            TLegend* leg = MakeLegend(0.68, 0.50, 0.90, 0.90);

            double pad = 0.1 * (ymax - ymin);
            bool first = true;
            for (size_t i = 0; i < hists.size(); ++i) {
                TH1D* h = hists[i];
                if (!h) continue;
                Color_t color = kPtColors[i % 9];
                TString title = TString(fam.title) + TString::Format(", #eta bin %d", j);
                StyleHist(h, title, kXTitleAlpha, fam.yTitle, color);
                if (first) {
                    h->GetYaxis()->SetRangeUser(ymin - pad, ymax + pad);
                    h->Draw("hist");
                    first = false;
                } else {
                    h->Draw("hist same");
                }
                leg->AddEntry(h, kPtBinPretty[i], "l");
            }
            leg->Draw();
            c.SaveAs(outDirReal + "/overlay_vsa_" + fam.canvasTag + "_" + etaTag + filetype);
        }
    }
}
// ---------------------------------------------------------------------------
// Driver
// ---------------------------------------------------------------------------
void DrawAll(const TString& inputFile = "/eos/cms/store/group/phys_heavyions/xirong/L3Residuals/0902BharadFolder//0902BharadInputs.root",
	       	const TString& outDir = "/afs/cern.ch/user/x/xirong/CMSSW_14_1_9_patch2/src/residualanalysis/L3Residual/PlotL3Intermediate/0902BharadPlots/") {
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(1);
    gStyle->SetTitleFontSize(0.05);
    gStyle->SetPadLeftMargin(0.14);
    EnsureDir(outDir);

    TFile* f = TFile::Open(inputFile);
    if (!f || f->IsZombie()) {
        std::cerr << "Could not open file: " << inputFile << std::endl;
        return;
    }

    DrawPtBinBalances(f, outDir);
    DrawPtBinBalancesOverlay(f, outDir);   // <-- add this line

    DrawVsAlpha(f, outDir);
    DrawVsAlphaOverlay(f, outDir);   // <-- add this
    
    DrawVsPhotonPt(f, outDir);
    DrawVsJetPt(f, outDir);
    Draw2DMaps(f, outDir);
    Draw3DObjects(f, outDir);

    std::cout << "Done. Plots written to: " << outDir << std::endl;
    f->Close();
}

// Allows `root -l -q draw_histograms.cpp` with default args, or
// `root -l -q 'draw_histograms.cpp("myfile.root","myplots")'`
void DrawL3ResidualBalance() { DrawAll(); }
void DrawL3ResidualBalance(const TString& inputFile, const TString& outDir = "plots") {
    DrawAll(inputFile, outDir);
}
