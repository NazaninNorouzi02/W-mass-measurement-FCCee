//root -l -b -q 'isr_truth_correlation.C("wmass_fit_pvalue_pulls_ISR_ecm365.root")'

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TTree.h"
#include "TH1D.h"
#include "TLegend.h"

struct CorrResult {
    Long64_t n = 0;
    double meanX = 0.0;
    double meanY = 0.0;
    double rmsX = 0.0;
    double rmsY = 0.0;
    double rho = 0.0;
};

CorrResult drawCorrelation(TTree *tree,
                           const char *xexpr,
                           const char *yexpr,
                           const TString &cut) {
    CorrResult r;
    const TString drawExpr = TString::Format("%s:%s", yexpr, xexpr);
    r.n = tree->Draw(drawExpr, cut, "goff");
    if (r.n <= 1) return r;

    const double *y = tree->GetV1();
    const double *x = tree->GetV2();

    double sx = 0.0, sy = 0.0;
    for (Long64_t i = 0; i < r.n; ++i) {
        sx += x[i];
        sy += y[i];
    }
    r.meanX = sx / r.n;
    r.meanY = sy / r.n;

    double vx = 0.0, vy = 0.0, cxy = 0.0;
    for (Long64_t i = 0; i < r.n; ++i) {
        const double dx = x[i] - r.meanX;
        const double dy = y[i] - r.meanY;
        vx += dx * dx;
        vy += dy * dy;
        cxy += dx * dy;
    }

    r.rmsX = std::sqrt(vx / r.n);
    r.rmsY = std::sqrt(vy / r.n);
    if (vx > 0.0 && vy > 0.0) {
        r.rho = cxy / std::sqrt(vx * vy);
    }
    return r;
}

void printCorrelation(const char *label, const CorrResult &r) {
    std::cout << label
              << ": N=" << r.n
              << " rho=" << r.rho
              << " mean_true=" << r.meanX
              << " mean_fit=" << r.meanY
              << " rms_true=" << r.rmsX
              << " rms_fit=" << r.rmsY
              << "\n";
}

void drawResidualSlices(TTree *tree,
                        const TString &baseCut,
                        const char *residualExpr,
                        const char *canvasName,
                        const char *title,
                        const char *outName,
                        double xmin,
                        double xmax) {
    const int colors[4] = {kBlue + 1, kGreen + 2, kOrange + 7, kRed + 1};

    const char *sliceLabels[4] = {
        "true ISR < 5 GeV",
        "5 <= true ISR < 20 GeV",
        "20 <= true ISR < 50 GeV",
        "true ISR >= 50 GeV"
    };

    const char *sliceCuts[4] = {
        "E_isr_true_collinear < 5",
        "E_isr_true_collinear >= 5 && E_isr_true_collinear < 20",
        "E_isr_true_collinear >= 20 && E_isr_true_collinear < 50",
        "E_isr_true_collinear >= 50"
    };

    TCanvas *c = new TCanvas(canvasName, title, 1000, 750);

    TLegend *leg = new TLegend(0.58, 0.68, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.030);

    double maxY = 0.0;
    TH1D *hist[4];

    for (int i = 0; i < 4; ++i) {
        TString hname = TString::Format("%s_h%d", canvasName, i);

        hist[i] = new TH1D(
            hname,
            title,
            80,
            xmin,
            xmax
        );

        TString cut = baseCut + " && " + sliceCuts[i];

        tree->Draw(
            TString::Format("(%s)>>%s", residualExpr, hname.Data()),
            cut,
            "goff"
        );

        const double n = hist[i]->GetEntries();

        hist[i]->SetLineColor(colors[i]);
        hist[i]->SetLineWidth(2);

        if (n > 0) {
            hist[i]->Scale(1.0 / hist[i]->Integral());
        }

        if (hist[i]->GetMaximum() > maxY) {
            maxY = hist[i]->GetMaximum();
        }

        leg->AddEntry(
            hist[i],
            TString::Format("%s, N = %.0f", sliceLabels[i], n),
            "l"
        );
    }

    for (int i = 0; i < 4; ++i) {
        hist[i]->SetMaximum(maxY * 1.25);
        hist[i]->GetXaxis()->SetTitle("fitted ISR - true ISR [GeV]");
        hist[i]->GetYaxis()->SetTitle("normalized entries");

        if (i == 0) {
            hist[i]->Draw("HIST");
        } else {
            hist[i]->Draw("HIST SAME");
        }
    }

    leg->Draw();

    c->SaveAs(TString::Format("%s.png", outName));
    c->SaveAs(TString::Format("%s.pdf", outName));
}


void isr_truth_correlation_fixed(const char *filename,
                                 const char *treename = "events") {
    gStyle->SetOptStat(0);

    TFile *file = TFile::Open(filename);
    if (!file || file->IsZombie()) {
        std::cerr << "ERROR: cannot open " << filename << "\n";
        return;
    }

    TTree *tree = (TTree*)file->Get(treename);
    if (!tree) {
        std::cerr << "ERROR: cannot find tree '" << treename << "'\n";
        return;
    }

    tree->SetAlias(
        "pz_isr_4C_fitted",
        "(y_isr_4C_fitted>=0 ? E_isr_4C_fitted : -E_isr_4C_fitted)"
    );
    tree->SetAlias(
        "pz_isr_5C_fitted",
        "(y_isr_5C_fitted>=0 ? E_isr_fitted : -E_isr_fitted)"
    );

    const TString cut4C =
        "standard_4C_valid>0.5 && isr_4C_fit_valid>0.5 && isr_4C_applied>0.5";
    const TString cut5C =
        "standard_5C_valid>0.5 && isr_5C_fit_valid>0.5 && isr_applied>0.5";

    std::cout << "\nISR truth-correlation diagnostic\n";
    std::cout << "File: " << filename << "\n";
    std::cout << "Interpretation: signed pz correlation is the main test of "
                 "event-by-event physical ISR recovery.\n\n";

    printCorrelation(
        "4C+ISR: fitted E vs true |pz_ISR|",
        drawCorrelation(tree, "E_isr_true_collinear", "E_isr_4C_fitted", cut4C)
    );
    printCorrelation(
        "4C+ISR: fitted signed pz vs true pz",
        drawCorrelation(tree, "pz_isr_true", "pz_isr_4C_fitted", cut4C)
    );
    printCorrelation(
        "5C+ISR: fitted E vs true |pz_ISR|",
        drawCorrelation(tree, "E_isr_true_collinear", "E_isr_fitted", cut5C)
    );
    printCorrelation(
        "5C+ISR: fitted signed pz vs true pz",
        drawCorrelation(tree, "pz_isr_true", "pz_isr_5C_fitted", cut5C)
    );

    TCanvas *c = new TCanvas(
        "c_isr_truth_correlation",
        "ISR truth correlation",
        1400,
        1000
    );
    c->Divide(2, 2);

    c->cd(1);
    tree->Draw(
        "E_isr_4C_fitted:E_isr_true_collinear>>h4cE(80,0,180,80,0,180)",
        cut4C,
        "COLZ"
    );
    ((TH2D*)gDirectory->Get("h4cE"))->SetTitle(
        "4C+ISR: fitted vs true ISR energy;true |p_{z,ISR}| [GeV];fitted E_{ISR} [GeV]"
    );

    c->cd(2);
    tree->Draw(
        "pz_isr_4C_fitted:pz_isr_true>>h4cpz(80,-180,180,80,-180,180)",
        cut4C,
        "COLZ"
    );
    ((TH2D*)gDirectory->Get("h4cpz"))->SetTitle(
        "4C+ISR: signed fitted vs true ISR p_{z};true p_{z,ISR} [GeV];fitted p_{z,ISR} [GeV]"
    );

    c->cd(3);
    tree->Draw(
        "E_isr_fitted:E_isr_true_collinear>>h5cE(80,0,180,80,0,180)",
        cut5C,
        "COLZ"
    );
    ((TH2D*)gDirectory->Get("h5cE"))->SetTitle(
        "5C+ISR: fitted vs true ISR energy;true |p_{z,ISR}| [GeV];fitted E_{ISR} [GeV]"
    );

    c->cd(4);
    tree->Draw(
        "pz_isr_5C_fitted:pz_isr_true>>h5cpz(80,-180,180,80,-180,180)",
        cut5C,
        "COLZ"
    );
    ((TH2D*)gDirectory->Get("h5cpz"))->SetTitle(
        "5C+ISR: signed fitted vs true ISR p_{z};true p_{z,ISR} [GeV];fitted p_{z,ISR} [GeV]"
    );

    c->SaveAs("isr_truth_correlation.png");
    c->SaveAs("isr_truth_correlation.pdf");
    
    const TString lowISR = "E_isr_true_collinear < 20";

    const TString cut4C_low =
        cut4C + " && " + lowISR;

    const TString cut5C_low =
        cut5C + " && " + lowISR;

    std::cout << "\nLow true-ISR diagnostic: E_isr_true_collinear < 20 GeV\n\n";

    printCorrelation(
        "4C+ISR true ISR < 20 GeV: fitted E vs true |pz_ISR|",
        drawCorrelation(tree, "E_isr_true_collinear", "E_isr_4C_fitted", cut4C_low)
    );

    printCorrelation(
        "4C+ISR true ISR < 20 GeV: fitted signed pz vs true pz",
        drawCorrelation(tree, "pz_isr_true", "pz_isr_4C_fitted", cut4C_low)
    );

    printCorrelation(
        "5C+ISR true ISR < 20 GeV: fitted E vs true |pz_ISR|",
        drawCorrelation(tree, "E_isr_true_collinear", "E_isr_fitted", cut5C_low)
    );

    printCorrelation(
        "5C+ISR true ISR < 20 GeV: fitted signed pz vs true pz",
        drawCorrelation(tree, "pz_isr_true", "pz_isr_5C_fitted", cut5C_low)
    );

    TCanvas *c_low = new TCanvas(
        "c_isr_truth_correlation_trueISR_lt20",
        "ISR truth correlation, true ISR below 20 GeV",
        1400,
        1000
    );

    c_low->Divide(2, 2);

    c_low->cd(1);
    tree->Draw(
        "E_isr_4C_fitted:E_isr_true_collinear>>h4cE_lt20(80,0,20,80,0,180)",
        cut4C_low,
        "COLZ"
    );
    ((TH2D*)gDirectory->Get("h4cE_lt20"))->SetTitle(
        "4C+ISR: true ISR < 20 GeV;true |p_{z,ISR}| [GeV];fitted E_{ISR} [GeV]"
    );

    c_low->cd(2);
    tree->Draw(
        "pz_isr_4C_fitted:pz_isr_true>>h4cpz_lt20(80,-20,20,80,-180,180)",
        cut4C_low,
        "COLZ"
    );
    ((TH2D*)gDirectory->Get("h4cpz_lt20"))->SetTitle(
        "4C+ISR: true ISR < 20 GeV;true p_{z,ISR} [GeV];fitted p_{z,ISR} [GeV]"
    );

    c_low->cd(3);
    tree->Draw(
        "E_isr_fitted:E_isr_true_collinear>>h5cE_lt20(80,0,20,80,0,180)",
        cut5C_low,
        "COLZ"
    );
    ((TH2D*)gDirectory->Get("h5cE_lt20"))->SetTitle(
        "5C+ISR: true ISR < 20 GeV;true |p_{z,ISR}| [GeV];fitted E_{ISR} [GeV]"
    );

    c_low->cd(4);
    tree->Draw(
        "pz_isr_5C_fitted:pz_isr_true>>h5cpz_lt20(80,-20,20,80,-180,180)",
        cut5C_low,
        "COLZ"
    );
    ((TH2D*)gDirectory->Get("h5cpz_lt20"))->SetTitle(
        "5C+ISR: true ISR < 20 GeV;true p_{z,ISR} [GeV];fitted p_{z,ISR} [GeV]"
    );

    c_low->SaveAs("isr_truth_correlation_trueISR_lt20.png");
    c_low->SaveAs("isr_truth_correlation_trueISR_lt20.pdf");

        std::cout << "\nResidual diagnostic by true ISR energy slice\n";
    std::cout << "Residual means fitted ISR minus true ISR.\n";
    std::cout << "Histograms are normalized, but the legend gives the event count in each slice.\n\n";

    drawResidualSlices(
        tree,
        cut4C,
        "E_isr_4C_fitted - E_isr_true_collinear",
        "c_residual_energy_4C",
        "4C+ISR: energy residual by true ISR slice;fitted E_{ISR} - true |p_{z,ISR}| [GeV];normalized entries",
        "isr_energy_residual_slices_4C",
        -80.0,
        160.0
    );

    drawResidualSlices(
        tree,
        cut5C,
        "E_isr_fitted - E_isr_true_collinear",
        "c_residual_energy_5C",
        "5C+ISR: energy residual by true ISR slice;fitted E_{ISR} - true |p_{z,ISR}| [GeV];normalized entries",
        "isr_energy_residual_slices_5C",
        -80.0,
        160.0
    );

    drawResidualSlices(
        tree,
        cut4C,
        "pz_isr_4C_fitted - pz_isr_true",
        "c_residual_pz_4C",
        "4C+ISR: signed p_{z} residual by true ISR slice;fitted p_{z,ISR} - true p_{z,ISR} [GeV];normalized entries",
        "isr_signed_pz_residual_slices_4C",
        -180.0,
        180.0
    );

    drawResidualSlices(
        tree,
        cut5C,
        "pz_isr_5C_fitted - pz_isr_true",
        "c_residual_pz_5C",
        "5C+ISR: signed p_{z} residual by true ISR slice;fitted p_{z,ISR} - true p_{z,ISR} [GeV];normalized entries",
        "isr_signed_pz_residual_slices_5C",
        -180.0,
        180.0
    );

    file->Close();
}

void isr_truth_correlation(const char *filename,
                           const char *treename = "events") {
    isr_truth_correlation_fixed(filename, treename);
}



//4C+ISR: fitted E vs true |pz_ISR|: N=20696 rho=0.580181 mean_true=22.1708 mean_fit=36.1308 rms_true=33.4227 rms_fit=40.5297
//4C+ISR: fitted signed pz vs true pz: N=20696 rho=0.678722 mean_true=-0.138269 mean_fit=0.20189 rms_true=40.1074 rms_fit=54.2959
//5C+ISR: fitted E vs true |pz_ISR|: N=40744 rho=0.595123 mean_true=12.9479 mean_fit=20.2724 rms_true=26.9927 rms_fit=31.9976
//5C+ISR: fitted signed pz vs true pz: N=40744 rho=0.5896 mean_true=-0.0213825 mean_fit=-0.0939875 rms_true=29.9375 rms_fit=37.8789
//Info in <TCanvas::Print>: png file isr_truth_correlation.png has been created
//Info in <TCanvas::Print>: pdf file isr_truth_correlation.pdf has been created

//Low true-ISR diagnostic: E_isr_true_collinear < 20 GeV

//4C+ISR true ISR < 20 GeV: fitted E vs true |pz_ISR|: N=14110 rho=-0.0668811 mean_true=3.67545 mean_fit=24.3378 rms_true=5.37706 rms_fit=36.1306
//4C+ISR true ISR < 20 GeV: fitted signed pz vs true pz: N=14110 rho=0.11893 mean_true=-0.0105033 mean_fit=0.438784 rms_true=6.51318 rms_fit=43.5609
//5C+ISR true ISR < 20 GeV: fitted E vs true |pz_ISR|: N=33219 rho=0.0880655 mean_true=2.23158 mean_fit=12.9252 rms_true=4.32008 rms_fit=24.3346
//5C+ISR true ISR < 20 GeV: fitted signed pz vs true pz: N=33219 rho=0.147685 mean_true=0.0166447 mean_fit=0.0504333 rms_true=4.86238 rms_fit=27.5541
