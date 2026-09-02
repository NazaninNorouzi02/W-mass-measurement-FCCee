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

void drawPanel(TTree *tree,
               int pad,
               const char *histName,
               const char *drawExpr,
               const TString &cut,
               const char *title) {
    gPad->cd(pad);
    tree->Draw(drawExpr, cut, "COLZ");
    TH2D *h = (TH2D*)gDirectory->Get(histName);
    if (h) h->SetTitle(title);
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

    file->Close();
}

void isr_truth_correlation(const char *filename,
                           const char *treename = "events") {
    isr_truth_correlation_fixed(filename, treename);
}