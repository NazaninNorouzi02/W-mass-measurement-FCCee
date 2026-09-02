// isr_truth_correlation_slices.C
#include <iostream>
#include <cmath>
#include "TFile.h"
#include "TTree.h"

struct R {
    Long64_t n;
    double rho;
    double meanTrue;
    double meanFit;
};

R corr(TTree* t, const char* xexpr, const char* yexpr, TString cut) {
    R r{0, 0, 0, 0};

    r.n = t->Draw(Form("%s:%s", yexpr, xexpr), cut, "goff");
    if (r.n <= 2) return r;

    const double* y = t->GetV1();
    const double* x = t->GetV2();

    double sx=0, sy=0;
    for (Long64_t i=0; i<r.n; i++) {
        sx += x[i];
        sy += y[i];
    }

    double mx = sx/r.n;
    double my = sy/r.n;
    r.meanTrue = mx;
    r.meanFit = my;

    double vx=0, vy=0, cxy=0;
    for (Long64_t i=0; i<r.n; i++) {
        double dx = x[i]-mx;
        double dy = y[i]-my;
        vx += dx*dx;
        vy += dy*dy;
        cxy += dx*dy;
    }

    if (vx > 0 && vy > 0) r.rho = cxy / std::sqrt(vx*vy);
    return r;
}

void print(const char* label, R r) {
    std::cout << label
              << " N=" << r.n
              << " rho=" << r.rho
              << " mean_true=" << r.meanTrue
              << " mean_fit=" << r.meanFit
              << "\n";
}

void isr_truth_correlation_slices(const char* filename,
                                  const char* treename="events") {
    TFile* f = TFile::Open(filename);
    TTree* t = (TTree*)f->Get(treename);

    t->SetAlias("pz_isr_4C_fitted",
                "(y_isr_4C_fitted>=0 ? E_isr_4C_fitted : -E_isr_4C_fitted)");
    t->SetAlias("pz_isr_5C_fitted",
                "(y_isr_5C_fitted>=0 ? E_isr_fitted : -E_isr_fitted)");

    TString base4C = "standard_4C_valid>0.5 && isr_4C_applied>0.5";
    TString base5C = "standard_5C_valid>0.5 && isr_applied>0.5";

    TString low  = "E_isr_true_collinear < 5";
    TString mid  = "E_isr_true_collinear >= 5 && E_isr_true_collinear < 20";
    TString high = "E_isr_true_collinear >= 20";
    TString very = "E_isr_true_collinear >= 50";

    std::cout << "\n=== 4C+ISR signed pz correlation by true ISR ===\n";
    print("true ISR < 5 GeV:    ", corr(t, "pz_isr_true", "pz_isr_4C_fitted", base4C+" && "+low));
    print("5 <= true ISR < 20: ", corr(t, "pz_isr_true", "pz_isr_4C_fitted", base4C+" && "+mid));
    print("true ISR >= 20 GeV: ", corr(t, "pz_isr_true", "pz_isr_4C_fitted", base4C+" && "+high));
    print("true ISR >= 50 GeV: ", corr(t, "pz_isr_true", "pz_isr_4C_fitted", base4C+" && "+very));

    std::cout << "\n=== 5C+ISR signed pz correlation by true ISR ===\n";
    print("true ISR < 5 GeV:    ", corr(t, "pz_isr_true", "pz_isr_5C_fitted", base5C+" && "+low));
    print("5 <= true ISR < 20: ", corr(t, "pz_isr_true", "pz_isr_5C_fitted", base5C+" && "+mid));
    print("true ISR >= 20 GeV: ", corr(t, "pz_isr_true", "pz_isr_5C_fitted", base5C+" && "+high));
    print("true ISR >= 50 GeV: ", corr(t, "pz_isr_true", "pz_isr_5C_fitted", base5C+" && "+very));
}