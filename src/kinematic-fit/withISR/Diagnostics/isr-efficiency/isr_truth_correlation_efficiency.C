// ISR truth-correlation and conditional-efficiency report.
// Kept separate from the general W-mass diagnostic by design.

//root -l -b -q 'isr_truth_correlation_efficiency.C("wmass_fit_pvalue_pulls_ISR_ecm162p5.root","162.5 GeV","isr_162p5")'

//root -l -b -q 'isr_truth_correlation_efficiency.C("wmass_fit_pvalue_pulls_ISR_ecm240.root","240 GeV","isr_240")'

//root -l -b -q 'isr_truth_correlation_efficiency.C("wmass_fit_pvalue_pulls_ISR_ecm365.root","365 GeV","isr_365")'

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "TFile.h"
#include "TString.h"
#include "TTree.h"

namespace ISRReport {

struct Correlation {
    Long64_t n = 0;
    double rho = std::numeric_limits<double>::quiet_NaN();
    double meanTrue = std::numeric_limits<double>::quiet_NaN();
    double meanFit = std::numeric_limits<double>::quiet_NaN();
    double bias = std::numeric_limits<double>::quiet_NaN();
    double residualRMS = std::numeric_limits<double>::quiet_NaN();
};

struct Row {
    std::string energy, fit, slice, sliceCut;
    Long64_t nStandard = 0, nEligible = 0, nAttempted = 0;
    Long64_t nFitValid = 0, nApplied = 0, nRecovered = 0;
    Correlation signedPz, energyCorr;
};

bool has(TTree *t, const char *name) { return t->GetBranch(name) != nullptr; }

double percent(Long64_t n, Long64_t d) {
    return d > 0 ? 100.0 * static_cast<double>(n) / d
                 : std::numeric_limits<double>::quiet_NaN();
}

std::string csvNumber(double x) {
    if (!std::isfinite(x)) return "";
    std::ostringstream s;
    s << std::setprecision(10) << x;
    return s.str();
}

std::string mdNumber(double x, int precision = 2) {
    if (!std::isfinite(x)) return "n/a";
    std::ostringstream s;
    s << std::fixed << std::setprecision(precision) << x;
    return s.str();
}

TString join(const TString &a, const TString &b) {
    if (a.IsNull() || a == "1") return b;
    if (b.IsNull() || b == "1") return a;
    return "(" + a + ") && (" + b + ")";
}

Correlation corr(TTree *t, const char *xexpr, const char *yexpr,
                 const TString &cut) {
    Correlation r;
    r.n = t->Draw(TString::Format("%s:%s", yexpr, xexpr), cut, "goff");
    if (r.n <= 0) return r;
    const double *y = t->GetV1();
    const double *x = t->GetV2();
    double sx = 0.0, sy = 0.0;
    for (Long64_t i = 0; i < r.n; ++i) { sx += x[i]; sy += y[i]; }
    r.meanTrue = sx / r.n;
    r.meanFit = sy / r.n;
    r.bias = r.meanFit - r.meanTrue;
    double vx = 0.0, vy = 0.0, cxy = 0.0, residual2 = 0.0;
    for (Long64_t i = 0; i < r.n; ++i) {
        const double dx = x[i] - r.meanTrue;
        const double dy = y[i] - r.meanFit;
        vx += dx * dx; vy += dy * dy; cxy += dx * dy;
        const double residual = y[i] - x[i];
        residual2 += residual * residual;
    }
    if (r.n > 1 && vx > 0.0 && vy > 0.0) r.rho = cxy/std::sqrt(vx*vy);
    r.residualRMS = std::sqrt(residual2/r.n);
    return r;
}

Row makeRow(TTree *t, const std::string &energy, const std::string &fit,
            const std::string &slice, const TString &sliceCut) {
    const bool is4C = fit == "4C+ISR";
    const TString standard = is4C ? "standard_4C_valid>0.5" : "standard_5C_valid>0.5";
    const TString probability = is4C ? "prob_4C" : "prob_5C";
    const TString attempted = is4C ? "isr_4C_attempted>0.5" : "isr_5C_attempted>0.5";
    const TString fitValid = is4C ? "isr_4C_fit_valid>0.5" : "isr_5C_fit_valid>0.5";
    const TString applied = is4C ? "isr_4C_applied>0.5" : "isr_applied>0.5";
    const TString recovered = is4C ? "isr_4C_recovered>0.5" : "isr_5C_recovered>0.5";
    const char *fitE = is4C ? "E_isr_4C_fitted" : "E_isr_fitted";
    const char *fitPz = is4C ? "pz_isr_4C_fitted" : "pz_isr_5C_fitted";

    Row r;
    r.energy = energy; r.fit = fit; r.slice = slice; r.sliceCut = sliceCut.Data();
    const TString standardSlice = join(standard, sliceCut);
    const TString eligible = join(standardSlice, probability + "<0.03");
    r.nStandard = t->GetEntries(standardSlice);
    r.nEligible = t->GetEntries(eligible);
    r.nAttempted = t->GetEntries(join(attempted, sliceCut));
    r.nFitValid = t->GetEntries(join(fitValid, sliceCut));
    r.nApplied = t->GetEntries(join(applied, sliceCut));
    r.nRecovered = t->GetEntries(join(recovered, sliceCut));
    const TString correlationCut = join(join(standard, applied), sliceCut);
    r.signedPz = corr(t, "pz_isr_true", fitPz, correlationCut);
    r.energyCorr = corr(t, "E_isr_true_collinear", fitE, correlationCut);
    return r;
}

void writeCSV(const std::vector<Row> &rows, const std::string &path) {
    std::ofstream out(path);
    out << "energy,fit,true_isr_slice,slice_cut,n_standard_valid,n_eligible_p_lt_003,"
           "n_attempted,n_fit_valid,n_applied,n_recovered,eligible_fraction_pct,"
           "attempted_over_eligible_pct,fit_valid_over_attempted_pct,"
           "applied_over_attempted_pct,recovered_over_applied_pct,"
           "applied_over_standard_pct,n_correlation,rho_signed_pz,rho_energy,"
           "mean_true_pz,mean_fitted_pz,pz_bias,pz_residual_rms,"
           "mean_true_energy,mean_fitted_energy,energy_bias,energy_residual_rms\n";
    for (const Row &r : rows) {
        out << r.energy << ',' << r.fit << ',' << r.slice << ",\"" << r.sliceCut << "\","
            << r.nStandard << ',' << r.nEligible << ',' << r.nAttempted << ','
            << r.nFitValid << ',' << r.nApplied << ',' << r.nRecovered << ','
            << csvNumber(percent(r.nEligible,r.nStandard)) << ','
            << csvNumber(percent(r.nAttempted,r.nEligible)) << ','
            << csvNumber(percent(r.nFitValid,r.nAttempted)) << ','
            << csvNumber(percent(r.nApplied,r.nAttempted)) << ','
            << csvNumber(percent(r.nRecovered,r.nApplied)) << ','
            << csvNumber(percent(r.nApplied,r.nStandard)) << ','
            << r.signedPz.n << ',' << csvNumber(r.signedPz.rho) << ','
            << csvNumber(r.energyCorr.rho) << ','
            << csvNumber(r.signedPz.meanTrue) << ',' << csvNumber(r.signedPz.meanFit) << ','
            << csvNumber(r.signedPz.bias) << ',' << csvNumber(r.signedPz.residualRMS) << ','
            << csvNumber(r.energyCorr.meanTrue) << ',' << csvNumber(r.energyCorr.meanFit) << ','
            << csvNumber(r.energyCorr.bias) << ',' << csvNumber(r.energyCorr.residualRMS) << '\n';
    }
}

void writeMarkdown(const std::vector<Row> &rows, const std::string &path,
                   const std::string &filename) {
    std::ofstream out(path);
    out << "# ISR truth-correlation and efficiency report\n\n"
        << "Input: `" << filename << "`\n\n"
        << "Correlations are calculated only for events on which the ISR fit was applied. "
           "Counts from different true-ISR slices are not compared as raw yields. Instead, "
           "each efficiency uses an explicit within-slice denominator.\n\n"
        << "- Eligible: standard fit valid and `0 < P < 0.03`.\n"
        << "- Attempt efficiency: attempted / eligible.\n"
        << "- Fit-valid efficiency: valid ISR solution / attempted.\n"
        << "- Recovery efficiency: final ISR probability above 3% / applied.\n"
        << "- Coverage: applied / standard-valid events in the same truth slice.\n\n"
        << "| Fit | True ISR slice | Standard valid | Eligible | Applied | Eligible fraction | Attempted/eligible | Fit-valid/attempted | Applied/attempted | Recovered/applied | Coverage | N corr. | rho signed pz | rho energy | pz bias [GeV] | pz residual RMS [GeV] |\n"
        << "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n";
    for (const Row &r : rows) {
        out << "| " << r.fit << " | " << r.slice << " | " << r.nStandard
            << " | " << r.nEligible << " | " << r.nApplied
            << " | " << mdNumber(percent(r.nEligible,r.nStandard)) << "%"
            << " | " << mdNumber(percent(r.nAttempted,r.nEligible)) << "%"
            << " | " << mdNumber(percent(r.nFitValid,r.nAttempted)) << "%"
            << " | " << mdNumber(percent(r.nApplied,r.nAttempted)) << "%"
            << " | " << mdNumber(percent(r.nRecovered,r.nApplied)) << "%"
            << " | " << mdNumber(percent(r.nApplied,r.nStandard)) << "%"
            << " | " << r.signedPz.n
            << " | " << mdNumber(r.signedPz.rho,3)
            << " | " << mdNumber(r.energyCorr.rho,3)
            << " | " << mdNumber(r.signedPz.bias)
            << " | " << mdNumber(r.signedPz.residualRMS) << " |\n";
    }
    out << "\nThe signed-pz coefficient is the main direction-sensitive test. The energy "
           "coefficient tests only the fitted magnitude. A high recovery efficiency does "
           "not by itself imply a high event-by-event truth correlation.\n";
}

} // namespace ISRReport

void isr_truth_correlation_efficiency(const char *filename,
                                      const char *energyLabel = "sample",
                                      const char *outputPrefix = "isr_truth_report",
                                      const char *treeName = "events") {
    using namespace ISRReport;
    TFile *file = TFile::Open(filename);
    if (!file || file->IsZombie()) {
        std::cerr << "ERROR: cannot open " << filename << '\n'; return;
    }
    TTree *tree = dynamic_cast<TTree*>(file->Get(treeName));
    if (!tree) {
        std::cerr << "ERROR: cannot find tree '" << treeName << "'\n"; return;
    }
    const std::vector<const char*> required = {
        "standard_4C_valid", "standard_5C_valid", "prob_4C", "prob_5C",
        "isr_4C_attempted", "isr_5C_attempted", "isr_4C_fit_valid",
        "isr_5C_fit_valid", "isr_4C_applied", "isr_applied",
        "isr_4C_recovered", "isr_5C_recovered", "E_isr_4C_fitted",
        "E_isr_fitted", "y_isr_4C_fitted", "y_isr_5C_fitted",
        "E_isr_true_collinear", "pz_isr_true"
    };
    for (const char *branch : required) {
        if (!has(tree, branch)) {
            std::cerr << "ERROR: missing branch " << branch << '\n'; return;
        }
    }
    tree->SetEstimate(tree->GetEntries()+1);
    tree->SetAlias("pz_isr_4C_fitted",
                   "(y_isr_4C_fitted>=0 ? E_isr_4C_fitted : -E_isr_4C_fitted)");
    tree->SetAlias("pz_isr_5C_fitted",
                   "(y_isr_5C_fitted>=0 ? E_isr_fitted : -E_isr_fitted)");

    const std::vector<std::pair<std::string,TString>> slices = {
        {"all", "1"},
        {"<5 GeV", "E_isr_true_collinear<5"},
        {"5-20 GeV", "E_isr_true_collinear>=5 && E_isr_true_collinear<20"},
        {"20-50 GeV", "E_isr_true_collinear>=20 && E_isr_true_collinear<50"},
        {">=50 GeV", "E_isr_true_collinear>=50"}
    };
    std::vector<Row> rows;
    for (const std::string fit : {std::string("4C+ISR"), std::string("5C+ISR")})
        for (const auto &slice : slices)
            rows.push_back(makeRow(tree, energyLabel, fit, slice.first, slice.second));

    const std::string prefix(outputPrefix);
    writeCSV(rows, prefix + ".csv");
    writeMarkdown(rows, prefix + ".md", filename);
    std::cout << "Wrote " << prefix << ".csv\n"
              << "Wrote " << prefix << ".md\n";
    file->Close();
}

// Backward-compatible short name.
void isr_truth_correlation_slices(const char *filename,
                                  const char *energyLabel = "sample",
                                  const char *outputPrefix = "isr_truth_report",
                                  const char *treeName = "events") {
    isr_truth_correlation_efficiency(filename, energyLabel, outputPrefix, treeName);
}
