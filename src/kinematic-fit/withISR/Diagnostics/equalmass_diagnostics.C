#include <TFile.h>
#include <TGraph.h>
#include <TMath.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

namespace {
bool requireBranches(TTree* tree, const std::vector<const char*>& names) {
    bool ok = true;
    for (const char* name : names) {
        if (!tree->GetBranch(name)) {
            std::cerr << "Missing branch: " << name << "\n";
            ok = false;
        }
    }
    if (!ok) {
        std::cerr << "This is a stale ROOT file. Re-run the revised producer.\n";
    }
    return ok;
}

void summarize(TTree* tree, const char* label, const char* expression,
               const TString& cut) {
    const Long64_t n = tree->Draw(expression, cut, "goff");
    if (n <= 0) {
        std::cout << label << ": N=0\n";
        return;
    }
    const double* values = tree->GetV1();
    std::cout << label
              << ": N=" << n
              << " mean=" << TMath::Mean(n, values)
              << " RMS=" << TMath::RMS(n, values)
              << " median=" << TMath::Median(n, values)
              << "\n";
}

void robustSummary(TTree* tree, const char* label, const TString& cut) {
    const Long64_t n = tree->Draw("z_delta_m_4C", cut, "goff");
    if (n <= 0) {
        std::cout << label << ": N=0\n";
        return;
    }
    std::vector<double> values(tree->GetV1(), tree->GetV1() + n);
    std::sort(values.begin(), values.end());
    const double q16 = values[static_cast<size_t>(0.16 * n)];
    const double q50 = values[static_cast<size_t>(0.50 * n)];
    const double q84 = values[static_cast<size_t>(0.84 * n)];
    std::cout << label
              << ": N=" << n
              << " q16=" << q16
              << " q50=" << q50
              << " q84=" << q84
              << " robust_sigma=" << 0.5 * (q84 - q16)
              << "\n";
}

void printTailFractions(TTree* tree, const char* label, const TString& cut) {
    const Long64_t n = tree->GetEntries(cut);
    if (n <= 0) {
        std::cout << label << ": N=0\n";
        return;
    }
    const double frac10 =
        100.0 * tree->GetEntries(cut + " && abs(z_delta_m_4C)>10") /
        static_cast<double>(n);
    const double frac50 =
        100.0 * tree->GetEntries(cut + " && abs(z_delta_m_4C)>50") /
        static_cast<double>(n);
    std::cout << label
              << ": frac |z|>10 " << frac10
              << "%  frac |z|>50 " << frac50 << "%\n";
}

void summarizeSigmaForZRegion(TTree* tree, const char* label,
                              const TString& cut) {
    const Long64_t n = tree->Draw("sigma_delta_m_4C", cut, "goff");
    if (n <= 0) {
        std::cout << label << ": N=0\n";
        return;
    }
    const double* values = tree->GetV1();
    std::cout << label
              << ": N=" << n
              << " median_sigma=" << TMath::Median(n, values)
              << " mean_sigma=" << TMath::Mean(n, values)
              << " RMS_sigma=" << TMath::RMS(n, values)
              << "\n";
}
}  // namespace

void equalmass_diagnostics(const char* filename,
                           const char* treeName = "events") {
    TFile file(filename, "READ");
    if (file.IsZombie()) {
        std::cerr << "Cannot open " << filename << "\n";
        return;
    }
    auto* tree = dynamic_cast<TTree*>(file.Get(treeName));
    if (!tree) {
        std::cerr << "Cannot find tree " << treeName << "\n";
        return;
    }

    const std::vector<const char*> required = {
        "standard_4C_valid", "standard_5C_valid",
        "selected_truth_pairing_correct", "E_isr_true",
        "E_isr_true_collinear", "prob_4C", "sigma_delta_m_4C",
        "pull_delta_m_4C", "z_delta_m_4C", "chi2_5C_minus_4C",
        "predicted_equalmass_chi2_4C", "equalmass_chi2_ratio",
        "E_side1_bias", "E_side2_bias", "truth_match_angular_cost",
        "truth_w_mass_1", "truth_w_mass_2"
    };
    if (!requireBranches(tree, required)) return;

    tree->SetEstimate(tree->GetEntries() + 1);

    const TString zBase =
        "standard_4C_valid>0.5"
        " && selected_truth_pairing_correct>0.5"
        " && E_isr_true<2.0 && E_isr_true_collinear<2.0"
        " && sigma_delta_m_4C>0"
        " && z_delta_m_4C==z_delta_m_4C";

    std::cout << "\n=== Option A: covariance calibration ===\n";
    summarize(tree, "z_dm_4C, all valid 4C", "z_delta_m_4C", zBase);
    summarize(tree, "z_dm_4C, prob_4C>0.03", "z_delta_m_4C",
              zBase + " && prob_4C>0.03");

    const Long64_t nZ = tree->GetEntries(zBase);
    if (nZ > 0) {
        const double within1 = tree->GetEntries(zBase + " && abs(z_delta_m_4C)<1") /
                               static_cast<double>(nZ);
        const double within2 = tree->GetEntries(zBase + " && abs(z_delta_m_4C)<2") /
                               static_cast<double>(nZ);
        std::cout << "coverage: |z|<1 " << 100.0 * within1
                  << "%  |z|<2 " << 100.0 * within2 << "%\n";
    }

    std::cout << "\n=== Robust tail diagnostic for z_dm_4C ===\n";
    robustSummary(tree, "central quantiles, all valid 4C", zBase);
    printTailFractions(tree, "tail fractions, all valid 4C", zBase);
    summarize(tree, "z_dm_4C with |z|<10 only", "z_delta_m_4C",
              zBase + " && abs(z_delta_m_4C)<10");
    summarizeSigmaForZRegion(tree, "sigma_delta_m_4C for |z|>10",
                             zBase + " && abs(z_delta_m_4C)>10");
    summarizeSigmaForZRegion(tree, "sigma_delta_m_4C for |z|<2",
                             zBase + " && abs(z_delta_m_4C)<2");

    const TString zNarrow =
        zBase + " && abs(truth_w_mass_1-truth_w_mass_2)<0.25";
    std::cout << "\n=== W-width cross-check: near-degenerate truth masses ===\n";
    summarize(tree, "z_dm_4C, |mW1_true-mW2_true|<0.25", "z_delta_m_4C",
              zNarrow);
    const Long64_t nZn = tree->GetEntries(zNarrow);
    if (nZn > 0) {
        const double within1n =
            tree->GetEntries(zNarrow + " && abs(z_delta_m_4C)<1") /
            static_cast<double>(nZn);
        const double within2n =
            tree->GetEntries(zNarrow + " && abs(z_delta_m_4C)<2") /
            static_cast<double>(nZn);
        std::cout << "coverage near truth: |z|<1 " << 100.0 * within1n
                  << "%  |z|<2 " << 100.0 * within2n << "%\n";
    }

    std::cout << "\n=== Robust tail diagnostic, near-degenerate truth masses ===\n";
    robustSummary(tree, "central quantiles, near truth", zNarrow);
    printTailFractions(tree, "tail fractions, near truth", zNarrow);
    summarize(tree, "near truth z_dm_4C with |z|<10 only", "z_delta_m_4C",
              zNarrow + " && abs(z_delta_m_4C)<10");
    summarizeSigmaForZRegion(tree, "near truth sigma_delta_m_4C for |z|>10",
                             zNarrow + " && abs(z_delta_m_4C)>10");
    summarizeSigmaForZRegion(tree, "near truth sigma_delta_m_4C for |z|<2",
                             zNarrow + " && abs(z_delta_m_4C)<2");

    const TString nested =
        zBase +
        " && standard_5C_valid>0.5"
        " && chi2_5C_minus_4C>=0"
        " && predicted_equalmass_chi2_4C>1e-12";
    std::cout << "\n=== Mechanical 4C-to-5C consistency ===\n";
    summarize(tree, "null pull", "pull_delta_m_4C", nested);
    summarize(tree, "dchi2/predicted", "equalmass_chi2_ratio", nested);

    std::cout << "\n=== Option C: upstream energy-transfer clue ===\n";
    const Long64_t nC = tree->Draw("-E_side2_bias:E_side1_bias", zBase, "goff");
    if (nC > 1) {
        TGraph graph(nC, tree->GetV2(), tree->GetV1());
        std::cout << "N=" << nC
                  << " corr(E_side1_bias,-E_side2_bias)="
                  << graph.GetCorrelationFactor() << "\n";
    } else {
        std::cout << "N=0\n";
    }
    std::cout << "Option C is suggestive only; it also contains showering, "
                 "hadronization and detector effects.\n";
}

void equalmass_covariance_test(const char* filename,
                               const char* treeName = "events") {
    equalmass_diagnostics(filename, treeName);
}