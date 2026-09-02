#include <TFile.h>
#include <TTree.h>
#include <TString.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <vector>

void check_dchi2_tail(
    const char* fileName =
        "wmass_fit_pvalue_pulls_ISR_ecm365.root"
) {
    TFile file(fileName, "READ");

    if (file.IsZombie()) {
        std::cerr << "ERROR: could not open " << fileName << std::endl;
        return;
    }

    TTree* events = nullptr;
    file.GetObject("events", events);

    if (!events) {
        std::cerr << "ERROR: tree 'events' was not found." << std::endl;
        return;
    }

    const TString quiet =
        "standard_4C_valid>0.5"
        " && standard_5C_valid>0.5"
        " && selected_truth_pairing_correct==1"
        " && E_isr_true<2.0"
        " && E_isr_true_collinear<2.0"
        " && abs(truth_w_mass_1-truth_w_mass_2)<0.25";

    // Ensure TTree::Draw stores every selected value.
    events->SetEstimate(events->GetEntries() + 1);

    const Long64_t nSelected = events->Draw(
        "chi2_5C-chi2_4C",
        quiet,
        "goff"
    );

    if (nSelected <= 0) {
        std::cerr << "ERROR: no events passed the quiet selection."
                  << std::endl;
        return;
    }

    const double* drawnValues = events->GetV1();

    std::vector<double> values;
    values.reserve(nSelected);

    for (Long64_t i = 0; i < nSelected; ++i) {
        if (std::isfinite(drawnValues[i])) {
            values.push_back(drawnValues[i]);
        }
    }

    if (values.empty()) {
        std::cerr << "ERROR: no finite delta-chi2 values." << std::endl;
        return;
    }

    std::sort(values.begin(), values.end());

    const auto quantile = [&](double probability) {
        const double position =
            probability * static_cast<double>(values.size() - 1);

        const std::size_t lower =
            static_cast<std::size_t>(std::floor(position));
        const std::size_t upper =
            static_cast<std::size_t>(std::ceil(position));

        if (lower == upper) {
            return values[lower];
        }

        const double fraction = position - lower;
        return values[lower] * (1.0 - fraction)
             + values[upper] * fraction;
    };

    const double sum =
        std::accumulate(values.begin(), values.end(), 0.0);
    const double mean = sum / values.size();

    double variance = 0.0;
    for (const double value : values) {
        variance += (value - mean) * (value - mean);
    }
    variance /= values.size();

    const auto fractionAbove = [&](double threshold) {
        const auto first = std::upper_bound(
            values.begin(), values.end(), threshold
        );

        return static_cast<double>(
            std::distance(first, values.end())
        ) / values.size();
    };

    const auto fractionBelow = [&](double threshold) {
        const auto first = std::lower_bound(
            values.begin(), values.end(), threshold
        );

        return static_cast<double>(
            std::distance(values.begin(), first)
        ) / values.size();
    };

    std::cout << std::fixed << std::setprecision(6);

    std::cout << "\nQuiet, truth-correct-pairing sample\n";
    std::cout << "N                  = " << values.size() << '\n';
    std::cout << "minimum            = " << values.front() << '\n';
    std::cout << "10% quantile       = " << quantile(0.10) << '\n';
    std::cout << "25% quantile       = " << quantile(0.25) << '\n';
    std::cout << "median             = " << quantile(0.50) << '\n';
    std::cout << "75% quantile       = " << quantile(0.75) << '\n';
    std::cout << "90% quantile       = " << quantile(0.90) << '\n';
    std::cout << "95% quantile       = " << quantile(0.95) << '\n';
    std::cout << "99% quantile       = " << quantile(0.99) << '\n';
    std::cout << "maximum            = " << values.back() << '\n';
    std::cout << "untruncated mean   = " << mean << '\n';
    std::cout << "untruncated sigma  = " << std::sqrt(variance) << '\n';

    std::cout << "\nFractions\n";
    std::cout << "dchi2 < 0          = "
              << fractionBelow(0.0) << '\n';
    std::cout << "dchi2 > 5          = "
              << fractionAbove(5.0) << '\n';
    std::cout << "dchi2 > 12.83      = "
              << fractionAbove(12.83) << '\n';
    std::cout << "dchi2 > 25         = "
              << fractionAbove(25.0) << '\n';
    std::cout << "dchi2 > 50         = "
              << fractionAbove(50.0) << '\n';
    std::cout << "dchi2 > 100        = "
              << fractionAbove(100.0) << '\n';
    std::cout << "dchi2 > 400        = "
              << fractionAbove(400.0) << '\n';
}