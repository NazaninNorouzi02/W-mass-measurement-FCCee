#ifndef EEWW_ANALYSIS_GETXRECO_V2_H
#define EEWW_ANALYSIS_GETXRECO_V2_H

#include <cmath>
#include <limits>

ROOT::VecOps::RVec<float> getXReco_v2(
    const ROOT::VecOps::RVec<fastjet::PseudoJet>& jets_reco,
    const ROOT::VecOps::RVec<int>& jet_match_indices
) {
    ROOT::VecOps::RVec<float> x_reco;
    x_reco.reserve(jet_match_indices.size());

    for (size_t i = 0; i < jet_match_indices.size(); ++i) {

        const int gen_idx = jet_match_indices[i];

        if (gen_idx >= 0 && i < jets_reco.size()) {

            const double p = jets_reco[i].modp();
            const double m = jets_reco[i].m();

            // CHANGED:
            // x is now log(p/m), not log(p/E).
            // Invalid or massless jets return NaN.
            if (
                std::isfinite(p) &&
                std::isfinite(m) &&
                p > 0.0 &&
                m > 0.0
            ) {
                x_reco.push_back(static_cast<float>(std::log(p / m)));
            } else {
                x_reco.push_back(std::numeric_limits<float>::quiet_NaN());
            }

        } else {
            x_reco.push_back(std::numeric_limits<float>::quiet_NaN());
        }
    }

    return x_reco;
}

#endif // EEWW_ANALYSIS_GETXRECO_V2_H