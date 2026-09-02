//
// getDeltaAlphaPartonFixed.h
//
// Fix vs original getDeltaAlphaParton.h:
//   Original silently skips push_back for unmatched jets (parton_idx < 0),
//   producing a vector shorter than jets_reco.size(). getElement then reads
//   out of bounds and returns garbage (-1000 or 0), corrupting all plots.
//
//   Fixed: always push_back for every jet. Unmatched or invalid -> -999.0f
//   sentinel. Output vector always has exactly jets_reco.size() elements.
//   Plot with cut variable > -998 to exclude sentinels without losing events.
//
#ifndef EEWW_ANALYSIS_GETDELTAALPHAPARTON_FIXED_H
#define EEWW_ANALYSIS_GETDELTAALPHAPARTON_FIXED_H

#include "ROOT/RVec.hxx"

ROOT::VecOps::RVec<float> getDeltaAlphaPartonFixed(
    const ROOT::VecOps::RVec<fastjet::PseudoJet>& jets_reco,
    const ROOT::VecOps::RVec<float>& parton_energies,
    const ROOT::VecOps::RVec<int>& jet_match_indices
) {
    ROOT::VecOps::RVec<float> delta_alpha;
    delta_alpha.reserve(jets_reco.size());

    for (size_t i = 0; i < jets_reco.size(); ++i) {
        int parton_idx = (i < jet_match_indices.size()) ? jet_match_indices[i] : -1;

        if (parton_idx >= 0 && parton_idx < (int)parton_energies.size()) {
            float E_parton = parton_energies[parton_idx];
            if (E_parton > 0.0f) {
                float E_reco = jets_reco[i].E();
                delta_alpha.push_back((E_reco - E_parton) / E_parton);
                continue;
            }
        }
        // unmatched or zero-energy parton -> sentinel
        delta_alpha.push_back(-999.0f);
    }

    return delta_alpha;
}

#endif // EEWW_ANALYSIS_GETDELTAALPHAPARTON_FIXED_H