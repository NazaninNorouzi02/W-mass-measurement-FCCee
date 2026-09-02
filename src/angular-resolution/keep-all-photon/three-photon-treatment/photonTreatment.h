#ifndef EEWW_ANALYSIS_PHOTONTREATMENT_H
#define EEWW_ANALYSIS_PHOTONTREATMENT_H

#include <cmath>

#include "ROOT/RVec.hxx"
#include "edm4hep/ReconstructedParticleData.h"
#include "fastjet/PseudoJet.hh"

// Select reconstructed photons that are considered ISR-like.
//
// Current candidate definition:
//   E_gamma > eMin
//   |cos(theta_gamma)| > absCosThetaMin
//
// This is a kinematic ISR candidate selection.
// It does not inspect generator ancestry.
inline ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>
selectISRLikePhotonsReco(
    const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& photons,
    float eMin,
    float absCosThetaMin
) {
    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> selected;
    selected.reserve(photons.size());

    for (const auto& photon : photons) {
        const float px = photon.momentum.x;
        const float py = photon.momentum.y;
        const float pz = photon.momentum.z;

        const float momentum =
            std::sqrt(px * px + py * py + pz * pz);

        if (momentum <= 0.0f) {
            continue;
        }

        const float absCosTheta =
            std::abs(pz / momentum);

        if (
            photon.energy > eMin &&
            absCosTheta > absCosThetaMin
        ) {
            selected.push_back(photon);
        }
    }

    return selected;
}

// Sum the energy of a clustered FastJet jet collection.
inline float sumJetEnergy(
    const ROOT::VecOps::RVec<fastjet::PseudoJet>& jets
) {
    float totalEnergy = 0.0f;

    for (const auto& jet : jets) {
        totalEnergy += jet.E();
    }

    return totalEnergy;
}

#endif // EEWW_ANALYSIS_PHOTONTREATMENT_H
