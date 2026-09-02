#ifndef EEWW_ANALYSIS_SELECTISRPHOTONS_H
#define EEWW_ANALYSIS_SELECTISRPHOTONS_H

#include <cmath>

#include "ROOT/RVec.hxx"
#include "edm4hep/ReconstructedParticleData.h"

ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> selectISRPhotons(
    const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& photons,
    float min_energy,
    float min_abs_cos_theta
) {
    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> result;

    for (const auto& photon : photons) {
        const float px = photon.momentum.x;
        const float py = photon.momentum.y;
        const float pz = photon.momentum.z;

        const float momentum =
            std::sqrt(px * px + py * py + pz * pz);

        if (momentum <= 0.0f) {
            continue;
        }

        const float abs_cos_theta =
            std::abs(pz / momentum);

        if (
            photon.energy > min_energy &&
            abs_cos_theta > min_abs_cos_theta
        ) {
            result.push_back(photon);
        }
    }

    return result;
}

#endif