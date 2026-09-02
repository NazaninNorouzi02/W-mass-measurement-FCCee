#ifndef EEWW_ANALYSIS_ENERGYDIAGNOSTICS_H
#define EEWW_ANALYSIS_ENERGYDIAGNOSTICS_H

#include <cmath>
#include <limits>
#include "ROOT/RVec.hxx"
#include "edm4hep/ReconstructedParticleData.h"
#include "edm4hep/MCParticleData.h"
#include "fastjet/PseudoJet.hh"

// ------------------------------------------------------------
// Safe scalar access
// ------------------------------------------------------------
float getElementOrNaN(const ROOT::VecOps::RVec<float>& values, int index) {
    if (index < 0) {
        return std::numeric_limits<float>::quiet_NaN();
    }

    if (static_cast<size_t>(index) >= values.size()) {
        return std::numeric_limits<float>::quiet_NaN();
    }

    return values[index];
}

// ------------------------------------------------------------
// MC energy helper
// edm4hep::MCParticleData does not necessarily contain p.energy.
// Therefore use E = sqrt(px^2 + py^2 + pz^2 + m^2).
// ------------------------------------------------------------
float getMCEnergy(const edm4hep::MCParticleData& p) {
    const float px = p.momentum.x;
    const float py = p.momentum.y;
    const float pz = p.momentum.z;
    const float m  = p.mass;

    const float e2 = px * px + py * py + pz * pz + m * m;

    if (e2 <= 0.0f) {
        return 0.0f;
    }

    return std::sqrt(e2);
}

// ------------------------------------------------------------
// Reconstructed-particle energy diagnostics
// ReconstructedParticleData has direct energy field.
// ------------------------------------------------------------
float sumRecoEnergy(
    const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& particles
) {
    float sum = 0.0f;

    for (const auto& p : particles) {
        sum += p.energy;
    }

    return sum;
}

float sumForwardPhotonEnergyReco(
    const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& photons,
    float eMin,
    float absCosThetaMin
) {
    float sum = 0.0f;

    for (const auto& p : photons) {
        const float px = p.momentum.x;
        const float py = p.momentum.y;
        const float pz = p.momentum.z;

        const float pAbs = std::sqrt(px * px + py * py + pz * pz);
        if (pAbs <= 0.0f) {
            continue;
        }

        const float absCosTheta = std::abs(pz / pAbs);

        if (p.energy > eMin && absCosTheta > absCosThetaMin) {
            sum += p.energy;
        }
    }

    return sum;
}

int countForwardPhotonsReco(
    const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& photons,
    float eMin,
    float absCosThetaMin
) {
    int count = 0;

    for (const auto& p : photons) {
        const float px = p.momentum.x;
        const float py = p.momentum.y;
        const float pz = p.momentum.z;

        const float pAbs = std::sqrt(px * px + py * py + pz * pz);
        if (pAbs <= 0.0f) {
            continue;
        }

        const float absCosTheta = std::abs(pz / pAbs);

        if (p.energy > eMin && absCosTheta > absCosThetaMin) {
            count++;
        }
    }

    return count;
}

// ------------------------------------------------------------
// MC photon diagnostics
// ------------------------------------------------------------
ROOT::VecOps::RVec<edm4hep::MCParticleData> selectMCFinalPhotons(
    const ROOT::VecOps::RVec<edm4hep::MCParticleData>& particles
) {
    ROOT::VecOps::RVec<edm4hep::MCParticleData> photons;

    for (const auto& p : particles) {
        if (p.PDG == 22 && p.generatorStatus == 1) {
            photons.push_back(p);
        }
    }

    return photons;
}

float sumMCEnergy(
    const ROOT::VecOps::RVec<edm4hep::MCParticleData>& particles
) {
    float sum = 0.0f;

    for (const auto& p : particles) {
        sum += getMCEnergy(p);
    }

    return sum;
}

float sumForwardPhotonEnergyMC(
    const ROOT::VecOps::RVec<edm4hep::MCParticleData>& photons,
    float eMin,
    float absCosThetaMin
) {
    float sum = 0.0f;

    for (const auto& p : photons) {
        const float px = p.momentum.x;
        const float py = p.momentum.y;
        const float pz = p.momentum.z;

        const float pAbs = std::sqrt(px * px + py * py + pz * pz);
        if (pAbs <= 0.0f) {
            continue;
        }

        const float absCosTheta = std::abs(pz / pAbs);
        const float energy = getMCEnergy(p);

        if (energy > eMin && absCosTheta > absCosThetaMin) {
            sum += energy;
        }
    }

    return sum;
}

int countForwardPhotonsMC(
    const ROOT::VecOps::RVec<edm4hep::MCParticleData>& photons,
    float eMin,
    float absCosThetaMin
) {
    int count = 0;

    for (const auto& p : photons) {
        const float px = p.momentum.x;
        const float py = p.momentum.y;
        const float pz = p.momentum.z;

        const float pAbs = std::sqrt(px * px + py * py + pz * pz);
        if (pAbs <= 0.0f) {
            continue;
        }

        const float absCosTheta = std::abs(pz / pAbs);
        const float energy = getMCEnergy(p);

        if (energy > eMin && absCosTheta > absCosThetaMin) {
            count++;
        }
    }

    return count;
}

// ------------------------------------------------------------
// FastJet diagnostics
// ------------------------------------------------------------
ROOT::VecOps::RVec<float> getPseudoJetEnergies(
    const ROOT::VecOps::RVec<fastjet::PseudoJet>& jets
) {
    ROOT::VecOps::RVec<float> energies;

    for (const auto& j : jets) {
        energies.push_back(j.E());
    }

    return energies;
}

float sumPseudoJetEnergy(
    const ROOT::VecOps::RVec<fastjet::PseudoJet>& jets
) {
    float sum = 0.0f;

    for (const auto& j : jets) {
        sum += j.E();
    }

    return sum;
}

#endif