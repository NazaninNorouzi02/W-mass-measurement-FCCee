//
// Created by kid-a on 2/12/26.
//

#ifndef EEWW_ANALYSIS_SELECTQUARKS_H
#define EEWW_ANALYSIS_SELECTQUARKS_H
#include <cmath>
#include <set>
#include "ROOT/RVec.hxx"
#include "edm4hep/MCParticleData.h"

ROOT::VecOps::RVec<edm4hep::MCParticleData> selectQuarks(
    const ROOT::VecOps::RVec<edm4hep::MCParticleData>& particles
) {
    ROOT::VecOps::RVec<edm4hep::MCParticleData> result;

    static const std::set<int> quark_pdgs = {1, -1, 2, -2, 3, -3, 4, -4, 5, -5, 6, -6};

    for (const auto& p : particles) {
        // Status 23: outgoing from hard subprocess in Pythia 8.
        // For ee->WW->qqqq this gives exactly the 4 hard-process quarks,
        // before parton shower — the correct truth reference for jet matching.
        if (quark_pdgs.count(p.PDG) > 0 && p.generatorStatus == 23) {
            result.push_back(p);
        }
    }

    return result;
}
#endif //EEWW_ANALYSIS_SELECTQUARKS_H
