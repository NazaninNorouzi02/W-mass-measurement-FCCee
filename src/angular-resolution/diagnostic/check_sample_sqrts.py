import ROOT

ROOT.gROOT.SetBatch(True)

# =========================================================
# Small helper functions
# =========================================================

ROOT.gInterpreter.Declare(r'''
#include "ROOT/RVec.hxx"
#include <cmath>
#include <algorithm>

using ROOT::VecOps::RVec;

float sumFloat(const RVec<float>& v) {
    float s = 0.0f;
    for (auto x : v) s += x;
    return s;
}

float invMassFromSums(float E, float px, float py, float pz) {
    float m2 = E*E - px*px - py*py - pz*pz;
    return std::sqrt(std::max(m2, 0.0f));
}
''')

# =========================================================
# FCCAnalyses configuration
# =========================================================

processList = {
    "p8_ee_WW_ecm160": {
        "fraction": 1.0e-5,
        "chunks": 1,
        "output": "check_sqrts_p8_ee_WW_ecm160"
    }
}

inputDir = "/eos/experiment/fcc/ee/generation/DelphesEvents/winter2023/IDEA/"
procDict = "FCCee_procDict_winter2023_IDEA.json"
outDir = "./outputs/check_sqrts"

nCPUS = 10
doTree = True

# =========================================================

class RDFanalysis:

    def analysers(df):
        df = (
            df
            .Define("RP_px", "ReconstructedParticle::get_px(ReconstructedParticles)")
            .Define("RP_py", "ReconstructedParticle::get_py(ReconstructedParticles)")
            .Define("RP_pz", "ReconstructedParticle::get_pz(ReconstructedParticles)")
            .Define("RP_e",  "ReconstructedParticle::get_e(ReconstructedParticles)")
            .Define("RP_m",  "ReconstructedParticle::get_mass(ReconstructedParticles)")

            # Full reconstructed visible system
            .Define("visibleE",  "sumFloat(RP_e)")
            .Define("visiblePx", "sumFloat(RP_px)")
            .Define("visiblePy", "sumFloat(RP_py)")
            .Define("visiblePz", "sumFloat(RP_pz)")
            .Define(
                "visibleP",
                "sqrt(visiblePx*visiblePx + visiblePy*visiblePy + visiblePz*visiblePz)"
            )
            .Define(
                "visibleM",
                "invMassFromSums(visibleE, visiblePx, visiblePy, visiblePz)"
            )
            .Define("nRecoParticles", "RP_e.size()")

            # Exclusive 4-jet clustering, same as your W-mass code
            .Define(
                "pseudo_jets",
                "JetClusteringUtils::set_pseudoJets_xyzm(RP_px, RP_py, RP_pz, RP_m)"
            )
            .Define(
                "Jets",
                "JetClustering::clustering_ee_kt(2, 4, 0, 0)(pseudo_jets)"
            )
            .Define("jets", "JetClusteringUtils::get_pseudoJets(Jets)")
            .Define("jet_px", "JetClusteringUtils::get_px(jets)")
            .Define("jet_py", "JetClusteringUtils::get_py(jets)")
            .Define("jet_pz", "JetClusteringUtils::get_pz(jets)")
            .Define("jet_e",  "JetClusteringUtils::get_e(jets)")
            .Define("nJets", "jet_e.size()")
            .Define("sumJetE", "sumFloat(jet_e)")
            .Define("sumJetPx", "sumFloat(jet_px)")
            .Define("sumJetPy", "sumFloat(jet_py)")
            .Define("sumJetPz", "sumFloat(jet_pz)")
            .Define(
                "sumJetP",
                "sqrt(sumJetPx*sumJetPx + sumJetPy*sumJetPy + sumJetPz*sumJetPz)"
            )
            .Define(
                "sumJetM",
                "invMassFromSums(sumJetE, sumJetPx, sumJetPy, sumJetPz)"
            )
        )

        return df

    def output():
        return [
            "visibleE", "visiblePx", "visiblePy", "visiblePz", "visibleP", "visibleM",
            "nRecoParticles",
            "nJets", "sumJetE", "sumJetPx", "sumJetPy", "sumJetPz", "sumJetP", "sumJetM",
        ]