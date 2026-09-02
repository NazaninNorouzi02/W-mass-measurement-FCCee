import ROOT, os
from glob import glob

ROOT.gROOT.SetBatch(True)

# =============================================================================
# matching_comparison.py
#
# Runs THREE parton-matching strategies on the same events and stores
# delta_alpha branches for each, so they can be compared directly.
#
# STRATEGY A — Original (non-exclusive):
#   jetPartonMatching: each reco-jet independently finds its nearest parton.
#   Two jets CAN claim the same parton. No filter on match quality.
#
# STRATEGY B — Unified greedy (exclusive, direct):
#   greedyJetMatching called with partons promoted to PseudoJet.
#   Same algorithm as reco<->gen matching. Strictly 1-to-1.
#   Filter: n_matched_partons_B == 4. Not implemented on this code for consistency. 
#   Although I believe this filter should be applied to resolve the issue of Strategy A.
#
# STRATEGY C — Transitive greedy (exclusive, two-step):
#   Step 1: reco -> gen  via greedyJetMatching (dR<0.1)
#   Step 2: gen  -> parton via greedyJetMatching (dR<0.4)
#   Compose the two maps: reco -> parton through gen.
#   Ensures the parton assigned to a reco-jet is the one whose shower
#   products are actually inside that gen-jet.
#
# EVENT TOPOLOGY VARIABLES (same for all strategies):
#   event_E_balance = sum(E_reco_jets) / sum(E_partons)
#   min_dr_partons  = minimum dR between any two partons
#   min_dr_jets     = minimum dR between any two reco-jets
#   softest_jet_E, hardest_jet_E
#
# Output: one ROOT file with all branches.
# Plot with: python3 plot_matching_comparison.py
# =============================================================================

# ── inline C++ ────────────────────────────────────────────────────────────────
ROOT.gInterpreter.Declare(r"""
#include <ROOT/RVec.hxx>
#include <fastjet/PseudoJet.hh>
#include <cmath>
#include <vector>
#include <algorithm>
#include <set>

// ── helpers ──────────────────────────────────────────────────────────────────

inline double dR_safe(double eta1, double phi1, double eta2, double phi2) {
    double deta = eta1 - eta2;
    double dphi = phi1 - phi2;
    while (dphi >  M_PI) dphi -= 2.0 * M_PI;
    while (dphi < -M_PI) dphi += 2.0 * M_PI;
    return std::sqrt(deta*deta + dphi*dphi);
}

// ── STRATEGY B: exclusive greedy matching of reco-jets to PseudoJet partons ─
// Same algorithm as greedyJetMatching in the header but written inline
// so it is self-contained and does not depend on the header's Pair struct.
ROOT::VecOps::RVec<int> greedyMatchToPseudoJets(
    const ROOT::VecOps::RVec<fastjet::PseudoJet>& jets_reco,
    const ROOT::VecOps::RVec<fastjet::PseudoJet>& targets,
    double dRThreshold = 0.4)
{
    struct P { int ri, ti; double dr; };
    std::vector<P> pairs;
    pairs.reserve(jets_reco.size() * targets.size());
    for (int i = 0; i < (int)jets_reco.size(); ++i)
        for (int j = 0; j < (int)targets.size(); ++j)
            pairs.push_back({i, j,
                dR_safe(jets_reco[i].eta(), jets_reco[i].phi(),
                        targets[j].eta(),  targets[j].phi())});
    std::sort(pairs.begin(), pairs.end(),
              [](const P& a, const P& b){ return a.dr < b.dr; });
    std::vector<bool> ru(jets_reco.size(), false), tu(targets.size(), false);
    ROOT::VecOps::RVec<int> result(jets_reco.size(), -1);
    for (const auto& p : pairs) {
        if (p.dr > dRThreshold) break;
        if (ru[p.ri] || tu[p.ti]) continue;
        result[p.ri] = p.ti;
        ru[p.ri] = tu[p.ti] = true;
    }
    return result;
}

// ── STRATEGY C helper: greedy match using eta/phi float arrays ────────────
ROOT::VecOps::RVec<int> greedyMatchEtaPhi(
    const ROOT::VecOps::RVec<float>& src_eta,
    const ROOT::VecOps::RVec<float>& src_phi,
    const ROOT::VecOps::RVec<float>& tgt_eta,
    const ROOT::VecOps::RVec<float>& tgt_phi,
    double dRThreshold = 0.4)
{
    struct P { int si, ti; double dr; };
    std::vector<P> pairs;
    for (int i = 0; i < (int)src_eta.size(); ++i)
        for (int j = 0; j < (int)tgt_eta.size(); ++j)
            pairs.push_back({i, j,
                dR_safe(src_eta[i], src_phi[i], tgt_eta[j], tgt_phi[j])});
    std::sort(pairs.begin(), pairs.end(),
              [](const P& a, const P& b){ return a.dr < b.dr; });
    std::vector<bool> su(src_eta.size(), false), tu(tgt_eta.size(), false);
    ROOT::VecOps::RVec<int> result(src_eta.size(), -1);
    for (const auto& p : pairs) {
        if (p.dr > dRThreshold) break;
        if (su[p.si] || tu[p.ti]) continue;
        result[p.si] = p.ti;
        su[p.si] = tu[p.ti] = true;
    }
    return result;
}

// ── STRATEGY C: compose reco->gen and gen->parton index maps ────────────────
ROOT::VecOps::RVec<int> transitiveMatch(
    const ROOT::VecOps::RVec<int>& reco_to_gen,
    const ROOT::VecOps::RVec<int>& gen_to_parton)
{
    ROOT::VecOps::RVec<int> result(reco_to_gen.size(), -1);
    for (int i = 0; i < (int)reco_to_gen.size(); ++i) {
        int g = reco_to_gen[i];
        if (g >= 0 && g < (int)gen_to_parton.size())
            result[i] = gen_to_parton[g];
    }
    return result;
}

// ── delta_alpha from any match index vector ──────────────────────────────────
// Always outputs a vector of jets_reco.size() elements.
// Unmatched or invalid -> -999 sentinel (excluded by all filters downstream).
ROOT::VecOps::RVec<float> computeDeltaAlpha(
    const ROOT::VecOps::RVec<fastjet::PseudoJet>& jets_reco,
    const ROOT::VecOps::RVec<float>& parton_energies,
    const ROOT::VecOps::RVec<int>& match_indices)
{
    ROOT::VecOps::RVec<float> out;
    out.reserve(jets_reco.size());
    for (int i = 0; i < (int)jets_reco.size(); ++i) {
        int j = (i < (int)match_indices.size()) ? match_indices[i] : -1;
        if (j >= 0 && j < (int)parton_energies.size() && parton_energies[j] > 0)
            out.push_back((jets_reco[i].E() - parton_energies[j]) / parton_energies[j]);
        else
            out.push_back(-999.0f);
    }
    return out;
}

// ── count valid (>=0) entries ────────────────────────────────────────────────
int countValid(const ROOT::VecOps::RVec<int>& v) {
    int n = 0; for (auto x : v) if (x >= 0) n++; return n;
}

// ── topology helpers ─────────────────────────────────────────────────────────
float minDR_pseudojets(const ROOT::VecOps::RVec<fastjet::PseudoJet>& jets) {
    float m = 999.0f;
    for (int i = 0; i < (int)jets.size(); ++i)
        for (int j = i+1; j < (int)jets.size(); ++j) {
            float dr = dR_safe(jets[i].eta(), jets[i].phi(),
                               jets[j].eta(), jets[j].phi());
            if (dr < m) m = dr;
        }
    return m;
}

float minDR_arrays(const ROOT::VecOps::RVec<float>& eta,
                   const ROOT::VecOps::RVec<float>& phi) {
    float m = 999.0f;
    for (int i = 0; i < (int)eta.size(); ++i)
        for (int j = i+1; j < (int)eta.size(); ++j) {
            float dr = dR_safe(eta[i], phi[i], eta[j], phi[j]);
            if (dr < m) m = dr;
        }
    return m;
}
""")

# ── configuration ─────────────────────────────────────────────────────────────
fractions = 1e-6

inputDir = "/eos/experiment/fcc/ee/generation/DelphesEvents/winter2023/IDEA/"

processList = {
    "p8_ee_WW_ecm160": {
        "fraction": fractions,
        "chunks": 4,
        "output": "matching_comparison_ecm160"
    }
}

outputDir = "outputs/matching_comparison/"
procDict  = "FCCee_procDict_winter2023_IDEA.json"
nCPUS     = 4
doTree    = True

includePaths = [
    "headers/greedyJetMatching.h",
    "headers/getDeltaTheta.h",
    "headers/getDeltaPhi.h",
    "headers/getDeltaEta.h",
    "headers/getDeltaMass.h",
    "headers/getXGen.h",
    "headers/getXReco.h",
    "headers/getElement.h",
    "headers/jetPartonMatching.h",   # original non-exclusive — Strategy A
    "headers/getDeltaAlphaParton.h",
    "headers/filterValues.h",
    "headers/selectQuarks.h"
]


class RDFanalysis:
    @staticmethod
    def analysers(df):

        # ── lepton removal ───────────────────────────────────────────────────
        df = df.Alias("Electron0", "Electron#0.index")
        df = df.Alias("Muon0",     "Muon#0.index")
        df = df.Alias("Photon0",   "Photon#0.index")
        df = df.Define("ele_all",
            "FCCAnalyses::ReconstructedParticle::get(Electron0, ReconstructedParticles)")
        df = df.Define("mu_all",
            "FCCAnalyses::ReconstructedParticle::get(Muon0, ReconstructedParticles)")
        df = df.Define("pho_all",
            "FCCAnalyses::ReconstructedParticle::get(Photon0, ReconstructedParticles)")
        df = df.Define("RP_noPho",
            "FCCAnalyses::ReconstructedParticle::remove(ReconstructedParticles, pho_all)")
        df = df.Define("RP_noEle",
            "FCCAnalyses::ReconstructedParticle::remove(RP_noPho, ele_all)")
        df = df.Define("reco_clean",
            "FCCAnalyses::ReconstructedParticle::remove(RP_noEle, mu_all)")

        # ── partons ──────────────────────────────────────────────────────────
        df = df.Define("partons_all", "selectQuarks(Particle)")
        df = df.Define("n_partons",   "(int)partons_all.size()")
        df = df.Filter("n_partons == 4", "Require exactly 4 partons")

        df = df.Define("parton_energies",
            "FCCAnalyses::MCParticle::get_e(partons_all)")
        df = df.Define("parton_eta",
            "FCCAnalyses::MCParticle::get_eta(partons_all)")
        df = df.Define("parton_phi",
            "FCCAnalyses::MCParticle::get_phi(partons_all)")

        # partons as PseudoJets for Strategy B
        df = df.Define("parton_px",
            "FCCAnalyses::MCParticle::get_px(partons_all)")
        df = df.Define("parton_py",
            "FCCAnalyses::MCParticle::get_py(partons_all)")
        df = df.Define("parton_pz",
            "FCCAnalyses::MCParticle::get_pz(partons_all)")
        df = df.Define("pseudo_jets_parton",
            "FCCAnalyses::JetClusteringUtils::set_pseudoJets("
            "parton_px, parton_py, parton_pz, parton_energies)")

        # ── gen jets ─────────────────────────────────────────────────────────
        df = df.Define("MC_final",
            "FCCAnalyses::MCParticle::sel_genStatus(1)(Particle)")
        df = df.Define("Particle_px",
            "FCCAnalyses::MCParticle::get_px(MC_final)")
        df = df.Define("Particle_py",
            "FCCAnalyses::MCParticle::get_py(MC_final)")
        df = df.Define("Particle_pz",
            "FCCAnalyses::MCParticle::get_pz(MC_final)")
        df = df.Define("Particle_e",
            "FCCAnalyses::MCParticle::get_e(MC_final)")
        df = df.Define("pseudo_jets_gen",
            "FCCAnalyses::JetClusteringUtils::set_pseudoJets("
            "Particle_px, Particle_py, Particle_pz, Particle_e)")
        df = df.Define("jets_gen_obj4",
            "JetClustering::clustering_ee_kt(2, 4, 0, 0)(pseudo_jets_gen)")
        df = df.Define("jets_gen4",
            "FCCAnalyses::JetClusteringUtils::get_pseudoJets(jets_gen_obj4)")

        # ── reco jets ────────────────────────────────────────────────────────
        df = df.Define("Reco_px",
            "FCCAnalyses::ReconstructedParticle::get_px(reco_clean)")
        df = df.Define("Reco_py",
            "FCCAnalyses::ReconstructedParticle::get_py(reco_clean)")
        df = df.Define("Reco_pz",
            "FCCAnalyses::ReconstructedParticle::get_pz(reco_clean)")
        df = df.Define("Reco_e",
            "FCCAnalyses::ReconstructedParticle::get_e(reco_clean)")
        df = df.Define("pseudo_jets_reco",
            "FCCAnalyses::JetClusteringUtils::set_pseudoJets("
            "Reco_px, Reco_py, Reco_pz, Reco_e)")
        df = df.Define("jets_reco_obj4",
            "JetClustering::clustering_ee_kt(2, 4, 0, 0)(pseudo_jets_reco)")
        df = df.Define("jets_reco4",
            "FCCAnalyses::JetClusteringUtils::get_pseudoJets(jets_reco_obj4)")

        df = df.Define("n_jets_gen",  "(int)jets_gen4.size()")
        df = df.Define("n_jets_reco", "(int)jets_reco4.size()")
        df = df.Filter("n_jets_gen == 4 && n_jets_reco == 4")

        # ── reco <-> gen matching (shared by all strategies) ─────────────────
        df = df.Define("jet_match_indices",
            "greedyJetMatching(jets_reco4, jets_gen4, 0.1)")
        df = df.Define("n_matched_jets",
            "countMatchedJets(jet_match_indices)")
        df = df.Filter("n_matched_jets == 4", "Require all 4 reco-gen jets matched")

        # angular resolutions (same for all strategies — use reco<->gen match)
        df = df.Define("delta_theta_matched",
            "getDeltaTheta(jets_reco4, jets_gen4, jet_match_indices)")
        df = df.Define("delta_phi_matched",
            "getDeltaPhi(jets_reco4, jets_gen4, jet_match_indices)")
        df = df.Define("delta_eta_matched",
            "getDeltaEta(jets_reco4, jets_gen4, jet_match_indices)")
        df = df.Define("x_gen_matched",
            "getXGen(jets_gen4, jet_match_indices)")
        df = df.Define("x_reco_matched",
            "getXReco(jets_reco4, jet_match_indices)")
        df = df.Define("delta_x_matched",
            "x_reco_matched - x_gen_matched")

        for jet_idx in range(1, 5):
            idx = jet_idx - 1
            df = df.Define(f"delta_theta_j{jet_idx}",
                f"getElement(delta_theta_matched, {idx})")
            df = df.Define(f"delta_phi_j{jet_idx}",
                f"getElement(delta_phi_matched, {idx})")
            df = df.Define(f"delta_x_j{jet_idx}",
                f"getElement(delta_x_matched, {idx})")

        # ── gen-jet eta/phi arrays for Strategy C ────────────────────────────
        df = df.Define("gen_eta",
            "FCCAnalyses::JetClusteringUtils::get_eta(jets_gen4)")
        df = df.Define("gen_phi",
            "FCCAnalyses::JetClusteringUtils::get_phi(jets_gen4)")
        df = df.Define("reco_eta",
            "FCCAnalyses::JetClusteringUtils::get_eta(jets_reco4)")
        df = df.Define("reco_phi",
            "FCCAnalyses::JetClusteringUtils::get_phi(jets_reco4)")

        # ════════════════════════════════════════════════════════════════════
        # STRATEGY A — original non-exclusive matching (jetPartonMatching.h)
        # Each reco-jet independently finds nearest parton. No uniqueness.
        # ════════════════════════════════════════════════════════════════════
        df = df.Define("match_A",
            "jetPartonMatching(jets_reco4, parton_eta, parton_phi)")
        df = df.Define("n_matched_A", "countValid(match_A)")
        df = df.Define("delta_alpha_A",
            "computeDeltaAlpha(jets_reco4, parton_energies, match_A)")
        for jet_idx in range(1, 5):
            df = df.Define(f"delta_alpha_A_j{jet_idx}",
                f"getElement(delta_alpha_A, {jet_idx-1})")

        # ════════════════════════════════════════════════════════════════════
        # STRATEGY B — unified greedy, direct reco->parton
        # Partons promoted to PseudoJet. Same greedyJetMatching algorithm
        # as reco<->gen. Strictly 1-to-1.
        # ════════════════════════════════════════════════════════════════════
        df = df.Define("match_B",
            "greedyMatchToPseudoJets(jets_reco4, pseudo_jets_parton, 0.4)")
        df = df.Define("n_matched_B", "countValid(match_B)")
        df = df.Define("delta_alpha_B",
            "computeDeltaAlpha(jets_reco4, parton_energies, match_B)")
        for jet_idx in range(1, 5):
            df = df.Define(f"delta_alpha_B_j{jet_idx}",
                f"getElement(delta_alpha_B, {jet_idx-1})")

        # ════════════════════════════════════════════════════════════════════
        # STRATEGY C — transitive greedy: reco -> gen -> parton
        # Step 1: reco->gen (already done above as jet_match_indices)
        # Step 2: gen->parton via greedy on eta/phi arrays
        # Compose: reco's parton = parton of its matched gen-jet
        # ════════════════════════════════════════════════════════════════════
        df = df.Define("gen_to_parton_C",
            "greedyMatchEtaPhi(gen_eta, gen_phi, parton_eta, parton_phi, 0.4)")
        df = df.Define("match_C",
            "transitiveMatch(jet_match_indices, gen_to_parton_C)")
        df = df.Define("n_matched_C", "countValid(match_C)")
        df = df.Define("delta_alpha_C",
            "computeDeltaAlpha(jets_reco4, parton_energies, match_C)")
        for jet_idx in range(1, 5):
            df = df.Define(f"delta_alpha_C_j{jet_idx}",
                f"getElement(delta_alpha_C, {jet_idx-1})")

        # ── event topology variables ─────────────────────────────────────────
        df = df.Define("reco_jet_energies",
            """ROOT::VecOps::RVec<float> e;
               for (auto& j : jets_reco4) e.push_back(j.E());
               return e;""")
        df = df.Define("total_reco_E",
            "ROOT::VecOps::Sum(reco_jet_energies)")
        df = df.Define("total_parton_E",
            "ROOT::VecOps::Sum(parton_energies)")
        df = df.Define("event_E_balance",
            "total_reco_E / total_parton_E")
        df = df.Define("softest_jet_E",
            "ROOT::VecOps::Min(reco_jet_energies)")
        df = df.Define("hardest_jet_E",
            "ROOT::VecOps::Max(reco_jet_energies)")
        df = df.Define("min_dr_jets",
            "minDR_pseudojets(jets_reco4)")
        df = df.Define("min_dr_partons",
            "minDR_arrays(parton_eta, parton_phi)")

        return df

    @staticmethod
    def output():
        branches = []

        # angular resolutions (strategy-independent)
        for v in ["delta_theta", "delta_phi", "delta_x"]:
            for j in range(1, 5):
                branches.append(f"{v}_j{j}")

        # match quality counters for all three strategies
        branches += ["n_matched_A", "n_matched_B", "n_matched_C"]

        # delta_alpha per strategy per jet
        for s in ["A", "B", "C"]:
            for j in range(1, 5):
                branches.append(f"delta_alpha_{s}_j{j}")

        # event topology
        branches += [
            "event_E_balance",
            "softest_jet_E",
            "hardest_jet_E",
            "min_dr_jets",
            "min_dr_partons",
            "total_reco_E",
            "total_parton_E",
            "parton_energies",
        ]

        return branches