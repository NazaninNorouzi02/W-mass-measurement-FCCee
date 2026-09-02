import ROOT, os
from glob import glob

ROOT.gROOT.SetBatch(True)

# ================================================
# following the thesis (equation 7.4):
#   alpha_i = E_rescaled_i / E_reco_i
#   where E_rescaled comes from solving the
#   4-momentum conservation linear system
#   with fixed jet velocities beta_i = p_i/E_i
#
# delta_alpha = alpha - 1 = (E_rescaled - E_reco)/E_reco
# Highly Biased result. Moving Back to Parton Matching.
# ================================================

fractions = 1e-6
inputDir  = "/eos/experiment/fcc/ee/generation/DelphesEvents/winter2023/IDEA/"

filter_config = {
    "delta_theta": {"min": -0.02, "max": 0.02},
    "delta_phi":   {"min": -0.02, "max": 0.02},
    "delta_eta":   {"min": -0.1,  "max": 0.1},
    "delta_x":     {"min": -0.6,  "max": 0.6},
    "delta_alpha": {"min": -0.25, "max": 0.25}
}

processList = {
    "p8_ee_WW_ecm160": {
        "fraction": fractions,
        "chunks": 4,
        "output": "angular_resolution_ecm160_v4"
    }
}

outputDir = "outputs/step1_angular_resolution_multiE/"
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
    "headers/filterValues.h",
    "headers/selectQuarks.h",
    # jetPartonMatching and getDeltaAlphaParton NOT included
 
]

# ------------------------------------------------
# alpha computation following thesis eq 7.3/7.4
# Solve: B * E_rescaled = [sqrt(s), 0, 0, 0]^T
# where B is the 4x4 matrix with rows:
#   row 0: [1, 1, 1, 1]                (energy)
#   row 1: [bx1, bx2, bx3, bx4]        (px/E)
#   row 2: [by1, by2, by3, by4]        (py/E)
#   row 3: [bz1, bz2, bz3, bz4]        (pz/E)
# alpha_i = E_rescaled_i / E_reco_i
# delta_alpha_i = alpha_i - 1
# ------------------------------------------------
ROOT.gInterpreter.Declare(r'''
#include <vector>
#include <cmath>
#include "ROOT/RVec.hxx"
#include "TMatrixD.h"
#include "TVectorD.h"

// Returns delta_alpha for each jet: (E_rescaled/E_reco) - 1
// Returns empty vector if system is singular or any E_rescaled <= 0
ROOT::VecOps::RVec<float>
computeDeltaAlpha(
    const ROOT::VecOps::RVec<fastjet::PseudoJet>& jets,
    float sqrts
){
    ROOT::VecOps::RVec<float> out;
    if(jets.size() != 4) return out;

    // Build 4x4 velocity matrix B
    // Column i = jet i, rows = [1, bx, by, bz]
    TMatrixD B(4, 4);
    TVectorD rhs(4);
    rhs[0] = sqrts;
    rhs[1] = 0.0;
    rhs[2] = 0.0;
    rhs[3] = 0.0;

    for(int i = 0; i < 4; i++){
        double E  = jets[i].E();
        if(E <= 0) return out;  // unphysical jet
        double px = jets[i].px();
        double py = jets[i].py();
        double pz = jets[i].pz();
        B[0][i] = 1.0;
        B[1][i] = px / E;   // bx
        B[2][i] = py / E;   // by
        B[3][i] = pz / E;   // bz
    }

    // Check determinant — discard singular events
    double det = 0.0;
    TMatrixD Bcopy(B);
    Bcopy.Invert(&det);
    if(std::abs(det) < 1e-10) return out;

    // Solve B * E_rescaled = rhs
    TVectorD E_rescaled = Bcopy * rhs;

    // Require all rescaled energies positive (thesis condition)
    for(int i = 0; i < 4; i++){
        if(E_rescaled[i] <= 0) return out;
    }

    // delta_alpha_i = E_rescaled_i / E_reco_i - 1
    for(int i = 0; i < 4; i++){
        double alpha = E_rescaled[i] / jets[i].E();
        out.push_back((float)(alpha - 1.0));
    }
    return out;
}
''')

# ------------------------------------------------
# RDFanalysis
# ------------------------------------------------
class RDFanalysis:
    @staticmethod
    def analysers(df):

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

        # parton selection 
        df = df.Define("partons_all", "selectQuarks(Particle)")
        df = df.Define("n_partons",   "(int)partons_all.size()")
        df = df.Filter("n_partons == 4", "Require exactly 4 partons")

        df = df.Define("parton_energies",
            "FCCAnalyses::MCParticle::get_e(partons_all)")
        df = df.Define("parton_eta",
            "FCCAnalyses::MCParticle::get_eta(partons_all)")
        df = df.Define("parton_phi",
            "FCCAnalyses::MCParticle::get_phi(partons_all)")
        df = df.Define("parton_y",
            "FCCAnalyses::MCParticle::get_y(partons_all)")

        # gen jets
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
            "FCCAnalyses::JetClusteringUtils::set_pseudoJets(Particle_px, Particle_py, Particle_pz, Particle_e)")
        df = df.Define("jets_gen_obj4",
            "JetClustering::clustering_ee_kt(2, 4, 0, 0)(pseudo_jets_gen)")
        df = df.Define("jets_gen4",
            "FCCAnalyses::JetClusteringUtils::get_pseudoJets(jets_gen_obj4)")

        # reco jets
        df = df.Define("Reco_px",
            "FCCAnalyses::ReconstructedParticle::get_px(reco_clean)")
        df = df.Define("Reco_py",
            "FCCAnalyses::ReconstructedParticle::get_py(reco_clean)")
        df = df.Define("Reco_pz",
            "FCCAnalyses::ReconstructedParticle::get_pz(reco_clean)")
        df = df.Define("Reco_e",
            "FCCAnalyses::ReconstructedParticle::get_e(reco_clean)")

        df = df.Define("pseudo_jets_reco",
            "FCCAnalyses::JetClusteringUtils::set_pseudoJets(Reco_px, Reco_py, Reco_pz, Reco_e)")
        df = df.Define("jets_reco_obj4",
            "JetClustering::clustering_ee_kt(2, 4, 0, 0)(pseudo_jets_reco)")
        df = df.Define("jets_reco4",
            "FCCAnalyses::JetClusteringUtils::get_pseudoJets(jets_reco_obj4)")

        df = df.Define("n_jets_gen",  "jets_gen4.size()")
        df = df.Define("n_jets_reco", "jets_reco4.size()")
        df = df.Filter("n_jets_gen == 4 && n_jets_reco == 4")

        # reco <-> gen matching 
        df = df.Define("jet_match_indices",
            "greedyJetMatching(jets_reco4, jets_gen4, 0.1)")
        df = df.Define("n_matched_jets",
            "countMatchedJets(jet_match_indices)")
        df = df.Filter("n_matched_jets == 4", "Require all 4 jets matched")

        # angular resolutions 
        df = df.Define("delta_theta_matched",
            "getDeltaTheta(jets_reco4, jets_gen4, jet_match_indices)")
        df = df.Define("delta_phi_matched",
            "getDeltaPhi(jets_reco4, jets_gen4, jet_match_indices)")
        df = df.Define("delta_eta_matched",
            "getDeltaEta(jets_reco4, jets_gen4, jet_match_indices)")
        df = df.Define("delta_mass_matched",
            "getDeltaMass(jets_reco4, jets_gen4, jet_match_indices)")
        df = df.Define("x_gen_matched",
            "getXGen(jets_gen4, jet_match_indices)")
        df = df.Define("x_reco_matched",
            "getXReco(jets_reco4, jet_match_indices)")
        df = df.Define("delta_x_matched",
            "x_reco_matched - x_gen_matched")


        # Compute Alpha
        df = df.Define("delta_alpha",
            "computeDeltaAlpha(jets_reco4, 160.0f)")

        # discard events where system was singular or
        # had unphysical (negative) rescaled energies
        df = df.Filter("delta_alpha.size() == 4",
            "Require valid kinematic rescaling solution")

        # individual jet variables
        for jet_idx in range(1, 5):
            idx = jet_idx - 1
            df = df.Define(f"delta_alpha_j{jet_idx}",
                f"getElement(delta_alpha, {idx})")
            df = df.Define(f"delta_theta_j{jet_idx}",
                f"getElement(delta_theta_matched, {idx})")
            df = df.Define(f"delta_phi_j{jet_idx}",
                f"getElement(delta_phi_matched, {idx})")
            df = df.Define(f"delta_eta_j{jet_idx}",
                f"getElement(delta_eta_matched, {idx})")
            df = df.Define(f"delta_x_j{jet_idx}",
                f"getElement(delta_x_matched, {idx})")

        df = df.Define("theta_gen_all",
            "FCCAnalyses::JetClusteringUtils::get_theta(jets_gen4)")
        df = df.Define("theta_reco_all",
            "FCCAnalyses::JetClusteringUtils::get_theta(jets_reco4)")
        df = df.Define("phi_gen_all",
            "FCCAnalyses::JetClusteringUtils::get_phi(jets_gen4)")
        df = df.Define("phi_reco_all",
            "FCCAnalyses::JetClusteringUtils::get_phi(jets_reco4)")

        # filtered versions
        df = df.Define("filtered_delta_theta_matched",
            f"filterValues(delta_theta_matched, {filter_config['delta_theta']['min']}, {filter_config['delta_theta']['max']})")
        df = df.Define("filtered_delta_phi_matched",
            f"filterValues(delta_phi_matched, {filter_config['delta_phi']['min']}, {filter_config['delta_phi']['max']})")
        df = df.Define("filtered_delta_eta_matched",
            f"filterValues(delta_eta_matched, {filter_config['delta_eta']['min']}, {filter_config['delta_eta']['max']})")
        df = df.Define("filtered_delta_x_matched",
            f"filterValues(delta_x_matched, {filter_config['delta_x']['min']}, {filter_config['delta_x']['max']})")
        df = df.Define("filtered_delta_alpha",
            f"filterValues(delta_alpha, {filter_config['delta_alpha']['min']}, {filter_config['delta_alpha']['max']})")

        for jet_idx in range(1, 5):
            df = df.Define(f"filtered_delta_theta_j{jet_idx}",
                f"(delta_theta_j{jet_idx} >= {filter_config['delta_theta']['min']} && delta_theta_j{jet_idx} <= {filter_config['delta_theta']['max']}) ? delta_theta_j{jet_idx} : -999.0f")
            df = df.Define(f"filtered_delta_phi_j{jet_idx}",
                f"(delta_phi_j{jet_idx} >= {filter_config['delta_phi']['min']} && delta_phi_j{jet_idx} <= {filter_config['delta_phi']['max']}) ? delta_phi_j{jet_idx} : -999.0f")
            df = df.Define(f"filtered_delta_x_j{jet_idx}",
                f"(delta_x_j{jet_idx} >= {filter_config['delta_x']['min']} && delta_x_j{jet_idx} <= {filter_config['delta_x']['max']}) ? delta_x_j{jet_idx} : -999.0f")
            df = df.Define(f"filtered_delta_alpha_j{jet_idx}",
                f"(delta_alpha_j{jet_idx} >= {filter_config['delta_alpha']['min']} && delta_alpha_j{jet_idx} <= {filter_config['delta_alpha']['max']}) ? delta_alpha_j{jet_idx} : -999.0f")

        return df

    @staticmethod
    def output():
        return [
            "delta_theta_j1", "delta_theta_j2", "delta_theta_j3", "delta_theta_j4",
            "delta_phi_j1",   "delta_phi_j2",   "delta_phi_j3",   "delta_phi_j4",
            "delta_eta_j1",   "delta_eta_j2",   "delta_eta_j3",   "delta_eta_j4",
            "delta_x_j1",     "delta_x_j2",     "delta_x_j3",     "delta_x_j4",
            "delta_alpha_j1", "delta_alpha_j2", "delta_alpha_j3", "delta_alpha_j4",

            "filtered_delta_theta_matched", "filtered_delta_phi_matched",
            "filtered_delta_eta_matched",   "filtered_delta_x_matched",
            "filtered_delta_alpha",

            "filtered_delta_theta_j1", "filtered_delta_theta_j2",
            "filtered_delta_theta_j3", "filtered_delta_theta_j4",
            "filtered_delta_phi_j1",   "filtered_delta_phi_j2",
            "filtered_delta_phi_j3",   "filtered_delta_phi_j4",
            "filtered_delta_x_j1",     "filtered_delta_x_j2",
            "filtered_delta_x_j3",     "filtered_delta_x_j4",
            "filtered_delta_alpha_j1", "filtered_delta_alpha_j2",
            "filtered_delta_alpha_j3", "filtered_delta_alpha_j4",

            "parton_energies",
            "parton_eta",
            "parton_phi",
            "parton_y",

            "n_partons",
            "n_matched_jets",
        ]
