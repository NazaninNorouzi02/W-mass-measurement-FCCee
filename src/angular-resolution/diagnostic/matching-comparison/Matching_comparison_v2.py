import ROOT, os
from glob import glob

ROOT.gROOT.SetBatch(True)

# ================================================
# CONFIGURATION: 160 GeV only
# ================================================

fractions = 1e-4

inputDir = "/eos/experiment/fcc/ee/generation/DelphesEvents/winter2023/IDEA/"

processList = {
    "p8_ee_WW_ecm160": {
        "fraction": fractions,
        "chunks": 4,
        "output": "angular_resolution_matching_test_ecm160"
    }
}

filter_config = {
    "delta_theta": {"min": -0.02, "max": 0.02},
    "delta_phi":   {"min": -0.02, "max": 0.02},
    "delta_eta":   {"min": -0.1,  "max": 0.1},
    "delta_x":     {"min": -0.6,  "max": 0.6},
    "delta_alpha": {"min": -0.25, "max": 0.25}
}

outputDir = "outputs/step1_matching_test_ecm160/"
procDict = "FCCee_procDict_winter2023_IDEA.json"
nCPUS = 10
doTree = True

includePaths = [
    "headers/greedyJetMatching.h",
    "headers/getDeltaTheta.h",
    "headers/getDeltaPhi.h",
    "headers/getDeltaEta.h",
    "headers/getDeltaMass.h",
    "headers/getXGen.h",
    "headers/getXReco.h",
    "headers/getElement.h",
    "headers/jetPartonMatching.h",
    "headers/getDeltaAlphaParton.h",
    "headers/filterValues.h",
    "headers/selectQuarks.h",
    "headers/matchingDiagnostics.h"
]

# ------------------------------------------------
# RDFanalysis
# ------------------------------------------------
class RDFanalysis:
    @staticmethod
    def analysers(df):
        # ---------------------------
        # Remove identified leptons
        # ---------------------------
        df = df.Alias("Electron0", "Electron#0.index")
        df = df.Alias("Muon0", "Muon#0.index")
        df = df.Alias("Photon0", "Photon#0.index")
        df = df.Define("ele_all", "FCCAnalyses::ReconstructedParticle::get(Electron0, ReconstructedParticles)")
        df = df.Define("mu_all", "FCCAnalyses::ReconstructedParticle::get(Muon0, ReconstructedParticles)")
        df = df.Define("pho_all", "FCCAnalyses::ReconstructedParticle::get(Photon0, ReconstructedParticles)")
        df = df.Define("RP_noPho", "FCCAnalyses::ReconstructedParticle::remove(ReconstructedParticles, pho_all)")
        df = df.Define("RP_noEle", "FCCAnalyses::ReconstructedParticle::remove(RP_noPho, ele_all)")
        df = df.Define("reco_clean", "FCCAnalyses::ReconstructedParticle::remove(RP_noEle, mu_all)")

        # ------------------------------------------------
        # PARTON GENERATION
        # ------------------------------------------------
        df = df.Define("partons_all", "selectQuarks(Particle)")
        df = df.Define("n_partons", "partons_all.size()")
        df = df.Filter("n_partons == 4", "Require exactly 4 quarks")

        df = df.Define("parton_px", "FCCAnalyses::MCParticle::get_px(partons_all)")
        df = df.Define("parton_py", "FCCAnalyses::MCParticle::get_py(partons_all)")
        df = df.Define("parton_pz", "FCCAnalyses::MCParticle::get_pz(partons_all)")
        df = df.Define("parton_energies", "FCCAnalyses::MCParticle::get_e(partons_all)")
        df = df.Define("parton_eta", "FCCAnalyses::MCParticle::get_eta(partons_all)")
        df = df.Define("parton_phi", "FCCAnalyses::MCParticle::get_phi(partons_all)")
        df = df.Define("parton_y", "FCCAnalyses::MCParticle::get_y(partons_all)")

        # Promote partons to pseudo-jets for Strategy B
        df = df.Define(
            "pseudo_jets_partons",
            "FCCAnalyses::JetClusteringUtils::set_pseudoJets(parton_px, parton_py, parton_pz, parton_energies)"
        )

        # ------------------------------------------------
        # GEN JETS
        # ------------------------------------------------
        df = df.Define("MC_final", "FCCAnalyses::MCParticle::sel_genStatus(1)(Particle)")
        df = df.Define("Particle_px", "FCCAnalyses::MCParticle::get_px(MC_final)")
        df = df.Define("Particle_py", "FCCAnalyses::MCParticle::get_py(MC_final)")
        df = df.Define("Particle_pz", "FCCAnalyses::MCParticle::get_pz(MC_final)")
        df = df.Define("Particle_e",  "FCCAnalyses::MCParticle::get_e(MC_final)")

        df = df.Define(
            "pseudo_jets_gen",
            "FCCAnalyses::JetClusteringUtils::set_pseudoJets(Particle_px, Particle_py, Particle_pz, Particle_e)"
        )
        df = df.Define(
            "jets_gen_obj4",
            "JetClustering::clustering_ee_kt(2, 4, 0, 0)(pseudo_jets_gen)"
        )
        df = df.Define(
            "jets_gen4",
            "FCCAnalyses::JetClusteringUtils::get_pseudoJets(jets_gen_obj4)"
        )

        # ------------------------------------------------
        # RECO JETS
        # ------------------------------------------------
        df = df.Define("Reco_px", "FCCAnalyses::ReconstructedParticle::get_px(reco_clean)")
        df = df.Define("Reco_py", "FCCAnalyses::ReconstructedParticle::get_py(reco_clean)")
        df = df.Define("Reco_pz", "FCCAnalyses::ReconstructedParticle::get_pz(reco_clean)")
        df = df.Define("Reco_e",  "FCCAnalyses::ReconstructedParticle::get_e(reco_clean)")

        df = df.Define(
            "pseudo_jets_reco",
            "FCCAnalyses::JetClusteringUtils::set_pseudoJets(Reco_px, Reco_py, Reco_pz, Reco_e)"
        )
        df = df.Define(
            "jets_reco_obj4",
            "JetClustering::clustering_ee_kt(2, 4, 0, 0)(pseudo_jets_reco)"
        )
        df = df.Define(
            "jets_reco4",
            "FCCAnalyses::JetClusteringUtils::get_pseudoJets(jets_reco_obj4)"
        )

        df = df.Define("n_jets_gen", "jets_gen4.size()")
        df = df.Define("n_jets_reco", "jets_reco4.size()")
        df = df.Filter("n_jets_gen == 4 && n_jets_reco == 4", "Require exactly 4 gen and 4 reco jets")

        # ------------------------------------------------
        # RECO -> GEN MATCHING (common preselection)
        # ------------------------------------------------
        df = df.Define("jet_match_indices", "greedyJetMatching(jets_reco4, jets_gen4, 0.1)")
        df = df.Define("n_matched_jets", "countMatchedJets(jet_match_indices)")
        df = df.Filter("n_matched_jets == 4", "Require all 4 reco jets matched to gen jets")

        # ------------------------------------------------
        # Angular resolutions based on reco->gen matching
        # ------------------------------------------------
        df = df.Define("delta_theta_matched", "getDeltaTheta(jets_reco4, jets_gen4, jet_match_indices)")
        df = df.Define("delta_phi_matched",   "getDeltaPhi(jets_reco4, jets_gen4, jet_match_indices)")
        df = df.Define("delta_eta_matched",   "getDeltaEta(jets_reco4, jets_gen4, jet_match_indices)")
        df = df.Define("delta_mass_matched",  "getDeltaMass(jets_reco4, jets_gen4, jet_match_indices)")

        df = df.Define("x_gen_matched",  "getXGen(jets_gen4, jet_match_indices)")
        df = df.Define("x_reco_matched", "getXReco(jets_reco4, jet_match_indices)")
        df = df.Define("delta_x_matched", "x_reco_matched - x_gen_matched")

        # ------------------------------------------------
        # PARTON MATCHING STRATEGIES
        # ------------------------------------------------

        # Strategy A: original non-exclusive direct jet->parton matching
        df = df.Define("parton_match_A", "jetPartonMatching(jets_reco4, parton_eta, parton_phi)")
        df = df.Define("n_matched_partons_A", "countValidMatches(parton_match_A)")
        df = df.Define("n_unique_partons_A", "countUniqueMatches(parton_match_A)")
        df = df.Define("has_duplicate_A", "hasDuplicateMatches(parton_match_A)")

        # Strategy B: direct greedy exclusive reco->parton matching
        df = df.Define("parton_match_B", "greedyJetMatching(jets_reco4, pseudo_jets_partons, 0.4)")
        df = df.Define("n_matched_partons_B", "countValidMatches(parton_match_B)")
        df = df.Define("n_unique_partons_B", "countUniqueMatches(parton_match_B)")

        # Strategy C: transitive reco->gen->parton
        df = df.Define("gen_to_parton_match", "greedyJetMatching(jets_gen4, pseudo_jets_partons, 0.4)")
        df = df.Define("parton_match_C", "transitiveMatch(jet_match_indices, gen_to_parton_match)")
        df = df.Define("n_matched_partons_C", "countValidMatches(parton_match_C)")
        df = df.Define("n_unique_partons_C", "countUniqueMatches(parton_match_C)")

        # ------------------------------------------------
        # Agreement diagnostics
        # ------------------------------------------------
        df = df.Define("agree_AB", "exactMatchAgreement(parton_match_A, parton_match_B)")
        df = df.Define("agree_AC", "exactMatchAgreement(parton_match_A, parton_match_C)")
        df = df.Define("agree_BC", "exactMatchAgreement(parton_match_B, parton_match_C)")

        df = df.Define("n_agree_AB", "countAssignmentAgreement(parton_match_A, parton_match_B)")
        df = df.Define("n_agree_AC", "countAssignmentAgreement(parton_match_A, parton_match_C)")
        df = df.Define("n_agree_BC", "countAssignmentAgreement(parton_match_B, parton_match_C)")

        # ------------------------------------------------
        # Delta alpha for each strategy
        # ------------------------------------------------
        df = df.Define("delta_alpha_A", "getDeltaAlphaParton(jets_reco4, parton_energies, parton_match_A)")
        df = df.Define("delta_alpha_B", "getDeltaAlphaParton(jets_reco4, parton_energies, parton_match_B)")
        df = df.Define("delta_alpha_C", "getDeltaAlphaParton(jets_reco4, parton_energies, parton_match_C)")

        # ------------------------------------------------
        # Individual jet variables
        # ------------------------------------------------
        for jet_idx in range(1, 5):
            idx = jet_idx - 1

            df = df.Define(f"delta_theta_j{jet_idx}", f"getElement(delta_theta_matched, {idx})")
            df = df.Define(f"delta_phi_j{jet_idx}",   f"getElement(delta_phi_matched, {idx})")
            df = df.Define(f"delta_eta_j{jet_idx}",   f"getElement(delta_eta_matched, {idx})")
            df = df.Define(f"delta_x_j{jet_idx}",     f"getElement(delta_x_matched, {idx})")

            df = df.Define(f"delta_alpha_A_j{jet_idx}", f"getElement(delta_alpha_A, {idx})")
            df = df.Define(f"delta_alpha_B_j{jet_idx}", f"getElement(delta_alpha_B, {idx})")
            df = df.Define(f"delta_alpha_C_j{jet_idx}", f"getElement(delta_alpha_C, {idx})")

        # ------------------------------------------------
        # Filtered versions
        # ------------------------------------------------
        df = df.Define(
            "filtered_delta_theta_matched",
            f"filterValues(delta_theta_matched, {filter_config['delta_theta']['min']}, {filter_config['delta_theta']['max']})"
        )
        df = df.Define(
            "filtered_delta_phi_matched",
            f"filterValues(delta_phi_matched, {filter_config['delta_phi']['min']}, {filter_config['delta_phi']['max']})"
        )
        df = df.Define(
            "filtered_delta_eta_matched",
            f"filterValues(delta_eta_matched, {filter_config['delta_eta']['min']}, {filter_config['delta_eta']['max']})"
        )
        df = df.Define(
            "filtered_delta_x_matched",
            f"filterValues(delta_x_matched, {filter_config['delta_x']['min']}, {filter_config['delta_x']['max']})"
        )

        df = df.Define(
            "filtered_delta_alpha_A",
            f"filterValues(delta_alpha_A, {filter_config['delta_alpha']['min']}, {filter_config['delta_alpha']['max']})"
        )
        df = df.Define(
            "filtered_delta_alpha_B",
            f"filterValues(delta_alpha_B, {filter_config['delta_alpha']['min']}, {filter_config['delta_alpha']['max']})"
        )
        df = df.Define(
            "filtered_delta_alpha_C",
            f"filterValues(delta_alpha_C, {filter_config['delta_alpha']['min']}, {filter_config['delta_alpha']['max']})"
        )

        for jet_idx in range(1, 5):
            df = df.Define(
                f"filtered_delta_theta_j{jet_idx}",
                f"(delta_theta_j{jet_idx} >= {filter_config['delta_theta']['min']} && delta_theta_j{jet_idx} <= {filter_config['delta_theta']['max']}) ? delta_theta_j{jet_idx} : -999.0f"
            )
            df = df.Define(
                f"filtered_delta_phi_j{jet_idx}",
                f"(delta_phi_j{jet_idx} >= {filter_config['delta_phi']['min']} && delta_phi_j{jet_idx} <= {filter_config['delta_phi']['max']}) ? delta_phi_j{jet_idx} : -999.0f"
            )
            df = df.Define(
                f"filtered_delta_x_j{jet_idx}",
                f"(delta_x_j{jet_idx} >= {filter_config['delta_x']['min']} && delta_x_j{jet_idx} <= {filter_config['delta_x']['max']}) ? delta_x_j{jet_idx} : -999.0f"
            )

            df = df.Define(
                f"filtered_delta_alpha_A_j{jet_idx}",
                f"(delta_alpha_A_j{jet_idx} >= {filter_config['delta_alpha']['min']} && delta_alpha_A_j{jet_idx} <= {filter_config['delta_alpha']['max']}) ? delta_alpha_A_j{jet_idx} : -999.0f"
            )
            df = df.Define(
                f"filtered_delta_alpha_B_j{jet_idx}",
                f"(delta_alpha_B_j{jet_idx} >= {filter_config['delta_alpha']['min']} && delta_alpha_B_j{jet_idx} <= {filter_config['delta_alpha']['max']}) ? delta_alpha_B_j{jet_idx} : -999.0f"
            )
            df = df.Define(
                f"filtered_delta_alpha_C_j{jet_idx}",
                f"(delta_alpha_C_j{jet_idx} >= {filter_config['delta_alpha']['min']} && delta_alpha_C_j{jet_idx} <= {filter_config['delta_alpha']['max']}) ? delta_alpha_C_j{jet_idx} : -999.0f"
            )

        return df

    @staticmethod
    def output():
        outputs = [
            # Angular / energy response
            "delta_theta_j1", "delta_theta_j2", "delta_theta_j3", "delta_theta_j4",
            "delta_phi_j1", "delta_phi_j2", "delta_phi_j3", "delta_phi_j4",
            "delta_eta_j1", "delta_eta_j2", "delta_eta_j3", "delta_eta_j4",
            "delta_x_j1", "delta_x_j2", "delta_x_j3", "delta_x_j4",

            # Delta alpha by strategy
            "delta_alpha_A_j1", "delta_alpha_A_j2", "delta_alpha_A_j3", "delta_alpha_A_j4",
            "delta_alpha_B_j1", "delta_alpha_B_j2", "delta_alpha_B_j3", "delta_alpha_B_j4",
            "delta_alpha_C_j1", "delta_alpha_C_j2", "delta_alpha_C_j3", "delta_alpha_C_j4",

            # Filtered vector outputs
            "filtered_delta_theta_matched",
            "filtered_delta_phi_matched",
            "filtered_delta_eta_matched",
            "filtered_delta_x_matched",
            "filtered_delta_alpha_A",
            "filtered_delta_alpha_B",
            "filtered_delta_alpha_C",

            # Filtered per-jet outputs
            "filtered_delta_theta_j1", "filtered_delta_theta_j2", "filtered_delta_theta_j3", "filtered_delta_theta_j4",
            "filtered_delta_phi_j1", "filtered_delta_phi_j2", "filtered_delta_phi_j3", "filtered_delta_phi_j4",
            "filtered_delta_x_j1", "filtered_delta_x_j2", "filtered_delta_x_j3", "filtered_delta_x_j4",

            "filtered_delta_alpha_A_j1", "filtered_delta_alpha_A_j2", "filtered_delta_alpha_A_j3", "filtered_delta_alpha_A_j4",
            "filtered_delta_alpha_B_j1", "filtered_delta_alpha_B_j2", "filtered_delta_alpha_B_j3", "filtered_delta_alpha_B_j4",
            "filtered_delta_alpha_C_j1", "filtered_delta_alpha_C_j2", "filtered_delta_alpha_C_j3", "filtered_delta_alpha_C_j4",

            # Matching diagnostics
            "jet_match_indices",
            "parton_match_A",
            "parton_match_B",
            "parton_match_C",
            "gen_to_parton_match",

            "n_matched_jets",
            "n_matched_partons_A",
            "n_matched_partons_B",
            "n_matched_partons_C",

            "n_unique_partons_A",
            "n_unique_partons_B",
            "n_unique_partons_C",

            "has_duplicate_A",
            "agree_AB",
            "agree_AC",
            "agree_BC",
            "n_agree_AB",
            "n_agree_AC",
            "n_agree_BC",

            # Parton kinematics
            "parton_energies",
            "parton_eta",
            "parton_phi",
            "parton_y"
        ]
        return outputs
