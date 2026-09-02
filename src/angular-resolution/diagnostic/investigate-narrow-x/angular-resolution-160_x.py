import ROOT, os

ROOT.gROOT.SetBatch(True)

# ================================================
# STAGE 3: ANGULAR RESOLUTION
# ================================================

fractions = 1.0

inputDir = "outputs/stage2_jetclustering/"

energy_configs = {
    "p8_ee_WW_ecm160": {"fraction": fractions, "output_suffix": "ecm160", "label": "√s = 160 GeV"},
#   "p8_ee_WW_ecm240": {"fraction": fractions, "output_suffix": "ecm240", "label": "√s = 240 GeV"},
#   "p8_ee_WW_ecm340": {"fraction": fractions, "output_suffix": "ecm340", "label": "√s = 340 GeV"},
#   "p8_ee_WW_ecm345": {"fraction": fractions, "output_suffix": "ecm345", "label": "√s = 345 GeV"},
#   "p8_ee_WW_ecm350": {"fraction": fractions, "output_suffix": "ecm350", "label": "√s = 350 GeV"},
#   "p8_ee_WW_ecm355": {"fraction": fractions, "output_suffix": "ecm355", "label": "√s = 355 GeV"},
#   "p8_ee_WW_ecm365": {"fraction": fractions, "output_suffix": "ecm365", "label": "√s = 365 GeV"},
}

processList = {}
 
for process_name, config in energy_configs.items():
    processList[f"jetclustering_{config['output_suffix']}"] = {
        "fraction": config["fraction"],
        "chunks": 1,
        "output": f"angular_resolution_{config['output_suffix']}"
    }

outputDir = "outputs/stage3_angular_resolution/"

procDict = "FCCee_procDict_winter2023_IDEA.json"

nCPUS = 10

doTree = True

filter_config = {
    "delta_theta": {"min": -0.02, "max": 0.02},
    "delta_phi": {"min": -0.02, "max": 0.02},
    "delta_eta": {"min": -0.1, "max": 0.1},
    "delta_x": {"min": -0.6, "max": 0.6},
    "delta_alpha": {"min": -0.25, "max": 0.25}
}

includePaths = [
    # Changed
    "headers/greedyJetMatching_v2.h",

    "headers/getDeltaTheta.h",
    "headers/getDeltaPhi.h",
    "headers/getDeltaEta.h",
    "headers/getDeltaMass.h",
    #The only changed header is getXGen_v2.h and getXReco_v2.h, 
    #which are used to calculate x = log(p/m) for gen and reco jets.
    "headers/getXGen_v2.h",
    "headers/getXReco_v2.h",
    #getElement_v2 is for replacing -999 with NaN for failed values, 
    #and is used in the per-jet filtered values.
    "headers/getElement_v2.h",
    "headers/jetPartonMatching.h",
    "headers/getDeltaAlphaParton.h",
    "headers/filterValues.h",
]


class RDFanalysis:
    @staticmethod
    def analysers(df):

        df = df.Define(
            "jets_gen4",
            "FCCAnalyses::JetClusteringUtils::set_pseudoJets(jet_gen_px, jet_gen_py, jet_gen_pz, jet_gen_e)"
        )

        df = df.Define(
            "jets_reco4",
            "FCCAnalyses::JetClusteringUtils::set_pseudoJets(jet_reco_px, jet_reco_py, jet_reco_pz, jet_reco_e)"
        )
        df = df.Define("jet_gen_p_all", "FCCAnalyses::JetClusteringUtils::get_p(jets_gen4)")
        df = df.Define("jet_reco_p_all", "FCCAnalyses::JetClusteringUtils::get_p(jets_reco4)")
        df = df.Define("jet_gen_m_all", "FCCAnalyses::JetClusteringUtils::get_m(jets_gen4)")
        df = df.Define("jet_reco_m_all", "FCCAnalyses::JetClusteringUtils::get_m(jets_reco4)")


        # ------------------------------------------------
        # JET MATCHING
        # ------------------------------------------------
        # Changed 
        df = df.Define("jet_match_indices", "greedyJetMatching(jets_reco4, jets_gen4)")
        df = df.Define("n_matched_jets", "countMatchedJets(jet_match_indices)")
        #df = df.Filter("n_matched_jets == 4", "Require all 4 jets matched")

        # Added
        df = df.Define("matched_deltaR", "getMatchedDeltaR(jets_reco4, jets_gen4, jet_match_indices)")

        df = df.Define("max_matched_deltaR", "getMaxMatchedDeltaR(matched_deltaR)")

        df = df.Define("n_matched_deltaR_gt_0p1", "countMatchedDeltaRAbove(matched_deltaR, 0.1)")

        df = df.Define("n_matched_deltaR_gt_1p0", "countMatchedDeltaRAbove(matched_deltaR, 1.0)")

        # ------------------------------------------------
        # ANGULAR RESOLUTIONS
        # ------------------------------------------------
        df = df.Define("delta_theta_matched", "getDeltaTheta(jets_reco4, jets_gen4, jet_match_indices)")
        df = df.Define("delta_phi_matched", "getDeltaPhi(jets_reco4, jets_gen4, jet_match_indices)")
        df = df.Define("delta_eta_matched", "getDeltaEta(jets_reco4, jets_gen4, jet_match_indices)")
        df = df.Define("delta_mass_matched", "getDeltaMass(jets_reco4, jets_gen4, jet_match_indices)")

        # ------------------------------------------------
        # x = log(p/m)
        # ------------------------------------------------
        df = df.Define("x_gen_matched", "getXGen_v2(jets_gen4, jet_match_indices)")
        df = df.Define("x_reco_matched", "getXReco_v2(jets_reco4, jet_match_indices)")
        df = df.Define("delta_x_matched", "x_reco_matched - x_gen_matched")

        # ------------------------------------------------
        # ALPHA AND PARTON ENERGY
        # ------------------------------------------------
        df = df.Define("parton_matched", "jetPartonMatching(jets_reco4, parton_eta, parton_phi)")
        df = df.Define("delta_alpha", "getDeltaAlphaParton(jets_reco4, parton_energies, parton_matched)")

        # ------------------------------------------------
        # PER-JET VALUES
        # ------------------------------------------------
        for jet_idx in range(1, 5):
            idx = jet_idx - 1
            df = df.Define(f"delta_alpha_j{jet_idx}", f"getElement_v2(delta_alpha, {idx})")
            df = df.Define(f"delta_theta_j{jet_idx}", f"getElement_v2(delta_theta_matched, {idx})")
            df = df.Define(f"delta_phi_j{jet_idx}", f"getElement_v2(delta_phi_matched, {idx})")
            df = df.Define(f"delta_eta_j{jet_idx}", f"getElement_v2(delta_eta_matched, {idx})")
            df = df.Define(f"delta_x_j{jet_idx}", f"getElement_v2(delta_x_matched, {idx})")

        # ------------------------------------------------
        # DEBUG ANGLES
        # ------------------------------------------------
        df = df.Define("theta_gen_all", "FCCAnalyses::JetClusteringUtils::get_theta(jets_gen4)")
        df = df.Define("theta_reco_all", "FCCAnalyses::JetClusteringUtils::get_theta(jets_reco4)")
        df = df.Define("phi_gen_all", "FCCAnalyses::JetClusteringUtils::get_phi(jets_gen4)")
        df = df.Define("phi_reco_all", "FCCAnalyses::JetClusteringUtils::get_phi(jets_reco4)")

        # ------------------------------------------------
        # FILTERED RVEC VALUES
        # ------------------------------------------------
        df = df.Define("filtered_delta_theta_matched", f"filterValues(delta_theta_matched, {filter_config['delta_theta']['min']}, {filter_config['delta_theta']['max']})")
        df = df.Define("filtered_delta_phi_matched", f"filterValues(delta_phi_matched, {filter_config['delta_phi']['min']}, {filter_config['delta_phi']['max']})")
        df = df.Define("filtered_delta_eta_matched", f"filterValues(delta_eta_matched, {filter_config['delta_eta']['min']}, {filter_config['delta_eta']['max']})")
        df = df.Define("filtered_delta_x_matched", f"filterValues(delta_x_matched, {filter_config['delta_x']['min']}, {filter_config['delta_x']['max']})")
        df = df.Define("filtered_delta_alpha", f"filterValues(delta_alpha, {filter_config['delta_alpha']['min']}, {filter_config['delta_alpha']['max']})")

        # ------------------------------------------------
        # FILTERED PER-JET VALUES
        # keep NaN for failed values
        # ------------------------------------------------
        for jet_idx in range(1, 5):
            df = df.Define(f"filtered_delta_theta_j{jet_idx}", f"(std::isfinite(delta_theta_j{jet_idx}) && delta_theta_j{jet_idx} >= {filter_config['delta_theta']['min']} && delta_theta_j{jet_idx} <= {filter_config['delta_theta']['max']}) ? delta_theta_j{jet_idx} : std::numeric_limits<float>::quiet_NaN()")
            df = df.Define(f"filtered_delta_phi_j{jet_idx}", f"(std::isfinite(delta_phi_j{jet_idx}) && delta_phi_j{jet_idx} >= {filter_config['delta_phi']['min']} && delta_phi_j{jet_idx} <= {filter_config['delta_phi']['max']}) ? delta_phi_j{jet_idx} : std::numeric_limits<float>::quiet_NaN()")
            df = df.Define(f"filtered_delta_eta_j{jet_idx}", f"(std::isfinite(delta_eta_j{jet_idx}) && delta_eta_j{jet_idx} >= {filter_config['delta_eta']['min']} && delta_eta_j{jet_idx} <= {filter_config['delta_eta']['max']}) ? delta_eta_j{jet_idx} : std::numeric_limits<float>::quiet_NaN()")
            df = df.Define(f"filtered_delta_x_j{jet_idx}", f"(std::isfinite(delta_x_j{jet_idx}) && delta_x_j{jet_idx} >= {filter_config['delta_x']['min']} && delta_x_j{jet_idx} <= {filter_config['delta_x']['max']}) ? delta_x_j{jet_idx} : std::numeric_limits<float>::quiet_NaN()")
            df = df.Define(f"filtered_delta_alpha_j{jet_idx}", f"(std::isfinite(delta_alpha_j{jet_idx}) && delta_alpha_j{jet_idx} >= {filter_config['delta_alpha']['min']} && delta_alpha_j{jet_idx} <= {filter_config['delta_alpha']['max']}) ? delta_alpha_j{jet_idx} : std::numeric_limits<float>::quiet_NaN()")
        # Added 
        for jet_idx in range(1, 5):
            idx = jet_idx - 1
            df = df.Define(f"jet_gen_p_j{jet_idx}", f"getElement_v2(jet_gen_p_all, {idx})")
            df = df.Define(f"jet_reco_p_j{jet_idx}", f"getElement_v2(jet_reco_p_all, {idx})")
            df = df.Define(f"jet_gen_m_j{jet_idx}", f"getElement_v2(jet_gen_m_all, {idx})")
            df = df.Define(f"jet_reco_m_j{jet_idx}", f"getElement_v2(jet_reco_m_all, {idx})")
            df = df.Define(f"delta_logp_j{jet_idx}", f"std::log(jet_reco_p_j{jet_idx} / jet_gen_p_j{jet_idx})")
            df = df.Define(f"delta_logm_j{jet_idx}", f"std::log(jet_reco_m_j{jet_idx} / jet_gen_m_j{jet_idx})")
        return df
        

    @staticmethod
    def output():
        outputs = [
            "delta_theta_j1", "delta_theta_j2", "delta_theta_j3", "delta_theta_j4",
            "delta_phi_j1", "delta_phi_j2", "delta_phi_j3", "delta_phi_j4",
            "delta_eta_j1", "delta_eta_j2", "delta_eta_j3", "delta_eta_j4",
            "delta_x_j1", "delta_x_j2", "delta_x_j3", "delta_x_j4",
            "delta_alpha_j1", "delta_alpha_j2", "delta_alpha_j3", "delta_alpha_j4",

            "filtered_delta_theta_matched",
            "filtered_delta_phi_matched",
            "filtered_delta_eta_matched",
            "filtered_delta_x_matched",
            "filtered_delta_alpha",

            "filtered_delta_theta_j1", "filtered_delta_theta_j2", "filtered_delta_theta_j3", "filtered_delta_theta_j4",
            "filtered_delta_phi_j1", "filtered_delta_phi_j2", "filtered_delta_phi_j3", "filtered_delta_phi_j4",
            "filtered_delta_x_j1", "filtered_delta_x_j2", "filtered_delta_x_j3", "filtered_delta_x_j4",
            "filtered_delta_alpha_j1", "filtered_delta_alpha_j2", "filtered_delta_alpha_j3", "filtered_delta_alpha_j4",

            "parton_energies", "parton_eta", "parton_phi", "parton_y",

            "n_matched_jets",
            # Added Plots
            "jet_gen_p_j1", "jet_gen_p_j2", "jet_gen_p_j3", "jet_gen_p_j4",
            "jet_reco_p_j1", "jet_reco_p_j2", "jet_reco_p_j3", "jet_reco_p_j4",
            "jet_gen_m_j1", "jet_gen_m_j2", "jet_gen_m_j3", "jet_gen_m_j4",
            "jet_reco_m_j1", "jet_reco_m_j2", "jet_reco_m_j3", "jet_reco_m_j4",
            "delta_logp_j1", "delta_logp_j2", "delta_logp_j3", "delta_logp_j4",
            "delta_logm_j1", "delta_logm_j2", "delta_logm_j3", "delta_logm_j4",
            "matched_deltaR",
            "max_matched_deltaR",
            "n_matched_deltaR_gt_0p1",
            "n_matched_deltaR_gt_1p0",
            ]
        return outputs