import ROOT, os
from glob import glob

ROOT.gROOT.SetBatch(True)

# ================================================
# CONFIGURATION: Multiple Energy Folders
# ================================================

fractions = 1e-4 # Changed number in entries is because of this, change to your preference 

inputDir = "/eos/experiment/fcc/ee/generation/DelphesEvents/winter2023/IDEA/"

energy_configs = {
    "p8_ee_WW_ecm160": {
        "fraction": fractions,
        "output_suffix": "ecm160",
        "label": "√s = 160 GeV"
    },
    "p8_ee_WW_ecm240": {
        "fraction": fractions,
        "output_suffix": "ecm240",
        "label": "√s = 240 GeV"
    },
    "p8_ee_WW_ecm340": {
        "fraction": fractions,
        "output_suffix": "ecm340",
        "label": "√s = 340 GeV"
    },
    "p8_ee_WW_ecm345": {
        "fraction": fractions,
        "output_suffix": "ecm345",
        "label": "√s = 345 GeV"
    },
    "p8_ee_WW_ecm350": {
        "fraction": fractions,
        "output_suffix": "ecm350",
        "label": "√s = 350 GeV"
    },
    "p8_ee_WW_ecm355": {
        "fraction": fractions,
        "output_suffix": "ecm355",
        "label": "√s = 355 GeV"
    },
    "p8_ee_WW_ecm365": {
        "fraction": fractions,
        "output_suffix": "ecm365",
        "label": "√s = 365 GeV"
    },
}

filter_config = {
    "delta_theta": {"min": -0.02, "max": 0.02},
    "delta_phi":   {"min": -0.02, "max": 0.02},
    "delta_eta":   {"min": -0.1,  "max": 0.1},
    "delta_x":     {"min": -0.6,  "max": 0.6},
    "delta_alpha": {"min": -0.25, "max": 0.25}
}

# ISR photon selection
isr_min_energy = 1.0
isr_min_abs_cos_theta = 0.95

# ================================================
# Build processList for FCCAnalyses
# ================================================

processList = {}
chunks = 4

for process_name, config in energy_configs.items():

    file_pattern = os.path.join(
        inputDir,
        process_name,
        "events_*.root"
    )

    all_files = glob(file_pattern)

    if all_files:
        processList[process_name] = {
            "fraction": config["fraction"],
            "chunks": chunks,
            "output": (
                f"angular_resolution_ISR_"
                f"{config['output_suffix']}"
            )
        }

        print(
            f"Added process: {process_name} "
            f"({config['label']}) - "
            f"{len(all_files)} files found"
        )

    else:
        print(f"WARNING: No files found for {process_name}")


outputDir = "outputs/step1_angular_resolution_ISR_multiE/"

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
    "headers/getElement_v2.h", # Changed to remove peak on -999
    "headers/jetPartonMatching.h",
    "headers/getDeltaAlphaParton.h",
    "headers/filterValues.h",
    "headers/selectQuarks.h",
    "headers/selectISRPhotons.h"
]


# ------------------------------------------------
# RDFanalysis
# ------------------------------------------------

class RDFanalysis:

    @staticmethod
    def analysers(df):

        # ---------------------------
        # Remove leptons and ISR photons
        # ---------------------------

        df = df.Alias(
            "Electron0",
            "Electron#0.index"
        )

        df = df.Alias(
            "Muon0",
            "Muon#0.index"
        )

        df = df.Alias(
            "Photon0",
            "Photon#0.index"
        )

        df = df.Define(
            "ele_all",
            "FCCAnalyses::ReconstructedParticle::get("
            "Electron0, ReconstructedParticles)"
        )

        df = df.Define(
            "mu_all",
            "FCCAnalyses::ReconstructedParticle::get("
            "Muon0, ReconstructedParticles)"
        )

        df = df.Define(
            "pho_all",
            "FCCAnalyses::ReconstructedParticle::get("
            "Photon0, ReconstructedParticles)"
        )

        df = df.Define(
            "isr_photons",
            f"selectISRPhotons("
            f"pho_all, "
            f"{isr_min_energy}f, "
            f"{isr_min_abs_cos_theta}f)"
        )

        df = df.Define(
            "RP_noISR",
            "FCCAnalyses::ReconstructedParticle::remove("
            "ReconstructedParticles, isr_photons)"
        )

        df = df.Define(
            "RP_noEle",
            "FCCAnalyses::ReconstructedParticle::remove("
            "RP_noISR, ele_all)"
        )

        df = df.Define(
            "reco_clean",
            "FCCAnalyses::ReconstructedParticle::remove("
            "RP_noEle, mu_all)"
        )

        # ------------------------------------------------
        # PARTON GENERATION
        # ------------------------------------------------

        df = df.Define(
            "partons_all",
            "selectQuarks(Particle)"
        )

        df = df.Define(
            "n_partons",
            "partons_all.size()"
        )

        df = df.Filter(
            "n_partons == 4",
            "Require only 4 partons"
        )

        df = df.Define(
            "parton_energies",
            "FCCAnalyses::MCParticle::get_e(partons_all)"
        )

        df = df.Define(
            "parton_eta",
            "FCCAnalyses::MCParticle::get_eta(partons_all)"
        )

        df = df.Define(
            "parton_phi",
            "FCCAnalyses::MCParticle::get_phi(partons_all)"
        )

        df = df.Define(
            "parton_y",
            "FCCAnalyses::MCParticle::get_y(partons_all)"
        )

        # ------------------------------------------------
        # GEN JETS
        # ------------------------------------------------

        df = df.Define(
            "MC_final",
            "FCCAnalyses::MCParticle::sel_genStatus(1)(Particle)"
        )

        df = df.Define(
            "Particle_px",
            "FCCAnalyses::MCParticle::get_px(MC_final)"
        )

        df = df.Define(
            "Particle_py",
            "FCCAnalyses::MCParticle::get_py(MC_final)"
        )

        df = df.Define(
            "Particle_pz",
            "FCCAnalyses::MCParticle::get_pz(MC_final)"
        )

        df = df.Define(
            "Particle_e",
            "FCCAnalyses::MCParticle::get_e(MC_final)"
        )

        df = df.Define(
            "pseudo_jets_gen",
            "FCCAnalyses::JetClusteringUtils::set_pseudoJets("
            "Particle_px, Particle_py, Particle_pz, Particle_e)"
        )

        df = df.Define(
            "jets_gen_obj4",
            "JetClustering::clustering_ee_kt("
            "2, 4, 0, 0)(pseudo_jets_gen)"
        )

        df = df.Define(
            "jets_gen4",
            "FCCAnalyses::JetClusteringUtils::get_pseudoJets("
            "jets_gen_obj4)"
        )

        # ------------------------------------------------
        # RECO JETS
        # ------------------------------------------------

        df = df.Define(
            "Reco_px",
            "FCCAnalyses::ReconstructedParticle::get_px("
            "reco_clean)"
        )

        df = df.Define(
            "Reco_py",
            "FCCAnalyses::ReconstructedParticle::get_py("
            "reco_clean)"
        )

        df = df.Define(
            "Reco_pz",
            "FCCAnalyses::ReconstructedParticle::get_pz("
            "reco_clean)"
        )

        df = df.Define(
            "Reco_e",
            "FCCAnalyses::ReconstructedParticle::get_e("
            "reco_clean)"
        )

        df = df.Define(
            "pseudo_jets_reco",
            "FCCAnalyses::JetClusteringUtils::set_pseudoJets("
            "Reco_px, Reco_py, Reco_pz, Reco_e)"
        )

        df = df.Define(
            "jets_reco_obj4",
            "JetClustering::clustering_ee_kt("
            "2, 4, 0, 0)(pseudo_jets_reco)"
        )

        df = df.Define(
            "jets_reco4",
            "FCCAnalyses::JetClusteringUtils::get_pseudoJets("
            "jets_reco_obj4)"
        )

        df = df.Define(
            "n_jets_gen",
            "jets_gen4.size()"
        )

        df = df.Define(
            "n_jets_reco",
            "jets_reco4.size()"
        )

        df = df.Filter(
            "n_jets_gen == 4 && n_jets_reco == 4"
        )

        # ------------------------------------------------
        # JET MATCHING
        # ------------------------------------------------

        df = df.Define(
            "jet_match_indices",
            "greedyJetMatching("
            "jets_reco4, jets_gen4, 0.1)"
        )

        df = df.Define(
            "n_matched_jets",
            "countMatchedJets(jet_match_indices)"
        )

        df = df.Filter(
            "n_matched_jets == 4",
            "Require all 4 jets matched"
        )

        # ------------------------------------------------
        # ANGULAR RESOLUTIONS
        # ------------------------------------------------

        df = df.Define(
            "delta_theta_matched",
            "getDeltaTheta("
            "jets_reco4, jets_gen4, jet_match_indices)"
        )

        df = df.Define(
            "delta_phi_matched",
            "getDeltaPhi("
            "jets_reco4, jets_gen4, jet_match_indices)"
        )

        df = df.Define(
            "delta_eta_matched",
            "getDeltaEta("
            "jets_reco4, jets_gen4, jet_match_indices)"
        )

        df = df.Define(
            "delta_mass_matched",
            "getDeltaMass("
            "jets_reco4, jets_gen4, jet_match_indices)"
        )

        # ------------------------------------------------
        # x = log(p/E)
        # ------------------------------------------------

        df = df.Define(
            "x_gen_matched",
            "getXGen(jets_gen4, jet_match_indices)"
        )

        df = df.Define(
            "x_reco_matched",
            "getXReco(jets_reco4, jet_match_indices)"
        )

        df = df.Define(
            "delta_x_matched",
            "x_reco_matched - x_gen_matched"
        )

        # ------------------------------------------------
        # ALPHA AND PARTON ENERGY
        # ------------------------------------------------

        df = df.Define(
            "parton_matched",
            "jetPartonMatching("
            "jets_reco4, parton_eta, parton_phi)"
        )

        df = df.Define(
            "delta_alpha",
            "getDeltaAlphaParton("
            "jets_reco4, parton_energies, parton_matched)"
        )

        # ------------------------------------------------
        # Individual jet variables
        # ------------------------------------------------

        for jet_idx in range(1, 5):

            idx = jet_idx - 1

            df = df.Define(
                f"delta_alpha_j{jet_idx}",
                f"getElement_v2(delta_alpha, {idx})"
            )

            df = df.Define(
                f"delta_theta_j{jet_idx}",
                f"getElement_v2(delta_theta_matched, {idx})"
            )

            df = df.Define(
                f"delta_phi_j{jet_idx}",
                f"getElement_v2(delta_phi_matched, {idx})"
            )

            df = df.Define(
                f"delta_eta_j{jet_idx}",
                f"getElement_v2(delta_eta_matched, {idx})"
            )

            df = df.Define(
                f"delta_x_j{jet_idx}",
                f"getElement_v2(delta_x_matched, {idx})"
            )

        df = df.Define(
            "theta_gen_all",
            "FCCAnalyses::JetClusteringUtils::get_theta("
            "jets_gen4)"
        )

        df = df.Define(
            "theta_reco_all",
            "FCCAnalyses::JetClusteringUtils::get_theta("
            "jets_reco4)"
        )

        df = df.Define(
            "phi_gen_all",
            "FCCAnalyses::JetClusteringUtils::get_phi("
            "jets_gen4)"
        )

        df = df.Define(
            "phi_reco_all",
            "FCCAnalyses::JetClusteringUtils::get_phi("
            "jets_reco4)"
        )

        # ------------------------------------------------
        # FILTERED VERSIONS
        # ------------------------------------------------

        df = df.Define(
            "filtered_delta_theta_matched",
            f"filterValues("
            f"delta_theta_matched, "
            f"{filter_config['delta_theta']['min']}, "
            f"{filter_config['delta_theta']['max']})"
        )

        df = df.Define(
            "filtered_delta_phi_matched",
            f"filterValues("
            f"delta_phi_matched, "
            f"{filter_config['delta_phi']['min']}, "
            f"{filter_config['delta_phi']['max']})"
        )

        df = df.Define(
            "filtered_delta_eta_matched",
            f"filterValues("
            f"delta_eta_matched, "
            f"{filter_config['delta_eta']['min']}, "
            f"{filter_config['delta_eta']['max']})"
        )

        df = df.Define(
            "filtered_delta_x_matched",
            f"filterValues("
            f"delta_x_matched, "
            f"{filter_config['delta_x']['min']}, "
            f"{filter_config['delta_x']['max']})"
        )

        df = df.Define(
            "filtered_delta_alpha",
            f"filterValues("
            f"delta_alpha, "
            f"{filter_config['delta_alpha']['min']}, "
            f"{filter_config['delta_alpha']['max']})"
        )
        for jet_idx in range(1, 5): 

# changed -999 to NaN to remove the peak on the root files
            df = df.Define(
                f"filtered_delta_theta_j{jet_idx}",
                f"(std::isfinite(delta_theta_j{jet_idx}) && "
                f"delta_theta_j{jet_idx} >= {filter_config['delta_theta']['min']} && "
                f"delta_theta_j{jet_idx} <= {filter_config['delta_theta']['max']}) "
                f"? delta_theta_j{jet_idx} : "
                f"std::numeric_limits<float>::quiet_NaN()"
            )

            df = df.Define(
                f"filtered_delta_phi_j{jet_idx}",
                f"(std::isfinite(delta_phi_j{jet_idx}) && "
                f"delta_phi_j{jet_idx} >= {filter_config['delta_phi']['min']} && "
                f"delta_phi_j{jet_idx} <= {filter_config['delta_phi']['max']}) "
                f"? delta_phi_j{jet_idx} : "
                f"std::numeric_limits<float>::quiet_NaN()"
            )

            df = df.Define(
                f"filtered_delta_eta_j{jet_idx}",
                f"(std::isfinite(delta_eta_j{jet_idx}) && "
                f"delta_eta_j{jet_idx} >= {filter_config['delta_eta']['min']} && "
                f"delta_eta_j{jet_idx} <= {filter_config['delta_eta']['max']}) "
                f"? delta_eta_j{jet_idx} : "
                f"std::numeric_limits<float>::quiet_NaN()"
            )

            df = df.Define(
                f"filtered_delta_x_j{jet_idx}",
                f"(std::isfinite(delta_x_j{jet_idx}) && "
                f"delta_x_j{jet_idx} >= {filter_config['delta_x']['min']} && "
                f"delta_x_j{jet_idx} <= {filter_config['delta_x']['max']}) "
                f"? delta_x_j{jet_idx} : "
                f"std::numeric_limits<float>::quiet_NaN()"
            )

            df = df.Define(
                f"filtered_delta_alpha_j{jet_idx}",
                f"(std::isfinite(delta_alpha_j{jet_idx}) && "
                f"delta_alpha_j{jet_idx} >= {filter_config['delta_alpha']['min']} && "
                f"delta_alpha_j{jet_idx} <= {filter_config['delta_alpha']['max']}) "
                f"? delta_alpha_j{jet_idx} : "
                f"std::numeric_limits<float>::quiet_NaN()"
            )

      
        return df

    @staticmethod
    def output():

        outputs = [
            "delta_theta_j1",
            "delta_theta_j2",
            "delta_theta_j3",
            "delta_theta_j4",

            "delta_phi_j1",
            "delta_phi_j2",
            "delta_phi_j3",
            "delta_phi_j4",

            "delta_eta_j1",
            "delta_eta_j2",
            "delta_eta_j3",
            "delta_eta_j4",

            "delta_x_j1",
            "delta_x_j2",
            "delta_x_j3",
            "delta_x_j4",

            "delta_alpha_j1",
            "delta_alpha_j2",
            "delta_alpha_j3",
            "delta_alpha_j4",

            "filtered_delta_theta_matched",
            "filtered_delta_phi_matched",
            "filtered_delta_eta_matched",
            "filtered_delta_x_matched",
            "filtered_delta_alpha",

            "filtered_delta_theta_j1",
            "filtered_delta_theta_j2",
            "filtered_delta_theta_j3",
            "filtered_delta_theta_j4",

            "filtered_delta_phi_j1",
            "filtered_delta_phi_j2",
            "filtered_delta_phi_j3",
            "filtered_delta_phi_j4",

            "filtered_delta_x_j1",
            "filtered_delta_x_j2",
            "filtered_delta_x_j3",
            "filtered_delta_x_j4",

            "filtered_delta_alpha_j1",
            "filtered_delta_alpha_j2",
            "filtered_delta_alpha_j3",
            "filtered_delta_alpha_j4",

            "parton_energies",
            "parton_eta",
            "parton_phi",
            "parton_y",

            "n_matched_jets"
        ]

        return outputs