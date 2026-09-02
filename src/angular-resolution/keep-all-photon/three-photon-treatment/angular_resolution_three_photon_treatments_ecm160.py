import ROOT
import os

ROOT.gROOT.SetBatch(True)

# ============================================================
# General configuration
# ============================================================

fractions = 1e-6

inputDir = (
    "/eos/experiment/fcc/ee/generation/"
    "DelphesEvents/winter2023/IDEA/"
)

processList = {
    "p8_ee_WW_ecm160": {
        "fraction": fractions,
        "chunks": 4,
        "output": (
            "angular_resolution_"
            "three_photon_treatments_ecm160"
        )
    }
}

outputDir = (
    "outputs/angular_resolution_"
    "three_photon_treatments/"
)

procDict = "FCCee_procDict_winter2023_IDEA.json"

nCPUS = 4
doTree = True

# ============================================================
# ISR-like reconstructed photon definition
# ============================================================

ISR_PHOTON_MIN_ENERGY = 1.0
ISR_PHOTON_MIN_ABS_COSTHETA = 0.95

# ============================================================
# Existing headers plus one new photon-treatment header
# ============================================================

includePaths = [
    "headers/greedyJetMatching.h",
    "headers/getDeltaTheta.h",
    "headers/getDeltaPhi.h",
    "headers/getXGen.h",
    "headers/getXReco.h",
    "headers/getElement.h",
    "headers/jetPartonMatching.h",
    "headers/getDeltaAlphaParton.h",
    "headers/selectQuarks.h",
    "headers/photonTreatment.h",
]


class RDFanalysis:

    @staticmethod
    def analysers(df):

        # ====================================================
        # Reconstructed particle collections
        # ====================================================

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

        # ====================================================
        # Photon treatment 1: remove all photons
        #
        # This reproduces the original reconstruction logic.
        # ====================================================

        df = df.Define(
            "RP_noPhotons",
            "FCCAnalyses::ReconstructedParticle::remove("
            "ReconstructedParticles, pho_all)"
        )

        df = df.Define(
            "RP_noPhotons_noElectrons",
            "FCCAnalyses::ReconstructedParticle::remove("
            "RP_noPhotons, ele_all)"
        )

        df = df.Define(
            "reco_noPhotons",
            "FCCAnalyses::ReconstructedParticle::remove("
            "RP_noPhotons_noElectrons, mu_all)"
        )

        # ====================================================
        # Photon treatment 2: keep all photons
        #
        # Only identified electrons and muons are removed.
        # ====================================================

        df = df.Define(
            "RP_allPhotons_noElectrons",
            "FCCAnalyses::ReconstructedParticle::remove("
            "ReconstructedParticles, ele_all)"
        )

        df = df.Define(
            "reco_allPhotons",
            "FCCAnalyses::ReconstructedParticle::remove("
            "RP_allPhotons_noElectrons, mu_all)"
        )

        # ====================================================
        # Photon treatment 3: remove only ISR-like photons
        # ====================================================

        df = df.Define(
            "isr_like_photons",
            "selectISRLikePhotonsReco("
            "pho_all, "
            f"{ISR_PHOTON_MIN_ENERGY}f, "
            f"{ISR_PHOTON_MIN_ABS_COSTHETA}f)"
        )

        df = df.Define(
            "RP_noISRPhotons",
            "FCCAnalyses::ReconstructedParticle::remove("
            "ReconstructedParticles, isr_like_photons)"
        )

        df = df.Define(
            "RP_noISRPhotons_noElectrons",
            "FCCAnalyses::ReconstructedParticle::remove("
            "RP_noISRPhotons, ele_all)"
        )

        df = df.Define(
            "reco_noISRPhotons",
            "FCCAnalyses::ReconstructedParticle::remove("
            "RP_noISRPhotons_noElectrons, mu_all)"
        )

        # ====================================================
        # Select the four hard-process quarks
        # ====================================================

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
            "Require exactly four selected partons"
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
            "sumE_partons",
            "ROOT::VecOps::Sum(parton_energies)"
        )

        # ====================================================
        # Generated jets
        #
        # This section remains equivalent to the original code.
        # ====================================================

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

        df = df.Define(
            "n_jets_gen",
            "jets_gen4.size()"
        )

        df = df.Filter(
            "n_jets_gen == 4",
            "Require four generated jets"
        )

        # ====================================================
        # Build reconstructed jets for all three treatments
        # ====================================================

        reco_collections = {
            "noPhotons": "reco_noPhotons",
            "allPhotons": "reco_allPhotons",
            "noISRPhotons": "reco_noISRPhotons",
        }

        for treatment, collection in reco_collections.items():

            df = df.Define(
                f"Reco_px_{treatment}",
                "FCCAnalyses::ReconstructedParticle::get_px("
                f"{collection})"
            )

            df = df.Define(
                f"Reco_py_{treatment}",
                "FCCAnalyses::ReconstructedParticle::get_py("
                f"{collection})"
            )

            df = df.Define(
                f"Reco_pz_{treatment}",
                "FCCAnalyses::ReconstructedParticle::get_pz("
                f"{collection})"
            )

            df = df.Define(
                f"Reco_e_{treatment}",
                "FCCAnalyses::ReconstructedParticle::get_e("
                f"{collection})"
            )

            df = df.Define(
                f"pseudo_jets_reco_{treatment}",
                "FCCAnalyses::JetClusteringUtils::set_pseudoJets("
                f"Reco_px_{treatment}, "
                f"Reco_py_{treatment}, "
                f"Reco_pz_{treatment}, "
                f"Reco_e_{treatment})"
            )

            df = df.Define(
                f"jets_reco_obj4_{treatment}",
                "JetClustering::clustering_ee_kt("
                f"2, 4, 0, 0)(pseudo_jets_reco_{treatment})"
            )

            df = df.Define(
                f"jets_reco4_{treatment}",
                "FCCAnalyses::JetClusteringUtils::get_pseudoJets("
                f"jets_reco_obj4_{treatment})"
            )

            df = df.Define(
                f"n_jets_reco_{treatment}",
                f"jets_reco4_{treatment}.size()"
            )

            df = df.Define(
                f"sumE_recoJets_{treatment}",
                f"sumJetEnergy(jets_reco4_{treatment})"
            )

        # Use the same events for the three comparisons.
        df = df.Filter(
            "n_jets_reco_noPhotons == 4 && "
            "n_jets_reco_allPhotons == 4 && "
            "n_jets_reco_noISRPhotons == 4",
            "Require four reconstructed jets in all treatments"
        )

        # ====================================================
        # Jet-to-generated-jet matching for each treatment
        # ====================================================

        for treatment in reco_collections:

            df = df.Define(
                f"jet_match_indices_{treatment}",
                "greedyJetMatching("
                f"jets_reco4_{treatment}, "
                "jets_gen4, "
                "0.1)"
            )

            df = df.Define(
                f"n_matched_jets_{treatment}",
                "countMatchedJets("
                f"jet_match_indices_{treatment})"
            )

        # Common matched event sample gives a direct comparison.
        df = df.Filter(
            "n_matched_jets_noPhotons == 4 && "
            "n_matched_jets_allPhotons == 4 && "
            "n_matched_jets_noISRPhotons == 4",
            "Require four matched jets in all treatments"
        )

        # ====================================================
        # Calculate the four fit variables
        #
        #   delta_theta
        #   delta_phi
        #   delta_x
        #   delta_alpha
        # ====================================================

        for treatment in reco_collections:

            # -------------------------------
            # Angular differences
            # -------------------------------

            df = df.Define(
                f"delta_theta_{treatment}",
                "getDeltaTheta("
                f"jets_reco4_{treatment}, "
                "jets_gen4, "
                f"jet_match_indices_{treatment})"
            )

            df = df.Define(
                f"delta_phi_{treatment}",
                "getDeltaPhi("
                f"jets_reco4_{treatment}, "
                "jets_gen4, "
                f"jet_match_indices_{treatment})"
            )

            # -------------------------------
            # x = log(p/E)
            # -------------------------------

            df = df.Define(
                f"x_gen_{treatment}",
                "getXGen("
                "jets_gen4, "
                f"jet_match_indices_{treatment})"
            )

            df = df.Define(
                f"x_reco_{treatment}",
                "getXReco("
                f"jets_reco4_{treatment}, "
                f"jet_match_indices_{treatment})"
            )

            df = df.Define(
                f"delta_x_{treatment}",
                f"x_reco_{treatment} - x_gen_{treatment}"
            )

            # -------------------------------
            # Jet-to-parton matching
            # -------------------------------

            df = df.Define(
                f"parton_matched_{treatment}",
                "jetPartonMatching("
                f"jets_reco4_{treatment}, "
                "parton_eta, "
                "parton_phi)"
            )

            # Fixed alpha function preserves one output per jet.
            df = df.Define(
                f"delta_alpha_{treatment}",
                "getDeltaAlphaPartonFixed("
                f"jets_reco4_{treatment}, "
                "parton_energies, "
                f"parton_matched_{treatment})"
            )

            # -------------------------------
            # Scalar branches for jets 1–4
            # -------------------------------

            for jet_number in range(1, 5):
                index = jet_number - 1

                df = df.Define(
                    f"delta_theta_{treatment}_j{jet_number}",
                    f"getElement(delta_theta_{treatment}, {index})"
                )

                df = df.Define(
                    f"delta_phi_{treatment}_j{jet_number}",
                    f"getElement(delta_phi_{treatment}, {index})"
                )

                df = df.Define(
                    f"delta_x_{treatment}_j{jet_number}",
                    f"getElement(delta_x_{treatment}, {index})"
                )

                df = df.Define(
                    f"delta_alpha_{treatment}_j{jet_number}",
                    f"getElement(delta_alpha_{treatment}, {index})"
                )

        return df

    @staticmethod
    def output():

        branches = [
            "sumE_partons",

            "sumE_recoJets_noPhotons",
            "sumE_recoJets_allPhotons",
            "sumE_recoJets_noISRPhotons",

            "n_matched_jets_noPhotons",
            "n_matched_jets_allPhotons",
            "n_matched_jets_noISRPhotons",
        ]

        treatments = [
            "noPhotons",
            "allPhotons",
            "noISRPhotons",
        ]

        variables = [
            "delta_theta",
            "delta_phi",
            "delta_x",
            "delta_alpha",
        ]

        for variable in variables:
            for treatment in treatments:
                for jet_number in range(1, 5):
                    branches.append(
                        f"{variable}_{treatment}_j{jet_number}"
                    )

        return branches
