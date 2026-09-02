import ROOT
import os

ROOT.gROOT.SetBatch(True)

# ------------------------------------------------------------
# Dataset / framework settings
# ------------------------------------------------------------
fractions = 1e-6

inputDir = "/eos/experiment/fcc/ee/generation/DelphesEvents/winter2023/IDEA/"

processList = {
    "p8_ee_WW_ecm160": {
        "fraction": fractions,
        "chunks": 4,
        "output": "ww160_deltaalpha_energy_photonISR_diagnostics"
    }
}

outputDir = "outputs/angular_resolution_energy_diagnostics/"
procDict = "FCCee_procDict_winter2023_IDEA.json"
nCPUS = 4
doTree = True

includePaths = [
    "headers/selectQuarks.h",
    "headers/jetPartonMatching.h",
    "headers/getDeltaAlphaParton.h",
    "headers/energyDiagnostics.h",
]


class RDFanalysis:
    @staticmethod
    def analysers(df):

        # ------------------------------------------------------------
        # Aliases
        # ------------------------------------------------------------
        df = df.Alias("Electron0", "Electron#0.index")
        df = df.Alias("Muon0", "Muon#0.index")
        df = df.Alias("Photon0", "Photon#0.index")

        # ------------------------------------------------------------
        # Reco object collections
        # ------------------------------------------------------------
        df = df.Define(
            "ele_all",
            "FCCAnalyses::ReconstructedParticle::get(Electron0, ReconstructedParticles)"
        )

        df = df.Define(
            "mu_all",
            "FCCAnalyses::ReconstructedParticle::get(Muon0, ReconstructedParticles)"
        )

        df = df.Define(
            "pho_all",
            "FCCAnalyses::ReconstructedParticle::get(Photon0, ReconstructedParticles)"
        )

        # Keep photons, remove only identified electrons/muons.
        df = df.Define(
            "RP_noEle",
            "FCCAnalyses::ReconstructedParticle::remove(ReconstructedParticles, ele_all)"
        )

        df = df.Define(
            "reco_withPhotons",
            "FCCAnalyses::ReconstructedParticle::remove(RP_noEle, mu_all)"
        )

        # Current problematic strategy: photons removed before clustering.
        df = df.Define(
            "reco_photonsRemoved",
            "FCCAnalyses::ReconstructedParticle::remove(reco_withPhotons, pho_all)"
        )

        # ------------------------------------------------------------
        # Basic reco photon / lepton diagnostics
        # ------------------------------------------------------------
        df = df.Define("n_reco_electrons", "ele_all.size()")
        df = df.Define("n_reco_muons", "mu_all.size()")
        df = df.Define("n_reco_photons", "pho_all.size()")

        df = df.Define("sumE_reco_all", "sumRecoEnergy(ReconstructedParticles)")
        df = df.Define("sumE_reco_electrons", "sumRecoEnergy(ele_all)")
        df = df.Define("sumE_reco_muons", "sumRecoEnergy(mu_all)")
        df = df.Define("sumE_reco_photons", "sumRecoEnergy(pho_all)")

        df = df.Define(
            "photon_energy_fraction",
            "sumE_reco_all > 0.0f ? sumE_reco_photons / sumE_reco_all : 0.0f"
        )

        # ISR-like reco photons:
        # forward/backward photon with |cos(theta)| > 0.95 and E > 1 GeV.
        df = df.Define(
            "sumE_ISR_candidate_reco",
            "sumForwardPhotonEnergyReco(pho_all, 1.0f, 0.95f)"
        )

        df = df.Define(
            "n_ISR_candidate_reco",
            "countForwardPhotonsReco(pho_all, 1.0f, 0.95f)"
        )

        df = df.Define(
            "has_ISR_candidate_reco",
            "sumE_ISR_candidate_reco > 1.0f"
        )

        # ------------------------------------------------------------
        # MC final-state photon diagnostics
        # ------------------------------------------------------------
        df = df.Define(
            "mc_final_photons",
            "selectMCFinalPhotons(Particle)"
        )

        df = df.Define("n_mc_photons", "mc_final_photons.size()")
        df = df.Define("sumE_mc_final_photons", "sumMCEnergy(mc_final_photons)")

        df = df.Define(
            "sumE_ISR_candidate_mc",
            "sumForwardPhotonEnergyMC(mc_final_photons, 1.0f, 0.95f)"
        )

        df = df.Define(
            "n_ISR_candidate_mc",
            "countForwardPhotonsMC(mc_final_photons, 1.0f, 0.95f)"
        )

        df = df.Define(
            "has_ISR_candidate_mc",
            "sumE_ISR_candidate_mc > 1.0f"
        )

        # ------------------------------------------------------------
        # Parton truth reference
        # No additional matched-parton filter is added.
        # ------------------------------------------------------------
        df = df.Define("partons_all", "selectQuarks(Particle)")
        df = df.Define("n_partons", "partons_all.size()")

        # Keep only events with four hard-process quarks.
        df = df.Filter("n_partons == 4", "Require exactly four selected partons")

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

        # ------------------------------------------------------------
        # Reco jet clustering: photons kept
        # ------------------------------------------------------------
        df = df.Define(
            "Reco_px_withPhotons",
            "FCCAnalyses::ReconstructedParticle::get_px(reco_withPhotons)"
        )

        df = df.Define(
            "Reco_py_withPhotons",
            "FCCAnalyses::ReconstructedParticle::get_py(reco_withPhotons)"
        )

        df = df.Define(
            "Reco_pz_withPhotons",
            "FCCAnalyses::ReconstructedParticle::get_pz(reco_withPhotons)"
        )

        df = df.Define(
            "Reco_e_withPhotons",
            "FCCAnalyses::ReconstructedParticle::get_e(reco_withPhotons)"
        )

        df = df.Define(
            "pseudo_jets_withPhotons",
            "FCCAnalyses::JetClusteringUtils::set_pseudoJets("
            "Reco_px_withPhotons, Reco_py_withPhotons, "
            "Reco_pz_withPhotons, Reco_e_withPhotons)"
        )

        df = df.Define(
            "jets_obj4_withPhotons",
            "JetClustering::clustering_ee_kt(2, 4, 0, 0)(pseudo_jets_withPhotons)"
        )

        df = df.Define(
            "jets_reco4_withPhotons",
            "FCCAnalyses::JetClusteringUtils::get_pseudoJets(jets_obj4_withPhotons)"
        )

        df = df.Define("n_jets_withPhotons", "jets_reco4_withPhotons.size()")

        # ------------------------------------------------------------
        # Reco jet clustering: photons removed
        # ------------------------------------------------------------
        df = df.Define(
            "Reco_px_photonsRemoved",
            "FCCAnalyses::ReconstructedParticle::get_px(reco_photonsRemoved)"
        )

        df = df.Define(
            "Reco_py_photonsRemoved",
            "FCCAnalyses::ReconstructedParticle::get_py(reco_photonsRemoved)"
        )

        df = df.Define(
            "Reco_pz_photonsRemoved",
            "FCCAnalyses::ReconstructedParticle::get_pz(reco_photonsRemoved)"
        )

        df = df.Define(
            "Reco_e_photonsRemoved",
            "FCCAnalyses::ReconstructedParticle::get_e(reco_photonsRemoved)"
        )

        df = df.Define(
            "pseudo_jets_photonsRemoved",
            "FCCAnalyses::JetClusteringUtils::set_pseudoJets("
            "Reco_px_photonsRemoved, Reco_py_photonsRemoved, "
            "Reco_pz_photonsRemoved, Reco_e_photonsRemoved)"
        )

        df = df.Define(
            "jets_obj4_photonsRemoved",
            "JetClustering::clustering_ee_kt(2, 4, 0, 0)(pseudo_jets_photonsRemoved)"
        )

        df = df.Define(
            "jets_reco4_photonsRemoved",
            "FCCAnalyses::JetClusteringUtils::get_pseudoJets(jets_obj4_photonsRemoved)"
        )

        df = df.Define("n_jets_photonsRemoved", "jets_reco4_photonsRemoved.size()")

        # Require that both diagnostic jet collections exist.
        df = df.Filter(
            "n_jets_withPhotons == 4 && n_jets_photonsRemoved == 4",
            "Require four reco jets in both photon treatments"
        )

        # ------------------------------------------------------------
        # Event-level energy closure
        # ------------------------------------------------------------
        df = df.Define(
            "jetE_withPhotons",
            "getPseudoJetEnergies(jets_reco4_withPhotons)"
        )

        df = df.Define(
            "jetE_photonsRemoved",
            "getPseudoJetEnergies(jets_reco4_photonsRemoved)"
        )

        df = df.Define(
            "sumE_recoJets_withPhotons",
            "sumPseudoJetEnergy(jets_reco4_withPhotons)"
        )

        df = df.Define(
            "sumE_recoJets_photonsRemoved",
            "sumPseudoJetEnergy(jets_reco4_photonsRemoved)"
        )

        df = df.Define(
            "response_withPhotons",
            "sumE_partons > 0.0f ? sumE_recoJets_withPhotons / sumE_partons : "
            "std::numeric_limits<float>::quiet_NaN()"
        )

        df = df.Define(
            "response_photonsRemoved",
            "sumE_partons > 0.0f ? sumE_recoJets_photonsRemoved / sumE_partons : "
            "std::numeric_limits<float>::quiet_NaN()"
        )

        df = df.Define(
            "energy_closure_withPhotons",
            "sumE_partons > 0.0f ? "
            "(sumE_recoJets_withPhotons - sumE_partons) / sumE_partons : "
            "std::numeric_limits<float>::quiet_NaN()"
        )

        df = df.Define(
            "energy_closure_photonsRemoved",
            "sumE_partons > 0.0f ? "
            "(sumE_recoJets_photonsRemoved - sumE_partons) / sumE_partons : "
            "std::numeric_limits<float>::quiet_NaN()"
        )

        # ------------------------------------------------------------
        # Delta-alpha: photons kept
        # ------------------------------------------------------------
        df = df.Define(
            "parton_match_withPhotons",
            "jetPartonMatching(jets_reco4_withPhotons, parton_eta, parton_phi)"
        )

        df = df.Define(
            "delta_alpha_withPhotons",
            "getDeltaAlphaParton("
            "jets_reco4_withPhotons, parton_energies, parton_match_withPhotons)"
        )

        # ------------------------------------------------------------
        # Delta-alpha: photons removed
        # ------------------------------------------------------------
        df = df.Define(
            "parton_match_photonsRemoved",
            "jetPartonMatching(jets_reco4_photonsRemoved, parton_eta, parton_phi)"
        )

        df = df.Define(
            "delta_alpha_photonsRemoved",
            "getDeltaAlphaParton("
            "jets_reco4_photonsRemoved, parton_energies, parton_match_photonsRemoved)"
        )

        # ------------------------------------------------------------
        # Per-jet branches and ISR/no-ISR split
        # ------------------------------------------------------------
        for jet_idx in range(1, 5):
            idx = jet_idx - 1

            df = df.Define(
                f"delta_alpha_withPhotons_j{jet_idx}",
                f"getElementOrNaN(delta_alpha_withPhotons, {idx})"
            )

            df = df.Define(
                f"delta_alpha_photonsRemoved_j{jet_idx}",
                f"getElementOrNaN(delta_alpha_photonsRemoved, {idx})"
            )

            df = df.Define(
                f"partonE_j{jet_idx}",
                f"getElementOrNaN(parton_energies, {idx})"
            )

            df = df.Define(
                f"jetE_withPhotons_j{jet_idx}",
                f"getElementOrNaN(jetE_withPhotons, {idx})"
            )

            df = df.Define(
                f"jetE_photonsRemoved_j{jet_idx}",
                f"getElementOrNaN(jetE_photonsRemoved, {idx})"
            )

            df = df.Define(
                f"delta_alpha_withPhotons_j{jet_idx}_ISR",
                f"has_ISR_candidate_reco ? delta_alpha_withPhotons_j{jet_idx} : "
                f"std::numeric_limits<float>::quiet_NaN()"
            )

            df = df.Define(
                f"delta_alpha_withPhotons_j{jet_idx}_noISR",
                f"!has_ISR_candidate_reco ? delta_alpha_withPhotons_j{jet_idx} : "
                f"std::numeric_limits<float>::quiet_NaN()"
            )

            df = df.Define(
                f"delta_alpha_photonsRemoved_j{jet_idx}_ISR",
                f"has_ISR_candidate_reco ? delta_alpha_photonsRemoved_j{jet_idx} : "
                f"std::numeric_limits<float>::quiet_NaN()"
            )

            df = df.Define(
                f"delta_alpha_photonsRemoved_j{jet_idx}_noISR",
                f"!has_ISR_candidate_reco ? delta_alpha_photonsRemoved_j{jet_idx} : "
                f"std::numeric_limits<float>::quiet_NaN()"
            )

        return df

    @staticmethod
    def output():
        branches = [
            # Counts
            "n_partons",
            "n_reco_electrons",
            "n_reco_muons",
            "n_reco_photons",
            "n_mc_photons",
            "n_ISR_candidate_reco",
            "n_ISR_candidate_mc",
            "has_ISR_candidate_reco",
            "has_ISR_candidate_mc",

            # Energy sums
            "sumE_reco_all",
            "sumE_reco_electrons",
            "sumE_reco_muons",
            "sumE_reco_photons",
            "sumE_mc_final_photons",
            "sumE_ISR_candidate_reco",
            "sumE_ISR_candidate_mc",
            "photon_energy_fraction",

            # Parton / reco jet energy closure
            "sumE_partons",
            "sumE_recoJets_withPhotons",
            "sumE_recoJets_photonsRemoved",
            "response_withPhotons",
            "response_photonsRemoved",
            "energy_closure_withPhotons",
            "energy_closure_photonsRemoved",

            # Vector outputs useful for debugging
            "parton_energies",
            "jetE_withPhotons",
            "jetE_photonsRemoved",
        ]

        for jet_idx in range(1, 5):
            branches += [
                f"partonE_j{jet_idx}",
                f"jetE_withPhotons_j{jet_idx}",
                f"jetE_photonsRemoved_j{jet_idx}",

                f"delta_alpha_withPhotons_j{jet_idx}",
                f"delta_alpha_photonsRemoved_j{jet_idx}",

                f"delta_alpha_withPhotons_j{jet_idx}_ISR",
                f"delta_alpha_withPhotons_j{jet_idx}_noISR",

                f"delta_alpha_photonsRemoved_j{jet_idx}_ISR",
                f"delta_alpha_photonsRemoved_j{jet_idx}_noISR",
            ]

        return branches
