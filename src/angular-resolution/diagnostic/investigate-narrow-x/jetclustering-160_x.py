import ROOT, os

ROOT.gROOT.SetBatch(True)

# ================================================
# STAGE 2: JET CLUSTERING
# ================================================
# This stage reads the stage-1 (clean-data.py) output and produces
# 4 gen jets and 4 reco jets per event.
#
# IMPORTANT TECHNICAL NOTE:
# fastjet::PseudoJet has no ROOT dictionary, so an RVec<PseudoJet>
# cannot be written to a TTree branch with Snapshot. Instead we
# snapshot the *kinematic components* of each jet (px, py, pz, E)
# as plain RVec<float> branches. Stage 3 (angular-resolution.py)
# rebuilds the PseudoJet objects from these floats using
# set_pseudoJets(), exactly as this stage does when first building
# pseudo_jets_gen / pseudo_jets_reco from MC_final_visible_had / reco_clean.
#
# GEN/RECO SYMMETRY (fix applied here):
# reco_clean (from stage 1) = ReconstructedParticles with identified
#   electrons and muons removed; photons are KEPT.
# Previously, MC_final (= all status-1 truth particles) had NO
# removal at all, so gen jets additionally contained neutrinos
# (invisible to any real detector) and any un-removed leptons from
# semileptonic heavy-flavor decays inside a jet cone -- an asymmetry
# against reco_clean that biases Delta_alpha / Delta_x in a
# flavor-dependent way. MC_final_visible_had now vetoes e+/-, mu+/-,
# and all three neutrino flavors (by absolute PDG ID, so both
# particle and antiparticle are removed), while explicitly keeping
# photons -- mirroring reco_clean's treatment on the truth side.
# ================================================

# NOTE: fraction here is relative to stage 1's OUTPUT, not the original
# raw Delphes dataset. Stage 1 already applied its own fraction to subsample
# the full dataset, so this stage should process all of what stage 1 produced.
fractions = 1.0

# Input directory: this points at stage 1's output, NOT the raw Delphes samples.
inputDir = "outputs/stage1_clean_data/"

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
    # Stage 1 output is named clean_data_<suffix>; point each process entry at it.
    processList[f"clean_data_{config['output_suffix']}"] = {
        "fraction": config["fraction"],
        "chunks": 1,
        "output": f"jetclustering_{config['output_suffix']}"
    }

outputDir = "outputs/stage2_jetclustering/"

procDict = "FCCee_procDict_winter2023_IDEA.json"

nCPUS = 10

doTree = True

# ------------------------------------------------
# Declare C++ functions needed for this stage
# ------------------------------------------------

# Veto e+/-, mu+/-, and nu_e/nu_mu/nu_tau (+ antiparticles) from the
# stable truth collection, keeping everything else (charged/neutral
# hadrons, photons). This mirrors reco_clean's removal of identified
# electrons and muons on the reconstructed side, while correctly
# excluding neutrinos -- which have no reco-side counterpart at all --
# from the gen-jet truth reference.
ROOT.gInterpreter.Declare(r"""
#include <cstdlib>

ROOT::VecOps::RVec<edm4hep::MCParticleData> selectVisibleHadronic(
    const ROOT::VecOps::RVec<edm4hep::MCParticleData>& mc_final
) {
    ROOT::VecOps::RVec<edm4hep::MCParticleData> result;
    result.reserve(mc_final.size());

    for (const auto& p : mc_final) {
        int abspdg = std::abs(p.PDG);
        // Veto: e- /e+ (11), mu-/mu+ (13), nu_e (12), nu_mu (14), nu_tau (16)
        if (abspdg == 11 || abspdg == 13 || abspdg == 12 || abspdg == 14 || abspdg == 16) {
            continue;
        }
        result.push_back(p);
    }
    return result;
}
""")


# ------------------------------------------------
# RDFanalysis
# ------------------------------------------------
class RDFanalysis:
    @staticmethod
    def analysers(df):

        # ------------------------------------------------
        # GEN JETS (4 jets)
        # Stable truth particles, with e/mu/neutrinos removed
        # (photons kept) -- symmetric with reco_clean below.
        # ------------------------------------------------
        df = df.Define("MC_final", "FCCAnalyses::MCParticle::sel_genStatus(1)(Particle)")
        df = df.Define("MC_final_visible_had", "selectVisibleHadronic(MC_final)")

        df = df.Define("Particle_px", "FCCAnalyses::MCParticle::get_px(MC_final_visible_had)")
        df = df.Define("Particle_py", "FCCAnalyses::MCParticle::get_py(MC_final_visible_had)")
        df = df.Define("Particle_pz", "FCCAnalyses::MCParticle::get_pz(MC_final_visible_had)")
        df = df.Define("Particle_e",  "FCCAnalyses::MCParticle::get_e(MC_final_visible_had)")

        df = df.Define("pseudo_jets_gen",
                       "FCCAnalyses::JetClusteringUtils::set_pseudoJets(Particle_px, Particle_py, Particle_pz, Particle_e)")

        df = df.Define("jets_gen_obj4",
                       "JetClustering::clustering_ee_kt(2, 4, 0, 0)(pseudo_jets_gen)")

        df = df.Define("jets_gen4",
                       "FCCAnalyses::JetClusteringUtils::get_pseudoJets(jets_gen_obj4)")

        # ------------------------------------------------
        # RECO JETS (4 jets)
        # reco_clean (from stage 1) = ReconstructedParticles with
        # identified electrons and muons removed; photons are kept.
        # ------------------------------------------------
        df = df.Define("Reco_px", "FCCAnalyses::ReconstructedParticle::get_px(reco_clean)")
        df = df.Define("Reco_py", "FCCAnalyses::ReconstructedParticle::get_py(reco_clean)")
        df = df.Define("Reco_pz", "FCCAnalyses::ReconstructedParticle::get_pz(reco_clean)")
        df = df.Define("Reco_e", "FCCAnalyses::ReconstructedParticle::get_e(reco_clean)")

        df = df.Define("pseudo_jets_reco",
                       "FCCAnalyses::JetClusteringUtils::set_pseudoJets(Reco_px, Reco_py, Reco_pz, Reco_e)")
        df = df.Define("jets_reco_obj4",
                       "JetClustering::clustering_ee_kt(2, 4, 0, 0)(pseudo_jets_reco)")
        df = df.Define("jets_reco4",
                       "FCCAnalyses::JetClusteringUtils::get_pseudoJets(jets_reco_obj4)")

        # ------------------------------------------------
        # Require exactly 4 jets on both sides
        # ------------------------------------------------
        df = df.Define("n_jets_gen", "jets_gen4.size()")
        df = df.Define("n_jets_reco", "jets_reco4.size()")
        df = df.Filter("n_jets_gen == 4 && n_jets_reco == 4", "Require four gen jets and four reco jets")

        # ------------------------------------------------
        # Convert PseudoJet kinematics to plain float branches
        # (PseudoJet itself is not snapshot-able to a TTree)
        # ------------------------------------------------
        df = df.Define("jet_gen_px", "FCCAnalyses::JetClusteringUtils::get_px(jets_gen4)")
        df = df.Define("jet_gen_py", "FCCAnalyses::JetClusteringUtils::get_py(jets_gen4)")
        df = df.Define("jet_gen_pz", "FCCAnalyses::JetClusteringUtils::get_pz(jets_gen4)")
        df = df.Define("jet_gen_e",  "FCCAnalyses::JetClusteringUtils::get_e(jets_gen4)")

        df = df.Define("jet_reco_px", "FCCAnalyses::JetClusteringUtils::get_px(jets_reco4)")
        df = df.Define("jet_reco_py", "FCCAnalyses::JetClusteringUtils::get_py(jets_reco4)")
        df = df.Define("jet_reco_pz", "FCCAnalyses::JetClusteringUtils::get_pz(jets_reco4)")
        df = df.Define("jet_reco_e",  "FCCAnalyses::JetClusteringUtils::get_e(jets_reco4)")

        return df

    @staticmethod
    def output():
        outputs = [
            # Gen jet kinematics (rebuild PseudoJets from these in stage 3)
            "jet_gen_px", "jet_gen_py", "jet_gen_pz", "jet_gen_e",

            # Reco jet kinematics
            "jet_reco_px", "jet_reco_py", "jet_reco_pz", "jet_reco_e",

            # Jet counts (sanity check / debugging)
            "n_jets_gen", "n_jets_reco",

            # Passthrough: parton kinematics needed for jet-parton matching in stage 3
            "parton_energies", "parton_eta", "parton_phi", "parton_y",
        ]
        return outputs
        