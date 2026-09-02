import ROOT, os
from glob import glob

ROOT.gROOT.SetBatch(True)

# ================================================
# STAGE 1: DATA CLEANING
# ================================================
# This stage reads the raw Delphes/EDM4hep samples and produces:
#   - reco_clean       : ReconstructedParticles with identified
#                         electrons, muons, and photons removed
#   - partons_all       : the 4 hard-process quarks (status 23)
#   - parton_energies/eta/phi/y : pre-computed parton kinematics
#   - Particle          : passthrough of the FULL truth collection,
#                         needed by stage 2 (jetclustering.py) to
#                         build MC_final (all genStatus==1 particles)
#                         for gen-jet clustering.
# ================================================

# This is the ONLY stage where `fraction` subsamples the raw dataset.
# Downstream stages (jetclustering.py, angular-resolution.py) read 100%
# of whatever this stage outputs (their own `fraction` is set to 1.0),
# since fraction is always relative to a stage's own input, not the
# original raw dataset.
fractions = 1e-6

inputDir = "/eos/experiment/fcc/ee/generation/DelphesEvents/winter2023/IDEA/"

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
    file_pattern = os.path.join(inputDir, process_name, "events_*.root")
    all_files = glob(file_pattern)

    if all_files:
        processList[f"{process_name}"] = {
            "fraction": config["fraction"],
            "chunks": 1,
            "output": f"clean_data_{config['output_suffix']}"
        }
        print(f"Added process: {process_name} ({config['label']}) - {len(all_files)} files found")
    else:
        print(f"WARNING: No files found for {process_name}")

# Output directory for stage 1 (stage 2 will read from here via inputDir)
outputDir = "outputs/stage1_clean_data/"

procDict = "FCCee_procDict_winter2023_IDEA.json"

nCPUS = 10

doTree = True

# Nazanin headers: only Stage 1 needs selectQuarks.h for parton selection.
# Do not also inline-declare selectQuarks here, otherwise ROOT/Cling can see
# a duplicate C++ definition if headers/selectQuarks.h is loaded.
includePaths = [
    "headers/selectQuarks.h"
]













# ------------------------------------------------
# RDFanalysis
# ------------------------------------------------
class RDFanalysis:
    @staticmethod
    def analysers(df):

        # ------------------------------------------------
        # LEPTON / PHOTON REMOVAL
        # ------------------------------------------------
        df = df.Alias("Electron0", "Electron#0.index")
        df = df.Alias("Muon0", "Muon#0.index")
        df = df.Alias("Photon0", "Photon#0.index")
        df = df.Define("ele_all", "FCCAnalyses::ReconstructedParticle::get(Electron0, ReconstructedParticles)")
        df = df.Define("mu_all", "FCCAnalyses::ReconstructedParticle::get(Muon0, ReconstructedParticles)")
        df = df.Define("pho_all", "FCCAnalyses::ReconstructedParticle::get(Photon0, ReconstructedParticles)")
        df = df.Define("RP_noPho", "ReconstructedParticles")
        df = df.Define("RP_noEle", "FCCAnalyses::ReconstructedParticle::remove(RP_noPho, ele_all)")
        df = df.Define("reco_clean", "FCCAnalyses::ReconstructedParticle::remove(RP_noEle, mu_all)")

        # ------------------------------------------------
        # PARTON SELECTION
        # ------------------------------------------------
        # Select quarks (PDG ID 1-6) with status 2 or 3 (partons after PS, before hadronization)
        df = df.Define("partons_all", "selectQuarks(Particle)")
        df = df.Define("n_partons", "partons_all.size()")

        # Get the 4 partons from WW -> 4 quarks
        df = df.Filter("n_partons == 4", "Require only 4 partons")

        # Calculate parton kinematics using FCCAnalyses::MCParticle functions
        df = df.Define("parton_energies", "FCCAnalyses::MCParticle::get_e(partons_all)")
        df = df.Define("parton_eta", "FCCAnalyses::MCParticle::get_eta(partons_all)")
        df = df.Define("parton_phi", "FCCAnalyses::MCParticle::get_phi(partons_all)")
        df = df.Define("parton_y", "FCCAnalyses::MCParticle::get_y(partons_all)")

        return df

    @staticmethod
    def output():
        outputs = [
            # Cleaned reco particle collection (full EDM4hep ReconstructedParticleData,
            # consumed directly by stage 2's jet clustering)
            "reco_clean",

            # Full truth collection, passed through untouched so stage 2 can
            # build MC_final (genStatus==1) for gen-jet clustering
            "Particle",

            # Parton kinematics
            "parton_energies", "parton_eta", "parton_phi", "parton_y",
            "n_partons",
        ]
        return outputs