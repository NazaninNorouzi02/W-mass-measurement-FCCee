import ROOT

ROOT.gROOT.SetBatch(True)


ROOT.gInterpreter.Declare(r'''
#include "ROOT/RVec.hxx"
#include "edm4hep/MCParticleData.h"
#include <cmath>

bool is_quark(const edm4hep::MCParticleData& p){
    int apdg = std::abs(p.PDG);
    if(apdg < 1 || apdg > 6) return false;
    return true;
}

ROOT::VecOps::RVec<edm4hep::MCParticleData>
select_W_quarks(
    const ROOT::VecOps::RVec<edm4hep::MCParticleData>& particles,
    const ROOT::VecOps::RVec<int>& parent_indices
){
    ROOT::VecOps::RVec<edm4hep::MCParticleData> out;
    for(const auto& p : particles){
        if(!is_quark(p)) continue;
        if(p.generatorStatus != 23) continue;
        for(unsigned int i = p.parents_begin; i < p.parents_end; i++){
            if(i >= parent_indices.size()) continue;
            int parent_index = parent_indices[i];
            if(parent_index < 0 || parent_index >= (int)particles.size()) continue;
            int parent_pdg = particles[parent_index].PDG;
            if(std::abs(parent_pdg) == 24){
                out.push_back(p);
                break;
            }
        }
    }
    return out;
}

ROOT::VecOps::RVec<int>
get_parent_pdg(
    const ROOT::VecOps::RVec<edm4hep::MCParticleData>& particles,
    const ROOT::VecOps::RVec<edm4hep::MCParticleData>& parts,
    const ROOT::VecOps::RVec<int>& parent_indices
){
    ROOT::VecOps::RVec<int> out;
    for(const auto& p : parts){
        bool found_parent = false;
        for(unsigned int i = p.parents_begin; i < p.parents_end; i++){
            if(i >= parent_indices.size()) continue;
            int parent_index = parent_indices[i];
            if(parent_index < 0 || parent_index >= (int)particles.size()) continue;
            out.push_back(particles[parent_index].PDG);
            found_parent = true;
            break;
        }
        if(!found_parent) out.push_back(0);
    }
    return out;
}

int get_int_at(const ROOT::VecOps::RVec<int>& values, int index){
    if(index < 0 || index >= (int)values.size()) return 0;
    return values[index];
}
''')

# ======================================================
# process configuration
# ======================================================

processList = {
    "p8_ee_WW_ecm160": {"fraction": 2e-6}
}

inputDir = "/eos/experiment/fcc/ee/generation/DelphesEvents/winter2023/IDEA/"
procDict = "FCCee_procDict_winter2023_IDEA.json"

outDir = "./outputs/sanity_ww160_gen"
outputFile = "sanity_ecm160_gen.root"

# ======================================================
# Analysis
# ======================================================

class RDFanalysis:

    def analysers(df):

    
        df = df.Alias("Particle_parent_index", "Particle#0.index")

        df = (
            df
            # ==================================================
            # GEN PARTICLES SELECTION (Final State Particles, Status 1)
            # ==================================================
         
            .Define("StableParticles", "MCParticle::sel_genStatus(1)(Particle)")
            
            .Define("GenP_px", "MCParticle::get_px(StableParticles)")
            .Define("GenP_py", "MCParticle::get_py(StableParticles)")
            .Define("GenP_pz", "MCParticle::get_pz(StableParticles)")
            .Define("GenP_e",  "MCParticle::get_e(StableParticles)")
            .Define("GenP_m",  "MCParticle::get_mass(StableParticles)")

            # ==================================================
            # TOTAL VISIBLE ENERGY (GEN LEVEL)
            # ==================================================
            .Define("visible_E", "Sum(GenP_e)")

            # ==================================================
            # TRUTH-LEVEL W-QUARK SELECTION (STATUS 23)
            # ==================================================
            .Define("HardQuarks_W", "select_W_quarks(Particle, Particle_parent_index)")
            .Define("hard_quark_pdg", "MCParticle::get_pdg(HardQuarks_W)")
            .Define("hard_quark_parent_pdg", "get_parent_pdg(Particle, HardQuarks_W, Particle_parent_index)")
            .Define("n_hard_quarks_W", "HardQuarks_W.size()")

         
            .Define("q1_pdg", "get_int_at(hard_quark_pdg, 0)")
            .Define("q2_pdg", "get_int_at(hard_quark_pdg, 1)")
            .Define("q3_pdg", "get_int_at(hard_quark_pdg, 2)")
            .Define("q4_pdg", "get_int_at(hard_quark_pdg, 3)")

            .Define("q1_parent_pdg", "get_int_at(hard_quark_parent_pdg, 0)")
            .Define("q2_parent_pdg", "get_int_at(hard_quark_parent_pdg, 1)")
            .Define("q3_parent_pdg", "get_int_at(hard_quark_parent_pdg, 2)")
            .Define("q4_parent_pdg", "get_int_at(hard_quark_parent_pdg, 3)")

            # ==================================================
            # GEN-LEVEL JET CLUSTERING
            # ==================================================
            .Define("pseudo_jets", "JetClusteringUtils::set_pseudoJets_xyzm(GenP_px, GenP_py, GenP_pz, GenP_m)")
            .Define("Jets", "JetClustering::clustering_ee_kt(2, 4, 0, 0)(pseudo_jets)")
            .Define("jets", "JetClusteringUtils::get_pseudoJets(Jets)")

            .Define("jet_px", "JetClusteringUtils::get_px(jets)")
            .Define("jet_py", "JetClusteringUtils::get_py(jets)")
            .Define("jet_pz", "JetClusteringUtils::get_pz(jets)")
            .Define("jet_e",  "JetClusteringUtils::get_e(jets)")

            # ==================================================
            # REQUIRE EXACTLY 4 JETS
            # ==================================================
            .Filter("jet_px.size() == 4")

            # ==================================================
            # GEN-LEVEL PHYSICS OBSERVABLES
            # ==================================================
            .Define("jet_E_sum", "Sum(jet_e)")
            .Define("jet_system_mass", """
                sqrt(
                    pow(Sum(jet_e), 2)
                    - pow(Sum(jet_px), 2)
                    - pow(Sum(jet_py), 2)
                    - pow(Sum(jet_pz), 2)
                )
                """)

            .Define("sum_px", "Sum(jet_px)")
            .Define("sum_py", "Sum(jet_py)")
            .Define("sum_pz", "Sum(jet_pz)")
            .Define("sum_p", "sqrt(sum_px*sum_px + sum_py*sum_py + sum_pz*sum_pz)")
        )

        return df

    def output():
        return [
            "hard_quark_pdg",
            "hard_quark_parent_pdg",
            "n_hard_quarks_W",
            "q1_pdg", "q2_pdg", "q3_pdg", "q4_pdg",
            "q1_parent_pdg", "q2_parent_pdg", "q3_parent_pdg", "q4_parent_pdg",
            "visible_E",
            "jet_E_sum",
            "jet_system_mass",
            "sum_px",
            "sum_py",
            "sum_pz",
            "sum_p"
        ]

