import ROOT, os
from glob import glob

ROOT.gROOT.SetBatch(True)

fractions = 1e-6

inputDir = "/eos/experiment/fcc/ee/generation/DelphesEvents/winter2023/IDEA/"

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
        "output": "angular_resolution_ecm160_v11_diagnostics"
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
    "headers/jetPartonMatching.h",
    "headers/getDeltaAlphaParton.h",
    "headers/getDeltaAlphaPartonFixed.h",
    "headers/filterValues.h",
    "headers/selectQuarks.h"
]


class RDFanalysis:

    @staticmethod
    def analysers(df):

        # ================================================================
        # Remove isolated photons/electrons/muons from reco particles
        # ================================================================
        df = df.Alias("Electron0", "Electron#0.index")
        df = df.Alias("Muon0",     "Muon#0.index")
        df = df.Alias("Photon0",   "Photon#0.index")

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

        df = df.Define(
            "RP_noPho",
            "FCCAnalyses::ReconstructedParticle::remove(ReconstructedParticles, pho_all)"
        )

        df = df.Define(
            "RP_noEle",
            "FCCAnalyses::ReconstructedParticle::remove(RP_noPho, ele_all)"
        )

        df = df.Define(
            "reco_clean",
            "FCCAnalyses::ReconstructedParticle::remove(RP_noEle, mu_all)"
        )

        # ================================================================
        # Select partons
        # ================================================================
        df = df.Define("partons_all", "selectQuarks(Particle)")
        df = df.Define("n_partons",   "partons_all.size()")
        df = df.Filter("n_partons == 4", "Require only 4 partons")

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

        # ================================================================
        # Build generator-level jets from final-state MC particles
        # ================================================================
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
            "JetClustering::clustering_ee_kt(2, 4, 0, 0)(pseudo_jets_gen)"
        )

        df = df.Define(
            "jets_gen4",
            "FCCAnalyses::JetClusteringUtils::get_pseudoJets(jets_gen_obj4)"
        )

        # ================================================================
        # Build reco-level jets from cleaned reconstructed particles
        # ================================================================
        df = df.Define(
            "Reco_px",
            "FCCAnalyses::ReconstructedParticle::get_px(reco_clean)"
        )

        df = df.Define(
            "Reco_py",
            "FCCAnalyses::ReconstructedParticle::get_py(reco_clean)"
        )

        df = df.Define(
            "Reco_pz",
            "FCCAnalyses::ReconstructedParticle::get_pz(reco_clean)"
        )

        df = df.Define(
            "Reco_e",
            "FCCAnalyses::ReconstructedParticle::get_e(reco_clean)"
        )

        df = df.Define(
            "pseudo_jets_reco",
            "FCCAnalyses::JetClusteringUtils::set_pseudoJets("
            "Reco_px, Reco_py, Reco_pz, Reco_e)"
        )

        df = df.Define(
            "jets_reco_obj4",
            "JetClustering::clustering_ee_kt(2, 4, 0, 0)(pseudo_jets_reco)"
        )

        df = df.Define(
            "jets_reco4",
            "FCCAnalyses::JetClusteringUtils::get_pseudoJets(jets_reco_obj4)"
        )

        # ================================================================
        # Optional y34-style diagnostics
        # NOTE:
        # This may depend on your FCCAnalyses version.
        # If it fails to compile, comment these 2 lines temporarily.
        # ================================================================
        df = df.Define(
            "y34_gen",
            "FCCAnalyses::JetClusteringUtils::get_exclusive_dmerge(jets_gen_obj4, 3)"
        )

        df = df.Define(
            "y34_reco",
            "FCCAnalyses::JetClusteringUtils::get_exclusive_dmerge(jets_reco_obj4, 3)"
        )

        # ================================================================
        # Require 4 gen jets and 4 reco jets
        # ================================================================
        df = df.Define("n_jets_gen",  "jets_gen4.size()")
        df = df.Define("n_jets_reco", "jets_reco4.size()")

        df = df.Filter(
            "n_jets_gen == 4 && n_jets_reco == 4",
            "Require 4 gen jets and 4 reco jets"
        )

        # ================================================================
        # Match reco jets to gen jets
        # ================================================================
        df = df.Define(
            "jet_match_indices",
            "greedyJetMatching(jets_reco4, jets_gen4, 0.1)"
        )

        df = df.Define(
            "n_matched_jets",
            "countMatchedJets(jet_match_indices)"
        )

        df = df.Filter(
            "n_matched_jets == 4",
            "Require all 4 jets matched"
        )

        # ================================================================
        # Angular and x resolutions, reco jet matched to gen jet
        # ================================================================
        df = df.Define(
            "delta_theta_matched",
            "getDeltaTheta(jets_reco4, jets_gen4, jet_match_indices)"
        )

        df = df.Define(
            "delta_phi_matched",
            "getDeltaPhi(jets_reco4, jets_gen4, jet_match_indices)"
        )

        df = df.Define(
            "delta_eta_matched",
            "getDeltaEta(jets_reco4, jets_gen4, jet_match_indices)"
        )

        df = df.Define(
            "delta_mass_matched",
            "getDeltaMass(jets_reco4, jets_gen4, jet_match_indices)"
        )

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

        # ================================================================
        # Match reco jets to partons
        # ================================================================
        df = df.Define(
            "parton_matched",
            "jetPartonMatching(jets_reco4, parton_eta, parton_phi)"
        )

        # ================================================================
        # Main alpha: reco jet vs matched parton
        # ================================================================
        df = df.Define(
            "delta_alpha",
            "getDeltaAlphaPartonFixed(jets_reco4, parton_energies, parton_matched)"
        )

        # ================================================================
        # alpha_gen = (E_reco - E_gen) / E_gen
        # detector/reco loss relative to gen jet
        # ================================================================
        df = df.Define(
            "delta_alpha_gen",
            """
            ROOT::VecOps::RVec<float> out;
            out.reserve(4);

            for (int i = 0; i < 4; ++i) {
                int g = jet_match_indices[i];

                if (g < 0 || g >= (int)jets_gen4.size()) {
                    out.push_back(-999.f);
                    continue;
                }

                float Eg = jets_gen4[g].E();
                float Er = jets_reco4[i].E();

                out.push_back(Eg > 0.f ? (Er - Eg) / Eg : -999.f);
            }

            return out;
            """
        )

        # ================================================================
        # alpha_gen_vs_parton = (E_gen_jet - E_parton) / E_parton
        # pure MC diagnostic
        # ================================================================
        df = df.Define(
            "delta_alpha_gen_vs_parton",
            """
            ROOT::VecOps::RVec<float> out;
            out.reserve(4);

            for (int i = 0; i < 4; ++i) {
                int g = jet_match_indices[i];
                int p = (i < (int)parton_matched.size()) ? parton_matched[i] : -1;

                if (g < 0 || p < 0 ||
                    g >= (int)jets_gen4.size() ||
                    p >= (int)parton_energies.size()) {
                    out.push_back(-999.f);
                    continue;
                }

                float Eg = jets_gen4[g].E();
                float Ep = parton_energies[p];

                out.push_back(Ep > 0.f ? (Eg - Ep) / Ep : -999.f);
            }

            return out;
            """
        )

        # ================================================================
        # absolute delta E = E_reco - E_parton
        # ================================================================
        df = df.Define(
            "delta_E_parton",
            """
            ROOT::VecOps::RVec<float> out;
            out.reserve(4);

            for (int i = 0; i < 4; ++i) {
                int p = (i < (int)parton_matched.size()) ? parton_matched[i] : -1;

                if (p < 0 || p >= (int)parton_energies.size()) {
                    out.push_back(-999.f);
                    continue;
                }

                out.push_back(jets_reco4[i].E() - parton_energies[p]);
            }

            return out;
            """
        )

        # ================================================================
        # is_perfect_reco = 1 if E_reco == E_gen within 0.1%
        # ================================================================
        df = df.Define(
            "is_perfect_reco",
            """
            ROOT::VecOps::RVec<int> out;
            out.reserve(4);

            for (int i = 0; i < 4; ++i) {
                int g = jet_match_indices[i];

                if (g < 0 || g >= (int)jets_gen4.size()) {
                    out.push_back(-1);
                    continue;
                }

                float Eg = jets_gen4[g].E();
                float Er = jets_reco4[i].E();

                float rel = Eg > 0.f ? std::abs(Er - Eg) / Eg : 999.f;

                out.push_back(rel < 0.001f ? 1 : 0);
            }

            return out;
            """
        )

        # ================================================================
        # Store gen/reco theta, phi, eta arrays
        # ================================================================
        df = df.Define(
            "theta_gen_all",
            "FCCAnalyses::JetClusteringUtils::get_theta(jets_gen4)"
        )

        df = df.Define(
            "theta_reco_all",
            "FCCAnalyses::JetClusteringUtils::get_theta(jets_reco4)"
        )

        df = df.Define(
            "phi_gen_all",
            "FCCAnalyses::JetClusteringUtils::get_phi(jets_gen4)"
        )

        df = df.Define(
            "phi_reco_all",
            "FCCAnalyses::JetClusteringUtils::get_phi(jets_reco4)"
        )

        df = df.Define(
            "eta_gen_all",
            "FCCAnalyses::JetClusteringUtils::get_eta(jets_gen4)"
        )

        df = df.Define(
            "eta_reco_all",
            "FCCAnalyses::JetClusteringUtils::get_eta(jets_reco4)"
        )

        # ================================================================
        # Reco and matched-gen jet energies as scalar branches
        # Also store scalar theta/phi/eta for easy ROOT plotting
        # ================================================================
        for jet_idx in range(1, 5):

            idx = jet_idx - 1

            df = df.Define(
                f"E_reco_j{jet_idx}",
                f"""
                if ((int)jets_reco4.size() <= {idx}) {{
                    return -999.0f;
                }}
                return (float)jets_reco4[{idx}].E();
                """
            )

            df = df.Define(
                f"E_gen_j{jet_idx}",
                f"""
                if ((int)jets_reco4.size() <= {idx}) {{
                    return -999.0f;
                }}

                int g = jet_match_indices[{idx}];

                if (g < 0 || g >= (int)jets_gen4.size()) {{
                    return -999.0f;
                }}

                return (float)jets_gen4[g].E();
                """
            )

            df = df.Define(
                f"theta_reco_j{jet_idx}",
                f"""
                if ((int)theta_reco_all.size() <= {idx}) return -999.0f;
                return (float)theta_reco_all[{idx}];
                """
            )

            df = df.Define(
                f"phi_reco_j{jet_idx}",
                f"""
                if ((int)phi_reco_all.size() <= {idx}) return -999.0f;
                return (float)phi_reco_all[{idx}];
                """
            )

            df = df.Define(
                f"eta_reco_j{jet_idx}",
                f"""
                if ((int)eta_reco_all.size() <= {idx}) return -999.0f;
                return (float)eta_reco_all[{idx}];
                """
            )

            df = df.Define(
                f"theta_gen_j{jet_idx}",
                f"""
                int g = jet_match_indices[{idx}];
                if (g < 0 || g >= (int)theta_gen_all.size()) return -999.0f;
                return (float)theta_gen_all[g];
                """
            )

            df = df.Define(
                f"phi_gen_j{jet_idx}",
                f"""
                int g = jet_match_indices[{idx}];
                if (g < 0 || g >= (int)phi_gen_all.size()) return -999.0f;
                return (float)phi_gen_all[g];
                """
            )

            df = df.Define(
                f"eta_gen_j{jet_idx}",
                f"""
                int g = jet_match_indices[{idx}];
                if (g < 0 || g >= (int)eta_gen_all.size()) return -999.0f;
                return (float)eta_gen_all[g];
                """
            )

        # ================================================================
        # Event-level total 4-jet energy diagnostics
        # ================================================================
        df = df.Define(
            "E_reco_total_4j",
            "E_reco_j1 + E_reco_j2 + E_reco_j3 + E_reco_j4"
        )

        df = df.Define(
            "E_gen_total_4j",
            "E_gen_j1 + E_gen_j2 + E_gen_j3 + E_gen_j4"
        )

        df = df.Define(
            "E_ratio_total_4j",
            """
            if (E_gen_total_4j <= 1e-6f) return -999.0f;
            return E_reco_total_4j / E_gen_total_4j;
            """
        )

        # ================================================================
        # Jet-pair geometry diagnostics for overlap studies
        # ================================================================
        df = df.Define(
            "dtheta_reco_j4j3",
            """
            if (theta_reco_j4 < -900.f || theta_reco_j3 < -900.f) return -999.0f;
            return theta_reco_j4 - theta_reco_j3;
            """
        )

        df = df.Define(
            "dphi_reco_j4j3",
            """
            if (phi_reco_j4 < -900.f || phi_reco_j3 < -900.f) return -999.0f;

            float dphi = phi_reco_j4 - phi_reco_j3;
            const float pi = std::acos(-1.0f);

            while (dphi > pi)  dphi -= 2.0f * pi;
            while (dphi < -pi) dphi += 2.0f * pi;

            return dphi;
            """
        )

        df = df.Define(
            "dR_reco_j4j3",
            """
            if (eta_reco_j4 < -900.f || eta_reco_j3 < -900.f ||
                phi_reco_j4 < -900.f || phi_reco_j3 < -900.f) return -999.0f;

            float dphi = phi_reco_j4 - phi_reco_j3;
            const float pi = std::acos(-1.0f);

            while (dphi > pi)  dphi -= 2.0f * pi;
            while (dphi < -pi) dphi += 2.0f * pi;

            float deta = eta_reco_j4 - eta_reco_j3;
            return std::sqrt(deta * deta + dphi * dphi);
            """
        )

        df = df.Define(
            "dtheta_reco_j4j1",
            """
            if (theta_reco_j4 < -900.f || theta_reco_j1 < -900.f) return -999.0f;
            return theta_reco_j4 - theta_reco_j1;
            """
        )

        df = df.Define(
            "dphi_reco_j4j1",
            """
            if (phi_reco_j4 < -900.f || phi_reco_j1 < -900.f) return -999.0f;

            float dphi = phi_reco_j4 - phi_reco_j1;
            const float pi = std::acos(-1.0f);

            while (dphi > pi)  dphi -= 2.0f * pi;
            while (dphi < -pi) dphi += 2.0f * pi;

            return dphi;
            """
        )

        df = df.Define(
            "dR_reco_j4j1",
            """
            if (eta_reco_j4 < -900.f || eta_reco_j1 < -900.f ||
                phi_reco_j4 < -900.f || phi_reco_j1 < -900.f) return -999.0f;

            float dphi = phi_reco_j4 - phi_reco_j1;
            const float pi = std::acos(-1.0f);

            while (dphi > pi)  dphi -= 2.0f * pi;
            while (dphi < -pi) dphi += 2.0f * pi;

            float deta = eta_reco_j4 - eta_reco_j1;
            return std::sqrt(deta * deta + dphi * dphi);
            """
        )

        df = df.Define(
            "min_dR_reco_j4_others",
            """
            if (eta_reco_j4 < -900.f || phi_reco_j4 < -900.f) return -999.0f;

            auto calcDR = [](float eta1, float phi1, float eta2, float phi2) {
                const float pi = std::acos(-1.0f);
                float dphi = phi1 - phi2;
                while (dphi > pi)  dphi -= 2.0f * pi;
                while (dphi < -pi) dphi += 2.0f * pi;
                float deta = eta1 - eta2;
                return std::sqrt(deta * deta + dphi * dphi);
            };

            float dr1 = (eta_reco_j1 < -900.f || phi_reco_j1 < -900.f) ? 999.f :
                        calcDR(eta_reco_j4, phi_reco_j4, eta_reco_j1, phi_reco_j1);

            float dr2 = (eta_reco_j2 < -900.f || phi_reco_j2 < -900.f) ? 999.f :
                        calcDR(eta_reco_j4, phi_reco_j4, eta_reco_j2, phi_reco_j2);

            float dr3 = (eta_reco_j3 < -900.f || phi_reco_j3 < -900.f) ? 999.f :
                        calcDR(eta_reco_j4, phi_reco_j4, eta_reco_j3, phi_reco_j3);

            float min12 = std::min(dr1, dr2);
            return std::min(min12, dr3);
            """
        )

        # ================================================================
        # Count reconstructed particles around each reco jet.
        # Counts reco_clean particles with DeltaR(particle, reco jet) < 0.4.
        # ================================================================
        for jet_idx in range(1, 5):

            idx = jet_idx - 1

            df = df.Define(
                f"n_reco_particles_j{jet_idx}",
                f"""
                int count = 0;

                if ((int)jets_reco4.size() <= {idx}) {{
                    return -1;
                }}

                float jet_eta = jets_reco4[{idx}].eta();
                float jet_phi = jets_reco4[{idx}].phi();

                auto px = FCCAnalyses::ReconstructedParticle::get_px(reco_clean);
                auto py = FCCAnalyses::ReconstructedParticle::get_py(reco_clean);
                auto pz = FCCAnalyses::ReconstructedParticle::get_pz(reco_clean);

                const float pi = std::acos(-1.0f);

                for (int k = 0; k < (int)px.size(); ++k) {{

                    float pt = std::sqrt(px[k] * px[k] + py[k] * py[k]);

                    if (pt <= 1e-6f) {{
                        continue;
                    }}

                    float eta = std::asinh(pz[k] / pt);
                    float phi = std::atan2(py[k], px[k]);

                    float dphi = phi - jet_phi;

                    while (dphi > pi) {{
                        dphi -= 2.0f * pi;
                    }}

                    while (dphi < -pi) {{
                        dphi += 2.0f * pi;
                    }}

                    float deta = eta - jet_eta;
                    float dR = std::sqrt(deta * deta + dphi * dphi);

                    if (dR < 0.4f) {{
                        count++;
                    }}
                }}

                return count;
                """
            )

        # ================================================================
        # Charged/neutral reconstructed energy around each reco jet.
        #
        # charged: abs(charge) > 0.5
        # neutral: abs(charge) < 0.5
        # cone: DeltaR(particle, reco jet) < 0.4
        # ================================================================
        for jet_idx in range(1, 5):

            idx = jet_idx - 1

            df = df.Define(
                f"E_charged_j{jet_idx}",
                f"""
                float E = 0.0f;

                if ((int)jets_reco4.size() <= {idx}) {{
                    return -999.0f;
                }}

                float jet_eta = jets_reco4[{idx}].eta();
                float jet_phi = jets_reco4[{idx}].phi();

                auto px = FCCAnalyses::ReconstructedParticle::get_px(reco_clean);
                auto py = FCCAnalyses::ReconstructedParticle::get_py(reco_clean);
                auto pz = FCCAnalyses::ReconstructedParticle::get_pz(reco_clean);
                auto en = FCCAnalyses::ReconstructedParticle::get_e(reco_clean);
                auto ch = FCCAnalyses::ReconstructedParticle::get_charge(reco_clean);

                const float pi = std::acos(-1.0f);

                for (int k = 0; k < (int)px.size(); ++k) {{

                    if (std::abs(ch[k]) < 0.5f) {{
                        continue;
                    }}

                    float pt = std::sqrt(px[k] * px[k] + py[k] * py[k]);

                    if (pt <= 1e-6f) {{
                        continue;
                    }}

                    float eta = std::asinh(pz[k] / pt);
                    float phi = std::atan2(py[k], px[k]);

                    float dphi = phi - jet_phi;

                    while (dphi > pi) {{
                        dphi -= 2.0f * pi;
                    }}

                    while (dphi < -pi) {{
                        dphi += 2.0f * pi;
                    }}

                    float deta = eta - jet_eta;
                    float dR = std::sqrt(deta * deta + dphi * dphi);

                    if (dR < 0.4f) {{
                        E += en[k];
                    }}
                }}

                return E;
                """
            )

            df = df.Define(
                f"E_neutral_j{jet_idx}",
                f"""
                float E = 0.0f;

                if ((int)jets_reco4.size() <= {idx}) {{
                    return -999.0f;
                }}

                float jet_eta = jets_reco4[{idx}].eta();
                float jet_phi = jets_reco4[{idx}].phi();

                auto px = FCCAnalyses::ReconstructedParticle::get_px(reco_clean);
                auto py = FCCAnalyses::ReconstructedParticle::get_py(reco_clean);
                auto pz = FCCAnalyses::ReconstructedParticle::get_pz(reco_clean);
                auto en = FCCAnalyses::ReconstructedParticle::get_e(reco_clean);
                auto ch = FCCAnalyses::ReconstructedParticle::get_charge(reco_clean);

                const float pi = std::acos(-1.0f);

                for (int k = 0; k < (int)px.size(); ++k) {{

                    if (std::abs(ch[k]) > 0.5f) {{
                        continue;
                    }}

                    float pt = std::sqrt(px[k] * px[k] + py[k] * py[k]);

                    if (pt <= 1e-6f) {{
                        continue;
                    }}

                    float eta = std::asinh(pz[k] / pt);
                    float phi = std::atan2(py[k], px[k]);

                    float dphi = phi - jet_phi;

                    while (dphi > pi) {{
                        dphi -= 2.0f * pi;
                    }}

                    while (dphi < -pi) {{
                        dphi += 2.0f * pi;
                    }}

                    float deta = eta - jet_eta;
                    float dR = std::sqrt(deta * deta + dphi * dphi);

                    if (dR < 0.4f) {{
                        E += en[k];
                    }}
                }}

                return E;
                """
            )

            df = df.Define(
                f"n_charged_j{jet_idx}",
                f"""
                int count = 0;

                if ((int)jets_reco4.size() <= {idx}) {{
                    return -1;
                }}

                float jet_eta = jets_reco4[{idx}].eta();
                float jet_phi = jets_reco4[{idx}].phi();

                auto px = FCCAnalyses::ReconstructedParticle::get_px(reco_clean);
                auto py = FCCAnalyses::ReconstructedParticle::get_py(reco_clean);
                auto pz = FCCAnalyses::ReconstructedParticle::get_pz(reco_clean);
                auto ch = FCCAnalyses::ReconstructedParticle::get_charge(reco_clean);

                const float pi = std::acos(-1.0f);

                for (int k = 0; k < (int)px.size(); ++k) {{

                    if (std::abs(ch[k]) < 0.5f) {{
                        continue;
                    }}

                    float pt = std::sqrt(px[k] * px[k] + py[k] * py[k]);

                    if (pt <= 1e-6f) {{
                        continue;
                    }}

                    float eta = std::asinh(pz[k] / pt);
                    float phi = std::atan2(py[k], px[k]);

                    float dphi = phi - jet_phi;

                    while (dphi > pi) {{
                        dphi -= 2.0f * pi;
                    }}

                    while (dphi < -pi) {{
                        dphi += 2.0f * pi;
                    }}

                    float deta = eta - jet_eta;
                    float dR = std::sqrt(deta * deta + dphi * dphi);

                    if (dR < 0.4f) {{
                        count++;
                    }}
                }}

                return count;
                """
            )

            df = df.Define(
                f"n_neutral_j{jet_idx}",
                f"""
                int count = 0;

                if ((int)jets_reco4.size() <= {idx}) {{
                    return -1;
                }}

                float jet_eta = jets_reco4[{idx}].eta();
                float jet_phi = jets_reco4[{idx}].phi();

                auto px = FCCAnalyses::ReconstructedParticle::get_px(reco_clean);
                auto py = FCCAnalyses::ReconstructedParticle::get_py(reco_clean);
                auto pz = FCCAnalyses::ReconstructedParticle::get_pz(reco_clean);
                auto ch = FCCAnalyses::ReconstructedParticle::get_charge(reco_clean);

                const float pi = std::acos(-1.0f);

                for (int k = 0; k < (int)px.size(); ++k) {{

                    if (std::abs(ch[k]) > 0.5f) {{
                        continue;
                    }}

                    float pt = std::sqrt(px[k] * px[k] + py[k] * py[k]);

                    if (pt <= 1e-6f) {{
                        continue;
                    }}

                    float eta = std::asinh(pz[k] / pt);
                    float phi = std::atan2(py[k], px[k]);

                    float dphi = phi - jet_phi;

                    while (dphi > pi) {{
                        dphi -= 2.0f * pi;
                    }}

                    while (dphi < -pi) {{
                        dphi += 2.0f * pi;
                    }}

                    float deta = eta - jet_eta;
                    float dR = std::sqrt(deta * deta + dphi * dphi);

                    if (dR < 0.4f) {{
                        count++;
                    }}
                }}

                return count;
                """
            )

            df = df.Define(
                f"neutral_fraction_clean_j{jet_idx}",
                f"""
                float Ech = E_charged_j{jet_idx};
                float Ene = E_neutral_j{jet_idx};

                if (Ech < -900.0f || Ene < -900.0f) {{
                    return -999.0f;
                }}

                float Etot = Ech + Ene;

                if (Etot <= 1e-6f) {{
                    return -999.0f;
                }}

                return Ene / Etot;
                """
            )

            df = df.Define(
                f"E_neutral_over_reco_j{jet_idx}",
                f"""
                float Ene = E_neutral_j{jet_idx};
                float Er  = E_reco_j{jet_idx};

                if (Ene < -900.0f || Er <= 1e-6f) {{
                    return -999.0f;
                }}

                return Ene / Er;
                """
            )

            df = df.Define(
                f"E_charged_over_reco_j{jet_idx}",
                f"""
                float Ech = E_charged_j{jet_idx};
                float Er  = E_reco_j{jet_idx};

                if (Ech < -900.0f || Er <= 1e-6f) {{
                    return -999.0f;
                }}

                return Ech / Er;
                """
            )

            df = df.Define(
                f"E_neutral_over_gen_j{jet_idx}",
                f"""
                float Ene = E_neutral_j{jet_idx};
                float Eg  = E_gen_j{jet_idx};

                if (Ene < -900.0f || Eg <= 1e-6f) {{
                    return -999.0f;
                }}

                return Ene / Eg;
                """
            )

            df = df.Define(
                f"E_charged_over_gen_j{jet_idx}",
                f"""
                float Ech = E_charged_j{jet_idx};
                float Eg  = E_gen_j{jet_idx};

                if (Ech < -900.0f || Eg <= 1e-6f) {{
                    return -999.0f;
                }}

                return Ech / Eg;
                """
            )

        # ================================================================
        # Extract per-jet scalar values from RVec branches
        # ================================================================
        for jet_idx in range(1, 5):

            idx = jet_idx - 1

            df = df.Define(
                f"delta_alpha_j{jet_idx}",
                f"getElement(delta_alpha, {idx})"
            )

            df = df.Define(
                f"delta_theta_j{jet_idx}",
                f"getElement(delta_theta_matched, {idx})"
            )

            df = df.Define(
                f"delta_phi_j{jet_idx}",
                f"getElement(delta_phi_matched, {idx})"
            )

            df = df.Define(
                f"delta_eta_j{jet_idx}",
                f"getElement(delta_eta_matched, {idx})"
            )

            df = df.Define(
                f"delta_x_j{jet_idx}",
                f"getElement(delta_x_matched, {idx})"
            )

            df = df.Define(
                f"delta_alpha_gen_j{jet_idx}",
                f"getElement(delta_alpha_gen, {idx})"
            )

            df = df.Define(
                f"delta_alpha_gen_vs_parton_j{jet_idx}",
                f"getElement(delta_alpha_gen_vs_parton, {idx})"
            )

            df = df.Define(
                f"delta_E_parton_j{jet_idx}",
                f"getElement(delta_E_parton, {idx})"
            )

            df = df.Define(
                f"is_perfect_reco_j{jet_idx}",
                f"getElement(is_perfect_reco, {idx})"
            )

        # ================================================================
        # Vector-level filtered quantities
        # ================================================================
        df = df.Define(
            "filtered_delta_theta_matched",
            f"filterValues(delta_theta_matched, "
            f"{filter_config['delta_theta']['min']}, "
            f"{filter_config['delta_theta']['max']})"
        )

        df = df.Define(
            "filtered_delta_phi_matched",
            f"filterValues(delta_phi_matched, "
            f"{filter_config['delta_phi']['min']}, "
            f"{filter_config['delta_phi']['max']})"
        )

        df = df.Define(
            "filtered_delta_eta_matched",
            f"filterValues(delta_eta_matched, "
            f"{filter_config['delta_eta']['min']}, "
            f"{filter_config['delta_eta']['max']})"
        )

        df = df.Define(
            "filtered_delta_x_matched",
            f"filterValues(delta_x_matched, "
            f"{filter_config['delta_x']['min']}, "
            f"{filter_config['delta_x']['max']})"
        )

        df = df.Define(
            "filtered_delta_alpha",
            f"filterValues(delta_alpha, "
            f"{filter_config['delta_alpha']['min']}, "
            f"{filter_config['delta_alpha']['max']})"
        )

        # ================================================================
        # Scalar filtered quantities per jet
        # ================================================================
        for jet_idx in range(1, 5):

            df = df.Define(
                f"filtered_delta_theta_j{jet_idx}",
                f"(delta_theta_j{jet_idx} >= {filter_config['delta_theta']['min']} && "
                f"delta_theta_j{jet_idx} <= {filter_config['delta_theta']['max']}) "
                f"? delta_theta_j{jet_idx} : -999.0f"
            )

            df = df.Define(
                f"filtered_delta_phi_j{jet_idx}",
                f"(delta_phi_j{jet_idx} >= {filter_config['delta_phi']['min']} && "
                f"delta_phi_j{jet_idx} <= {filter_config['delta_phi']['max']}) "
                f"? delta_phi_j{jet_idx} : -999.0f"
            )

            df = df.Define(
                f"filtered_delta_x_j{jet_idx}",
                f"(delta_x_j{jet_idx} >= {filter_config['delta_x']['min']} && "
                f"delta_x_j{jet_idx} <= {filter_config['delta_x']['max']}) "
                f"? delta_x_j{jet_idx} : -999.0f"
            )

            df = df.Define(
                f"filtered_delta_alpha_j{jet_idx}",
                f"(delta_alpha_j{jet_idx} >= {filter_config['delta_alpha']['min']} && "
                f"delta_alpha_j{jet_idx} <= {filter_config['delta_alpha']['max']}) "
                f"? delta_alpha_j{jet_idx} : -999.0f"
            )

        return df

    @staticmethod
    def output():

        branches = [

            # ============================================================
            # Angular and x resolutions
            # ============================================================
            "delta_theta_j1", "delta_theta_j2", "delta_theta_j3", "delta_theta_j4",
            "delta_phi_j1",   "delta_phi_j2",   "delta_phi_j3",   "delta_phi_j4",
            "delta_eta_j1",   "delta_eta_j2",   "delta_eta_j3",   "delta_eta_j4",
            "delta_x_j1",     "delta_x_j2",     "delta_x_j3",     "delta_x_j4",

            # ============================================================
            # Main alpha: reco jet vs parton
            # ============================================================
            "delta_alpha_j1", "delta_alpha_j2", "delta_alpha_j3", "delta_alpha_j4",

            # ============================================================
            # Filtered vector branches
            # ============================================================
            "filtered_delta_theta_matched",
            "filtered_delta_phi_matched",
            "filtered_delta_eta_matched",
            "filtered_delta_x_matched",
            "filtered_delta_alpha",

            # ============================================================
            # Filtered scalar branches
            # ============================================================
            "filtered_delta_theta_j1", "filtered_delta_theta_j2",
            "filtered_delta_theta_j3", "filtered_delta_theta_j4",

            "filtered_delta_phi_j1",   "filtered_delta_phi_j2",
            "filtered_delta_phi_j3",   "filtered_delta_phi_j4",

            "filtered_delta_x_j1",     "filtered_delta_x_j2",
            "filtered_delta_x_j3",     "filtered_delta_x_j4",

            "filtered_delta_alpha_j1", "filtered_delta_alpha_j2",
            "filtered_delta_alpha_j3", "filtered_delta_alpha_j4",

            # ============================================================
            # Parton kinematics
            # ============================================================
            "parton_energies",
            "parton_eta",
            "parton_phi",
            "parton_y",

            # ============================================================
            # Matching info
            # ============================================================
            "n_matched_jets",

            # ============================================================
            # Stored jet angle arrays
            # ============================================================
            "theta_gen_all",
            "theta_reco_all",
            "phi_gen_all",
            "phi_reco_all",
            "eta_gen_all",
            "eta_reco_all",

            # ============================================================
            # y34-style diagnostics
            # ============================================================
            "y34_gen",
            "y34_reco",

            # ============================================================
            # Event-level energy diagnostics
            # ============================================================
            "E_reco_total_4j",
            "E_gen_total_4j",
            "E_ratio_total_4j",

            # ============================================================
            # Jet overlap / geometry diagnostics
            # ============================================================
            "dtheta_reco_j4j3",
            "dphi_reco_j4j3",
            "dR_reco_j4j3",
            "dtheta_reco_j4j1",
            "dphi_reco_j4j1",
            "dR_reco_j4j1",
            "min_dR_reco_j4_others",
        ]

        for j in range(1, 5):

            branches += [

                # Reco jet vs gen jet
                f"delta_alpha_gen_j{j}",

                # Gen jet vs parton
                f"delta_alpha_gen_vs_parton_j{j}",

                # Absolute reco-parton energy difference
                f"delta_E_parton_j{j}",

                # Flag: reco energy approximately equal to gen energy
                f"is_perfect_reco_j{j}",

                # Reco and matched gen jet energy
                f"E_reco_j{j}",
                f"E_gen_j{j}",

                # Scalar reco/gen angle branches
                f"theta_reco_j{j}",
                f"phi_reco_j{j}",
                f"eta_reco_j{j}",
                f"theta_gen_j{j}",
                f"phi_gen_j{j}",
                f"eta_gen_j{j}",

                # Number of reco_clean particles within dR < 0.4 of reco jet
                f"n_reco_particles_j{j}",

                # Charged/neutral reco energy inside dR < 0.4 of reco jet
                f"E_charged_j{j}",
                f"E_neutral_j{j}",

                # Clean neutral fraction
                f"neutral_fraction_clean_j{j}",

                # Charged/neutral energy normalized by reco jet energy
                f"E_charged_over_reco_j{j}",
                f"E_neutral_over_reco_j{j}",

                # Charged/neutral energy normalized by matched gen jet energy
                f"E_charged_over_gen_j{j}",
                f"E_neutral_over_gen_j{j}",

                # Charged/neutral multiplicity inside dR < 0.4 of reco jet
                f"n_charged_j{j}",
                f"n_neutral_j{j}",
            ]

        return branches
