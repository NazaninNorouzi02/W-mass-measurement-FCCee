#!/usr/bin/env python3

import os
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(1110)
ROOT.gStyle.SetPadRightMargin(0.15)

root_file = (
    "outputs/angular_resolution_energy_diagnostics/"
    "ww160_deltaalpha_energy_photonISR_diagnostics.root"
)

out_dir = "outputs/angular_resolution_energy_diagnostics/plots/"
os.makedirs(out_dir, exist_ok=True)


def safe_name(name):
    out = name
    for ch in [" ", "/", "(", ")", "[", "]", "{", "}", "*", "+", "-", ".", ",", ":", ";"]:
        out = out.replace(ch, "_")
    return out


def finite_cut(branch):
    return f"({branch} == {branch})"


def combine_cuts(*cuts):
    valid = [c for c in cuts if c and c.strip()]
    if not valid:
        return ""
    return " && ".join([f"({c})" for c in valid])


def save_canvas(canvas, name):
    png = os.path.join(out_dir, f"{safe_name(name)}.png")
    pdf = os.path.join(out_dir, f"{safe_name(name)}.pdf")
    canvas.SaveAs(png)
    canvas.SaveAs(pdf)
    print(f"saved: {png}")
    print(f"saved: {pdf}")


def draw_hist(tree, branch, bins, xmin, xmax, title, xlabel, cut=""):
    cname = f"c_{safe_name(branch)}"
    hname = f"h_{safe_name(branch)}"

    canvas = ROOT.TCanvas(cname, cname, 1100, 800)
    hist = ROOT.TH1F(hname, title, bins, xmin, xmax)
    hist.Sumw2()

    selection = combine_cuts(finite_cut(branch), cut)
    tree.Draw(f"{branch}>>{hname}", selection, "goff")

    if hist.GetEntries() == 0:
        print(f"[WARN] no entries for {branch}")
        return

    hist.SetLineWidth(2)
    hist.GetXaxis().SetTitle(xlabel)
    hist.GetYaxis().SetTitle("N_{entries}")
    hist.Draw("HIST")

    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.030)
    latex.DrawLatex(0.18, 0.84, f"Entries = {hist.GetEntries():.0f}")
    latex.DrawLatex(0.18, 0.79, f"Mean = {hist.GetMean():.5f}")
    latex.DrawLatex(0.18, 0.74, f"RMS = {hist.GetRMS():.5f}")

    save_canvas(canvas, branch)


def draw_overlay(tree, branch_a, branch_b, bins, xmin, xmax, title, xlabel, label_a, label_b):
    cname = f"c_overlay_{safe_name(branch_a)}_{safe_name(branch_b)}"
    hname_a = f"h_{safe_name(branch_a)}"
    hname_b = f"h_{safe_name(branch_b)}"

    canvas = ROOT.TCanvas(cname, cname, 1100, 800)

    hist_a = ROOT.TH1F(hname_a, title, bins, xmin, xmax)
    hist_b = ROOT.TH1F(hname_b, title, bins, xmin, xmax)

    hist_a.Sumw2()
    hist_b.Sumw2()

    tree.Draw(f"{branch_a}>>{hname_a}", finite_cut(branch_a), "goff")
    tree.Draw(f"{branch_b}>>{hname_b}", finite_cut(branch_b), "goff")

    if hist_a.GetEntries() == 0 and hist_b.GetEntries() == 0:
        print(f"[WARN] no entries for overlay {branch_a}, {branch_b}")
        return

    hist_a.SetLineColor(ROOT.kBlue + 1)
    hist_b.SetLineColor(ROOT.kRed + 1)
    hist_a.SetLineWidth(2)
    hist_b.SetLineWidth(2)

    max_y = max(hist_a.GetMaximum(), hist_b.GetMaximum())
    hist_a.SetMaximum(1.25 * max_y)

    hist_a.GetXaxis().SetTitle(xlabel)
    hist_a.GetYaxis().SetTitle("N_{entries}")

    hist_a.Draw("HIST")
    hist_b.Draw("HIST SAME")

    legend = ROOT.TLegend(0.58, 0.72, 0.88, 0.88)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.AddEntry(hist_a, label_a, "l")
    legend.AddEntry(hist_b, label_b, "l")
    legend.Draw()

    save_canvas(canvas, f"overlay_{branch_a}_vs_{branch_b}")


def draw_2d(tree, xbranch, ybranch, xbins, xmin, xmax, ybins, ymin, ymax, title, xlabel, ylabel):
    cname = f"c2_{safe_name(ybranch)}_vs_{safe_name(xbranch)}"
    hname = f"h2_{safe_name(ybranch)}_vs_{safe_name(xbranch)}"

    canvas = ROOT.TCanvas(cname, cname, 1200, 900)
    hist = ROOT.TH2F(hname, title, xbins, xmin, xmax, ybins, ymin, ymax)

    selection = combine_cuts(finite_cut(xbranch), finite_cut(ybranch))

    tree.Draw(f"{ybranch}:{xbranch}>>{hname}", selection, "COLZ")

    if hist.GetEntries() == 0:
        print(f"[WARN] no entries for 2D plot {ybranch} vs {xbranch}")
        return

    hist.GetXaxis().SetTitle(xlabel)
    hist.GetYaxis().SetTitle(ylabel)
    hist.Draw("COLZ")

    save_canvas(canvas, f"{ybranch}_vs_{xbranch}")


def print_summary(tree):
    total = tree.GetEntries()
    reco_photons = tree.GetEntries("n_reco_photons > 0")
    mc_photons = tree.GetEntries("n_mc_photons > 0")
    isr_reco = tree.GetEntries("has_ISR_candidate_reco")
    isr_mc = tree.GetEntries("has_ISR_candidate_mc")

    print("\n" + "=" * 70)
    print("DATASET SUMMARY")
    print("=" * 70)
    print(f"Total selected events:              {total}")
    print(f"Events with reco photons:           {reco_photons}")
    print(f"Events with MC final photons:       {mc_photons}")
    print(f"Events with reco ISR candidate:     {isr_reco}")
    print(f"Events with MC ISR candidate:       {isr_mc}")

    if total > 0:
        print(f"Reco photon fraction:               {reco_photons / total:.4f}")
        print(f"MC photon fraction:                 {mc_photons / total:.4f}")
        print(f"Reco ISR-candidate fraction:        {isr_reco / total:.4f}")
        print(f"MC ISR-candidate fraction:          {isr_mc / total:.4f}")

    print("=" * 70 + "\n")


def main():
    if not os.path.exists(root_file):
        raise SystemExit(f"ERROR: ROOT file not found: {root_file}")

    f = ROOT.TFile.Open(root_file, "READ")
    if not f or f.IsZombie():
        raise SystemExit(f"ERROR: cannot open ROOT file: {root_file}")

    tree = f.Get("events")
    if not tree:
        raise SystemExit("ERROR: tree 'events' not found")

    print_summary(tree)

    # ------------------------------------------------------------
    # Photon and ISR diagnostics
    # ------------------------------------------------------------
    draw_hist(
        tree,
        "n_reco_photons",
        80, 0, 160,
        "Reco photon multiplicity",
        "N_{reco photons}"
    )

    draw_hist(
        tree,
        "n_mc_photons",
        80, 0, 160,
        "MC final-state photon multiplicity",
        "N_{MC final photons}"
    )

    draw_hist(
        tree,
        "sumE_reco_photons",
        100, 0, 160,
        "Reco photon energy sum",
        "#Sigma E_{reco photons} [GeV]"
    )

    draw_hist(
        tree,
        "sumE_mc_final_photons",
        100, 0, 160,
        "MC final-state photon energy sum",
        "#Sigma E_{MC final photons} [GeV]"
    )

    draw_hist(
        tree,
        "photon_energy_fraction",
        100, 0, 1.0,
        "Reco photon energy fraction",
        "#Sigma E_{reco photons} / #Sigma E_{all reco}"
    )

    draw_hist(
        tree,
        "sumE_ISR_candidate_reco",
        100, 0, 100,
        "Reco ISR-candidate photon energy",
        "#Sigma E_{ISR candidate, reco} [GeV]"
    )

    draw_hist(
        tree,
        "sumE_ISR_candidate_mc",
        100, 0, 100,
        "MC ISR-candidate photon energy",
        "#Sigma E_{ISR candidate, MC} [GeV]"
    )

    # ------------------------------------------------------------
    # Event-level energy closure
    # ------------------------------------------------------------
    draw_hist(
        tree,
        "sumE_partons",
        100, 0, 200,
        "Parton energy sum",
        "#Sigma E_{partons} [GeV]"
    )

    draw_hist(
        tree,
        "sumE_recoJets_withPhotons",
        100, 0, 200,
        "Reco jet energy sum: photons kept",
        "#Sigma E_{reco jets} [GeV]"
    )

    draw_hist(
        tree,
        "sumE_recoJets_photonsRemoved",
        100, 0, 200,
        "Reco jet energy sum: photons removed",
        "#Sigma E_{reco jets} [GeV]"
    )

    draw_hist(
        tree,
        "response_withPhotons",
        100, 0, 1.5,
        "Event energy response: photons kept",
        "#Sigma E_{reco jets} / #Sigma E_{partons}"
    )

    draw_hist(
        tree,
        "response_photonsRemoved",
        100, 0, 1.5,
        "Event energy response: photons removed",
        "#Sigma E_{reco jets} / #Sigma E_{partons}"
    )

    draw_hist(
        tree,
        "energy_closure_withPhotons",
        100, -1.0, 0.5,
        "Energy closure: photons kept",
        "(#Sigma E_{reco jets} - #Sigma E_{partons}) / #Sigma E_{partons}"
    )

    draw_hist(
        tree,
        "energy_closure_photonsRemoved",
        100, -1.0, 0.5,
        "Energy closure: photons removed",
        "(#Sigma E_{reco jets} - #Sigma E_{partons}) / #Sigma E_{partons}"
    )

    draw_overlay(
        tree,
        "energy_closure_withPhotons",
        "energy_closure_photonsRemoved",
        100, -1.0, 0.5,
        "Energy closure comparison",
        "(#Sigma E_{reco jets} - #Sigma E_{partons}) / #Sigma E_{partons}",
        "Photons kept",
        "Photons removed"
    )

    draw_overlay(
        tree,
        "response_withPhotons",
        "response_photonsRemoved",
        100, 0, 1.5,
        "Event energy response comparison",
        "#Sigma E_{reco jets} / #Sigma E_{partons}",
        "Photons kept",
        "Photons removed"
    )

    # ------------------------------------------------------------
    # Delta-alpha per jet
    # ------------------------------------------------------------
    for jet_idx in range(1, 5):
        draw_hist(
            tree,
            f"delta_alpha_withPhotons_j{jet_idx}",
            120, -1.0, 0.5,
            f"Delta alpha j{jet_idx}: photons kept",
            f"(E_{{reco}} - E_{{parton}}) / E_{{parton}} [j{jet_idx}]"
        )

        draw_hist(
            tree,
            f"delta_alpha_photonsRemoved_j{jet_idx}",
            120, -1.0, 0.5,
            f"Delta alpha j{jet_idx}: photons removed",
            f"(E_{{reco}} - E_{{parton}}) / E_{{parton}} [j{jet_idx}]"
        )

        draw_overlay(
            tree,
            f"delta_alpha_withPhotons_j{jet_idx}",
            f"delta_alpha_photonsRemoved_j{jet_idx}",
            120, -1.0, 0.5,
            f"Delta alpha comparison j{jet_idx}",
            f"(E_{{reco}} - E_{{parton}}) / E_{{parton}} [j{jet_idx}]",
            "Photons kept",
            "Photons removed"
        )

        draw_overlay(
            tree,
            f"delta_alpha_withPhotons_j{jet_idx}_ISR",
            f"delta_alpha_withPhotons_j{jet_idx}_noISR",
            120, -1.0, 0.5,
            f"Delta alpha ISR split j{jet_idx}: photons kept",
            f"(E_{{reco}} - E_{{parton}}) / E_{{parton}} [j{jet_idx}]",
            "Reco ISR candidate",
            "No reco ISR candidate"
        )

        draw_overlay(
            tree,
            f"delta_alpha_photonsRemoved_j{jet_idx}_ISR",
            f"delta_alpha_photonsRemoved_j{jet_idx}_noISR",
            120, -1.0, 0.5,
            f"Delta alpha ISR split j{jet_idx}: photons removed",
            f"(E_{{reco}} - E_{{parton}}) / E_{{parton}} [j{jet_idx}]",
            "Reco ISR candidate",
            "No reco ISR candidate"
        )

    # ------------------------------------------------------------
    # 2D diagnostic plots
    # ------------------------------------------------------------
    draw_2d(
        tree,
        "sumE_partons",
        "sumE_recoJets_withPhotons",
        100, 0, 200,
        100, 0, 200,
        "Reco jet energy vs parton energy: photons kept",
        "#Sigma E_{partons} [GeV]",
        "#Sigma E_{reco jets} [GeV]"
    )

    draw_2d(
        tree,
        "sumE_partons",
        "sumE_recoJets_photonsRemoved",
        100, 0, 200,
        100, 0, 200,
        "Reco jet energy vs parton energy: photons removed",
        "#Sigma E_{partons} [GeV]",
        "#Sigma E_{reco jets} [GeV]"
    )

    draw_2d(
        tree,
        "photon_energy_fraction",
        "energy_closure_withPhotons",
        100, 0, 1.0,
        100, -1.0, 0.5,
        "Energy closure vs photon fraction: photons kept",
        "#Sigma E_{reco photons} / #Sigma E_{all reco}",
        "Energy closure"
    )

    draw_2d(
        tree,
        "photon_energy_fraction",
        "energy_closure_photonsRemoved",
        100, 0, 1.0,
        100, -1.0, 0.5,
        "Energy closure vs photon fraction: photons removed",
        "#Sigma E_{reco photons} / #Sigma E_{all reco}",
        "Energy closure"
    )

    draw_2d(
        tree,
        "sumE_reco_photons",
        "delta_alpha_withPhotons_j4",
        100, 0, 160,
        100, -1.0, 0.5,
        "Delta alpha j4 vs reco photon energy: photons kept",
        "#Sigma E_{reco photons} [GeV]",
        "#delta_{#alpha} j4"
    )

    draw_2d(
        tree,
        "sumE_reco_photons",
        "delta_alpha_photonsRemoved_j4",
        100, 0, 160,
        100, -1.0, 0.5,
        "Delta alpha j4 vs reco photon energy: photons removed",
        "#Sigma E_{reco photons} [GeV]",
        "#delta_{#alpha} j4"
    )

    draw_2d(
        tree,
        "photon_energy_fraction",
        "delta_alpha_withPhotons_j4",
        100, 0, 1.0,
        100, -1.0, 0.5,
        "Delta alpha j4 vs photon fraction: photons kept",
        "#Sigma E_{reco photons} / #Sigma E_{all reco}",
        "#delta_{#alpha} j4"
    )

    draw_2d(
        tree,
        "photon_energy_fraction",
        "delta_alpha_photonsRemoved_j4",
        100, 0, 1.0,
        100, -1.0, 0.5,
        "Delta alpha j4 vs photon fraction: photons removed",
        "#Sigma E_{reco photons} / #Sigma E_{all reco}",
        "#delta_{#alpha} j4"
    )

    draw_2d(
        tree,
        "sumE_ISR_candidate_reco",
        "delta_alpha_withPhotons_j4",
        100, 0, 100,
        100, -1.0, 0.5,
        "Delta alpha j4 vs reco ISR-candidate energy: photons kept",
        "#Sigma E_{ISR candidate, reco} [GeV]",
        "#delta_{#alpha} j4"
    )

    draw_2d(
        tree,
        "sumE_ISR_candidate_reco",
        "delta_alpha_photonsRemoved_j4",
        100, 0, 100,
        100, -1.0, 0.5,
        "Delta alpha j4 vs reco ISR-candidate energy: photons removed",
        "#Sigma E_{ISR candidate, reco} [GeV]",
        "#delta_{#alpha} j4"
    )

    f.Close()


if __name__ == "__main__":
    main()