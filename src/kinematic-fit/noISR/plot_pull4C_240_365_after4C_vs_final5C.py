#!/usr/bin/env python3
import os
import csv
import argparse
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)
ROOT.gStyle.SetTitleFont(42, "XYZ")
ROOT.gStyle.SetLabelFont(42, "XYZ")
ROOT.gStyle.SetTitleSize(0.045, "XYZ")
ROOT.gStyle.SetLabelSize(0.040, "XYZ")

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
TREE_NAME = "events"
VARIABLES = ["alpha", "theta", "phi", "x"]
NBINS, XMIN, XMAX = 100, -5.0, 5.0
FIT_MIN, FIT_MAX = -3.0, 3.0

CONFIGS = {
    "240": {"ecm": 240.0, "label": "#sqrt{s} = 240 GeV", "file": "wmass_fit_pull4c_ecm240.root"},
    "365": {"ecm": 365.0, "label": "#sqrt{s} = 365 GeV", "file": "wmass_fit_pull4c_ecm365.root"},
}

CUT_4C = "prob_4C > 0.03"
CUT_FINAL_5C = (
    "((isr_applied > 0.5 && prob_5C_isr > 0.03) || "
    "(isr_applied < 0.5 && prob_5C > 0.03))"
)
SELECTIONS = [
    ("after4C", "P_{4C} > 0.03", CUT_4C, ROOT.kRed + 1, ROOT.kRed + 2, 1),
    ("afterFinal5C", "Final 5C candidate: P > 0.03", CUT_FINAL_5C, ROOT.kBlue + 1, ROOT.kBlue + 2, 2),
]

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("energy", choices=["240", "365"])
    return p.parse_args()

def find_input(name):
    candidates = [os.path.join(SCRIPT_DIR, "outputs", "wmass", name), os.path.join(SCRIPT_DIR, name)]
    return next((p for p in candidates if os.path.exists(p)), None), candidates

def branch_exists(tree, branch):
    return bool(tree.GetBranch(branch))

def variable_latex(var):
    return {"alpha": "#alpha", "theta": "#theta", "phi": "#phi", "x": "x"}[var]

def check_required(tree):
    required = ["prob_4C", "prob_5C", "prob_5C_isr", "isr_applied"]
    required += [f"pull4C_{v}_j{j}" for j in range(1, 5) for v in VARIABLES]
    missing = [b for b in required if not branch_exists(tree, b)]
    if missing:
        raise SystemExit("ERROR: missing branch(es):\n" + "\n".join(missing))

def make_hist(tree, branch, name, cut):
    old = ROOT.gDirectory.Get(name)
    if old: old.Delete()
    h = ROOT.TH1F(name, "", NBINS, XMIN, XMAX); h.Sumw2(); h.SetStats(0)
    tree.Draw(f"{branch}>>{name}", f"TMath::Finite({branch}) && ({cut})", "goff")
    h = ROOT.gDirectory.Get(name)
    if not h: return None
    h.SetDirectory(0); h.SetStats(0)
    h._entries, h._mean, h._rms = h.GetEntries(), h.GetMean(), h.GetStdDev()
    integral = h.Integral(1, h.GetNbinsX())
    if integral > 0: h.Scale(1.0 / integral)
    return h

def fit_hist(h, name, color):
    f = ROOT.TF1(name, "gaus", FIT_MIN, FIT_MAX); f.SetLineColor(color); f.SetLineWidth(2); h.Fit(f, "QRM"); return f

def stat_box(x1, y1, x2, y2, color, title, h, fit):
    b = ROOT.TPaveText(x1, y1, x2, y2, "NDC")
    b.SetFillColor(0); b.SetBorderSize(1); b.SetLineColor(color); b.SetTextColor(color)
    b.SetTextFont(42); b.SetTextAlign(12); b.SetTextSize(0.028)
    b.AddText(title); b.AddText(f"#mu = {fit.GetParameter(1):+.3f}")
    b.AddText(f"#sigma = {fit.GetParameter(2):.3f}"); b.AddText(f"N = {h._entries:.0f}")
    b.Draw(); return b

def plot_jet(tree, jet, cfg, out_dir, energy):
    c = ROOT.TCanvas(f"c_jet{jet}_{energy}", "", 3200, 750); c.Divide(4, 1)
    rows, keep = [], []
    for pad_idx, var in enumerate(VARIABLES, 1):
        branch = f"pull4C_{var}_j{jet}"
        c.cd(pad_idx)
        p = ROOT.gPad; p.SetLeftMargin(0.14); p.SetRightMargin(0.24); p.SetBottomMargin(0.16); p.SetTopMargin(0.10); p.SetTicks(1, 1)
        hs, fs = [], []
        for key, label, cut, color, fit_color, line_style in SELECTIONS:
            h = make_hist(tree, branch, f"h_{key}_{branch}_{energy}", cut)
            if not h or h._entries == 0: raise RuntimeError(f"Empty histogram: {branch}, {key}")
            h.SetLineColor(color); h.SetLineWidth(2); h.SetLineStyle(line_style); h.SetFillStyle(0)
            fit = fit_hist(h, f"fit_{key}_{branch}_{energy}", fit_color)
            hs.append(h); fs.append(fit)
            rows.append({
                "energy_GeV": cfg["ecm"], "jet": jet, "variable": var, "branch": branch,
                "selection": key, "cut": cut, "entries": int(h._entries),
                "hist_mean": h._mean, "hist_rms": h._rms,
                "fit_mu": fit.GetParameter(1), "fit_sigma": fit.GetParameter(2),
                "fit_chi2": fit.GetChisquare(), "fit_ndf": fit.GetNDF(),
            })
        ymax = max(h.GetMaximum() for h in hs); label = variable_latex(var)
        hs[0].SetTitle(f"Jet {jet}: {label}"); hs[0].SetMinimum(0); hs[0].SetMaximum(1.35 * ymax)
        hs[0].GetXaxis().SetTitle(f"pull_{{4C}}({label})"); hs[0].GetYaxis().SetTitle("Normalized entries")
        hs[0].GetXaxis().CenterTitle(); hs[0].GetYaxis().CenterTitle()
        hs[0].GetXaxis().SetTitleSize(0.055); hs[0].GetYaxis().SetTitleSize(0.055)
        hs[0].GetXaxis().SetLabelSize(0.045); hs[0].GetYaxis().SetLabelSize(0.045)
        hs[0].Draw("HIST"); hs[1].Draw("HIST SAME"); fs[0].Draw("SAME"); fs[1].Draw("SAME")
        zero = ROOT.TLine(0, 0, 0, 1.15 * ymax); zero.SetLineColor(ROOT.kGray + 1); zero.SetLineStyle(3); zero.SetLineWidth(2); zero.Draw("SAME")
        leg = ROOT.TLegend(0.18, 0.69, 0.73, 0.82); leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextFont(42); leg.SetTextSize(0.029)
        leg.AddEntry(hs[0], SELECTIONS[0][1], "l"); leg.AddEntry(hs[1], SELECTIONS[1][1], "l"); leg.Draw()
        b0 = stat_box(0.77, 0.69, 0.98, 0.89, SELECTIONS[0][3], "After 4C cut", hs[0], fs[0])
        b1 = stat_box(0.77, 0.45, 0.98, 0.65, SELECTIONS[1][3], "After final 5C", hs[1], fs[1])
        latex = ROOT.TLatex(); latex.SetNDC(); latex.SetTextFont(42); latex.SetTextSize(0.040); latex.DrawLatex(0.19, 0.88, cfg["label"])
        keep.extend(hs + fs + [zero, leg, b0, b1, latex])
    png = os.path.join(out_dir, f"pull4C_jet{jet}_after4C_vs_final5C_ecm{energy}.png")
    pdf = os.path.join(out_dir, f"pull4C_jet{jet}_after4C_vs_final5C_ecm{energy}.pdf")
    c.Update(); c.SaveAs(png); c.SaveAs(pdf); c.Close()
    for r in rows: r["png"], r["pdf"] = png, pdf
    print(f"Saved: {png}"); print(f"Saved: {pdf}")
    return rows

def main():
    args = parse_args(); cfg = CONFIGS[args.energy]
    input_root, candidates = find_input(cfg["file"])
    if input_root is None:
        raise SystemExit("ERROR: cannot find input ROOT file. Checked:\n" + "\n".join(candidates))
    out_dir = os.path.join(SCRIPT_DIR, f"plots_pull4c_after4C_vs_final5C_ecm{args.energy}")
    os.makedirs(out_dir, exist_ok=True)
    summary = os.path.join(out_dir, f"pull4C_after4C_vs_final5C_ecm{args.energy}_summary.csv")
    print(f"Input ROOT: {input_root}")
    print("IMPORTANT: these are 4C pulls under two event selections, not separately calculated 5C pulls.")
    f = ROOT.TFile.Open(input_root, "READ")
    if not f or f.IsZombie(): raise SystemExit(f"ERROR: cannot open ROOT file: {input_root}")
    tree = f.Get(TREE_NAME)
    if not tree: raise SystemExit(f"ERROR: cannot find tree '{TREE_NAME}'")
    check_required(tree)
    print(f"All events: {tree.GetEntries()}")
    print(f"After P4C > 0.03: {tree.GetEntries(CUT_4C)}")
    print(f"After final 5C P > 0.03: {tree.GetEntries(CUT_FINAL_5C)}")
    print(f"Passing both: {tree.GetEntries(f'({CUT_4C}) && ({CUT_FINAL_5C})')}")
    rows = []
    for jet in range(1, 5): rows.extend(plot_jet(tree, jet, cfg, out_dir, args.energy))
    f.Close()
    with open(summary, "w", newline="") as fp:
        w = csv.DictWriter(fp, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    print(f"Summary CSV saved: {summary}")

if __name__ == "__main__":
    main()
