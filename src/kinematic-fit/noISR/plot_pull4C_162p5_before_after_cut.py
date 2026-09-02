#!/usr/bin/env python3
import os
import csv
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)
ROOT.gStyle.SetTitleFont(42, "XYZ")
ROOT.gStyle.SetLabelFont(42, "XYZ")
ROOT.gStyle.SetTitleSize(0.045, "XYZ")
ROOT.gStyle.SetLabelSize(0.040, "XYZ")

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
INPUT_CANDIDATES = [
    os.path.join(SCRIPT_DIR, "outputs", "wmass", "wmass_fit_pull4c_ecm160.root"),
    os.path.join(SCRIPT_DIR, "wmass_fit_pull4c_ecm160.root"),
]
INPUT_ROOT = next((p for p in INPUT_CANDIDATES if os.path.exists(p)), None)
if INPUT_ROOT is None:
    raise SystemExit("ERROR: cannot find input ROOT file. Checked:\n" + "\n".join(INPUT_CANDIDATES))

TREE_NAME = "events"
OUT_DIR = os.path.join(SCRIPT_DIR, "plots_pull4c_before_after_cut_ecm162p5")
os.makedirs(OUT_DIR, exist_ok=True)
SUMMARY_CSV = os.path.join(OUT_DIR, "pull4C_before_after_cut_summary.csv")

ECM_LABEL = "#sqrt{s} = 162.5 GeV"
VARIABLES = ["alpha", "theta", "phi", "x"]
NBINS, XMIN, XMAX = 100, -5.0, 5.0
FIT_MIN, FIT_MAX = -3.0, 3.0

SELECTIONS = [
    ("before", "Before P_{4C} cut", "", ROOT.kBlack, ROOT.kGray + 2, 2),
    ("after4C", "P_{4C} > 0.03", "prob_4C > 0.03", ROOT.kBlue + 1, ROOT.kRed + 1, 1),
]

def branch_exists(tree, branch):
    return bool(tree.GetBranch(branch))

def variable_latex(var):
    return {"alpha": "#alpha", "theta": "#theta", "phi": "#phi", "x": "x"}[var]

def check_required_branches(tree):
    required = ["prob_4C"] + [f"pull4C_{v}_j{j}" for j in range(1, 5) for v in VARIABLES]
    missing = [b for b in required if not branch_exists(tree, b)]
    if missing:
        raise SystemExit("ERROR: missing branch(es):\n" + "\n".join(missing))

def make_hist(tree, branch, name, cut):
    old = ROOT.gDirectory.Get(name)
    if old:
        old.Delete()
    h = ROOT.TH1F(name, "", NBINS, XMIN, XMAX)
    h.Sumw2(); h.SetStats(0)
    selection = f"TMath::Finite({branch})"
    if cut:
        selection = f"({selection}) && ({cut})"
    tree.Draw(f"{branch}>>{name}", selection, "goff")
    h = ROOT.gDirectory.Get(name)
    if not h:
        return None
    h.SetDirectory(0); h.SetStats(0)
    h._entries = h.GetEntries(); h._mean = h.GetMean(); h._rms = h.GetStdDev()
    integral = h.Integral(1, h.GetNbinsX())
    if integral > 0:
        h.Scale(1.0 / integral)
    return h

def fit_hist(h, name, color):
    f = ROOT.TF1(name, "gaus", FIT_MIN, FIT_MAX)
    f.SetLineColor(color); f.SetLineWidth(2)
    h.Fit(f, "QRM")
    return f

def stat_box(x1, y1, x2, y2, color, title, h, fit):
    b = ROOT.TPaveText(x1, y1, x2, y2, "NDC")
    b.SetFillColor(0); b.SetBorderSize(1); b.SetLineColor(color); b.SetTextColor(color)
    b.SetTextFont(42); b.SetTextAlign(12); b.SetTextSize(0.028)
    b.AddText(title); b.AddText(f"#mu = {fit.GetParameter(1):+.3f}")
    b.AddText(f"#sigma = {fit.GetParameter(2):.3f}"); b.AddText(f"N = {h._entries:.0f}")
    b.Draw(); return b

def plot_jet(tree, jet):
    c = ROOT.TCanvas(f"c_jet{jet}", "", 3200, 750)
    c.Divide(4, 1)
    rows, keep = [], []
    for pad_idx, var in enumerate(VARIABLES, 1):
        branch = f"pull4C_{var}_j{jet}"
        c.cd(pad_idx)
        p = ROOT.gPad
        p.SetLeftMargin(0.14); p.SetRightMargin(0.24); p.SetBottomMargin(0.16); p.SetTopMargin(0.10); p.SetTicks(1, 1)
        hs, fs = [], []
        for key, label, cut, color, fit_color, line_style in SELECTIONS:
            h = make_hist(tree, branch, f"h_{key}_{branch}", cut)
            if not h or h._entries == 0:
                raise RuntimeError(f"Empty histogram: {branch}, {key}")
            h.SetLineColor(color); h.SetLineWidth(2); h.SetLineStyle(line_style); h.SetFillStyle(0)
            fit = fit_hist(h, f"fit_{key}_{branch}", fit_color)
            hs.append(h); fs.append(fit)
            rows.append({
                "energy_GeV": 162.5, "jet": jet, "variable": var, "branch": branch,
                "selection": key, "cut": cut or "none", "entries": int(h._entries),
                "hist_mean": h._mean, "hist_rms": h._rms,
                "fit_mu": fit.GetParameter(1), "fit_sigma": fit.GetParameter(2),
                "fit_chi2": fit.GetChisquare(), "fit_ndf": fit.GetNDF(),
            })
        ymax = max(h.GetMaximum() for h in hs)
        label = variable_latex(var)
        hs[0].SetTitle(f"Jet {jet}: {label}"); hs[0].SetMinimum(0); hs[0].SetMaximum(1.35 * ymax)
        hs[0].GetXaxis().SetTitle(f"pull_{{4C}}({label})"); hs[0].GetYaxis().SetTitle("Normalized entries")
        hs[0].GetXaxis().CenterTitle(); hs[0].GetYaxis().CenterTitle()
        hs[0].GetXaxis().SetTitleSize(0.055); hs[0].GetYaxis().SetTitleSize(0.055)
        hs[0].GetXaxis().SetLabelSize(0.045); hs[0].GetYaxis().SetLabelSize(0.045)
        hs[0].Draw("HIST"); hs[1].Draw("HIST SAME"); fs[0].Draw("SAME"); fs[1].Draw("SAME")
        zero = ROOT.TLine(0, 0, 0, 1.15 * ymax); zero.SetLineColor(ROOT.kGray + 1); zero.SetLineStyle(3); zero.SetLineWidth(2); zero.Draw("SAME")
        leg = ROOT.TLegend(0.18, 0.69, 0.72, 0.82); leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextFont(42); leg.SetTextSize(0.029)
        leg.AddEntry(hs[0], SELECTIONS[0][1], "l"); leg.AddEntry(hs[1], SELECTIONS[1][1], "l"); leg.Draw()
        b0 = stat_box(0.77, 0.69, 0.98, 0.89, SELECTIONS[0][3], "Before cut", hs[0], fs[0])
        b1 = stat_box(0.77, 0.45, 0.98, 0.65, SELECTIONS[1][3], "After 4C cut", hs[1], fs[1])
        latex = ROOT.TLatex(); latex.SetNDC(); latex.SetTextFont(42); latex.SetTextSize(0.040); latex.DrawLatex(0.19, 0.88, ECM_LABEL)
        keep.extend(hs + fs + [zero, leg, b0, b1, latex])
    png = os.path.join(OUT_DIR, f"pull4C_jet{jet}_before_after_P4C_cut.png")
    pdf = os.path.join(OUT_DIR, f"pull4C_jet{jet}_before_after_P4C_cut.pdf")
    c.Update(); c.SaveAs(png); c.SaveAs(pdf); c.Close()
    for r in rows:
        r["png"], r["pdf"] = png, pdf
    print(f"Saved: {png}"); print(f"Saved: {pdf}")
    return rows

def main():
    print(f"Input ROOT: {INPUT_ROOT}")
    f = ROOT.TFile.Open(INPUT_ROOT, "READ")
    if not f or f.IsZombie():
        raise SystemExit(f"ERROR: cannot open ROOT file: {INPUT_ROOT}")
    tree = f.Get(TREE_NAME)
    if not tree:
        raise SystemExit(f"ERROR: cannot find tree '{TREE_NAME}'")
    check_required_branches(tree)
    print(f"All events: {tree.GetEntries()}")
    print(f"P4C > 0.03: {tree.GetEntries('prob_4C > 0.03')}")
    rows = []
    for jet in range(1, 5):
        rows.extend(plot_jet(tree, jet))
    f.Close()
    with open(SUMMARY_CSV, "w", newline="") as fp:
        w = csv.DictWriter(fp, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    print(f"Summary CSV saved: {SUMMARY_CSV}")

if __name__ == "__main__":
    main()
