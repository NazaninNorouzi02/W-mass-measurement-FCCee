#!/usr/bin/env python3
import os
import csv
import ROOT

ROOT.gROOT.SetBatch(True)

# =========================================================
# ROOT style
# =========================================================

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)
ROOT.gStyle.SetCanvasColor(0)
ROOT.gStyle.SetPadColor(0)
ROOT.gStyle.SetFrameFillColor(0)
ROOT.gStyle.SetFrameBorderMode(0)
ROOT.gStyle.SetTitleFont(42, "XYZ")
ROOT.gStyle.SetLabelFont(42, "XYZ")
ROOT.gStyle.SetTitleSize(0.045, "XYZ")
ROOT.gStyle.SetLabelSize(0.040, "XYZ")
ROOT.gStyle.SetLegendFont(42)
ROOT.gStyle.SetLineScalePS(1)

# =========================================================
# Input / output
# =========================================================

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

INPUT_CANDIDATES = [
    os.path.join(SCRIPT_DIR, "outputs", "wmass", "wmass_fit_pvalue_ecm240.root"),
    os.path.join(SCRIPT_DIR, "wmass_fit_pvalue_ecm240.root"),
]

INPUT_ROOT = next((p for p in INPUT_CANDIDATES if os.path.exists(p)), None)

if INPUT_ROOT is None:
    raise SystemExit(
        "ERROR: cannot find input ROOT file. Checked:\n"
        + "\n".join(INPUT_CANDIDATES)
    )

TREE_NAME = "events"

OUT_DIR = os.path.join(SCRIPT_DIR, "plots_240_thesis_style_clean")
os.makedirs(OUT_DIR, exist_ok=True)

SUMMARY_CSV = os.path.join(OUT_DIR, "summary_240_thesis_style_clean.csv")

# =========================================================
# Labels
# =========================================================

ECM_TEXT = "240.0 GeV"
CHANNEL_TEXT = "WW #rightarrow qqqq"

RAW_LABEL = "Raw Mass"
FIT4C_LABEL = "4C kinematic fit"
FIT5C_LABEL = "5C candidate"

FINAL_5C_BRANCH = "mW_5C_isr"

# =========================================================
# Histogram settings
# =========================================================

NBINS = 160
XMIN = 40.0
XMAX = 120.0

# =========================================================
# Method-specific cuts
# =========================================================
#
# Raw mass:
#   no p-value cut
#
# 4C mass:
#   prob_4C > 0.03
#
# 5C candidate:
#   if ISR candidate used:     prob_5C_isr > 0.03
#   if ISR candidate not used: prob_5C     > 0.03
#
# All curves are also restricted to the visible mass range 40-120 GeV.

CUT_RAW = ""

CUT_4C = "prob_4C > 0.03"

CUT_5C_CANDIDATE = (
    "((isr_applied > 0.5 && prob_5C_isr > 0.03) || "
    "(isr_applied < 0.5 && prob_5C > 0.03))"
)

# =========================================================
# Helpers
# =========================================================

def branch_exists(tree, branch):
    return bool(tree.GetBranch(branch))


def check_branches(tree, branches):
    missing = [br for br in branches if not branch_exists(tree, br)]
    if missing:
        raise SystemExit(
            "ERROR: missing branch(es):\n"
            + "\n".join(missing)
            + "\n\nCheck your producer output list."
        )


def finite(branch):
    return f"TMath::Finite({branch})"


def mass_window(branch):
    return f"({branch} >= {XMIN} && {branch} < {XMAX})"


def combine_selection(branch, method_cut):
    sel = f"({finite(branch)}) && ({mass_window(branch)})"

    if method_cut and method_cut.strip():
        sel = f"({sel}) && ({method_cut})"

    return sel


def make_hist(tree, branch, hist_name, method_cut):
    old = ROOT.gDirectory.Get(hist_name)
    if old:
        old.Delete()

    h = ROOT.TH1F(hist_name, "", NBINS, XMIN, XMAX)
    h.Sumw2()
    h.SetStats(0)

    selection = combine_selection(branch, method_cut)
    tree.Draw(f"{branch}>>{hist_name}", selection, "goff")

    h = ROOT.gDirectory.Get(hist_name)
    if not h:
        return None

    h.SetDirectory(0)
    h.SetStats(0)

    h._entries_in_range = h.GetEntries()
    h._mean_in_range = h.GetMean()
    h._rms_in_range = h.GetStdDev()

    integral = h.Integral(1, h.GetNbinsX())
    if integral > 0:
        h.Scale(1.0 / integral)

    return h


def style_hist(h, color):
    h.SetLineColor(color)
    h.SetLineWidth(2)
    h.SetFillStyle(0)


def draw_stat_box(x1, y1, x2, y2, color, label, h):
    box = ROOT.TPaveText(x1, y1, x2, y2, "NDC")
    box.SetFillColor(0)
    box.SetBorderSize(1)
    box.SetLineColor(color)
    box.SetTextColor(color)
    box.SetTextFont(42)
    box.SetTextAlign(12)
    box.SetTextSize(0.027)

    box.AddText(label)
    box.AddText(f"Mean        {h._mean_in_range:.3g}")
    box.AddText(f"Std Dev     {h._rms_in_range:.3g}")
    box.AddText(f"N           {h._entries_in_range:.0f}")

    box.Draw()
    return box


def draw_energy_label():
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextFont(42)
    latex.SetTextSize(0.038)
    latex.SetTextAlign(22)
    latex.DrawLatex(0.73, 0.245, ECM_TEXT)
    latex.DrawLatex(0.73, 0.205, CHANNEL_TEXT)
    return latex


def plot_mass_overlay(tree, raw_branch, fit4c_branch, outname, xtitle):
    h_raw = make_hist(tree, raw_branch, f"h_raw_{outname}", CUT_RAW)
    h_4c = make_hist(tree, fit4c_branch, f"h_4c_{outname}", CUT_4C)
    h_5c = make_hist(tree, FINAL_5C_BRANCH, f"h_5c_{outname}", CUT_5C_CANDIDATE)

    if not h_raw or h_raw.GetEntries() == 0:
        print(f"[SKIP] Empty raw histogram for {outname}")
        return None

    if not h_4c or h_4c.GetEntries() == 0:
        print(f"[SKIP] Empty 4C histogram for {outname}")
        return None

    if not h_5c or h_5c.GetEntries() == 0:
        print(f"[SKIP] Empty 5C candidate histogram for {outname}")
        return None

    c = ROOT.TCanvas(f"c_{outname}", "", 1200, 850)
    c.SetLeftMargin(0.12)
    c.SetRightMargin(0.05)
    c.SetBottomMargin(0.12)
    c.SetTopMargin(0.05)
    c.SetTicks(1, 1)

    style_hist(h_raw, ROOT.kBlack)
    style_hist(h_4c, ROOT.kRed + 1)
    style_hist(h_5c, ROOT.kBlue + 1)

    ymax = max(h_raw.GetMaximum(), h_4c.GetMaximum(), h_5c.GetMaximum())

    h_raw.SetMinimum(0.0)
    h_raw.SetMaximum(ymax * 1.18)
    h_raw.SetTitle("")
    h_raw.GetXaxis().SetTitle(xtitle)
    h_raw.GetYaxis().SetTitle("Normalized entries")
    h_raw.GetXaxis().CenterTitle()
    h_raw.GetYaxis().CenterTitle()
    h_raw.GetXaxis().SetTitleSize(0.045)
    h_raw.GetYaxis().SetTitleSize(0.045)
    h_raw.GetXaxis().SetLabelSize(0.038)
    h_raw.GetYaxis().SetLabelSize(0.038)
    h_raw.GetXaxis().SetTitleOffset(1.12)
    h_raw.GetYaxis().SetTitleOffset(1.22)
    h_raw.GetXaxis().SetNdivisions(510)
    h_raw.GetYaxis().SetNdivisions(510)

    h_raw.Draw("HIST")
    h_4c.Draw("HIST SAME")
    h_5c.Draw("HIST SAME")

    leg = ROOT.TLegend(0.15, 0.80, 0.42, 0.93)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.032)
    leg.AddEntry(h_5c, FIT5C_LABEL, "l")
    leg.AddEntry(h_4c, FIT4C_LABEL, "l")
    leg.AddEntry(h_raw, RAW_LABEL, "l")
    leg.Draw()

    box_5c = draw_stat_box(0.78, 0.785, 0.95, 0.93, ROOT.kBlue + 1, "5C", h_5c)
    box_4c = draw_stat_box(0.78, 0.625, 0.95, 0.77, ROOT.kRed + 1, "4C", h_4c)
    box_raw = draw_stat_box(0.78, 0.465, 0.95, 0.61, ROOT.kBlack, "Raw", h_raw)

    label = draw_energy_label()

    keepalive = [h_raw, h_4c, h_5c, leg, box_5c, box_4c, box_raw, label]

    out_png = os.path.join(OUT_DIR, f"{outname}.png")
    out_pdf = os.path.join(OUT_DIR, f"{outname}.pdf")

    c.SaveAs(out_png)
    c.SaveAs(out_pdf)
    c.Close()

    row = {
        "plot": outname,

        "raw_branch": raw_branch,
        "raw_cut": "none",
        "raw_entries_40_120": int(h_raw._entries_in_range),
        "raw_mean_40_120": h_raw._mean_in_range,
        "raw_rms_40_120": h_raw._rms_in_range,

        "fit4c_branch": fit4c_branch,
        "fit4c_cut": CUT_4C,
        "fit4c_entries_40_120": int(h_4c._entries_in_range),
        "fit4c_mean_40_120": h_4c._mean_in_range,
        "fit4c_rms_40_120": h_4c._rms_in_range,

        "fit5c_candidate_branch": FINAL_5C_BRANCH,
        "fit5c_candidate_cut": CUT_5C_CANDIDATE,
        "fit5c_candidate_entries_40_120": int(h_5c._entries_in_range),
        "fit5c_candidate_mean_40_120": h_5c._mean_in_range,
        "fit5c_candidate_rms_40_120": h_5c._rms_in_range,

        "png": out_png,
        "pdf": out_pdf,
    }

    print("")
    print(f"Saved: {out_png}")
    print(f"Saved: {out_pdf}")
    print(f"{outname}:")
    print(f"  Raw:          N={row['raw_entries_40_120']}, mean={row['raw_mean_40_120']:.4f}, RMS={row['raw_rms_40_120']:.4f}")
    print(f"  4C:           N={row['fit4c_entries_40_120']}, mean={row['fit4c_mean_40_120']:.4f}, RMS={row['fit4c_rms_40_120']:.4f}")
    print(f"  5C candidate: N={row['fit5c_candidate_entries_40_120']}, mean={row['fit5c_candidate_mean_40_120']:.4f}, RMS={row['fit5c_candidate_rms_40_120']:.4f}")

    return row


def print_cut_counts(tree):
    n_all = tree.GetEntries()
    n_4c = tree.GetEntries(CUT_4C)
    n_5c = tree.GetEntries(CUT_5C_CANDIDATE)
    n_both = tree.GetEntries(f"({CUT_4C}) && ({CUT_5C_CANDIDATE})")

    print("=" * 80)
    print("Input tree and cut counts")
    print("=" * 80)
    print(f"Input tree entries:        {n_all}")
    print(f"Pass P4C > 0.03:           {n_4c}")
    print(f"Pass P5C candidate > 0.03: {n_5c}")
    print(f"Pass both cuts:            {n_both}")
    print("=" * 80)


def write_summary(rows):
    if not rows:
        return

    with open(SUMMARY_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    print(f"Summary CSV saved: {SUMMARY_CSV}")


def main():
    print(f"Input ROOT: {INPUT_ROOT}")
    print(f"Output dir: {OUT_DIR}")

    f = ROOT.TFile.Open(INPUT_ROOT, "READ")
    if not f or f.IsZombie():
        raise SystemExit(f"ERROR: cannot open ROOT file: {INPUT_ROOT}")

    tree = f.Get(TREE_NAME)
    if not tree:
        raise SystemExit(f"ERROR: cannot find tree '{TREE_NAME}' in {INPUT_ROOT}")

    required = [
        "m_small_raw", "m_large_raw",
        "m_small_4C", "m_large_4C",
        FINAL_5C_BRANCH,
        "prob_4C", "prob_5C", "prob_5C_isr",
        "isr_applied",
    ]
    check_branches(tree, required)

    print_cut_counts(tree)

    rows = []

    rows.append(
        plot_mass_overlay(
            tree,
            raw_branch="m_small_raw",
            fit4c_branch="m_small_4C",
            outname="mass_small_240_thesis_style_clean",
            xtitle="M_{W} (smaller dijet mass) [GeV]",
        )
    )

    rows.append(
        plot_mass_overlay(
            tree,
            raw_branch="m_large_raw",
            fit4c_branch="m_large_4C",
            outname="mass_large_240_thesis_style_clean",
            xtitle="M_{W} (larger dijet mass) [GeV]",
        )
    )

    rows = [r for r in rows if r is not None]
    write_summary(rows)

    f.Close()
    print("Done.")


if __name__ == "__main__":
    main()