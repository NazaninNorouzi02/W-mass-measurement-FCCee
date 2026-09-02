#!/usr/bin/env python3
"""Create thesis Figure 7.8 from the three energy-specific ROOT files.

Run without an energy argument:

    python3 plot_wmass_fig7-8.py

Branch mapping
--------------
* 162.5 GeV uses the 4C+ISR result: E_isr_4C_fitted/isr_4C_applied.
* 240 GeV uses the 5C+ISR result: E_isr_fitted/isr_applied.
* 365 GeV uses the 5C+ISR result: E_isr_fitted/isr_applied.

The MC curve is E_isr_true_collinear = |pz_isr_true| because this is the
truth quantity corresponding to the fitted one-effective-collinear-photon
model.  At each energy, truth and fitted curves use exactly the same events
on which that energy's ISR refit was selected.  The curves are event counts,
not independently normalized, matching the presentation of Figure 7.8.
"""

import argparse
import csv
import os

import ROOT


ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetCanvasColor(ROOT.kWhite)
ROOT.gStyle.SetPadColor(ROOT.kWhite)
ROOT.gStyle.SetFrameFillColor(ROOT.kWhite)
ROOT.gStyle.SetFrameBorderMode(0)
ROOT.gStyle.SetTitleFont(42, "XYZ")
ROOT.gStyle.SetLabelFont(42, "XYZ")
ROOT.gStyle.SetTitleSize(0.046, "XYZ")
ROOT.gStyle.SetLabelSize(0.040, "XYZ")
ROOT.gStyle.SetLegendFont(42)
ROOT.gStyle.SetLineScalePS(2.0)


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
TREE_NAME = "events"
NBINS, XMIN, XMAX = 180, 0.0, 180.0

CONFIGS = {
    "162.5": {
        "energy": 162.5,
        "label": "162.5 GeV",
        "tag": "162p5",
        "fitted": "E_isr_4C_fitted",
        "applied": "isr_4C_applied",
        "fit": "4C+ISR",
        "color": "#1D3557",
    },
    "240": {
        "energy": 240.0,
        "label": "240 GeV",
        "tag": "240",
        "fitted": "E_isr_fitted",
        "applied": "isr_applied",
        "fit": "5C+ISR",
        "color": "#8B1E3F",
    },
    "365": {
        "energy": 365.0,
        "label": "365 GeV",
        "tag": "365",
        "fitted": "E_isr_fitted",
        "applied": "isr_applied",
        "fit": "5C+ISR",
        "color": "#24513A",
    },
}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Create thesis Figure 7.8 using 162.5, 240 and 365 GeV."
    )
    parser.add_argument("--input-162p5", default=None)
    parser.add_argument("--input-240", default=None)
    parser.add_argument("--input-365", default=None)
    parser.add_argument("--outdir", default=None)
    return parser.parse_args()


def find_input(config, user_path=None):
    if user_path:
        path = os.path.abspath(os.path.expanduser(user_path))
        if not os.path.isfile(path):
            raise SystemExit(f"ERROR: ROOT file does not exist:\n  {path}")
        return path

    filename = f"wmass_fit_pvalue_pulls_ISR_ecm{config['tag']}.root"
    directories = [
        os.path.join(SCRIPT_DIR, "outputs", "wmass"),
        SCRIPT_DIR,
        os.path.join(os.getcwd(), "outputs", "wmass"),
        os.getcwd(),
    ]
    checked, matches, seen = [], [], set()
    for directory in directories:
        path = os.path.abspath(os.path.join(directory, filename))
        if path in seen:
            continue
        seen.add(path)
        checked.append(path)
        if os.path.isfile(path):
            matches.append(path)

    if len(matches) == 1:
        return matches[0]
    if len(matches) > 1:
        raise SystemExit(
            f"ERROR: multiple {config['label']} input files were found:\n"
            + "\n".join(f"  {path}" for path in matches)
            + "\nUse the corresponding --input-* option explicitly."
        )
    raise SystemExit(
        f"ERROR: cannot find the {config['label']} ROOT file. Checked:\n"
        + "\n".join(f"  {path}" for path in checked)
    )


def check_branches(tree, config):
    required = [
        "E_isr_true_collinear",
        config["fitted"],
        config["applied"],
    ]
    missing = [branch for branch in required if not tree.GetBranch(branch)]
    if missing:
        raise SystemExit(
            f"ERROR: {config['label']} file is missing branch(es):\n"
            + "\n".join(f"  {branch}" for branch in missing)
        )


def delete_existing(name):
    old = ROOT.gDirectory.Get(name)
    if old:
        old.Delete()


def make_histogram(tree, branch, name, selection):
    delete_existing(name)
    histogram = ROOT.TH1D(name, "", NBINS, XMIN, XMAX)
    histogram.Sumw2()
    histogram.SetStats(0)
    finite_selection = f"({selection}) && TMath::Finite({branch})"
    visible_selection = (
        f"({finite_selection}) && ({branch} >= {XMIN:.17g}) && "
        f"({branch} < {XMAX:.17g})"
    )
    tree.Draw(f"{branch}>>{name}", visible_selection, "goff")
    histogram.SetDirectory(0)
    histogram._selected = int(tree.GetEntries(finite_selection))
    histogram._visible = int(histogram.GetEntries())
    histogram._mean = float(histogram.GetMean())
    histogram._rms = float(histogram.GetStdDev())
    return histogram


def configure_axes(histogram):
    histogram.SetTitle("")
    histogram.GetXaxis().SetTitle("ISR energy [GeV]")
    histogram.GetYaxis().SetTitle("Events / 1 GeV")
    histogram.GetXaxis().CenterTitle()
    histogram.GetYaxis().CenterTitle()
    histogram.GetXaxis().SetTitleSize(0.047)
    histogram.GetYaxis().SetTitleSize(0.047)
    histogram.GetXaxis().SetLabelSize(0.040)
    histogram.GetYaxis().SetLabelSize(0.040)
    histogram.GetXaxis().SetTitleOffset(1.17)
    histogram.GetYaxis().SetTitleOffset(1.30)
    histogram.GetXaxis().SetNdivisions(510)
    histogram.GetYaxis().SetNdivisions(510)


def main():
    args = parse_args()
    overrides = {
        "162.5": args.input_162p5,
        "240": args.input_240,
        "365": args.input_365,
    }
    output_directory = (
        os.path.abspath(os.path.expanduser(args.outdir))
        if args.outdir
        else os.path.join(SCRIPT_DIR, "wmass_figure_7_8_all_energies")
    )
    os.makedirs(output_directory, exist_ok=True)

    root_files = {}
    products = []
    rows = []
    for energy_key in ("162.5", "240", "365"):
        config = CONFIGS[energy_key]
        input_path = find_input(config, overrides[energy_key])
        root_file = ROOT.TFile.Open(input_path, "READ")
        if not root_file or root_file.IsZombie():
            raise SystemExit(f"ERROR: cannot open ROOT file:\n  {input_path}")
        tree = root_file.Get(TREE_NAME)
        if not tree:
            root_file.Close()
            raise SystemExit(
                f"ERROR: tree '{TREE_NAME}' not found in:\n  {input_path}"
            )
        check_branches(tree, config)
        root_files[energy_key] = root_file

        selection = f"{config['applied']} > 0.5"
        h_truth = make_histogram(
            tree,
            "E_isr_true_collinear",
            f"h78_truth_{config['tag']}",
            selection,
        )
        h_fitted = make_histogram(
            tree,
            config["fitted"],
            f"h78_fitted_{config['tag']}",
            selection,
        )
        if h_truth._visible <= 0 or h_fitted._visible <= 0:
            raise RuntimeError(
                f"Empty Figure 7.8 histogram at {config['label']}"
            )

        color = ROOT.TColor.GetColor(config["color"])
        h_truth.SetLineColor(color)
        h_truth.SetLineWidth(2)
        h_truth.SetLineStyle(2)
        h_truth.SetFillColor(color)
        h_truth.SetFillStyle(3004)
        h_fitted.SetLineColor(color)
        h_fitted.SetLineWidth(3)
        h_fitted.SetLineStyle(1)
        h_fitted.SetFillStyle(0)
        products.append((config, h_truth, h_fitted))

        for curve, branch, histogram in [
            ("MC effective collinear ISR", "E_isr_true_collinear", h_truth),
            ("fitted ISR", config["fitted"], h_fitted),
        ]:
            rows.append(
                {
                    "energy_GeV": config["energy"],
                    "fit": config["fit"],
                    "curve": curve,
                    "branch": branch,
                    "selection": selection,
                    "selected_entries": histogram._selected,
                    "visible_entries": histogram._visible,
                    "mean_GeV": histogram._mean,
                    "stddev_GeV": histogram._rms,
                }
            )
        print(
            f"{config['label']:>9}: {config['fit']}, "
            f"selected={h_fitted._selected}, input={input_path}"
        )

    canvas = ROOT.TCanvas("c_figure_7_8_all_energies", "", 1250, 830)
    canvas.SetLeftMargin(0.12)
    canvas.SetRightMargin(0.045)
    canvas.SetBottomMargin(0.13)
    canvas.SetTopMargin(0.055)
    canvas.SetTicks(1, 1)
    canvas.SetLogy(True)

    first = products[0][1]
    configure_axes(first)
    maximum = max(
        histogram.GetMaximum()
        for _, h_truth, h_fitted in products
        for histogram in (h_truth, h_fitted)
    )
    positive_bins = [
        histogram.GetBinContent(i)
        for _, h_truth, h_fitted in products
        for histogram in (h_truth, h_fitted)
        for i in range(1, histogram.GetNbinsX() + 1)
        if histogram.GetBinContent(i) > 0.0
    ]
    first.SetMinimum(max(0.5, min(positive_bins) * 0.5))
    first.SetMaximum(maximum * 3.0)
    first.Draw("HIST")
    products[0][2].Draw("HIST SAME")
    for _, h_truth, h_fitted in products[1:]:
        h_truth.Draw("HIST SAME")
        h_fitted.Draw("HIST SAME")

    legend = ROOT.TLegend(0.45, 0.74, 0.94, 0.93)
    legend.SetBorderSize(0)
    legend.SetFillColor(ROOT.kWhite)
    legend.SetFillStyle(1001)
    legend.SetTextFont(42)
    legend.SetTextSize(0.031)
    legend.SetNColumns(2)
    for config, h_truth, h_fitted in products:
        legend.AddEntry(h_truth, f"MC ISR, {config['label']}", "f")
        legend.AddEntry(h_fitted, f"Fitted ISR, {config['label']}", "l")
    legend.Draw()
    canvas.RedrawAxis()

    png_path = os.path.join(
        output_directory,
        "figure_7_8_ISR_truth_vs_fitted_all_energies.png",
    )
    csv_path = os.path.join(
        output_directory,
        "figure_7_8_ISR_truth_vs_fitted_summary.csv",
    )
    canvas.Update()
    canvas.SaveAs(png_path)
    canvas.Close()
    with open(csv_path, "w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    for root_file in root_files.values():
        root_file.Close()

    print(f"Saved PNG : {png_path}")
    print(f"Saved CSV : {csv_path}")


if __name__ == "__main__":
    main()