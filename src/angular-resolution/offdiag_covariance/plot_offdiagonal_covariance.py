#!/usr/bin/env python3
import argparse
import os
from glob import glob

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPalette(ROOT.kBird)


VARIABLES = {
    "alpha": {
        "raw": "delta_alpha_j{j}",
        "filtered": "filtered_delta_alpha_j{j}",
        "label": "#Delta#alpha",
        "range": (-0.25, 0.25),
    },
    "x": {
        "raw": "delta_x_j{j}",
        "filtered": "filtered_delta_x_j{j}",
        "label": "#Deltax = #Delta log(p/m)",
        "range": (-0.6, 0.6),
    },
    "theta": {
        "raw": "delta_theta_j{j}",
        "filtered": "filtered_delta_theta_j{j}",
        "label": "#Delta#theta [rad]",
        "range": (-0.02, 0.02),
    },
    "phi": {
        "raw": "delta_phi_j{j}",
        "filtered": "filtered_delta_phi_j{j}",
        "label": "#Delta#phi [rad]",
        "range": (-0.02, 0.02),
    },
}


PAIRS = [
    ("x", "alpha"),
    ("theta", "alpha"),
    ("phi", "alpha"),
    ("theta", "x"),
    ("phi", "x"),
    ("theta", "phi"),
]


def resolve_input(path_or_pattern):
    matches = sorted(glob(path_or_pattern))
    if matches:
        return matches
    if os.path.exists(path_or_pattern):
        return [path_or_pattern]
    raise FileNotFoundError(f"No ROOT file matches: {path_or_pattern}")


def get_chain(files, tree_name):
    chain = ROOT.TChain(tree_name)
    for path in files:
        chain.Add(path)
    if chain.GetEntries() <= 0:
        raise RuntimeError(f"No entries found in tree '{tree_name}'")
    return chain


def branch_exists(chain, branch_name):
    return chain.GetBranch(branch_name) is not None


def fill_hist(chain, hist, x_var, y_var, jets, use_filtered):
    mode = "filtered" if use_filtered else "raw"
    total = 0

    for jet in jets:
        x_branch = VARIABLES[x_var][mode].format(j=jet)
        y_branch = VARIABLES[y_var][mode].format(j=jet)

        if not branch_exists(chain, x_branch):
            raise RuntimeError(f"Missing branch: {x_branch}")
        if not branch_exists(chain, y_branch):
            raise RuntimeError(f"Missing branch: {y_branch}")

        cuts = []
        if use_filtered:
            cuts += [f"{x_branch}>-998", f"{y_branch}>-998"]
        cut = " && ".join(cuts)

        draw_expr = f"{y_branch}:{x_branch}>>+{hist.GetName()}"
        total += chain.Draw(draw_expr, cut, "goff")

    return total


def draw_info_box(hist, n_filled):
    corr = hist.GetCorrelationFactor()
    cov = hist.GetCovariance()

    box = ROOT.TPaveText(0.58, 0.72, 0.88, 0.88, "NDC")
    box.SetFillColor(0)
    box.SetBorderSize(1)
    box.SetTextAlign(12)
    box.SetTextSize(0.030)
    box.AddText(f"Entries = {n_filled}")
    box.AddText(f"#rho = {corr:.4g}")
    box.AddText(f"Cov = {cov:.4g}")
    box.AddText(f"#sigma_x = {hist.GetStdDev(1):.4g}")
    box.AddText(f"#sigma_y = {hist.GetStdDev(2):.4g}")
    box.Draw()
    return box


def make_plot(args):
    files = resolve_input(args.input)
    chain = get_chain(files, args.tree)

    jets = [1, 2, 3, 4] if args.jet == "all" else [int(args.jet)]
    mode_text = "filtered" if args.filtered else "raw"

    canvas = ROOT.TCanvas("c_offdiag", "off diagonal covariance", 1600, 2100)
    canvas.Divide(2, 3)

    boxes = []
    hists = []

    print("\nOff-diagonal covariance/correlation summary")
    print(f"Input files: {len(files)}")
    print(f"Tree: {args.tree}")
    print(f"Jets: {jets}")
    print(f"Mode: {mode_text}")
    print("pair,entries,corr,cov,std_x,std_y")

    for idx, (x_var, y_var) in enumerate(PAIRS, start=1):
        canvas.cd(idx)
        ROOT.gPad.SetRightMargin(0.16)
        ROOT.gPad.SetLeftMargin(0.12)
        ROOT.gPad.SetBottomMargin(0.11)

        xmin, xmax = VARIABLES[x_var]["range"]
        ymin, ymax = VARIABLES[y_var]["range"]
        hist = ROOT.TH2F(
            f"h_{y_var}_vs_{x_var}",
            "",
            args.bins,
            xmin,
            xmax,
            args.bins,
            ymin,
            ymax,
        )
        hist.GetXaxis().SetTitle(VARIABLES[x_var]["label"])
        hist.GetYaxis().SetTitle(VARIABLES[y_var]["label"])
        hist.GetZaxis().SetTitle("Events")
        hist.GetXaxis().SetTitleSize(0.045)
        hist.GetYaxis().SetTitleSize(0.045)
        hist.GetXaxis().SetLabelSize(0.035)
        hist.GetYaxis().SetLabelSize(0.035)
        hist.GetZaxis().SetLabelSize(0.030)

        n_filled = fill_hist(chain, hist, x_var, y_var, jets, args.filtered)
        hist.Draw("COLZ")
        box = draw_info_box(hist, n_filled)

        print(
            f"{y_var}_vs_{x_var},"
            f"{n_filled},"
            f"{hist.GetCorrelationFactor():.8g},"
            f"{hist.GetCovariance():.8g},"
            f"{hist.GetStdDev(1):.8g},"
            f"{hist.GetStdDev(2):.8g}"
        )

        hists.append(hist)
        boxes.append(box)

    canvas.cd()
    title = ROOT.TLatex()
    title.SetNDC(True)
    title.SetTextAlign(22)
    title.SetTextSize(0.025)
    title.DrawLatex(
        0.5,
        0.985,
        f"Off-diagonal residual covariance terms ({args.label}, {mode_text}, jets {args.jet})",
    )

    os.makedirs(args.output_dir, exist_ok=True)
    out_base = os.path.join(
        args.output_dir,
        f"offdiag_covariance_{args.label}_{mode_text}_jet{args.jet}",
    )
    canvas.SaveAs(out_base + ".png")
    canvas.SaveAs(out_base + ".pdf")
    print(f"\nSaved:\n  {out_base}.png\n  {out_base}.pdf")


def main():
    parser = argparse.ArgumentParser(
        description="Draw thesis-style off-diagonal residual covariance plots."
    )
    parser.add_argument(
        "--input",
        default="outputs/stage3_angular_resolution/angular_resolution_ecm160.root",
        help="Stage-3 ROOT file or glob pattern.",
    )
    parser.add_argument("--tree", default="events", help="TTree name.")
    parser.add_argument(
        "--label",
        default="ecm160",
        help="Label used in output filenames and plot title.",
    )
    parser.add_argument(
        "--jet",
        default="all",
        choices=["all", "1", "2", "3", "4"],
        help="Use all jets together or one specific jet.",
    )
    parser.add_argument(
        "--filtered",
        action="store_true",
        help="Use filtered_delta_*_jN branches and reject -999 sentinels.",
    )
    parser.add_argument("--bins", type=int, default=120)
    parser.add_argument("--output-dir", default="plots/offdiag_covariance")
    args = parser.parse_args()

    make_plot(args)


if __name__ == "__main__":
    main()
