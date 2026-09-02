#!/usr/bin/env python3
"""Create publication-ready thesis-style W-mass plots for one energy.

Examples
--------
    python3 plot_wmass.py 365
    python3 plot_wmass.py 240
    python3 plot_wmass.py 162.5

Outputs
-------
* Figure 7.9 (without ISR treatment versus with ISR treatment).
* Figure 7.10 (raw, ISR-treated 4C and, above threshold, ISR-treated 5C).
* Appendix Figure B.4 (fits without ISR treatment, before/after P > 0.03).
* Four 4C pull figures without ISR treatment, one per jet, with
  alpha/theta/phi/x panels.

At 162.5 GeV no 5C branch is read or plotted.  Figure 7.9 uses 4C+ISR,
Figure 7.10 contains raw and ISR-treated 4C only, and Appendix B.4 contains only
the two 4C panels.

The mass/probability definitions intentionally follow the thesis workflow.
Figures 7.9 and 7.10 show the ISR-treated result before the 3% event
rejection.  Appendix B.4 and the pull figures compare fits without explicit
ISR treatment before and after their own probability requirement.

PNG outputs are written.  As in the thesis, mass and
probability boxes show only the mean and standard deviation.  Only the pull
boxes show Gaussian-fit parameter uncertainties.  These pull-fit errors are
not the final W-mass measurement uncertainty.
"""

import argparse
import csv
import os

import ROOT


ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)
ROOT.gStyle.SetCanvasColor(ROOT.kWhite)
ROOT.gStyle.SetPadColor(ROOT.kWhite)
ROOT.gStyle.SetFrameFillColor(ROOT.kWhite)
ROOT.gStyle.SetFrameBorderMode(0)
ROOT.gStyle.SetTitleFont(42, "XYZ")
ROOT.gStyle.SetLabelFont(42, "XYZ")
ROOT.gStyle.SetTitleSize(0.052, "XYZ")
ROOT.gStyle.SetLabelSize(0.043, "XYZ")
ROOT.gStyle.SetLegendFont(42)
ROOT.gStyle.SetLineScalePS(2.0)
ROOT.gStyle.SetEndErrorSize(2.0)
ROOT.gStyle.SetErrorX(0.0)
ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetPadTickY(1)


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
TREE_NAME = "events"
P_CUT = 0.03

MASS_NBINS, MASS_XMIN, MASS_XMAX = 160, 40.0, 120.0
PROB_NBINS, PROB_XMIN, PROB_XMAX = 100, 0.0, 1.0
PULL_NBINS, PULL_XMIN, PULL_XMAX = 100, -5.0, 5.0
PULL_FIT_MIN, PULL_FIT_MAX = -3.0, 3.0

PULL_VARIABLES = ["alpha", "theta", "phi", "x"]
PULL_LABELS = {
    "alpha": "#alpha",
    "theta": "#theta",
    "phi": "#phi",
    "x": "x = ln(p/m)",
}

COLOR_RAW = ROOT.TColor.GetColor("#202020")
COLOR_4C = ROOT.TColor.GetColor("#C43B3B")
COLOR_5C = ROOT.TColor.GetColor("#2457A6")
COLOR_WITHOUT_ISR = ROOT.TColor.GetColor("#303030")
COLOR_WITH_ISR = ROOT.TColor.GetColor("#2457A6")
COLOR_BEFORE = ROOT.TColor.GetColor("#D1495B")
COLOR_AFTER = ROOT.TColor.GetColor("#1F5AA6")


CONFIGS = {
    "162.5": {
        "energy": 162.5,
        "label": "162.5 GeV",
        "tag": "162p5",
        "has_5c": False,
    },
    "240": {
        "energy": 240.0,
        "label": "240 GeV",
        "tag": "240",
        "has_5c": True,
    },
    "365": {
        "energy": 365.0,
        "label": "365 GeV",
        "tag": "365",
        "has_5c": True,
    },
}


def normalize_energy(value):
    aliases = {
        "162.5": "162.5",
        "162p5": "162.5",
        "162.6": "162.5",
        "162p6": "162.5",
        "240": "240",
        "240.0": "240",
        "365": "365",
        "365.0": "365",
    }
    key = str(value).strip().lower()
    if key not in aliases:
        raise argparse.ArgumentTypeError("energy must be 162.5, 240, or 365")
    return aliases[key]


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Create Figures 7.9, 7.10, Appendix B.4 and "
            "pull plots for one centre-of-mass energy."
        )
    )
    parser.add_argument("energy", type=normalize_energy)
    parser.add_argument(
        "--input",
        default=None,
        help="optional explicit ROOT input path",
    )
    parser.add_argument(
        "--outdir",
        default=None,
        help="optional output directory",
    )
    return parser.parse_args()


def find_input(config, user_path=None):
    if user_path:
        path = os.path.abspath(os.path.expanduser(user_path))
        if not os.path.isfile(path):
            raise SystemExit(f"ERROR: input ROOT file does not exist:\n  {path}")
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
            "ERROR: multiple copies of the exact input filename were found:\n"
            + "\n".join(f"  {path}" for path in matches)
            + "\nUse --input to choose the intended file explicitly."
        )
    raise SystemExit(
        "ERROR: input ROOT file was not found. Checked:\n"
        + "\n".join(f"  {path}" for path in checked)
    )


def branch_exists(tree, branch):
    return bool(tree.GetBranch(branch))


def check_required_branches(tree, config):
    required = [
        "m_small_raw",
        "m_large_raw",
        "m_small_4C",
        "m_large_4C",
        "m_small_4C_final",
        "m_large_4C_final",
        "prob_4C",
        "prob_final_4C",
        "standard_4C_valid",
        "pull4C_standard_valid",
    ]
    for jet in range(1, 5):
        for variable in PULL_VARIABLES:
            required.append(f"pull4C_standard_{variable}_j{jet}")

    # UPDATE: Added required branches for the new pull plot categories
    required.extend(
        [
            "pull4C_final_valid",
            "pull5C_standard_valid",
            "pull5C_final_valid",
        ]
    )
    for jet in range(1, 5):
        for variable in PULL_VARIABLES:
            required.append(f"pull4C_final_{variable}_j{jet}")
            required.append(f"pull5C_standard_{variable}_j{jet}")
            required.append(f"pull5C_final_{variable}_j{jet}")

    if config["has_5c"]:
        required.extend(
            [
                "mW_5C",
                "mW_5C_final",
                "prob_5C",
                "prob_final_5C",
                "standard_5C_valid",
            ]
        )

    missing = [branch for branch in required if not branch_exists(tree, branch)]
    if missing:
        special = (
            "\nAt 162.5 GeV no 5C branch is required."
            if not config["has_5c"]
            else ""
        )
        raise SystemExit(
            "ERROR: input tree is missing required branch(es):\n"
            + "\n".join(f"  {branch}" for branch in missing)
            + special
            + "\nUse the ROOT file produced by the supplied ISR/pull producer."
        )


def delete_existing(name):
    old = ROOT.gDirectory.Get(name)
    if old:
        old.Delete()


def finite_cut(expression):
    return f"TMath::Finite(({expression}))"


def make_histogram(
    tree,
    expressions,
    name,
    nbins,
    xmin,
    xmax,
    selection,
    normalize=True,
):
    if isinstance(expressions, str):
        expressions = [expressions]

    delete_existing(name)
    histogram = ROOT.TH1D(name, "", nbins, xmin, xmax)
    histogram.Sumw2()
    histogram.SetStats(0)

    finite_fills = 0
    for expression in expressions:
        finite_selection = f"({selection}) && ({finite_cut(expression)})"
        visible_selection = (
            f"({finite_selection}) && "
            f"(({expression}) >= {xmin:.17g}) && "
            f"(({expression}) < {xmax:.17g})"
        )
        finite_fills += int(tree.GetEntries(finite_selection))
        tree.Draw(f"({expression})>>+{name}", visible_selection, "goff")

    histogram.SetDirectory(0)
    histogram._expressions = list(expressions)
    histogram._selection = selection
    histogram._selected_events = int(tree.GetEntries(selection))
    histogram._finite_fills = finite_fills
    histogram._visible_fills = int(histogram.GetEntries())
    histogram._mean = float(histogram.GetMean())
    histogram._rms = float(histogram.GetStdDev())
    histogram._mean_error = float(histogram.GetMeanError())
    histogram._rms_error = float(histogram.GetStdDevError())

    integral = histogram.Integral(1, histogram.GetNbinsX())
    if normalize and integral > 0.0:
        histogram.Scale(1.0 / integral)
    return histogram


def require_nonempty(histogram, description):
    if not histogram or histogram._visible_fills <= 0:
        raise RuntimeError(f"Empty histogram: {description}")


def style_histogram(histogram, color, line_style=1, width=3, fill_style=0):
    histogram.SetLineColor(color)
    histogram.SetLineStyle(line_style)
    histogram.SetLineWidth(width)
    histogram.SetFillColor(color)
    histogram.SetFillStyle(fill_style)


def configure_axes(histogram, x_title):
    histogram.SetTitle("")
    histogram.GetXaxis().SetTitle(x_title)
    # The thesis figures show the y-axis scale without a y-axis title.
    histogram.GetYaxis().SetTitle("")
    histogram.GetXaxis().CenterTitle()
    histogram.GetXaxis().SetTitleSize(0.056)
    histogram.GetXaxis().SetLabelSize(0.047)
    histogram.GetYaxis().SetLabelSize(0.047)
    histogram.GetXaxis().SetTitleOffset(1.08)
    histogram.GetXaxis().SetNdivisions(510)
    histogram.GetYaxis().SetNdivisions(508)


def configure_pad(pad, right=0.035, bottom=0.145, top=0.055, left=0.135):
    pad.SetLeftMargin(left)
    pad.SetRightMargin(right)
    pad.SetBottomMargin(bottom)
    pad.SetTopMargin(top)
    pad.SetTicks(1, 1)


def plot_corner_labels(
    config,
    second_line="WW #rightarrow qqqq",
    x_right=0.935,
    y_top=0.315,
):
    """Draw the thesis-style energy and channel labels at lower right."""
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextFont(62)
    latex.SetTextSize(0.038)
    latex.SetTextAlign(33)
    latex.DrawLatex(x_right, y_top, f"#sqrt{{s}} = {config['label']}")
    latex.SetTextSize(0.036)
    latex.DrawLatex(x_right, y_top - 0.055, second_line)
    return latex


def thesis_stat_box(
    x1,
    y1,
    x2,
    y2,
    histogram,
    color,
    precision=2,
    text_size=0.032,
    show_errors=False,
    mean=None,
    mean_error=None,
    sigma=None,
    sigma_error=None,
):
    """Draw one compact thesis-style result box for one colored curve."""
    box = ROOT.TPaveText(x1, y1, x2, y2, "NDC")
    box.SetFillColor(ROOT.kWhite)
    box.SetFillStyle(1001)
    box.SetBorderSize(1)
    box.SetLineColor(color)
    box.SetLineWidth(1)
    box.SetTextFont(42)
    box.SetTextAlign(12)
    box.SetTextSize(text_size)
    box.SetTextColor(color)

    mean_value = histogram._mean if mean is None else mean
    sigma_value = histogram._rms if sigma is None else sigma
    if show_errors:
        mean_uncertainty = histogram._mean_error if mean_error is None else mean_error
        sigma_uncertainty = histogram._rms_error if sigma_error is None else sigma_error
        box.AddText(
            f"Mean  {mean_value:+.{precision}f} #pm {mean_uncertainty:.{precision}f}"
        )
        box.AddText(
            f"Sigma  {abs(sigma_value):.{precision}f} #pm {sigma_uncertainty:.{precision}f}"
        )
    else:
        box.AddText(f"Mean  {mean_value:.{precision}f}")
        box.AddText(f"Sigma  {abs(sigma_value):.{precision}f}")
    box.Draw()
    return box


def save_canvas(canvas, output_directory, base_name):
    png_path = os.path.join(output_directory, base_name + ".png")
    canvas.Update()
    canvas.SaveAs(png_path)
    print(f"Saved PNG : {png_path}")
    return png_path


def summary_row(config, figure, panel, curve, histogram):
    return {
        "energy_GeV": config["energy"],
        "figure": figure,
        "panel": panel,
        "curve": curve,
        "expressions": " + ".join(histogram._expressions),
        "selection": histogram._selection,
        "selected_events": histogram._selected_events,
        "finite_fills": histogram._finite_fills,
        "visible_fills": histogram._visible_fills,
        "mean": histogram._mean,
        "mean_error": histogram._mean_error,
        "stddev": histogram._rms,
        "stddev_error": histogram._rms_error,
    }


def draw_figure_7_9(tree, config, output_directory):
    if config["has_5c"]:
        fit_name = "5C"
        valid_cut = "standard_5C_valid > 0.5"
        mass_without, mass_with = ["mW_5C"], ["mW_5C_final"]
        prob_without, prob_with = ["prob_5C"], ["prob_final_5C"]
        mass_title = "M_{W} (5C) [GeV]"
    else:
        fit_name = "4C"
        valid_cut = "standard_4C_valid > 0.5"
        mass_without = ["m_small_4C", "m_large_4C"]
        mass_with = ["m_small_4C_final", "m_large_4C_final"]
        prob_without, prob_with = ["prob_4C"], ["prob_final_4C"]
        mass_title = "M_{W} candidates (4C) [GeV]"

    tag = config["tag"]
    h_mass_without = make_histogram(
        tree, mass_without, f"h79_mass_without_{tag}",
        MASS_NBINS, MASS_XMIN, MASS_XMAX, valid_cut,
    )
    h_mass_with = make_histogram(
        tree, mass_with, f"h79_mass_with_isr_{tag}",
        MASS_NBINS, MASS_XMIN, MASS_XMAX, valid_cut,
    )
    h_prob_without = make_histogram(
        tree, prob_without, f"h79_prob_without_{tag}",
        PROB_NBINS, PROB_XMIN, PROB_XMAX, valid_cut,
    )
    h_prob_with = make_histogram(
        tree, prob_with, f"h79_prob_with_isr_{tag}",
        PROB_NBINS, PROB_XMIN, PROB_XMAX, valid_cut,
    )
    for label, histogram in [
        ("mass without ISR", h_mass_without),
        ("mass with ISR treatment", h_mass_with),
        ("probability without ISR", h_prob_without),
        ("probability with ISR treatment", h_prob_with),
    ]:
        require_nonempty(histogram, f"Figure 7.9 {label}")

    for histogram in (h_mass_without, h_prob_without):
        style_histogram(histogram, COLOR_WITHOUT_ISR, line_style=2, width=3)
    for histogram in (h_mass_with, h_prob_with):
        style_histogram(histogram, COLOR_WITH_ISR, width=3)

    canvas = ROOT.TCanvas(f"c_figure_7_9_{tag}", "", 1800, 820)
    canvas.Divide(2, 1, 0.010, 0.0)
    keep = []

    canvas.cd(1)
    configure_pad(ROOT.gPad)
    configure_axes(h_mass_without, mass_title)
    h_mass_without.SetMinimum(0.0)
    h_mass_without.SetMaximum(
        1.20 * max(h_mass_without.GetMaximum(), h_mass_with.GetMaximum())
    )
    h_mass_without.Draw("HIST")
    h_mass_with.Draw("HIST SAME")
    legend_mass = ROOT.TLegend(0.16, 0.775, 0.50, 0.915)
    legend_mass.SetBorderSize(0)
    legend_mass.SetFillStyle(0)
    legend_mass.SetTextSize(0.034)
    legend_mass.AddEntry(h_mass_with, "With ISR treatment", "l")
    legend_mass.AddEntry(h_mass_without, "Without ISR treatment", "l")
    legend_mass.Draw()
    box_mass_with = thesis_stat_box(
        0.730, 0.830, 0.940, 0.925,
        h_mass_with, COLOR_WITH_ISR,
        precision=2, text_size=0.031,
    )
    box_mass_without = thesis_stat_box(
        0.730, 0.715, 0.940, 0.810,
        h_mass_without, COLOR_WITHOUT_ISR,
        precision=2, text_size=0.031,
    )
    labels_m = plot_corner_labels(config)
    keep += [
        h_mass_without,
        h_mass_with,
        legend_mass,
        box_mass_with,
        box_mass_without,
        labels_m,
    ]

    canvas.cd(2)
    configure_pad(ROOT.gPad)
    ROOT.gPad.SetLogy(True)
    configure_axes(h_prob_without, "#chi^{2} probability")
    positive = [
        h.GetBinContent(i)
        for h in (h_prob_without, h_prob_with)
        for i in range(1, h.GetNbinsX() + 1)
        if h.GetBinContent(i) > 0.0
    ]
    h_prob_without.SetMinimum(max(min(positive) * 0.5, 1.0e-7))
    h_prob_without.SetMaximum(
        3.0 * max(h_prob_without.GetMaximum(), h_prob_with.GetMaximum())
    )
    h_prob_without.Draw("HIST")
    h_prob_with.Draw("HIST SAME")
    legend_prob = ROOT.TLegend(0.16, 0.775, 0.50, 0.915)
    legend_prob.SetBorderSize(0)
    legend_prob.SetFillStyle(0)
    legend_prob.SetTextSize(0.034)
    legend_prob.AddEntry(h_prob_with, "With ISR treatment", "l")
    legend_prob.AddEntry(h_prob_without, "Without ISR treatment", "l")
    legend_prob.Draw()
    box_prob_with = thesis_stat_box(
        0.730, 0.830, 0.940, 0.925,
        h_prob_with, COLOR_WITH_ISR,
        precision=3, text_size=0.031,
    )
    box_prob_without = thesis_stat_box(
        0.730, 0.715, 0.940, 0.810,
        h_prob_without, COLOR_WITHOUT_ISR,
        precision=3, text_size=0.031,
    )
    # The probability tail rises near x=1, so place the physics label higher.
    labels_p = plot_corner_labels(config, y_top=0.420)
    keep += [
        h_prob_without,
        h_prob_with,
        legend_prob,
        box_prob_with,
        box_prob_without,
        labels_p,
    ]

    save_canvas(
        canvas,
        output_directory,
        f"figure_7_9_{fit_name}_with_without_ISR_treatment_ecm{tag}",
    )
    canvas.Close()
    return [
        summary_row(config, "7.9", "mass", "without ISR treatment", h_mass_without),
        summary_row(config, "7.9", "mass", "with ISR treatment", h_mass_with),
        summary_row(config, "7.9", "probability", "without ISR treatment", h_prob_without),
        summary_row(config, "7.9", "probability", "with ISR treatment", h_prob_with),
    ]


def draw_figure_7_10_panel(canvas, pad_number, tree, config, panel):
    canvas.cd(pad_number)
    configure_pad(ROOT.gPad)
    tag = config["tag"]
    is_small = panel == "smaller"
    raw_branch = "m_small_raw" if is_small else "m_large_raw"
    four_branch = "m_small_4C_final" if is_small else "m_large_4C_final"
    x_title = f"M_{{W}} ({panel} dijet mass) [GeV]"

    histograms = {
        "raw": make_histogram(
            tree, raw_branch, f"h710_raw_{panel}_{tag}",
            MASS_NBINS, MASS_XMIN, MASS_XMAX, "1",
        ),
        "4C": make_histogram(
            tree, four_branch, f"h710_4c_{panel}_{tag}",
            MASS_NBINS, MASS_XMIN, MASS_XMAX, "standard_4C_valid > 0.5",
        ),
    }
    if config["has_5c"]:
        histograms["5C"] = make_histogram(
            tree, "mW_5C_final", f"h710_5c_{panel}_{tag}",
            MASS_NBINS, MASS_XMIN, MASS_XMAX, "standard_5C_valid > 0.5",
        )
    for key, histogram in histograms.items():
        require_nonempty(histogram, f"Figure 7.10 {panel} {key}")

    style_histogram(histograms["raw"], COLOR_RAW, line_style=2, width=3)
    style_histogram(histograms["4C"], COLOR_4C, width=3)
    if config["has_5c"]:
        style_histogram(histograms["5C"], COLOR_5C, width=3)

    maximum = max(h.GetMaximum() for h in histograms.values())
    configure_axes(histograms["raw"], x_title)
    histograms["raw"].SetMinimum(0.0)
    histograms["raw"].SetMaximum(1.20 * maximum)
    histograms["raw"].Draw("HIST")
    histograms["4C"].Draw("HIST SAME")
    if config["has_5c"]:
        histograms["5C"].Draw("HIST SAME")

    legend = ROOT.TLegend(0.16, 0.750, 0.49, 0.920)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.034)
    legend.AddEntry(histograms["raw"], "Raw dijet mass", "l")
    legend.AddEntry(histograms["4C"], "4C with ISR treatment", "l")
    if config["has_5c"]:
        legend.AddEntry(histograms["5C"], "5C with ISR treatment", "l")
    legend.Draw()

    stats = []
    if config["has_5c"]:
        stats = [
            thesis_stat_box(
                0.745, 0.835, 0.940, 0.925,
                histograms["5C"], COLOR_5C,
                precision=2, text_size=0.030,
            ),
            thesis_stat_box(
                0.745, 0.725, 0.940, 0.815,
                histograms["4C"], COLOR_4C,
                precision=2, text_size=0.030,
            ),
            thesis_stat_box(
                0.745, 0.615, 0.940, 0.705,
                histograms["raw"], COLOR_RAW,
                precision=2, text_size=0.030,
            ),
        ]
    else:
        stats = [
            thesis_stat_box(
                0.745, 0.835, 0.940, 0.925,
                histograms["4C"], COLOR_4C,
                precision=2, text_size=0.030,
            ),
            thesis_stat_box(
                0.745, 0.725, 0.940, 0.815,
                histograms["raw"], COLOR_RAW,
                precision=2, text_size=0.030,
            ),
        ]
    labels = plot_corner_labels(config)
    ROOT.gPad.RedrawAxis()

    rows = [
        summary_row(config, "7.10", panel, "raw dijet mass", histograms["raw"]),
        summary_row(config, "7.10", panel, "4C with ISR treatment", histograms["4C"]),
    ]
    if config["has_5c"]:
        rows.append(summary_row(config, "7.10", panel, "5C with ISR treatment", histograms["5C"]))
    return rows, list(histograms.values()) + [legend, labels] + stats


def draw_figure_7_10(tree, config, output_directory):
    canvas = ROOT.TCanvas(f"c_figure_7_10_{config['tag']}", "", 1800, 820)
    canvas.Divide(2, 1, 0.002, 0.002)
    rows, keep = [], []
    for pad_number, panel in enumerate(("smaller", "larger"), start=1):
        panel_rows, panel_keep = draw_figure_7_10_panel(
            canvas, pad_number, tree, config, panel
        )
        rows += panel_rows
        keep += panel_keep
    save_canvas(canvas, output_directory, f"figure_7_10_ecm{config['tag']}")
    canvas.Close()
    return rows


def draw_before_after_mass_panel(
    pad,
    tree,
    config,
    fit_name,
    panel,
    mass_expression,
    probability_expression,
    valid_selection,
    x_title,
):
    pad.cd()
    configure_pad(pad)
    tag = config["tag"]
    after_selection = (
        f"({valid_selection}) && ({probability_expression} > {P_CUT:.17g})"
    )
    h_before = make_histogram(
        tree, mass_expression, f"hb4_{fit_name}_{panel}_before_{tag}",
        MASS_NBINS, MASS_XMIN, MASS_XMAX, valid_selection,
    )
    h_after = make_histogram(
        tree, mass_expression, f"hb4_{fit_name}_{panel}_after_{tag}",
        MASS_NBINS, MASS_XMIN, MASS_XMAX, after_selection,
    )
    require_nonempty(h_before, f"Appendix B.4 {fit_name} {panel} before")
    require_nonempty(h_after, f"Appendix B.4 {fit_name} {panel} after")
    style_histogram(h_before, COLOR_BEFORE, line_style=1, width=3)
    style_histogram(h_after, COLOR_AFTER, width=3)
    configure_axes(h_before, x_title)
    h_before.SetMinimum(0.0)
    h_before.SetMaximum(1.20 * max(h_before.GetMaximum(), h_after.GetMaximum()))
    h_before.Draw("HIST")
    h_after.Draw("HIST SAME")

    legend = ROOT.TLegend(0.16, 0.785, 0.49, 0.915)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.034)
    legend.AddEntry(h_before, "Without P cut", "l")
    legend.AddEntry(h_after, f"With P_{{{fit_name}}} > {P_CUT:.2f}", "l")
    legend.Draw()
    stats_before = thesis_stat_box(
        0.745, 0.835, 0.940, 0.925,
        h_before, COLOR_BEFORE,
        precision=2, text_size=0.030,
    )
    stats_after = thesis_stat_box(
        0.745, 0.725, 0.940, 0.815,
        h_after, COLOR_AFTER,
        precision=2, text_size=0.030,
    )
    labels = plot_corner_labels(config)
    pad.RedrawAxis()
    rows = [
        summary_row(
            config,
            "Appendix B.4",
            panel,
            f"{fit_name} without ISR treatment; no probability requirement",
            h_before,
        ),
        summary_row(
            config,
            "Appendix B.4",
            panel,
            f"{fit_name} without ISR treatment; P_{fit_name}>0.03",
            h_after,
        ),
    ]
    return rows, [h_before, h_after, legend, stats_before, stats_after, labels]


def draw_appendix_b4(tree, config, output_directory):
    tag = config["tag"]
    rows, keep, pads = [], [], []
    if config["has_5c"]:
        canvas = ROOT.TCanvas(f"c_appendix_B4_{tag}", "", 1700, 1180)
        pads = [
            ROOT.TPad(f"pb4_small_{tag}", "", 0.00, 0.53, 0.50, 1.00),
            ROOT.TPad(f"pb4_large_{tag}", "", 0.50, 0.53, 1.00, 1.00),
            ROOT.TPad(f"pb4_5c_{tag}", "", 0.25, 0.00, 0.75, 0.47),
        ]
    else:
        canvas = ROOT.TCanvas(f"c_appendix_B4_{tag}", "", 1700, 760)
        pads = [
            ROOT.TPad(f"pb4_small_{tag}", "", 0.00, 0.00, 0.50, 1.00),
            ROOT.TPad(f"pb4_large_{tag}", "", 0.50, 0.00, 1.00, 1.00),
        ]
    for pad in pads:
        pad.Draw()

    panel_specs = [
        (pads[0], "4C", "smaller", "m_small_4C", "prob_4C", "standard_4C_valid > 0.5", "M_{W} (smaller 4C dijet mass) [GeV]"),
        (pads[1], "4C", "larger", "m_large_4C", "prob_4C", "standard_4C_valid > 0.5", "M_{W} (larger 4C dijet mass) [GeV]"),
    ]
    if config["has_5c"]:
        panel_specs.append(
            (pads[2], "5C", "equal mass", "mW_5C", "prob_5C", "standard_5C_valid > 0.5", "M_{W} (5C equal dijet mass) [GeV]")
        )
    for (
        pad,
        fit_name,
        panel,
        mass_expression,
        probability_expression,
        valid_selection,
        x_title,
    ) in panel_specs:
        panel_rows, panel_keep = draw_before_after_mass_panel(
            pad,
            tree,
            config,
            fit_name,
            panel,
            mass_expression,
            probability_expression,
            valid_selection,
            x_title,
        )
        rows += panel_rows
        keep += panel_keep

    save_canvas(
        canvas,
        output_directory,
        f"appendix_B4_without_ISR_before_after_P003_ecm{tag}",
    )
    canvas.Close()
    return rows


def fit_pull_gaussian(histogram, name, color, line_style):
    gaussian = ROOT.TF1(name, "gaus", PULL_FIT_MIN, PULL_FIT_MAX)
    gaussian.SetParameters(
        histogram.GetMaximum(),
        histogram._mean,
        max(min(histogram._rms, 2.5), 0.25),
    )
    gaussian.SetLineColor(color)
    gaussian.SetLineWidth(3)
    gaussian.SetLineStyle(line_style)
    gaussian.SetNpx(500)
    fit_result = histogram.Fit(gaussian, "Q0RMS")
    gaussian._fit_status = int(fit_result)
    return gaussian


# UPDATE: Generalized pull plot function that works for all 4 categories
def draw_pull_jet_general(
    tree,
    config,
    jet,
    output_directory,
    pull_prefix,
    fit_label,
    valid_branch,
    probability_branch,
    output_suffix,
):
    """Draw pull plots for one jet and one pull category.

    Parameters
    ----------
    pull_prefix : str
        Branch prefix, e.g. "pull4C_standard", "pull4C_final",
        "pull5C_standard", or "pull5C_final".
    fit_label : str
        Fit name used in labels, e.g. "4C", "5C".
    valid_branch : str
        Branch name of the validity flag, e.g. "pull4C_standard_valid".
    probability_branch : str
        Branch name of the probability, e.g. "prob_4C", "prob_final_4C".
    output_suffix : str
        Filename suffix for the output files.
    """
    tag = config["tag"]
    canvas = ROOT.TCanvas(
        f"c_pull_{pull_prefix}_jet{jet}_{tag}", "", 1800, 1350
    )
    canvas.Divide(2, 2, 0.012, 0.012)

    before_selection = f"{valid_branch} > 0.5"
    after_selection = (
        f"({valid_branch} > 0.5) && ({probability_branch} > {P_CUT:.17g})"
    )
    rows, keep = [], []

    for pad_number, variable in enumerate(PULL_VARIABLES, start=1):
        canvas.cd(pad_number)
        configure_pad(
            ROOT.gPad,
            right=0.035,
            bottom=0.145,
            top=0.055,
            left=0.135,
        )
        branch = f"{pull_prefix}_{variable}_j{jet}"
        h_before = make_histogram(
            tree, branch, f"h{pull_prefix}_before_{variable}_j{jet}_{tag}",
            PULL_NBINS, PULL_XMIN, PULL_XMAX, before_selection,
        )
        h_after = make_histogram(
            tree, branch, f"h{pull_prefix}_after_{variable}_j{jet}_{tag}",
            PULL_NBINS, PULL_XMIN, PULL_XMAX, after_selection,
        )
        require_nonempty(h_before, f"{pull_prefix} pull {branch} before cut")
        require_nonempty(h_after, f"{pull_prefix} pull {branch} after cut")
        bin_width = h_before.GetBinWidth(1)
        if bin_width > 0.0:
            h_before.Scale(1.0 / bin_width)
            h_after.Scale(1.0 / bin_width)
        style_histogram(h_before, COLOR_BEFORE, line_style=1, width=3)
        style_histogram(h_after, COLOR_AFTER, width=3)
        g_before = fit_pull_gaussian(
            h_before, f"g{pull_prefix}_before_{variable}_j{jet}_{tag}",
            COLOR_BEFORE, 1,
        )
        g_after = fit_pull_gaussian(
            h_after, f"g{pull_prefix}_after_{variable}_j{jet}_{tag}",
            COLOR_AFTER, 1,
        )

        label = PULL_LABELS[variable]
        h_before.SetTitle("")
        configure_axes(h_before, f"Pull_{{{fit_label}}}({label})")
        h_before.SetMinimum(0.0)
        h_before.SetMaximum(
            1.20
            * max(
                h_before.GetMaximum(),
                h_after.GetMaximum(),
                g_before.GetMaximum(PULL_FIT_MIN, PULL_FIT_MAX),
                g_after.GetMaximum(PULL_FIT_MIN, PULL_FIT_MAX),
            )
        )
        h_before.Draw("HIST")
        h_after.Draw("HIST SAME")
        g_before.Draw("SAME")
        g_after.Draw("SAME")

        panel_title = ROOT.TLatex()
        panel_title.SetNDC()
        panel_title.SetTextFont(62)
        panel_title.SetTextSize(0.038)
        panel_title.SetTextAlign(13)
        panel_title.DrawLatex(0.160, 0.920, f"Jet {jet}: {label}")

        legend = ROOT.TLegend(0.16, 0.760, 0.48, 0.875)
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        legend.SetTextSize(0.032)
        legend.AddEntry(h_before, "Without P cut", "l")
        legend.AddEntry(
            h_after, f"With P_{{{fit_label}}} > 0.03", "l"
        )
        legend.Draw()

        stats_before = thesis_stat_box(
            0.660, 0.815, 0.940, 0.925,
            h_before, COLOR_BEFORE,
            precision=3, text_size=0.025, show_errors=True,
            mean=float(g_before.GetParameter(1)),
            mean_error=float(g_before.GetParError(1)),
            sigma=abs(float(g_before.GetParameter(2))),
            sigma_error=float(g_before.GetParError(2)),
        )
        stats_after = thesis_stat_box(
            0.660, 0.680, 0.940, 0.790,
            h_after, COLOR_AFTER,
            precision=3, text_size=0.025, show_errors=True,
            mean=float(g_after.GetParameter(1)),
            mean_error=float(g_after.GetParError(1)),
            sigma=abs(float(g_after.GetParameter(2))),
            sigma_error=float(g_after.GetParError(2)),
        )
        labels = plot_corner_labels(config)
        ROOT.gPad.RedrawAxis()
        keep += [
            h_before,
            h_after,
            g_before,
            g_after,
            panel_title,
            stats_before,
            stats_after,
            labels,
            legend,
        ]

        for stage, selection, histogram, gaussian in [
            ("before", before_selection, h_before, g_before),
            ("after", after_selection, h_after, g_after),
        ]:
            ndf = int(gaussian.GetNDF())
            rows.append(
                {
                    "energy_GeV": config["energy"],
                    "fit": fit_label,
                    "pull_category": pull_prefix,
                    "jet": jet,
                    "variable": variable,
                    "branch": branch,
                    "stage": stage,
                    "selection": selection,
                    "selected_events": histogram._selected_events,
                    "visible_fills": histogram._visible_fills,
                    "hist_mean": histogram._mean,
                    "hist_rms": histogram._rms,
                    "gaussian_mu": float(gaussian.GetParameter(1)),
                    "gaussian_mu_error": float(gaussian.GetParError(1)),
                    "gaussian_sigma": abs(float(gaussian.GetParameter(2))),
                    "gaussian_sigma_error": float(gaussian.GetParError(2)),
                    "gaussian_fit_status": gaussian._fit_status,
                    "gaussian_chi2_ndf": (
                        float(gaussian.GetChisquare()) / ndf if ndf > 0 else float("nan")
                    ),
                }
            )

    save_canvas(
        canvas,
        output_directory,
        f"figure_7_3_style_{output_suffix}_jet{jet}_ecm{tag}",
    )
    canvas.Close()
    return rows


def draw_4c_pull_jet(tree, config, jet, output_directory):
    """Backward-compatible wrapper for the original 4C standard pulls."""
    return draw_pull_jet_general(
        tree,
        config,
        jet,
        output_directory,
        pull_prefix="pull4C_standard",
        fit_label="4C",
        valid_branch="pull4C_standard_valid",
        probability_branch="prob_4C",
        output_suffix="4C_pulls_standard",
    )


# UPDATE: New functions for the additional pull plot categories
def draw_4c_final_pull_jet(tree, config, jet, output_directory):
    """Pull plots for final 4C (with ISR treatment where applied)."""
    return draw_pull_jet_general(
        tree,
        config,
        jet,
        output_directory,
        pull_prefix="pull4C_final",
        fit_label="4C",
        valid_branch="pull4C_final_valid",
        probability_branch="prob_final_4C",
        output_suffix="4C_pulls_final",
    )


def draw_5c_standard_pull_jet(tree, config, jet, output_directory):
    """Pull plots for standard 5C (without ISR treatment)."""
    return draw_pull_jet_general(
        tree,
        config,
        jet,
        output_directory,
        pull_prefix="pull5C_standard",
        fit_label="5C",
        valid_branch="pull5C_standard_valid",
        probability_branch="prob_5C",
        output_suffix="5C_pulls_standard",
    )


def draw_5c_final_pull_jet(tree, config, jet, output_directory):
    """Pull plots for final 5C (with ISR treatment where applied)."""
    return draw_pull_jet_general(
        tree,
        config,
        jet,
        output_directory,
        pull_prefix="pull5C_final",
        fit_label="5C",
        valid_branch="pull5C_final_valid",
        probability_branch="prob_final_5C",
        output_suffix="5C_pulls_final",
    )


def write_csv(path, rows):
    if not rows:
        return
    with open(path, "w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"Saved CSV : {path}")


def main():
    args = parse_args()
    config = CONFIGS[args.energy]
    input_path = find_input(config, args.input)
    output_directory = (
        os.path.abspath(os.path.expanduser(args.outdir))
        if args.outdir
        else os.path.join(SCRIPT_DIR, f"wmass_thesis_plots_ecm{config['tag']}")
    )
    os.makedirs(output_directory, exist_ok=True)

    root_file = ROOT.TFile.Open(input_path, "READ")
    if not root_file or root_file.IsZombie():
        raise SystemExit(f"ERROR: cannot open ROOT file:\n  {input_path}")
    tree = root_file.Get(TREE_NAME)
    if not tree:
        root_file.Close()
        raise SystemExit(f"ERROR: tree '{TREE_NAME}' not found in:\n  {input_path}")
    check_required_branches(tree, config)

    print("=" * 78)
    print(f"Input ROOT file : {input_path}")
    print(f"Energy          : {config['label']}")
    print(f"Tree entries    : {tree.GetEntries()}")
    print(f"5C enabled      : {'yes' if config['has_5c'] else 'no'}")
    print(f"Output directory: {output_directory}")
    print("=" * 78)

    figure_rows = []
    figure_rows += draw_figure_7_9(tree, config, output_directory)
    figure_rows += draw_figure_7_10(tree, config, output_directory)
    figure_rows += draw_appendix_b4(tree, config, output_directory)

    # UPDATE: Collect pull rows from all four pull categories
    pull_rows = []
    for jet in range(1, 5):
        # Standard 4C pulls (original)
        pull_rows += draw_4c_pull_jet(tree, config, jet, output_directory)
        # Final 4C pulls (new)
        pull_rows += draw_4c_final_pull_jet(tree, config, jet, output_directory)
        
        # 5C pulls only when 5C is enabled
        if config["has_5c"]:
            # Standard 5C pulls (new)
            pull_rows += draw_5c_standard_pull_jet(tree, config, jet, output_directory)
            # Final 5C pulls (new)
            pull_rows += draw_5c_final_pull_jet(tree, config, jet, output_directory)

    write_csv(
        os.path.join(output_directory, f"wmass_thesis_summary_ecm{config['tag']}.csv"),
        figure_rows,
    )
    write_csv(
        os.path.join(output_directory, f"pull_summary_ecm{config['tag']}.csv"),
        pull_rows,
    )
    root_file.Close()

    print("=" * 78)
    print("Created: Figure 7.9, Figure 7.10, Appendix B.4")
    print("Created: 4C standard pull plots (4 images, one per jet)")
    print("Created: 4C final pull plots (4 images, one per jet)")
    if config["has_5c"]:
        print("Created: 5C standard pull plots (4 images, one per jet)")
        print("Created: 5C final pull plots (4 images, one per jet)")
        print("Created: 16 pull plot images in total")
    else:
        print("Created: 8 pull plot images in total")
        print("162.5-GeV rule: no 5C branches were read or plotted")
    print("=" * 78)


if __name__ == "__main__":
    main()