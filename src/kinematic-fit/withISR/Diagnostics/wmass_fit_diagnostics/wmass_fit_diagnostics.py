#!/usr/bin/env python3
"""Final W-mass diagnostics focused on the 5C and ISR questions.

This script is intentionally separate from the reconstruction producer and the
plotting script. It reads the existing ROOT trees and writes a compact README
plus CSV tables with only the numbers needed to discuss:

1. plot-style means and sigmas, computed with the same visible mass window as
   the plotting code;
2. whether the 5C problem follows the pre-5C mass splitting;
3. whether the ISR refit recovers low-probability events and tracks truth.

The 5C fit is skipped as a physics result at 162.5 GeV, matching the thesis
logic and the plotting script's `has_5c = False` setting.

python3 wmass_fit_diagnostics.py \
  --input 162p5=wmass_fit_pvalue_pulls_ISR_ecm162p5.root \
  --input 240=wmass_fit_pvalue_pulls_ISR_ecm240.root \
  --input 365=wmass_fit_pvalue_pulls_ISR_ecm365.root \
  --outdir wmass_final_diagnostics
"""

from __future__ import annotations

import argparse
import csv
import math
import statistics
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import ROOT

ROOT.gROOT.SetBatch(True)

MASS_XMIN = 40.0
MASS_XMAX = 120.0
P_CUT = 0.03
MW_REF = 80.385
TREE_NAME = "events"


def energy_key(label: str) -> str:
    text = label.strip().lower()
    if text in {"162.5", "162p5", "162p6", "162.6"}:
        return "162p5"
    if text in {"240", "240.0"}:
        return "240"
    if text in {"365", "365.0"}:
        return "365"
    return text


def has_5c_physics(label: str) -> bool:
    return energy_key(label) in {"240", "365"}


def ecm_value(label: str) -> float | None:
    key = energy_key(label)
    if key == "162p5":
        return 162.5
    if key == "240":
        return 240.0
    if key == "365":
        return 365.0
    try:
        return float(label)
    except ValueError:
        return None


def finite_expr(expr: str) -> str:
    return f"(({expr})==({expr}) && abs({expr})<1.0e20)"


def and_cut(*parts: str) -> str:
    clean = [f"({p})" for p in parts if p and p.strip() and p.strip() != "1"]
    return " && ".join(clean) if clean else "1"


def fmt(x: float | int | str, digits: int = 3) -> str:
    try:
        value = float(x)
    except (TypeError, ValueError):
        return str(x)
    if not math.isfinite(value):
        return "n/a"
    return f"{value:.{digits}f}"


def pct(n: float, d: float) -> float:
    return 100.0 * n / d if d else float("nan")


def quantile(values: Sequence[float], q: float) -> float:
    if not values:
        return float("nan")
    ordered = sorted(values)
    pos = (len(ordered) - 1) * q
    lo = math.floor(pos)
    hi = math.ceil(pos)
    if lo == hi:
        return ordered[int(lo)]
    return ordered[int(lo)] * (hi - pos) + ordered[int(hi)] * (pos - lo)


def mean(values: Sequence[float]) -> float:
    return statistics.fmean(values) if values else float("nan")


def stddev(values: Sequence[float]) -> float:
    if not values:
        return float("nan")
    m = statistics.fmean(values)
    return math.sqrt(sum((x - m) ** 2 for x in values) / len(values))


def corr_stats(x: Sequence[float], y: Sequence[float]) -> tuple[float, float, float]:
    if len(x) != len(y) or len(x) < 2:
        return float("nan"), float("nan"), float("nan")
    mx = statistics.fmean(x)
    my = statistics.fmean(y)
    vx = sum((v - mx) ** 2 for v in x)
    vy = sum((v - my) ** 2 for v in y)
    cov = sum((a - mx) * (b - my) for a, b in zip(x, y))
    rho = cov / math.sqrt(vx * vy) if vx > 0 and vy > 0 else float("nan")
    bias = my - mx
    residual = math.sqrt(sum((b - a) ** 2 for a, b in zip(x, y)) / len(x))
    return rho, bias, residual


@dataclass
class Sample:
    label: str
    path: Path
    file: ROOT.TFile
    tree: ROOT.TTree

    @property
    def branches(self) -> set[str]:
        return {b.GetName() for b in self.tree.GetListOfBranches()}

    def has(self, *names: str) -> bool:
        b = self.branches
        return all(name in b for name in names)

    def count(self, cut: str = "1") -> int:
        return int(self.tree.GetEntries(cut))

    def values(self, expr: str, cut: str = "1") -> list[float]:
        selection = and_cut(cut, finite_expr(expr))
        n = int(self.tree.Draw(expr, selection, "goff"))
        return [float(self.tree.GetV1()[i]) for i in range(n)]

    def values2(self, xexpr: str, yexpr: str, cut: str = "1") -> tuple[list[float], list[float]]:
        selection = and_cut(cut, finite_expr(xexpr), finite_expr(yexpr))
        n = int(self.tree.Draw(f"{yexpr}:{xexpr}", selection, "goff"))
        y = [float(self.tree.GetV1()[i]) for i in range(n)]
        x = [float(self.tree.GetV2()[i]) for i in range(n)]
        return x, y


def plot_stats(sample: Sample, expressions: Sequence[str], cut: str) -> dict[str, float]:
    visible: list[float] = []
    selected = sample.count(cut)
    finite = 0
    for expr in expressions:
        vals = sample.values(expr, cut)
        finite += len(vals)
        visible.extend(v for v in vals if MASS_XMIN <= v < MASS_XMAX)
    return {
        "selected_events": selected,
        "finite_fills": finite,
        "visible_fills": len(visible),
        "mean": mean(visible),
        "sigma": stddev(visible),
        "outside_visible_pct": pct(finite - len(visible), finite),
    }


def add_plot_rows(sample: Sample) -> list[dict]:
    rows: list[dict] = []
    h5 = has_5c_physics(sample.label)
    specs = [
        ("Fig 7.9", "mass", "without ISR treatment", ["mW_5C"] if h5 else ["m_small_4C", "m_large_4C"],
         "standard_5C_valid>0.5" if h5 else "standard_4C_valid>0.5"),
        ("Fig 7.9", "mass", "with ISR treatment", ["mW_5C_final"] if h5 else ["m_small_4C_final", "m_large_4C_final"],
         "standard_5C_valid>0.5" if h5 else "standard_4C_valid>0.5"),
        ("Fig 7.10", "smaller", "raw", ["m_small_raw"], "1"),
        ("Fig 7.10", "larger", "raw", ["m_large_raw"], "1"),
        ("Fig 7.10", "smaller", "4C with ISR treatment", ["m_small_4C_final"], "standard_4C_valid>0.5"),
        ("Fig 7.10", "larger", "4C with ISR treatment", ["m_large_4C_final"], "standard_4C_valid>0.5"),
    ]
    if h5:
        specs += [
            ("Fig 7.10", "smaller", "5C with ISR treatment", ["mW_5C_final"], "standard_5C_valid>0.5"),
            ("Fig 7.10", "larger", "5C with ISR treatment", ["mW_5C_final"], "standard_5C_valid>0.5"),
            ("Appendix B.4", "5C equal mass", "without P cut", ["mW_5C"], "standard_5C_valid>0.5"),
            ("Appendix B.4", "5C equal mass", "with P>0.03", ["mW_5C"], "standard_5C_valid>0.5 && prob_5C>0.03"),
            ("Final selected", "5C equal mass", "with final P>0.03", ["mW_5C_final"], "pass_final_5C_p03>0.5"),
        ]
    specs += [
        ("Appendix B.4", "4C smaller", "without P cut", ["m_small_4C"], "standard_4C_valid>0.5"),
        ("Appendix B.4", "4C larger", "without P cut", ["m_large_4C"], "standard_4C_valid>0.5"),
        ("Appendix B.4", "4C smaller", "with P>0.03", ["m_small_4C"], "standard_4C_valid>0.5 && prob_4C>0.03"),
        ("Appendix B.4", "4C larger", "with P>0.03", ["m_large_4C"], "standard_4C_valid>0.5 && prob_4C>0.03"),
    ]
    for figure, panel, curve, exprs, cut in specs:
        missing = [expr for expr in exprs if not sample.has(expr)]
        if missing:
            continue
        stats = plot_stats(sample, exprs, cut)
        rows.append({"energy": sample.label, "figure": figure, "panel": panel,
                     "curve": curve, "expressions": ";".join(exprs),
                     "selection": cut, **stats})
    return rows


def add_cutflow_rows(sample: Sample) -> list[dict]:
    rows = []
    stored = sample.count()
    valid4 = sample.count("standard_4C_valid>0.5")
    final4 = sample.count("pass_final_4C_p03>0.5") if sample.has("pass_final_4C_p03") else 0
    rows += [
        {"energy": sample.label, "quantity": "stored events", "events": stored, "denominator": stored, "percent": 100.0},
        {"energy": sample.label, "quantity": "standard 4C P>0.03", "events": sample.count("standard_4C_valid>0.5 && prob_4C>0.03"), "denominator": valid4, "percent": pct(sample.count("standard_4C_valid>0.5 && prob_4C>0.03"), valid4)},
        {"energy": sample.label, "quantity": "final 4C accepted", "events": final4, "denominator": stored, "percent": pct(final4, stored)},
    ]
    if has_5c_physics(sample.label):
        valid5 = sample.count("standard_5C_valid>0.5")
        pass5 = sample.count("standard_5C_valid>0.5 && prob_5C>0.03")
        final5 = sample.count("pass_final_5C_p03>0.5") if sample.has("pass_final_5C_p03") else 0
        rows += [
            {"energy": sample.label, "quantity": "standard 5C P>0.03", "events": pass5, "denominator": valid5, "percent": pct(pass5, valid5)},
            {"energy": sample.label, "quantity": "final 5C accepted", "events": final5, "denominator": stored, "percent": pct(final5, stored)},
        ]
    return rows


def add_mass_split_rows(sample: Sample) -> list[dict]:
    if not has_5c_physics(sample.label) or not sample.has("delta_m_4C_pre5C"):
        return []
    rows = []
    bins = [
        ("0-2", "delta_m_4C_pre5C>=0 && delta_m_4C_pre5C<2"),
        ("2-5", "delta_m_4C_pre5C>=2 && delta_m_4C_pre5C<5"),
        ("5-10", "delta_m_4C_pre5C>=5 && delta_m_4C_pre5C<10"),
        ("10-20", "delta_m_4C_pre5C>=10 && delta_m_4C_pre5C<20"),
        (">=20", "delta_m_4C_pre5C>=20"),
    ]
    base = "standard_4C_valid>0.5 && standard_5C_valid>0.5"
    total = sample.count(base)
    for label, region in bins:
        cut = and_cut(base, region)
        n = sample.count(cut)
        if n == 0:
            continue
        prob5 = sample.values("prob_5C", cut)
        probf = sample.values("prob_final_5C", cut) if sample.has("prob_final_5C") else []
        delta = sample.values("delta_m_4C_pre5C", cut)
        m5 = plot_stats(sample, ["mW_5C"], cut)
        m5f = plot_stats(sample, ["mW_5C_final"], cut) if sample.has("mW_5C_final") else {}
        applied = sample.count(and_cut(cut, "isr_applied>0.5")) if sample.has("isr_applied") else 0
        recovered = sample.count(and_cut(cut, "isr_5C_recovered>0.5")) if sample.has("isr_5C_recovered") else 0
        low_before = sum(p < P_CUT for p in prob5)
        low_after = sum(p < P_CUT for p in probf) if probf else 0
        rows.append({
            "energy": sample.label,
            "delta_m_4C_bin_GeV": label,
            "events": n,
            "fraction_of_valid_4C5C_pct": pct(n, total),
            "mean_delta_m_4C_GeV": mean(delta),
            "p5c_lt_003_pct": pct(low_before, len(prob5)),
            "pfinal5c_lt_003_pct": pct(low_after, len(probf)) if probf else float("nan"),
            "mean_prob_5C": mean(prob5),
            "mean_prob_final_5C": mean(probf),
            "sigma_mW_5C_plot_GeV": m5["sigma"],
            "sigma_mW_5C_final_plot_GeV": m5f.get("sigma", float("nan")),
            "isr_applied": applied,
            "isr_recovered": recovered,
            "recovered_over_applied_pct": pct(recovered, applied),
        })
    return rows


def add_isr_rows(sample: Sample) -> list[dict]:
    rows = []
    if sample.has("pz_isr_4C_fitted"):
        pass
    else:
        sample.tree.SetAlias("pz_isr_4C_fitted",
                             "(y_isr_4C_fitted>=0 ? E_isr_4C_fitted : -E_isr_4C_fitted)")
    if sample.has("pz_isr_5C_fitted"):
        pass
    else:
        sample.tree.SetAlias("pz_isr_5C_fitted",
                             "(y_isr_5C_fitted>=0 ? E_isr_fitted : -E_isr_fitted)")
    slices = [
        ("all", "1"),
        ("<5", "E_isr_true_collinear<5"),
        ("5-20", "E_isr_true_collinear>=5 && E_isr_true_collinear<20"),
        ("20-50", "E_isr_true_collinear>=20 && E_isr_true_collinear<50"),
        (">=50", "E_isr_true_collinear>=50"),
    ]
    fits = [
        ("4C+ISR", "standard_4C_valid", "prob_4C", "isr_4C_applied", "isr_4C_recovered",
         "pz_isr_4C_fitted", "E_isr_4C_fitted"),
    ]
    if has_5c_physics(sample.label):
        fits.append(("5C+ISR", "standard_5C_valid", "prob_5C", "isr_applied", "isr_5C_recovered",
                     "pz_isr_5C_fitted", "E_isr_fitted"))
    for fit, valid, prob, applied, recovered, fit_pz, fit_e in fits:
        if not sample.has(valid, prob, applied, recovered, "pz_isr_true", "E_isr_true_collinear"):
            continue
        for slabel, scut in slices:
            standard = and_cut(f"{valid}>0.5", scut)
            eligible = and_cut(standard, f"{prob}>0", f"{prob}<0.03")
            appcut = and_cut(standard, f"{applied}>0.5")
            n_standard = sample.count(standard)
            n_eligible = sample.count(eligible)
            n_applied = sample.count(appcut)
            n_recovered = sample.count(and_cut(standard, f"{recovered}>0.5"))
            true_pz, fitted_pz = sample.values2("pz_isr_true", fit_pz, appcut)
            rho_pz, bias_pz, rms_pz = corr_stats(true_pz, fitted_pz)
            true_e, fitted_e = sample.values2("E_isr_true_collinear", fit_e, appcut)
            rho_e, bias_e, rms_e = corr_stats(true_e, fitted_e)
            rows.append({
                "energy": sample.label,
                "fit": fit,
                "true_isr_slice_GeV": slabel,
                "standard_valid": n_standard,
                "eligible_p_lt_003": n_eligible,
                "applied": n_applied,
                "recovered": n_recovered,
                "eligible_fraction_pct": pct(n_eligible, n_standard),
                "recovered_over_applied_pct": pct(n_recovered, n_applied),
                "rho_signed_pz": rho_pz,
                "pz_bias_GeV": bias_pz,
                "pz_residual_rms_GeV": rms_pz,
                "rho_energy": rho_e,
                "energy_residual_rms_GeV": rms_e,
            })
    return rows


def add_emax_rows(sample: Sample) -> list[dict]:
    if not sample.has("Emax_isr"):
        return []
    ecm = ecm_value(sample.label)
    values = sample.values("Emax_isr", "1")
    s = ecm * ecm if ecm else float("nan")
    single_w = 0.5 * ecm * (1.0 - MW_REF * MW_REF / s) if ecm else float("nan")
    ww_pair = max(0.0, 0.5 * ecm * (1.0 - 4.0 * MW_REF * MW_REF / s)) if ecm else float("nan")
    return [{
        "energy": sample.label,
        "stored_Emax_mean_GeV": mean(values),
        "stored_Emax_median_GeV": quantile(values, 0.5),
        "formula_single_W_mass_GeV": single_w,
        "formula_WW_pair_threshold_GeV": ww_pair,
        "warning": "If mW_ref is one W mass, the WW-pair threshold formula should be used."
    }]


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    fields = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def table(lines: list[str], headers: Sequence[str], rows: Iterable[Sequence[str]]) -> None:
    lines.append("| " + " | ".join(headers) + " |")
    lines.append("|" + "|".join(["---"] * len(headers)) + "|")
    for row in rows:
        lines.append("| " + " | ".join(row) + " |")
    lines.append("")


def generate_readme(outdir: Path, cutflow: list[dict], plot_rows: list[dict],
                    split_rows: list[dict], isr_rows: list[dict], emax_rows: list[dict]) -> None:
    lines: list[str] = []
    lines += [
        "# Final W-Mass Diagnostic Report",
        "",
        "This report keeps only the numbers needed to understand the current 5C and ISR problems.",
        "",
        "Important conventions:",
        "",
        "- The 162.5 GeV sample is treated as 4C-only for physics conclusions. Any 5C branches in the ROOT tree are ignored here.",
        "- Plot means and sigmas are computed from the same visible mass window as the plotting code: 40 <= M < 120 GeV.",
        "- `with ISR treatment` means the final hybrid variable before applying the final 3% rejection, matching Figures 7.9 and 7.10.",
        "- The final accepted rows use `P_final > 0.03`.",
        "- Full RMS values are intentionally omitted because they overemphasize extreme tails and do not match the plotted sigma boxes.",
        "",
        "Generated CSV files:",
        "",
        "- `plot_style_numbers.csv`: numbers matching the plot definitions.",
        "- `cutflow_summary.csv`: compact acceptance and rejection summary.",
        "- `mass_split_stress_test.csv`: whether the 5C problem follows the pre-5C 4C mass difference.",
        "- `isr_truth_recovery.csv`: ISR recovery and truth-correlation checks.",
        "- `emax_check.csv`: check of the stored ISR photon energy limit.",
        "",
    ]

    lines += ["## Executive Summary", ""]
    table(lines, ["Energy", "standard 5C P>0.03", "final 5C accepted", "note"], (
        [r["energy"], f"{fmt(r['percent'], 2)}%", f"{int(r['events'])}/{int(r['denominator'])}",
         "not used at threshold" if not has_5c_physics(r["energy"]) else "5C physics sample"]
        for r in cutflow if r["quantity"] in {"standard 5C P>0.03", "final 5C accepted"}
    ))

    fig710 = [r for r in plot_rows if r["figure"] == "Fig 7.10" and r["curve"] in {"raw", "4C with ISR treatment", "5C with ISR treatment"}]
    table(lines, ["Energy", "Panel", "Curve", "Mean [GeV]", "Plot Sigma [GeV]", "Visible entries"], (
        [r["energy"], r["panel"], r["curve"], fmt(r["mean"], 2), fmt(r["sigma"], 2), str(int(r["visible_fills"]))]
        for r in fig710
    ))

    lines += [
        "## 5C Stress Test",
        "",
        "This is the main diagnostic for the current problem. If the low 5C probability and large 5C sigma increase with `delta_m_4C_pre5C`, the equal-mass constraint is being forced onto events whose two W candidates are already incompatible after 4C.",
        "",
    ]
    table(lines, ["Energy", "Delta m bin [GeV]", "Events", "P5C<0.03", "Pfinal<0.03", "Sigma 5C final [GeV]", "ISR recovered/applied"], (
        [r["energy"], r["delta_m_4C_bin_GeV"], str(int(r["events"])),
         f"{fmt(r['p5c_lt_003_pct'], 1)}%", f"{fmt(r['pfinal5c_lt_003_pct'], 1)}%",
         fmt(r["sigma_mW_5C_final_plot_GeV"], 2), f"{fmt(r['recovered_over_applied_pct'], 1)}%"]
        for r in split_rows
    ))

    lines += [
        "## ISR Truth And Recovery",
        "",
        "The most useful ISR numbers are the eligible fraction, recovery fraction, signed-pz correlation, and pz residual RMS. The signed-pz correlation tests whether the fitted ISR photon follows the true ISR direction event by event.",
        "",
    ]
    compact_isr = [r for r in isr_rows if r["true_isr_slice_GeV"] in {"all", "<5", ">=50"}]
    table(lines, ["Energy", "Fit", "True ISR", "Eligible", "Recovered/applied", "rho signed pz", "pz RMS [GeV]"], (
        [r["energy"], r["fit"], r["true_isr_slice_GeV"], f"{fmt(r['eligible_fraction_pct'], 1)}%",
         f"{fmt(r['recovered_over_applied_pct'], 1)}%", fmt(r["rho_signed_pz"], 3),
         fmt(r["pz_residual_rms_GeV"], 2)]
        for r in compact_isr
    ))

    lines += [
        "## ISR Emax Check",
        "",
        "This table is included because a wrong ISR energy limit can make the ISR refit look mathematically valid while being physically too flexible, especially near threshold.",
        "",
    ]
    table(lines, ["Energy", "Stored Emax [GeV]", "single-W formula [GeV]", "WW-threshold formula [GeV]"], (
        [r["energy"], fmt(r["stored_Emax_median_GeV"], 3),
         fmt(r["formula_single_W_mass_GeV"], 3), fmt(r["formula_WW_pair_threshold_GeV"], 3)]
        for r in emax_rows
    ))

    lines += [
        "## How To Read This",
        "",
        "Use `plot_style_numbers.csv` when comparing to the plot boxes. Use `mass_split_stress_test.csv` when discussing why 5C fails. Use `isr_truth_recovery.csv` when discussing whether the ISR photon fit is physically meaningful.",
        "",
        "The central question is not only whether ISR improves the mass plot. The sharper question is whether the equal-mass 5C constraint is being applied to events whose two 4C masses are already too different, and whether the ISR parameter is then absorbing that tension.",
        "",
    ]

    (outdir / "README.md").write_text("\n".join(lines), encoding="utf-8")


def parse_input(text: str) -> tuple[str, Path]:
    if "=" not in text:
        raise argparse.ArgumentTypeError("use LABEL=/path/to/file.root")
    label, filename = text.split("=", 1)
    return label, Path(filename).expanduser().resolve()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", action="append", required=True, type=parse_input,
                        metavar="LABEL=FILE", help="repeat, e.g. 365=/path/file.root")
    parser.add_argument("--tree", default=TREE_NAME)
    parser.add_argument("--outdir", default="wmass_final_diagnostics")
    args = parser.parse_args()

    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    samples: list[Sample] = []
    for label, path in args.input:
        root_file = ROOT.TFile.Open(str(path))
        if not root_file or root_file.IsZombie():
            raise RuntimeError(f"cannot open {path}")
        tree = root_file.Get(args.tree)
        if not tree:
            raise RuntimeError(f"cannot find tree '{args.tree}' in {path}")
        tree.SetEstimate(max(int(tree.GetEntries()) + 1, 1_000_000))
        samples.append(Sample(label=energy_key(label), path=path, file=root_file, tree=tree))

    cutflow: list[dict] = []
    plot_rows: list[dict] = []
    split_rows: list[dict] = []
    isr_rows: list[dict] = []
    emax_rows: list[dict] = []
    for sample in samples:
        cutflow.extend(add_cutflow_rows(sample))
        plot_rows.extend(add_plot_rows(sample))
        split_rows.extend(add_mass_split_rows(sample))
        isr_rows.extend(add_isr_rows(sample))
        emax_rows.extend(add_emax_rows(sample))

    write_csv(outdir / "cutflow_summary.csv", cutflow)
    write_csv(outdir / "plot_style_numbers.csv", plot_rows)
    write_csv(outdir / "mass_split_stress_test.csv", split_rows)
    write_csv(outdir / "isr_truth_recovery.csv", isr_rows)
    write_csv(outdir / "emax_check.csv", emax_rows)
    generate_readme(outdir, cutflow, plot_rows, split_rows, isr_rows, emax_rows)

    for sample in samples:
        sample.file.Close()

    print(f"Wrote {outdir / 'README.md'}")
    print(f"Wrote {outdir / 'plot_style_numbers.csv'}")
    print(f"Wrote {outdir / 'mass_split_stress_test.csv'}")
    print(f"Wrote {outdir / 'isr_truth_recovery.csv'}")
    print(f"Wrote {outdir / 'emax_check.csv'}")


if __name__ == "__main__":
    main()
