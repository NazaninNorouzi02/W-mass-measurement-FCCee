#!/usr/bin/env python3
import os
import ROOT

ROOT.gROOT.SetBatch(True)

# -----------------------------------------------------------------------------
# Input / output
# -----------------------------------------------------------------------------
INPUT_FILE = "outputs/matching_comparison/matching_comparison_ecm160.root"
TREE_NAME  = "events"
OUTDIR     = "outputs/matching_comparison/plots1/"

os.makedirs(OUTDIR, exist_ok=True)

# -----------------------------------------------------------------------------
# Style
# -----------------------------------------------------------------------------
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetTitleBorderSize(0)
ROOT.gStyle.SetTitleFillColor(0)
ROOT.gStyle.SetPadTopMargin(0.12)
ROOT.gStyle.SetPadBottomMargin(0.13)
ROOT.gStyle.SetPadLeftMargin(0.14)
ROOT.gStyle.SetPadRightMargin(0.05)
ROOT.gStyle.SetTitleSize(0.045, "XY")
ROOT.gStyle.SetLabelSize(0.04, "XY")
ROOT.gStyle.SetTitleOffset(1.25, "Y")
ROOT.gStyle.SetLegendBorderSize(0)

# -----------------------------------------------------------------------------
# Open file / tree
# -----------------------------------------------------------------------------
f = ROOT.TFile.Open(INPUT_FILE)
if not f or f.IsZombie():
    raise RuntimeError(f"Could not open file: {INPUT_FILE}")

tree = f.Get(TREE_NAME)
if not tree:
    raise RuntimeError(f"Could not find tree '{TREE_NAME}' in {INPUT_FILE}")

print(f"[INFO] Opened file: {INPUT_FILE}")
print(f"[INFO] Tree: {TREE_NAME}")
print(f"[INFO] Entries: {tree.GetEntries()}")

# -----------------------------------------------------------------------------
# Strategy style map
# -----------------------------------------------------------------------------
styles = {
    "A": {
        "label": "A: non-exclusive",
        "color": ROOT.kRed + 1,
        "marker": 20,
        "line": 1
    },
    "B": {
        "label": "B: unified greedy",
        "color": ROOT.kBlue + 1,
        "marker": 21,
        "line": 2
    },
    "C": {
        "label": "C: transitive",
        "color": ROOT.kGreen + 2,
        "marker": 22,
        "line": 3
    },
}

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
def save_canvas(c, name):
    png = os.path.join(OUTDIR, f"{name}.png")
    pdf = os.path.join(OUTDIR, f"{name}.pdf")
    c.SaveAs(png)
    c.SaveAs(pdf)
    print(f"[SAVED] {png}")
    print(f"[SAVED] {pdf}")

def setup_hist(h, color, marker, line_style=1, line_width=2):
    h.SetLineColor(color)
    h.SetMarkerColor(color)
    h.SetMarkerStyle(marker)
    h.SetMarkerSize(1.0)
    h.SetLineStyle(line_style)
    h.SetLineWidth(line_width)

def make_latex():
    lat = ROOT.TLatex()
    lat.SetNDC()
    lat.SetTextFont(42)
    lat.SetTextSize(0.035)
    return lat

def make_legend(x1=0.58, y1=0.68, x2=0.88, y2=0.88):
    leg = ROOT.TLegend(x1, y1, x2, y2)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.032)
    return leg

def get_hist_mean_rms(branch, cut="", hname="htmp", nbins=200, xmin=-2, xmax=2):
    h = ROOT.TH1F(hname, "", nbins, xmin, xmax)
    draw_expr = f"{branch}>>{hname}"
    tree.Draw(draw_expr, cut, "goff")
    mean = h.GetMean()
    rms  = h.GetRMS()
    entries = h.GetEntries()
    return mean, rms, entries

# -----------------------------------------------------------------------------
# Plot 1: delta_alpha distributions per jet, comparing A/B/C
# -----------------------------------------------------------------------------
for jidx in range(1, 5):
    c = ROOT.TCanvas(f"c_delta_alpha_j{jidx}", "", 900, 700)

    leg = make_legend(0.57, 0.67, 0.88, 0.87)
    hists = []
    maxy = 0.0

    for s in ["A", "B", "C"]:
        branch = f"delta_alpha_{s}_j{jidx}"
        hname = f"h_delta_alpha_{s}_j{jidx}"
        h = ROOT.TH1F(hname, "", 60, -1.5, 1.5)

        cut = f"{branch} > -998"
        tree.Draw(f"{branch}>>{hname}", cut, "goff")

        if h.Integral() > 0:
            h.Scale(1.0 / h.Integral())

        setup_hist(
            h,
            styles[s]["color"],
            styles[s]["marker"],
            styles[s]["line"],
            3
        )

        h.SetTitle("")
        h.GetXaxis().SetTitle(f"#Delta#alpha (jet {jidx})")
        h.GetYaxis().SetTitle("Normalized entries")

        if h.GetMaximum() > maxy:
            maxy = h.GetMaximum()

        hists.append((s, h))

    for i, (s, h) in enumerate(hists):
        h.SetMaximum(maxy * 1.25 if maxy > 0 else 1.0)
        drawopt = "HIST" if i == 0 else "HIST SAME"
        h.Draw(drawopt)
        leg.AddEntry(h, styles[s]["label"], "l")

    leg.Draw()

    lat = make_latex()
    lat.DrawLatex(0.15, 0.94, "#sqrt{s} = 160 GeV, WW #rightarrow 4q")
    lat.DrawLatex(0.15, 0.89, f"Jet {jidx}: #Delta#alpha comparison")

    c.RedrawAxis()
    save_canvas(c, f"delta_alpha_compare_jet{jidx}")

# -----------------------------------------------------------------------------
# Plot 2: matched-parton multiplicity comparison
# FIXED: branch names are n_matched_A/B/C, not n_partons_matched_A/B/C
# -----------------------------------------------------------------------------
c_mult = ROOT.TCanvas("c_mult", "", 900, 700)
leg_mult = make_legend(0.57, 0.68, 0.88, 0.88)

mult_hists = []
maxy = 0.0

for s in ["A", "B", "C"]:
    branch = f"n_matched_{s}"
    hname = f"h_nmatched_{s}"
    h = ROOT.TH1F(hname, "", 6, -0.5, 5.5)

    tree.Draw(f"{branch}>>{hname}", "", "goff")

    if h.Integral() > 0:
        h.Scale(1.0 / h.Integral())

    setup_hist(
        h,
        styles[s]["color"],
        styles[s]["marker"],
        styles[s]["line"],
        3
    )

    h.SetTitle("")
    h.GetXaxis().SetTitle("Number of matched partons")
    h.GetYaxis().SetTitle("Fraction of events")

    if h.GetMaximum() > maxy:
        maxy = h.GetMaximum()

    mult_hists.append((s, h))

for i, (s, h) in enumerate(mult_hists):
    h.SetMaximum(maxy * 1.25 if maxy > 0 else 1.0)
    drawopt = "HIST" if i == 0 else "HIST SAME"
    h.Draw(drawopt)
    leg_mult.AddEntry(h, styles[s]["label"], "l")

leg_mult.Draw()
lat = make_latex()
lat.DrawLatex(0.15, 0.94, "#sqrt{s} = 160 GeV, WW #rightarrow 4q")
lat.DrawLatex(0.15, 0.89, "Matched-parton multiplicity comparison")

c_mult.RedrawAxis()
save_canvas(c_mult, "matched_parton_multiplicity_comparison")

# -----------------------------------------------------------------------------
# Plot 3: summary graph of mean(delta_alpha) and RMS(delta_alpha) vs jet index
# -----------------------------------------------------------------------------
jet_x = [1.0, 2.0, 3.0, 4.0]

graphs_mean = {}
graphs_rms  = {}

for s in ["A", "B", "C"]:
    g_mean = ROOT.TGraph(4)
    g_rms  = ROOT.TGraph(4)

    g_mean.SetName(f"g_mean_{s}")
    g_rms.SetName(f"g_rms_{s}")

    g_mean.SetLineColor(styles[s]["color"])
    g_mean.SetMarkerColor(styles[s]["color"])
    g_mean.SetMarkerStyle(styles[s]["marker"])
    g_mean.SetMarkerSize(1.2)
    g_mean.SetLineStyle(styles[s]["line"])
    g_mean.SetLineWidth(3)

    g_rms.SetLineColor(styles[s]["color"])
    g_rms.SetMarkerColor(styles[s]["color"])
    g_rms.SetMarkerStyle(styles[s]["marker"])
    g_rms.SetMarkerSize(1.2)
    g_rms.SetLineStyle(styles[s]["line"])
    g_rms.SetLineWidth(3)

    for i, jidx in enumerate(range(1, 5)):
        branch = f"delta_alpha_{s}_j{jidx}"
        mean, rms, _ = get_hist_mean_rms(
            branch,
            cut=f"{branch} > -998",
            hname=f"h_tmp_{s}_{jidx}",
            nbins=200,
            xmin=-2,
            xmax=2
        )
        g_mean.SetPoint(i, jet_x[i], mean)
        g_rms.SetPoint(i, jet_x[i], rms)

    graphs_mean[s] = g_mean
    graphs_rms[s]  = g_rms

# --- mean plot
c_mean = ROOT.TCanvas("c_mean", "", 900, 700)
frame_mean = ROOT.TH1F("frame_mean", "", 4, 0.5, 4.5)
frame_mean.SetMinimum(-0.2)
frame_mean.SetMaximum(0.2)
frame_mean.SetTitle("")
frame_mean.GetXaxis().SetTitle("Jet index")
frame_mean.GetYaxis().SetTitle("Mean(#Delta#alpha)")
frame_mean.Draw()

leg_mean = make_legend(0.56, 0.67, 0.88, 0.87)

for s in ["A", "B", "C"]:
    graphs_mean[s].Draw("PL SAME")
    leg_mean.AddEntry(graphs_mean[s], styles[s]["label"], "lp")

leg_mean.Draw()
lat = make_latex()
lat.DrawLatex(0.15, 0.94, "#sqrt{s} = 160 GeV, WW #rightarrow 4q")
lat.DrawLatex(0.15, 0.89, "Mean of #Delta#alpha per jet")

c_mean.RedrawAxis()
save_canvas(c_mean, "delta_alpha_mean_summary")

# --- RMS plot
c_rms = ROOT.TCanvas("c_rms", "", 900, 700)
frame_rms = ROOT.TH1F("frame_rms", "", 4, 0.5, 4.5)
frame_rms.SetMinimum(0.0)
frame_rms.SetMaximum(1.0)
frame_rms.SetTitle("")
frame_rms.GetXaxis().SetTitle("Jet index")
frame_rms.GetYaxis().SetTitle("RMS(#Delta#alpha)")
frame_rms.Draw()

leg_rms = make_legend(0.56, 0.67, 0.88, 0.87)

for s in ["A", "B", "C"]:
    graphs_rms[s].Draw("PL SAME")
    leg_rms.AddEntry(graphs_rms[s], styles[s]["label"], "lp")

leg_rms.Draw()
lat = make_latex()
lat.DrawLatex(0.15, 0.94, "#sqrt{s} = 160 GeV, WW #rightarrow 4q")
lat.DrawLatex(0.15, 0.89, "RMS of #Delta#alpha per jet")

c_rms.RedrawAxis()
save_canvas(c_rms, "delta_alpha_rms_summary")

# -----------------------------------------------------------------------------
# Plot 4: event energy balance
# -----------------------------------------------------------------------------
c_bal = ROOT.TCanvas("c_bal", "", 900, 700)
h_bal = ROOT.TH1F("h_bal", "", 60, 0.0, 2.0)

tree.Draw("event_E_balance>>h_bal", "", "goff")

if h_bal.Integral() > 0:
    h_bal.Scale(1.0 / h_bal.Integral())

setup_hist(h_bal, ROOT.kMagenta + 1, 20, 1, 3)
h_bal.SetTitle("")
h_bal.GetXaxis().SetTitle("E(reco jets) / E(partons)")
h_bal.GetYaxis().SetTitle("Normalized entries")
h_bal.Draw("HIST")

mean_bal = h_bal.GetMean()
rms_bal  = h_bal.GetRMS()
nent_bal = int(h_bal.GetEntries())

lat = make_latex()
lat.DrawLatex(0.15, 0.94, "#sqrt{s} = 160 GeV, WW #rightarrow 4q")
lat.DrawLatex(0.15, 0.89, "Event energy balance")
lat.DrawLatex(0.60, 0.84, f"Mean = {mean_bal:.3f}")
lat.DrawLatex(0.60, 0.79, f"RMS = {rms_bal:.3f}")
lat.DrawLatex(0.60, 0.74, f"N = {nent_bal}")

c_bal.RedrawAxis()
save_canvas(c_bal, "event_energy_balance")

# -----------------------------------------------------------------------------
# Text summary table
# -----------------------------------------------------------------------------
print()
print("="*78)
print("SUMMARY TABLE — delta_alpha mean and RMS per strategy per jet")
print("="*78)
print(f"{'Jet':<8} {'Strategy A':>22} {'Strategy B':>22} {'Strategy C':>22}")
print("-"*78)

for jidx in range(1, 5):
    row = f"{jidx:<8}"
    for s in ["A", "B", "C"]:
        branch = f"delta_alpha_{s}_j{jidx}"
        mean, rms, entries = get_hist_mean_rms(
            branch,
            cut=f"{branch} > -998",
            hname=f"h_summary_{s}_{jidx}",
            nbins=200,
            xmin=-2,
            xmax=2
        )
        row += f"{mean:+.3f} ± {rms:.3f}".rjust(22)
    print(row)

print("-"*78)

h_balance_summary = ROOT.TH1F("h_balance_summary", "", 200, 0, 2)
tree.Draw("event_E_balance>>h_balance_summary", "", "goff")
bal_mean = h_balance_summary.GetMean()

print(f"event_E_balance mean = {bal_mean:.3f}")
print(f"-> reco jets recover about {bal_mean*100:.1f}% of total parton energy")
print("="*78)

# -----------------------------------------------------------------------------
# Matching fractions for n=0..5
# FIXED: use n_matched_A/B/C
# -----------------------------------------------------------------------------
print()
print("="*78)
print("MATCHED PARTON MULTIPLICITY FRACTIONS")
print("="*78)

for s in ["A", "B", "C"]:
    htmp = ROOT.TH1F(f"h_frac_{s}", "", 6, -0.5, 5.5)
    tree.Draw(f"n_matched_{s}>>h_frac_{s}", "", "goff")
    total = htmp.Integral()

    print(f"\n{styles[s]['label']}")
    if total > 0:
        for ibin in range(1, htmp.GetNbinsX() + 1):
            nmatch = int(htmp.GetBinCenter(ibin))
            frac = htmp.GetBinContent(ibin) / total
            print(f"  n={nmatch}: {frac:.4f}")
    else:
        print("  no entries")

print("\n[INFO] Done.")
