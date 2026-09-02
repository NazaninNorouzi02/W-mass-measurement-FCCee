"""
diagnostic_alpha.py
====================
Diagnoses the delta_alpha distribution for WW->qqqq at FCC-ee.

Requires: outputs/step1_angular_resolution_multiE/angular_resolution_ecm160.root
Produced by: angular_resolution.py

Two tests:

PART 1 — j1 alpha in events where j4 alpha is bad (< -0.3)
  If j1 stays narrow and centred while j4 is terrible, the problem is
  soft-jet out-of-cone energy loss, not a code bug.

PART 2 — Pearson correlation between delta_alpha and delta_x per jet
  If uncorrelated: alpha and x carry independent info, keep both.

PART 3 — All four alpha distributions overlaid for visual comparison.

Results at sqrt(s) = 160 GeV:
  j1 alpha mean (all):      -0.018
  j1 alpha mean (j4<-0.3):  -0.024  -> unchanged -> out-of-cone physics confirmed
  Pearson corr j1: 0.133, j4: 0.351 -> independent -> keep both
"""

import ROOT
import os

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

INPUT  = "outputs/step1_angular_resolution_multiE/angular_resolution_ecm160.root"
OUTDIR = "outputs/diagnostic_alpha/"
os.makedirs(OUTDIR, exist_ok=True)

f = ROOT.TFile.Open(INPUT)
if not f or f.IsZombie():
    raise RuntimeError(f"Cannot open {INPUT}")

tree = f.Get("events")
if not tree:
    key = f.GetListOfKeys()[0].GetName()
    tree = f.Get(key).Get("events")
if not tree:
    raise RuntimeError("Cannot find 'events' tree")

print(f"Tree entries: {tree.GetEntries()}")
c = ROOT.TCanvas("c", "c", 1200, 900)

# ── PART 1 ─────────────────────────────────────────────────────────────────
print("\n── PART 1: j1 alpha in events where j4 alpha < -0.3 ──")

h_j1_bad_j4 = ROOT.TH1F("h_j1_bad_j4", "j1 alpha when j4<-0.3;(E_{jet}-E_{parton})/E_{parton};Entries (norm)", 80, -0.5, 0.5)
h_j1_all    = ROOT.TH1F("h_j1_all",    "j1 alpha all events;(E_{jet}-E_{parton})/E_{parton};Entries (norm)",   80, -0.5, 0.5)
h_j1_bad_j4.SetLineColor(ROOT.kRed);  h_j1_bad_j4.SetLineWidth(2)
h_j1_all.SetLineColor(ROOT.kBlue);    h_j1_all.SetLineWidth(2)

tree.Draw("delta_alpha_j1>>h_j1_bad_j4", "delta_alpha_j4 < -0.3 && delta_alpha_j4 > -998", "goff")
tree.Draw("delta_alpha_j1>>h_j1_all",    "delta_alpha_j1 > -998", "goff")

c.Clear()
h_j1_all.DrawNormalized()
h_j1_bad_j4.DrawNormalized("same")

leg1 = ROOT.TLegend(0.50, 0.72, 0.92, 0.88)
leg1.SetBorderSize(0)
leg1.AddEntry(h_j1_all,    f"j1 all events (N={int(h_j1_all.GetEntries())})", "l")
leg1.AddEntry(h_j1_bad_j4, f"j1 when j4<-0.3 (N={int(h_j1_bad_j4.GetEntries())})", "l")
leg1.Draw()

latex = ROOT.TLatex(); latex.SetNDC(); latex.SetTextSize(0.032)
latex.DrawLatex(0.15, 0.88, f"j1 mean (all):        {h_j1_all.GetMean():.4f}")
latex.DrawLatex(0.15, 0.83, f"j1 mean (j4 < -0.3):  {h_j1_bad_j4.GetMean():.4f}")
c.SaveAs(f"{OUTDIR}j1_alpha_when_j4_bad.png")

m_all = h_j1_all.GetMean()
m_bad = h_j1_bad_j4.GetMean()
print(f"  j1 mean (all):        {m_all:.4f}")
print(f"  j1 mean (j4 < -0.3):  {m_bad:.4f}")
conclusion1 = "OUT-OF-CONE PHYSICS — not a code bug" if abs(m_all - m_bad) < 0.02 else "INVESTIGATE FURTHER"
print(f"  RESULT: {conclusion1}")

# ── PART 2 ─────────────────────────────────────────────────────────────────
print("\n── PART 2: alpha vs x correlation per jet ──")
ROOT.gStyle.SetPalette(ROOT.kBird)
correlations = {}

for jidx in [1, 2, 3, 4]:
    c.Clear()
    hname = f"h2_j{jidx}"
    h2 = ROOT.TH2F(hname,
        f"Jet {jidx}: #delta#alpha vs #deltax;#deltax;#delta#alpha",
        60, -0.6, 0.6, 60, -0.6, 0.6)
    tree.Draw(f"delta_alpha_j{jidx}:delta_x_j{jidx}>>{hname}",
              f"delta_alpha_j{jidx} > -998 && delta_x_j{jidx} > -998", "goff")
    h2.Draw("COLZ")
    corr = h2.GetCorrelationFactor()
    correlations[jidx] = corr
    interp = "independent" if abs(corr) < 0.5 else "correlated"
    lat2 = ROOT.TLatex(); lat2.SetNDC(); lat2.SetTextSize(0.038)
    lat2.DrawLatex(0.15, 0.88, f"Pearson corr = {corr:.3f}  ({interp})")
    c.SaveAs(f"{OUTDIR}alpha_vs_x_j{jidx}.png")
    print(f"  j{jidx}: corr = {corr:.3f}  [{interp}]")

# ── PART 3 ─────────────────────────────────────────────────────────────────
print("\n── PART 3: all four alpha distributions overlaid ──")
c.Clear()
colors_jet = [ROOT.kBlue, ROOT.kRed, ROOT.kGreen+2, ROOT.kOrange+1]
leg3 = ROOT.TLegend(0.12, 0.62, 0.55, 0.88)
leg3.SetBorderSize(0)
first = True
jet_stats = {}

for jidx in [1, 2, 3, 4]:
    hname = f"h_ov_j{jidx}"
    h = ROOT.TH1F(hname, "#delta#alpha all jets;(E_{jet}-E_{parton})/E_{parton};Entries (norm)", 80, -0.7, 0.4)
    h.SetLineColor(colors_jet[jidx-1]); h.SetLineWidth(2)
    tree.Draw(f"delta_alpha_j{jidx}>>{hname}", f"delta_alpha_j{jidx} > -998", "goff")
    h.DrawNormalized("" if first else "same")
    leg3.AddEntry(h, f"j{jidx}  #mu={h.GetMean():.3f}  RMS={h.GetRMS():.3f}", "l")
    jet_stats[jidx] = (h.GetMean(), h.GetRMS(), int(h.GetEntries()))
    first = False
    print(f"  j{jidx}: mean={h.GetMean():.4f}  RMS={h.GetRMS():.4f}  N={int(h.GetEntries())}")

leg3.Draw()
c.SaveAs(f"{OUTDIR}all_jets_alpha_overlay.png")

# ── SUMMARY ────────────────────────────────────────────────────────────────
print(f"""
{'='*60}
DIAGNOSTIC SUMMARY — delta_alpha at sqrt(s)=160 GeV
{'='*60}

PART 1: Is the j4 plateau a bug or physics?
  j1 mean (all):       {m_all:.4f}
  j1 mean (j4<-0.3):  {m_bad:.4f}
  -> {conclusion1}

PART 2: Are alpha and x independent?
  j1: {correlations.get(1,0):.3f}   j2: {correlations.get(2,0):.3f}
  j3: {correlations.get(3,0):.3f}   j4: {correlations.get(4,0):.3f}
  -> {'ALL INDEPENDENT — keep both in covariance matrix' if all(abs(v)<0.5 for v in correlations.values()) else 'SOME CORRELATED — review'}

PART 3: Distribution means
  j1: {jet_stats[1][0]:.4f}   j2: {jet_stats[2][0]:.4f}
  j3: {jet_stats[3][0]:.4f}   j4: {jet_stats[4][0]:.4f}
  Note: systematic leftward shift from j1->j4.
  Fixed absolute out-of-cone loss -> larger fractional bias for soft jets.

{'='*60}
""")

f.Close()
print("Plots saved to", OUTDIR)