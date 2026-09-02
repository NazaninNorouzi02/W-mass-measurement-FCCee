#!/usr/bin/env python3
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(1110)

input_file = "wmass_fit_pvalue_pulls_ISR_ecm365.root"
tree_name = "events"
branch = "delta_m_4C_pre5C"

df = ROOT.RDataFrame(tree_name, input_file)
df = df.Filter(f"std::isfinite({branch})")

h = df.Histo1D(
    (
        "h_delta_m_4C_pre5C",
        ";#Delta M_{4C}^{pre-5C} = M_{large,4C} - M_{small,4C} [GeV];Events",
        100,
        0.0,
        150.0,
    ),
    branch,
)

c = ROOT.TCanvas("c", "c", 900, 700)
h.SetLineColor(ROOT.kBlue + 1)
h.SetLineWidth(2)
h.Draw("HIST")

latex = ROOT.TLatex()
latex.SetNDC()
latex.SetTextSize(0.035)
latex.DrawLatex(0.62, 0.86, f"Entries = {int(h.GetEntries())}")
latex.DrawLatex(0.62, 0.81, f"Mean = {h.GetMean():.3f} GeV")
latex.DrawLatex(0.62, 0.76, f"RMS = {h.GetRMS():.3f} GeV")

c.SaveAs("delta_m_4C_pre5C.png")
c.SaveAs("delta_m_4C_pre5C.pdf")