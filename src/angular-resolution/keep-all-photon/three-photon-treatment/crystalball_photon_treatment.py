#!/usr/bin/env python3
import os
import ROOT
import time

ROOT.gROOT.SetBatch(True)


class CrystalBallFit:
    @staticmethod
    def evaluate(x, par):
        """
          Evaluate the double-sided Crystal Ball function
          par[0] = normalization (amplitude)
          par[1] = mu (mean)
          par[2] = sigma (width)
          par[3] = a1 (left tail parameter)
          par[4] = n1 (left tail power)
          par[5] = a2 (right tail parameter)
          par[6] = n2 (right tail power)
        """
        t = x[0]
        N = par[0]
        mu = par[1]
        sigma = par[2]
        a1 = abs(par[3])  # Ensure positive
        n1 = abs(par[4])  # Ensure positive
        a2 = abs(par[5])  # Ensure positive
        n2 = abs(par[6])  # Ensure positive

        # If these are zero we would have errors in computation!
        assert a1 != 0 and a2 != 0 and n1 != 0 and n2 != 0

        if sigma == 0:
            return 0
        # Normalized variable
        u = (t - mu) / sigma

        # Gaussian core: -a1 <= u <= a2
        if u >= -a1 and u <= a2:
            result = ROOT.TMath.Exp(-0.5 * u * u)

        # Left tail: u < -a1
        elif u < -a1:
            A1 = ROOT.TMath.Exp(-0.5 * a1 * a1)
            B1 = (n1 / a1) - a1 - u
            if B1 <= 0:
                return 0
            result = A1 * ROOT.TMath.Power((a1 / n1) * B1, -n1)

        # Right tail: u > a2
        else:  # u > a2
            A2 = ROOT.TMath.Exp(-0.5 * a2 * a2)
            B2 = (n2 / a2) - a2 + u
            if B2 <= 0:
                return 0
            result = A2 * ROOT.TMath.Power((a2 / n2) * B2, -n2)

        return N * result


if __name__ == "__main__":

    ### INITIALIZATION:
    # Input
    ecms = ["160"]
    for ecm in ecms:
        root_file = f"./outputs/angular_resolution_three_photon_treatments/angular_resolution_three_photon_treatments_ecm{ecm}.root"
        outDir = f"./outputs/resolutions_three_photon_treatments/{time.localtime().tm_year}-{time.localtime().tm_mon}-{time.localtime().tm_mday}-{time.localtime().tm_hour}-{time.localtime().tm_min}/{ecm}"
        os.makedirs(outDir, exist_ok=True)
        currenttime = time.time()

        plots_config = []

        # Branch format produced by the three-photon-treatment analysis:
        #   <variable>_<treatment>_j<jet number>
        #
        # Treatments:
        #   noPhotons     = all reconstructed photons removed
        #   allPhotons    = all reconstructed photons retained
        #   noISRPhotons  = only ISR-like reconstructed photons removed
        treatments = [
            ("noPhotons", "All photons removed"),
            ("allPhotons", "All photons retained"),
            ("noISRPhotons", "ISR-like photons removed"),
        ]

        variables = {
            "delta_alpha": {
                "label": "#sigma_{#alpha}",
                "title": "(E_{reco jet} - E_{parton})/E_{parton}",
                "bins": 100,
                "range": (-0.25, 0.25),
            },
            "delta_x": {
                "label": "#sigma_{x}",
                "title": "log(p/E)_{reco} - log(p/E)_{gen}",
                "bins": 100,
                "range": (-0.6, 0.6),
            },
            "delta_theta": {
                "label": "#sigma_{#theta}",
                "title": "#theta_{reco} - #theta_{gen}",
                "bins": 100,
                "range": (-0.02, 0.02),
            },
            "delta_phi": {
                "label": "#sigma_{#phi}",
                "title": "#phi_{reco} - #phi_{gen}",
                "bins": 100,
                "range": (-0.02, 0.02),
            },
        }

        for treatment, treatment_label in treatments:
            for jetnumber in range(1, 5):
                for variable_name, variable_config in variables.items():
                    plots_config.append({
                        "branch": f"{variable_name}_{treatment}_j{jetnumber}",
                        "label": variable_config["label"],
                        "title": variable_config["title"],
                        "bins": variable_config["bins"],
                        "range": variable_config["range"],
                        "treatment": treatment,
                        "treatment_label": treatment_label,
                        "jetnumber": jetnumber,
                    })

        ### CANVAS INITIALIZATIONS:
        # Set ROOT style to match reference plots
        ROOT.gStyle.SetOptStat(0)  # Remove stat box
        ROOT.gStyle.SetOptFit(0)  # Remove fit box
        ROOT.gStyle.SetPadLeftMargin(0.14)
        ROOT.gStyle.SetPadBottomMargin(0.12)
        ROOT.gStyle.SetPadTopMargin(0.05)
        ROOT.gStyle.SetPadRightMargin(0.05)
        ROOT.gStyle.SetTitleSize(0.05, "XY")
        ROOT.gStyle.SetLabelSize(0.045, "XY")
        ROOT.gStyle.SetLineWidth(1)
        ROOT.gStyle.SetFrameLineWidth(1)

        ### OPENNING THE FILE:
        if not os.path.exists(root_file):
            raise SystemExit(f"ERROR: root file '{root_file}' not found.")

        f = ROOT.TFile.Open(root_file, "READ")
        if f is None or f.IsZombie():
            raise SystemExit(f"ERROR: cannot open {root_file}")

        tree = f.Get("events")
        if not tree:
            raise SystemExit(f"ERROR: cannot find 'events' tree in {root_file}")

        print(f"Found tree with {tree.GetEntries()} entries")

        results = {}
        for i, config in enumerate(plots_config, 1):
            # Create separate canvas for each plot
            canvas = ROOT.TCanvas(f"c{i}", "Resolution Fits", 1920, 1080)
            print(f"\nProcessing {config['branch']}...")
            histogramName = f"h_{config['branch']}"
            tree.Draw(
                f"{config['branch']}>>{histogramName}({config['bins']}, {config['range'][0]}, {config['range'][1]})",
                f"({config['branch']} == {config['branch']})",
                "goff")
            histogram = ROOT.gDirectory.Get(histogramName)

            if not histogram or histogram.GetEntries() == 0:
                print(f"[WARN] No data for {config['branch']}")
                continue
            print(f"  Entries: {histogram.GetEntries()}")
            print(f"  Mean: {histogram.GetMean():.6f}, RMS: {histogram.GetRMS():.6f}")

            xmin = config['range'][0]
            xmax = config['range'][1]

            ## Creating a Fitting Function
            dscb = CrystalBallFit()
            fitFunction = ROOT.TF1("dscb", dscb.evaluate, xmin, xmax, 7)

            ## Setting Parameter Names:
            fitFunction.SetParName(0, "Amplitude")
            fitFunction.SetParName(1, "#mu")
            fitFunction.SetParName(2, "#sigma")
            fitFunction.SetParName(3, "a_{1}")
            fitFunction.SetParName(4, "n_{1}")
            fitFunction.SetParName(5, "a_{2}")
            fitFunction.SetParName(6, "n_{2}")

            ## Setting initial parameter values
            # CRITICAL: Calculate mean and RMS only within the fit range!
            # Find the bin corresponding to xmin and xmax
            bin_min = histogram.FindBin(xmin)
            bin_max = histogram.FindBin(xmax)

            # Get statistics within the fit range
            maximumBin = histogram.GetMaximum()
            meanGuess = histogram.GetMean()  # This will be recalculated below
            sigmaGuess = histogram.GetRMS()  # This will be recalculated below

            # Calculate mean and RMS only within fit range
            sum_weights = 0.0
            sum_wx = 0.0
            sum_wx2 = 0.0

            for i in range(bin_min, bin_max + 1):
                bin_center = histogram.GetBinCenter(i)
                bin_content = histogram.GetBinContent(i)

                sum_weights += bin_content
                sum_wx += bin_content * bin_center
                sum_wx2 += bin_content * bin_center * bin_center

            if sum_weights > 0:
                meanGuess = sum_wx / sum_weights
                variance = (sum_wx2 / sum_weights) - (meanGuess * meanGuess)
                sigmaGuess = ROOT.TMath.Sqrt(variance) if variance > 0 else sigmaGuess

            print(f"  Fit range mean: {meanGuess:.6f}, RMS: {sigmaGuess:.6f}")

            fitFunction.SetParameter(0, maximumBin)
            fitFunction.SetParameter(1, meanGuess)
            fitFunction.SetParameter(2, sigmaGuess)
            fitFunction.SetParameter(3, 1.0)
            fitFunction.SetParameter(4, 2.0)
            fitFunction.SetParameter(5, 1.0)
            fitFunction.SetParameter(6, 2.0)

            ## Setting parameter limits:
            fitFunction.SetParLimits(0, 0, 2 * maximumBin)
            fitFunction.SetParLimits(1, meanGuess - 3 * sigmaGuess, meanGuess + 3 * sigmaGuess)  # Constrain mu
            fitFunction.SetParLimits(2, 0.01 * sigmaGuess, 5 * sigmaGuess)  # Reasonable sigma range
            fitFunction.SetParLimits(3, 0.1, 3.0)  # a1
            fitFunction.SetParLimits(4, 0.1, 10.0)  # n1
            fitFunction.SetParLimits(5, 0.1, 3.0)  # a2
            fitFunction.SetParLimits(6, 0.1, 10.0)  # n2

            ## PERFORMING FIT:
            print("\n" + "=" * 60)
            print("Performing Double-Sided Crystal Ball Fit...")
            print("=" * 60)

            # fitResult = histogram.Fit("dscb", "SRML")
            # Fit
            ROOT.gErrorIgnoreLevel = ROOT.kInfo
            # histogram.Fit(fitFunction, "QSR")
            histogram.Fit(fitFunction, "SMRL")  # Q=quiet, S=save result, M=improve, R=range

            # Extract fit parameters
            fitParameters = {
                'amplitude': fitFunction.GetParameter(0),
                'amplitude_err': fitFunction.GetParError(0),
                'mu': fitFunction.GetParameter(1),
                'mu_err': fitFunction.GetParError(1),
                'sigma': fitFunction.GetParameter(2),
                'sigma_err': fitFunction.GetParError(2),
                'a1': fitFunction.GetParameter(3),
                'a1_err': fitFunction.GetParError(3),
                'n1': fitFunction.GetParameter(4),
                'n1_err': fitFunction.GetParError(4),
                'a2': fitFunction.GetParameter(5),
                'a2_err': fitFunction.GetParError(5),
                'n2': fitFunction.GetParameter(6),
                'n2_err': fitFunction.GetParError(6),
                'chi2': fitFunction.GetChisquare(),
                'ndf': fitFunction.GetNDF(),
                'chi2_ndf': fitFunction.GetChisquare() / fitFunction.GetNDF() if fitFunction.GetNDF() > 0 else 0
            }

            # Store results
            results[config['branch']] = fitParameters

            # Print results
            print("\n" + "=" * 60)
            print("FIT RESULTS")
            print("=" * 60)
            print(f"Amplitude:  {fitParameters['amplitude']:.2f} ± {fitParameters['amplitude_err']:.2f}")
            print(f"μ (mean):   {fitParameters['mu']:.4f} ± {fitParameters['mu_err']:.4f}")
            print(f"σ (width):  {fitParameters['sigma']:.4f} ± {fitParameters['sigma_err']:.4f}")
            print(f"a₁ (left):  {fitParameters['a1']:.4f} ± {fitParameters['a1_err']:.4f}")
            print(f"n₁ (left):  {fitParameters['n1']:.4f} ± {fitParameters['n1_err']:.4f}")
            print(f"a₂ (right): {fitParameters['a2']:.4f} ± {fitParameters['a2_err']:.4f}")
            print(f"n₂ (right): {fitParameters['n2']:.4f} ± {fitParameters['n2_err']:.4f}")
            print(
                f"\nχ²/ndf:     {fitParameters['chi2']:.2f}/{fitParameters['ndf']:.0f} = {fitParameters['chi2_ndf']:.2f}")
            print("=" * 60 + "\n")

            histogram.SetLineColor(ROOT.kBlue)
            histogram.SetLineWidth(1)
            histogram.SetFillColor(0)

            histogram.SetTitle(f"")
            histogram.GetXaxis().SetTitle(config['title'])
            histogram.GetYaxis().SetTitle("N_{entries}")
            histogram.GetXaxis().CenterTitle()
            histogram.GetYaxis().CenterTitle()
            histogram.GetXaxis().SetTitleSize(0.05)
            histogram.GetYaxis().SetTitleSize(0.05)
            histogram.GetXaxis().SetLabelSize(0.045)
            histogram.GetYaxis().SetLabelSize(0.045)

            fitFunction.SetLineColor(ROOT.kBlack)
            fitFunction.SetLineWidth(2)

            # Set up the main pad for histogram (left side, ~70% width)
            canvas.cd()
            pad1 = ROOT.TPad("pad1", "pad1", 0.0, 0.0, 0.7, 1.0)
            pad1.SetLeftMargin(0.14)
            pad1.SetBottomMargin(0.12)
            pad1.SetTopMargin(0.08)
            pad1.SetRightMargin(0.02)
            pad1.Draw()
            pad1.cd()

            histogram.Draw("HIST")

            fitFunction.Draw("same")

            # Add legend on the plot
            legend = ROOT.TLegend(0.18, 0.70, 0.45, 0.88)
            legend.SetBorderSize(0)
            legend.SetFillStyle(0)
            legend.AddEntry(histogram, "Data", "l")
            legend.AddEntry(fitFunction, "Double-Sided Crystal Ball", "l")
            legend.Draw()

            # Add sigma label on the plot
            sigmaLabel = ROOT.TLatex()
            sigmaLabel.SetNDC()
            sigmaLabel.SetTextSize(0.045)
            sigmaLabel.SetTextFont(42)
            sigmaLabel.DrawLatex(0.18, 0.62,
                                 f"{config['label']} = {fitParameters['sigma']:.5f} #pm {fitParameters['sigma_err']:.5f}")

            # Set up the text pad (right side, ~30% width)
            canvas.cd()
            pad2 = ROOT.TPad("pad2", "pad2", 0.7, 0.0, 1.0, 1.0)
            pad2.SetLeftMargin(0.05)
            pad2.SetRightMargin(0.05)
            pad2.SetTopMargin(0.08)
            pad2.SetBottomMargin(0.12)
            pad2.Draw()
            pad2.cd()

            # Add fit parameters text box (beside the graph)
            text = ROOT.TPaveText(0.05, 0.20, 0.95, 0.92, "NDC")
            text.SetBorderSize(1)
            text.SetFillColor(0)
            text.SetTextAlign(12)
            text.SetTextSize(0.045)
            text.AddText("Fit Parameters:")
            text.AddText("")
            text.AddText(f"#mu = {fitParameters['mu']:.5f}")
            text.AddText(f"      #pm {fitParameters['mu_err']:.5f}")
            text.AddText("")
            text.AddText(f"#sigma = {fitParameters['sigma']:.5f}")
            text.AddText(f"         #pm {fitParameters['sigma_err']:.5f}")
            text.AddText("")
            text.AddText(f"a_{{1}} = {fitParameters['a1']:.4f}")
            text.AddText(f"       #pm {fitParameters['a1_err']:.4f}")
            text.AddText("")
            text.AddText(f"n_{{1}} = {fitParameters['n1']:.4f}")
            text.AddText(f"       #pm {fitParameters['n1_err']:.4f}")
            text.AddText("")
            text.AddText(f"a_{{2}} = {fitParameters['a2']:.4f}")
            text.AddText(f"       #pm {fitParameters['a2_err']:.4f}")
            text.AddText("")
            text.AddText(f"n_{{2}} = {fitParameters['n2']:.4f}")
            text.AddText(f"       #pm {fitParameters['n2_err']:.4f}")
            text.AddText("")
            text.AddText(f"#chi^{{2}}/ndf = {fitParameters['chi2_ndf']:.3f}")
            text.Draw()

            # Update and save
            canvas.cd()
            canvas.Update()

            # Save only PNG for this plot
            png_path = os.path.join(outDir, f"{config['branch']}.png")
            canvas.SaveAs(png_path)
            print(f"  Saved: {png_path}")

        print("\n" + "=" * 60)
        print("RESULTS:")
        for branch, res in results.items():
            print(f"{branch}: σ = {res['sigma']:.5f} ± {res['sigma_err']:.5f}")
        print("=" * 60)
        print(f"Execution time: {time.time() - currenttime:.2f} seconds")

        f.Close()

