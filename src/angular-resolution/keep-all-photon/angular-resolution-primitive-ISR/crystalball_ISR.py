#!/usr/bin/env python3

import os
import json
import time
import ROOT

ROOT.gROOT.SetBatch(True)


class CrystalBallFit:
    @staticmethod
    def evaluate(x, par):
        """
        Double-sided Crystal Ball function

        par[0] = amplitude
        par[1] = mean
        par[2] = sigma
        par[3] = left transition
        par[4] = left tail power
        par[5] = right transition
        par[6] = right tail power
        """

        t = x[0]

        amplitude = par[0]
        mean = par[1]
        sigma = par[2]

        a1 = abs(par[3])
        n1 = abs(par[4])
        a2 = abs(par[5])
        n2 = abs(par[6])

        if sigma <= 0:
            return 0.0

        if a1 <= 0 or n1 <= 0 or a2 <= 0 or n2 <= 0:
            return 0.0

        u = (t - mean) / sigma

        if -a1 <= u <= a2:
            result = ROOT.TMath.Exp(-0.5 * u * u)

        elif u < -a1:
            A1 = ROOT.TMath.Exp(-0.5 * a1 * a1)
            B1 = (n1 / a1) - a1 - u

            if B1 <= 0:
                return 0.0

            result = A1 * ROOT.TMath.Power(
                (a1 / n1) * B1,
                -n1
            )

        else:
            A2 = ROOT.TMath.Exp(-0.5 * a2 * a2)
            B2 = (n2 / a2) - a2 + u

            if B2 <= 0:
                return 0.0

            result = A2 * ROOT.TMath.Power(
                (a2 / n2) * B2,
                -n2
            )

        return amplitude * result


if __name__ == "__main__":

    ecms = [
        "160",
        "240",
        "340",
        "345",
        "350",
        "355",
        "365"
    ]

    input_directory = (
        "./outputs/"
        "step1_angular_resolution_ISR_multiE"
    )

    output_directory = (
        "./outputs/"
        "resolutions_ISR"
    )

    current_time = time.localtime()

    run_name = (
        f"{current_time.tm_year}-"
        f"{current_time.tm_mon:02d}-"
        f"{current_time.tm_mday:02d}-"
        f"{current_time.tm_hour:02d}-"
        f"{current_time.tm_min:02d}"
    )

    run_output_directory = os.path.join(
        output_directory,
        run_name
    )

    os.makedirs(
        run_output_directory,
        exist_ok=True
    )

    all_energy_results = {}

    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptFit(0)

    ROOT.gStyle.SetPadLeftMargin(0.14)
    ROOT.gStyle.SetPadBottomMargin(0.12)
    ROOT.gStyle.SetPadTopMargin(0.05)
    ROOT.gStyle.SetPadRightMargin(0.05)

    ROOT.gStyle.SetTitleSize(0.05, "XY")
    ROOT.gStyle.SetLabelSize(0.045, "XY")

    ROOT.gStyle.SetLineWidth(1)
    ROOT.gStyle.SetFrameLineWidth(1)

    for ecm in ecms:

        current_time_ecm = time.time()

        root_file = os.path.join(
            input_directory,
            f"angular_resolution_ISR_ecm{ecm}.root"
        )

        energy_output_directory = os.path.join(
            run_output_directory,
            ecm
        )

        os.makedirs(
            energy_output_directory,
            exist_ok=True
        )

        print("\n" + "=" * 70)
        print(f"Processing ECM = {ecm} GeV")
        print(f"Input file: {root_file}")
        print("=" * 70)

        if not os.path.exists(root_file):
            print(
                f"[WARNING] ROOT file does not exist: "
                f"{root_file}"
            )
            continue

        root_input = ROOT.TFile.Open(
            root_file,
            "READ"
        )

        if (
            root_input is None or
            root_input.IsZombie()
        ):
            print(
                f"[WARNING] Cannot open ROOT file: "
                f"{root_file}"
            )
            continue

        tree = root_input.Get("events")

        if not tree:
            print(
                f"[WARNING] Cannot find events tree in: "
                f"{root_file}"
            )
            root_input.Close()
            continue

        print(
            f"Found events tree with "
            f"{tree.GetEntries()} entries"
        )

        plots_config = []

        for jetnumber in range(1, 5):

            plots_config.append({
                "branch": f"delta_alpha_j{jetnumber}",
                "label": "#sigma_{#alpha}",
                "title": (
                    "(E_{jet} - E_{parton})"
                    "/E_{parton}"
                ),
                "bins": 100,
                "range": (-0.25, 0.25)
            })

            plots_config.append({
                "branch": f"delta_x_j{jetnumber}",
                "label": "#sigma_{log(#beta)}",
                "title": (
                    "log(#beta)_{gen} - "
                    "log(#beta)_{reco}"
                ),
                "bins": 100,
                "range": (-0.6, 0.6)
            })

            plots_config.append({
                "branch": f"delta_theta_j{jetnumber}",
                "label": "#sigma_{#theta}",
                "title": (
                    "#theta_{gen} - "
                    "#theta_{reco}"
                ),
                "bins": 100,
                "range": (-0.02, 0.02)
            })

            plots_config.append({
                "branch": f"delta_phi_j{jetnumber}",
                "label": "#sigma_{#phi}",
                "title": (
                    "#phi_{gen} - "
                    "#phi_{reco}"
                ),
                "bins": 100,
                "range": (-0.02, 0.02)
            })

            plots_config.append({
                "branch": (
                    f"filtered_delta_alpha_j"
                    f"{jetnumber}"
                ),
                "label": "#sigma_{#alpha}",
                "title": (
                    "(E_{jet} - E_{parton})"
                    "/E_{parton}"
                ),
                "bins": 100,
                "range": (-0.25, 0.25)
            })

            plots_config.append({
                "branch": (
                    f"filtered_delta_x_j"
                    f"{jetnumber}"
                ),
                "label": "#sigma_{log(#beta)}",
                "title": (
                    "log(#beta)_{gen} - "
                    "log(#beta)_{reco}"
                ),
                "bins": 100,
                "range": (-0.6, 0.6)
            })

            plots_config.append({
                "branch": (
                    f"filtered_delta_theta_j"
                    f"{jetnumber}"
                ),
                "label": "#sigma_{#theta}",
                "title": (
                    "#theta_{gen} - "
                    "#theta_{reco}"
                ),
                "bins": 100,
                "range": (-0.02, 0.02)
            })

            plots_config.append({
                "branch": (
                    f"filtered_delta_phi_j"
                    f"{jetnumber}"
                ),
                "label": "#sigma_{#phi}",
                "title": (
                    "#phi_{gen} - "
                    "#phi_{reco}"
                ),
                "bins": 100,
                "range": (-0.02, 0.02)
            })

        energy_results = {}

        for plot_number, config in enumerate(
            plots_config,
            start=1
        ):

            branch = config["branch"]
            xmin = config["range"][0]
            xmax = config["range"][1]

            print(
                f"\nProcessing {branch}..."
            )

            if not tree.GetBranch(branch):
                print(
                    f"[WARNING] Branch not found: "
                    f"{branch}"
                )
                continue

            canvas_name = (
                f"canvas_{ecm}_{plot_number}"
            )

            canvas = ROOT.TCanvas(
                canvas_name,
                "Resolution Fits",
                1920,
                1080
            )

            histogram_name = (
                f"histogram_{ecm}_{branch}"
            )

            draw_expression = (
                f"{branch}>>"
                f"{histogram_name}("
                f"{config['bins']},"
                f"{xmin},"
                f"{xmax})"
            )

            selection = (
                f"TMath::Finite({branch}) && "  # changed from  f"({branch} == {branch}) 
                f"({branch} >= {xmin}) && "
                f"({branch} <= {xmax})"
            )

            tree.Draw(
                draw_expression,
                selection,
                "goff"
            )

            histogram = ROOT.gDirectory.Get(
                histogram_name
            )

            if (
                histogram is None or
                histogram.GetEntries() == 0
            ):
                print(
                    f"[WARNING] No data for {branch}"
                )
                canvas.Close()
                continue

            entries = histogram.GetEntries()
            histogram_mean = histogram.GetMean()
            histogram_rms = histogram.GetRMS()

            print(f"  Entries: {entries}")
            print(
                f"  Mean: {histogram_mean:.6f}"
            )
            print(
                f"  RMS: {histogram_rms:.6f}"
            )

            bin_min = histogram.FindBin(xmin)
            bin_max = histogram.FindBin(xmax)

            maximum_bin = histogram.GetMaximum()

            sum_weights = 0.0
            sum_wx = 0.0
            sum_wx2 = 0.0

            for bin_index in range(
                bin_min,
                bin_max + 1
            ):
                bin_center = histogram.GetBinCenter(
                    bin_index
                )

                bin_content = histogram.GetBinContent(
                    bin_index
                )

                sum_weights += bin_content
                sum_wx += bin_content * bin_center
                sum_wx2 += (
                    bin_content *
                    bin_center *
                    bin_center
                )

            mean_guess = histogram_mean
            sigma_guess = histogram_rms

            if sum_weights > 0:
                mean_guess = (
                    sum_wx / sum_weights
                )

                variance = (
                    sum_wx2 / sum_weights
                ) - (
                    mean_guess * mean_guess
                )

                if variance > 0:
                    sigma_guess = (
                        ROOT.TMath.Sqrt(variance)
                    )

            if sigma_guess <= 0:
                print(
                    f"[WARNING] Invalid sigma guess "
                    f"for {branch}"
                )
                canvas.Close()
                continue

            print(
                f"  Fit range mean: "
                f"{mean_guess:.6f}"
            )

            print(
                f"  Fit range RMS: "
                f"{sigma_guess:.6f}"
            )

            fit_object = CrystalBallFit()

            fit_function_name = (
                f"dscb_{ecm}_{plot_number}"
            )

            fit_function = ROOT.TF1(
                fit_function_name,
                fit_object.evaluate,
                xmin,
                xmax,
                7
            )

            fit_function.SetParName(
                0,
                "Amplitude"
            )

            fit_function.SetParName(
                1,
                "#mu"
            )

            fit_function.SetParName(
                2,
                "#sigma"
            )

            fit_function.SetParName(
                3,
                "a_{1}"
            )

            fit_function.SetParName(
                4,
                "n_{1}"
            )

            fit_function.SetParName(
                5,
                "a_{2}"
            )

            fit_function.SetParName(
                6,
                "n_{2}"
            )

            fit_function.SetParameter(
                0,
                maximum_bin
            )

            fit_function.SetParameter(
                1,
                mean_guess
            )

            fit_function.SetParameter(
                2,
                sigma_guess
            )

            fit_function.SetParameter(
                3,
                1.0
            )

            fit_function.SetParameter(
                4,
                2.0
            )

            fit_function.SetParameter(
                5,
                1.0
            )

            fit_function.SetParameter(
                6,
                2.0
            )

            fit_function.SetParLimits(
                0,
                0.0,
                2.0 * maximum_bin
            )

            fit_function.SetParLimits(
                1,
                mean_guess - 3.0 * sigma_guess,
                mean_guess + 3.0 * sigma_guess
            )

            fit_function.SetParLimits(
                2,
                0.01 * sigma_guess,
                5.0 * sigma_guess
            )

            fit_function.SetParLimits(
                3,
                0.1,
                3.0
            )

            fit_function.SetParLimits(
                4,
                0.1,
                10.0
            )

            fit_function.SetParLimits(
                5,
                0.1,
                3.0
            )

            fit_function.SetParLimits(
                6,
                0.1,
                10.0
            )

            print("\n" + "=" * 60)
            print(
                "Performing Double-Sided "
                "Crystal Ball Fit..."
            )
            print("=" * 60)

            ROOT.gErrorIgnoreLevel = ROOT.kInfo

            fit_result = histogram.Fit(
                fit_function,
                "SMRL"
            )

            fit_status = int(fit_result)

            chi2 = fit_function.GetChisquare()
            ndf = fit_function.GetNDF()

            if ndf > 0:
                chi2_ndf = chi2 / ndf
            else:
                chi2_ndf = 0.0

            fit_parameters = {
                "ecm": int(ecm),
                "branch": branch,
                "entries": int(entries),
                "histogram_mean": float(
                    histogram_mean
                ),
                "histogram_rms": float(
                    histogram_rms
                ),
                "fit_status": fit_status,

                "amplitude": float(
                    fit_function.GetParameter(0)
                ),
                "amplitude_err": float(
                    fit_function.GetParError(0)
                ),

                "mu": float(
                    fit_function.GetParameter(1)
                ),
                "mu_err": float(
                    fit_function.GetParError(1)
                ),

                "sigma": float(
                    fit_function.GetParameter(2)
                ),
                "sigma_err": float(
                    fit_function.GetParError(2)
                ),

                "a1": float(
                    fit_function.GetParameter(3)
                ),
                "a1_err": float(
                    fit_function.GetParError(3)
                ),

                "n1": float(
                    fit_function.GetParameter(4)
                ),
                "n1_err": float(
                    fit_function.GetParError(4)
                ),

                "a2": float(
                    fit_function.GetParameter(5)
                ),
                "a2_err": float(
                    fit_function.GetParError(5)
                ),

                "n2": float(
                    fit_function.GetParameter(6)
                ),
                "n2_err": float(
                    fit_function.GetParError(6)
                ),

                "chi2": float(chi2),
                "ndf": int(ndf),
                "chi2_ndf": float(chi2_ndf)
            }

            energy_results[branch] = fit_parameters

            print("\n" + "=" * 60)
            print("FIT RESULTS")
            print("=" * 60)

            print(
                f"Status:     "
                f"{fit_parameters['fit_status']}"
            )

            print(
                f"Amplitude:  "
                f"{fit_parameters['amplitude']:.2f} "
                f"± "
                f"{fit_parameters['amplitude_err']:.2f}"
            )

            print(
                f"μ:          "
                f"{fit_parameters['mu']:.6f} "
                f"± "
                f"{fit_parameters['mu_err']:.6f}"
            )

            print(
                f"σ:          "
                f"{fit_parameters['sigma']:.6f} "
                f"± "
                f"{fit_parameters['sigma_err']:.6f}"
            )

            print(
                f"a1:         "
                f"{fit_parameters['a1']:.4f} "
                f"± "
                f"{fit_parameters['a1_err']:.4f}"
            )

            print(
                f"n1:         "
                f"{fit_parameters['n1']:.4f} "
                f"± "
                f"{fit_parameters['n1_err']:.4f}"
            )

            print(
                f"a2:         "
                f"{fit_parameters['a2']:.4f} "
                f"± "
                f"{fit_parameters['a2_err']:.4f}"
            )

            print(
                f"n2:         "
                f"{fit_parameters['n2']:.4f} "
                f"± "
                f"{fit_parameters['n2_err']:.4f}"
            )

            print(
                f"χ²/ndf:     "
                f"{fit_parameters['chi2']:.2f}/"
                f"{fit_parameters['ndf']} = "
                f"{fit_parameters['chi2_ndf']:.3f}"
            )

            print("=" * 60 + "\n")

            histogram.SetLineColor(
                ROOT.kBlue
            )

            histogram.SetLineWidth(1)
            histogram.SetFillColor(0)
            histogram.SetTitle("")

            histogram.GetXaxis().SetTitle(
                config["title"]
            )

            histogram.GetYaxis().SetTitle(
                "N_{entries}"
            )

            histogram.GetXaxis().CenterTitle()
            histogram.GetYaxis().CenterTitle()

            histogram.GetXaxis().SetTitleSize(
                0.05
            )

            histogram.GetYaxis().SetTitleSize(
                0.05
            )

            histogram.GetXaxis().SetLabelSize(
                0.045
            )

            histogram.GetYaxis().SetLabelSize(
                0.045
            )

            fit_function.SetLineColor(
                ROOT.kBlack
            )

            fit_function.SetLineWidth(2)

            canvas.cd()

            pad1_name = (
                f"pad1_{ecm}_{plot_number}"
            )

            pad1 = ROOT.TPad(
                pad1_name,
                pad1_name,
                0.0,
                0.0,
                0.7,
                1.0
            )

            pad1.SetLeftMargin(0.14)
            pad1.SetBottomMargin(0.12)
            pad1.SetTopMargin(0.08)
            pad1.SetRightMargin(0.02)

            pad1.Draw()
            pad1.cd()

            histogram.Draw("HIST")
            fit_function.Draw("SAME")

            legend = ROOT.TLegend(
                0.18,
                0.70,
                0.48,
                0.88
            )

            legend.SetBorderSize(0)
            legend.SetFillStyle(0)

            legend.AddEntry(
                histogram,
                f"Data, #sqrt{{s}} = {ecm} GeV",
                "l"
            )

            legend.AddEntry(
                fit_function,
                "Double-Sided Crystal Ball",
                "l"
            )

            legend.Draw()

            sigma_label = ROOT.TLatex()

            sigma_label.SetNDC()
            sigma_label.SetTextSize(0.045)
            sigma_label.SetTextFont(42)

            sigma_label.DrawLatex(
                0.18,
                0.62,
                (
                    f"{config['label']} = "
                    f"{fit_parameters['sigma']:.5f} "
                    f"#pm "
                    f"{fit_parameters['sigma_err']:.5f}"
                )
            )

            canvas.cd()

            pad2_name = (
                f"pad2_{ecm}_{plot_number}"
            )

            pad2 = ROOT.TPad(
                pad2_name,
                pad2_name,
                0.7,
                0.0,
                1.0,
                1.0
            )

            pad2.SetLeftMargin(0.05)
            pad2.SetRightMargin(0.05)
            pad2.SetTopMargin(0.08)
            pad2.SetBottomMargin(0.12)

            pad2.Draw()
            pad2.cd()

            text = ROOT.TPaveText(
                0.05,
                0.15,
                0.95,
                0.92,
                "NDC"
            )

            text.SetBorderSize(1)
            text.SetFillColor(0)
            text.SetTextAlign(12)
            text.SetTextSize(0.042)

            text.AddText(
                f"#sqrt{{s}} = {ecm} GeV"
            )

            text.AddText(
                f"Entries = {int(entries)}"
            )

            text.AddText("")

            text.AddText(
                f"#mu = "
                f"{fit_parameters['mu']:.5f}"
            )

            text.AddText(
                f"     #pm "
                f"{fit_parameters['mu_err']:.5f}"
            )

            text.AddText("")

            text.AddText(
                f"#sigma = "
                f"{fit_parameters['sigma']:.5f}"
            )

            text.AddText(
                f"        #pm "
                f"{fit_parameters['sigma_err']:.5f}"
            )

            text.AddText("")

            text.AddText(
                f"a_{{1}} = "
                f"{fit_parameters['a1']:.4f}"
            )

            text.AddText(
                f"n_{{1}} = "
                f"{fit_parameters['n1']:.4f}"
            )

            text.AddText("")

            text.AddText(
                f"a_{{2}} = "
                f"{fit_parameters['a2']:.4f}"
            )

            text.AddText(
                f"n_{{2}} = "
                f"{fit_parameters['n2']:.4f}"
            )

            text.AddText("")

            text.AddText(
                f"#chi^{{2}}/ndf = "
                f"{fit_parameters['chi2_ndf']:.3f}"
            )

            text.AddText(
                f"Fit status = "
                f"{fit_parameters['fit_status']}"
            )

            text.Draw()

            canvas.cd()
            canvas.Update()

            png_path = os.path.join(
                energy_output_directory,
                f"{branch}.png"
            )

            canvas.SaveAs(png_path)

            print(
                f"Saved: {png_path}"
            )

            canvas.Close()

        energy_json_path = os.path.join(
            energy_output_directory,
            "all_results.json"
        )

        with open(
            energy_json_path,
            "w",
            encoding="utf-8"
        ) as json_file:
            json.dump(
                energy_results,
                json_file,
                indent=2
            )

        all_energy_results[ecm] = energy_results

        print("\n" + "=" * 60)
        print(f"RESULTS FOR ECM {ecm} GeV")
        print("=" * 60)

        for branch, result in energy_results.items():
            print(
                f"{branch}: "
                f"σ = {result['sigma']:.5f} "
                f"± {result['sigma_err']:.5f}"
            )

        print("=" * 60)

        print(
            f"JSON saved: {energy_json_path}"
        )

        print(
            f"Execution time for {ecm} GeV: "
            f"{time.time() - current_time_ecm:.2f} seconds"
        )

        root_input.Close()

    combined_json_path = os.path.join(
        run_output_directory,
        "all_energies_results.json"
    )

    with open(
        combined_json_path,
        "w",
        encoding="utf-8"
    ) as json_file:
        json.dump(
            all_energy_results,
            json_file,
            indent=2
        )

    print("\n" + "=" * 70)
    print("ALL ENERGY LEVELS COMPLETED")
    print("=" * 70)

    print(
        f"Combined JSON saved: "
        f"{combined_json_path}"
    )

    print(
        f"Output directory: "
        f"{run_output_directory}"
    )