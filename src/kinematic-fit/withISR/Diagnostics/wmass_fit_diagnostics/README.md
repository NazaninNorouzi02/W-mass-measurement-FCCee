# Final W-Mass Diagnostic Report

This report keeps only the numbers needed to understand the current 5C and ISR problems.

Important conventions:

- The 162.5 GeV sample is treated as 4C-only for physics conclusions. Any 5C branches in the ROOT tree are ignored here.
- Plot means and sigmas are computed from the same visible mass window as the plotting code: 40 <= M < 120 GeV.
- `with ISR treatment` means the final hybrid variable before applying the final 3% rejection, matching Figures 7.9 and 7.10.
- The final accepted rows use `P_final > 0.03`.
- Full RMS values are intentionally omitted because they overemphasize extreme tails and do not match the plotted sigma boxes.

Generated CSV files:

- `plot_style_numbers.csv`: numbers matching the plot definitions.
- `cutflow_summary.csv`: compact acceptance and rejection summary.
- `mass_split_stress_test.csv`: whether the 5C problem follows the pre-5C 4C mass difference.
- `isr_truth_recovery.csv`: ISR recovery and truth-correlation checks.
- `emax_check.csv`: check of the stored ISR photon energy limit.

## Executive Summary

| Energy | standard 5C P>0.03 | final 5C accepted | note |
|---|---|---|---|
| 240 | 8.01% | 3666/45771 | 5C physics sample |
| 240 | 22.34% | 10237/45830 | 5C physics sample |
| 365 | 9.82% | 4438/45204 | 5C physics sample |
| 365 | 24.10% | 10986/45583 | 5C physics sample |

| Energy | Panel | Curve | Mean [GeV] | Plot Sigma [GeV] | Visible entries |
|---|---|---|---|---|---|
| 162p5 | smaller | raw | 66.17 | 9.59 | 43980 |
| 162p5 | larger | raw | 80.05 | 7.67 | 45841 |
| 162p5 | smaller | 4C with ISR treatment | 67.80 | 9.84 | 44203 |
| 162p5 | larger | 4C with ISR treatment | 82.33 | 7.80 | 45840 |
| 240 | smaller | raw | 73.01 | 10.24 | 44385 |
| 240 | larger | raw | 85.30 | 10.03 | 45029 |
| 240 | smaller | 4C with ISR treatment | 73.80 | 10.24 | 44328 |
| 240 | larger | 4C with ISR treatment | 86.03 | 10.01 | 44944 |
| 240 | smaller | 5C with ISR treatment | 78.28 | 9.77 | 44804 |
| 240 | larger | 5C with ISR treatment | 78.28 | 9.77 | 44804 |
| 365 | smaller | raw | 74.30 | 10.52 | 42493 |
| 365 | larger | raw | 84.48 | 9.52 | 40897 |
| 365 | smaller | 4C with ISR treatment | 75.00 | 10.41 | 42390 |
| 365 | larger | 4C with ISR treatment | 85.24 | 9.53 | 40791 |
| 365 | smaller | 5C with ISR treatment | 78.86 | 10.63 | 43356 |
| 365 | larger | 5C with ISR treatment | 78.86 | 10.63 | 43356 |

## 5C Stress Test

This is the main diagnostic for the current problem. If the low 5C probability and large 5C sigma increase with `delta_m_4C_pre5C`, the equal-mass constraint is being forced onto events whose two W candidates are already incompatible after 4C.

| Energy | Delta m bin [GeV] | Events | P5C<0.03 | Pfinal<0.03 | Sigma 5C final [GeV] | ISR recovered/applied |
|---|---|---|---|---|---|---|
| 240 | 0-2 | 8402 | 59.3% | 42.2% | 4.08 | 28.8% |
| 240 | 2-5 | 8896 | 97.4% | 69.2% | 5.16 | 28.9% |
| 240 | 5-10 | 8089 | 99.9% | 81.6% | 6.91 | 18.3% |
| 240 | 10-20 | 7977 | 100.0% | 92.4% | 10.26 | 7.6% |
| 240 | >=20 | 12396 | 100.0% | 95.5% | 14.59 | 4.5% |
| 365 | 0-2 | 9354 | 55.6% | 39.3% | 3.44 | 29.3% |
| 365 | 2-5 | 8742 | 96.8% | 69.8% | 4.45 | 27.9% |
| 365 | 5-10 | 6958 | 99.9% | 83.4% | 6.48 | 16.6% |
| 365 | 10-20 | 5905 | 100.0% | 90.4% | 10.03 | 9.6% |
| 365 | >=20 | 14177 | 100.0% | 93.4% | 17.01 | 6.6% |

## ISR Truth And Recovery

The most useful ISR numbers are the eligible fraction, recovery fraction, signed-pz correlation, and pz residual RMS. The signed-pz correlation tests whether the fitted ISR photon follows the true ISR direction event by event.

| Energy | Fit | True ISR | Eligible | Recovered/applied | rho signed pz | pz RMS [GeV] |
|---|---|---|---|---|---|---|
| 162p5 | 4C+ISR | all | 29.5% | 52.3% | 0.401 | 6.73 |
| 162p5 | 4C+ISR | <5 | 28.7% | 50.5% | 0.140 | 6.80 |
| 162p5 | 4C+ISR | >=50 | 0.0% | 0.0% | n/a | n/a |
| 240 | 4C+ISR | all | 30.7% | 60.2% | 0.734 | 15.29 |
| 240 | 4C+ISR | <5 | 24.2% | 47.9% | 0.094 | 17.12 |
| 240 | 4C+ISR | >=50 | 5.2% | 69.0% | 0.965 | 14.95 |
| 240 | 5C+ISR | all | 55.5% | 15.6% | 0.474 | 15.79 |
| 240 | 5C+ISR | <5 | 62.0% | 15.5% | 0.094 | 13.64 |
| 240 | 5C+ISR | >=50 | 2.9% | 5.2% | 0.845 | 30.54 |
| 365 | 4C+ISR | all | 32.3% | 53.9% | 0.679 | 40.01 |
| 365 | 4C+ISR | <5 | 26.3% | 40.8% | 0.011 | 47.02 |
| 365 | 4C+ISR | >=50 | 14.3% | 62.9% | 0.937 | 32.31 |
| 365 | 5C+ISR | all | 51.4% | 16.1% | 0.590 | 31.53 |
| 365 | 5C+ISR | <5 | 60.3% | 15.3% | 0.050 | 27.00 |
| 365 | 5C+ISR | >=50 | 3.4% | 11.6% | 0.811 | 53.91 |

## ISR Emax Check

This table is included because a wrong ISR energy limit can make the ISR refit look mathematically valid while being physically too flexible, especially near threshold.

| Energy | Stored Emax [GeV] | single-W formula [GeV] | WW-threshold formula [GeV] |
|---|---|---|---|
| 162p5 | 61.368 | 61.368 | 1.721 |
| 240 | 106.538 | 106.538 | 66.152 |
| 365 | 173.648 | 173.648 | 147.093 |

## How To Read This

Use `plot_style_numbers.csv` when comparing to the plot boxes. Use `mass_split_stress_test.csv` when discussing why 5C fails. Use `isr_truth_recovery.csv` when discussing whether the ISR photon fit is physically meaningful.

The central question is not only whether ISR improves the mass plot. The sharper question is whether the equal-mass 5C constraint is being applied to events whose two 4C masses are already too different, and whether the ISR parameter is then absorbing that tension.
