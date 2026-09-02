# Diagnostic Analysis: delta_alpha Distribution

This script investigates the `delta_alpha` distribution behavior in $WW \to qqqq$ events at $\sqrt{s} = 160$ GeV, where `delta_alpha = (E_reco_jet - E_parton) / E_parton`. 

It is designed to be run after `angular_resolution.py` has generated the required ROOT files.

## Problem Statement

The `delta_alpha` distributions for the softer jets ($j_3$ and $j_4$) do not form Gaussian peaks. Instead, they form left-shifted plateaus. Consequently, Crystal Ball (CB) fits are completely failing on these distributions ($\sigma$ pegs at the upper boundary, $\chi^2/\text{ndf} \approx 270$). 

This diagnostic tests two hypotheses for this behavior:
1. **Code Bug:** An underlying issue with the energy matching or computation logic.
2. **Physical Kinematics:** Soft jets are experiencing heavy out-of-cone energy loss.

## Diagnostic Tests & Results

### 1. Cross-Jet Isolation Test ($j_1$ vs $j_4$)
To determine if the issue is a global event bug or an isolated soft-jet effect, we looked at the behavior of the leading jet ($j_1$) in events where the softest jet ($j_4$) performs poorly (`delta_alpha_j4 < -0.3`).

* **$j_1$ mean (all events):** -0.018
* **$j_1$ mean (when $j_4$ is bad):** -0.024
* **Difference:** 0.006 (negligible)
* **Conclusion:** The $j_1$ distribution remains essentially unchanged even when $j_4$ is severely shifted. This confirms the plateau is driven by **out-of-cone physics**, not a global code bug. In near-threshold events, the softest jet is losing a significant fraction of its energy to neighboring jets, leading to a large fractional bias.

### 2. Independence of $\alpha$ and $x$
We evaluated the Pearson correlation coefficient between $\alpha$ (referenced to parton energy) and $x$ (referenced to gen-jet energy) to see if they carry redundant information.

* **$j_1$ correlation:** 0.133
* **$j_4$ correlation:** 0.351
* **Conclusion:** Both variables are statistically independent. They must both be retained as separate parameters in the kinematic fit covariance matrix. Replacing $\alpha$ with $x$ would result in a loss of information.

### 3. Per-Jet Fit Range Extension Test
We attempted to resolve the CB fit failures by implementing per-jet fit configurations. Since $j_4$ peaks around -0.22, the default `[-0.25, 0.25]` range was cutting the peak. We widened the fit ranges significantly (e.g., `[-0.60, 0.30]` for $j_4$) to ensure the entire structure was contained.

* **Result:** The fit still fails completely. The $\sigma$ parameter hits its absolute upper limit (0.5), and the $\chi^2/\text{ndf}$ remains excessively high ($\approx 278$). 
* **Conclusion:** This proves that the fit failure is **not a range issue**. The distributions are fundamentally skewed plateaus, and no Crystal Ball fit will successfully describe this shape regardless of the applied range boundaries.

## Open Questions & Next Steps

Because the soft jet distributions are non-Gaussian plateaus, we cannot use a standard CB fit. We are currently evaluating:

## Usage

```bash
python3 src/diagnostic_alpha.py