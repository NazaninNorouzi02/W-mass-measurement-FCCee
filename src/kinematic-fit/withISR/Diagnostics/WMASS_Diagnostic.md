# Diagnostic Report on the 5C Fit and ISR Reconstruction

## Purpose of the Diagnostic Study

After the ISR reconstruction was implemented, the 4C pulls were reasonably well behaved, but the 5C pulls were extremely broad and a large fraction of events had P<sub>5C</sub> < 0.03. A temporary diagnostic version of the code was therefore created to determine whether this came from:

- An incorrect derivative of the fifth constraint
- An unstable matrix inversion in the 5C pull calculation
- Numerical failure of the minimisation
- Incompatibility of the events with the additional equal-mass constraint

The diagnostic code did not change the fitted masses, probabilities, ISR decisions, jet pairing, or event selection. It only appended extra quantities after the existing outputs. Therefore, returning afterwards to the main ISR code removed the diagnostic branches but did not invalidate the conclusions obtained from them.

## Diagnostic Information Added to the Code

For every valid fit, the temporary version saved:

- The full 5×5 matrix **AVA<sup>T</sup>** used in the 5C pull denominator
- Its eigenvalues and condition number
- Whether the matrix was finite and positive definite
- The numerical derivatives of the fifth constraint with respect to every jet's θ and ϕ
- An independent analytic calculation of the same derivatives
- The difference between the analytic and numerical derivatives
- The difference **Δχ<sup>2</sup> = χ<sup>2</sup><sub>5C</sub> − χ<sup>2</sup><sub>4C</sub>**
- The dijet-mass difference before imposing the fifth constraint, **ΔM<sub>4C</sub> = M<sub>1,4C</sub> − M<sub>2,4C</sub>**
- The largest 5C angular pull in each event
- The selected pairing, active parameter boundaries, and event index

A separate diagnostic program then examined the worst-pull events and tested whether the pull magnitude was associated with the derivative error, matrix condition number, or pre-5C mass difference.

## What the 365 GeV Audit Demonstrated

The audit was performed on the preliminary processed subset containing 45,583 events. This output was produced using only 10<sup>−6</sup> of each selected process sample, rather than the full production samples, so the numerical values remain preliminary.

| Result in processed subset | 4C | 5C |
|---|---|---|
| Valid ordinary fits | 45,581 | 45,562 |
| Already satisfying P > 0.03 | 24,819 | 4,499 |
| Events triggering ISR treatment | 20,762 | 41,063 |
| Recovered above P = 0.03 by ISR | 11,051 | 6,757 |
| Total satisfying P > 0.03 after recovery | 35,870 | 11,256 |

Thus:

- Approximately 54.5% of valid 4C fits passed directly, compared with only 9.9% of valid 5C fits
- ISR recovered approximately 53.2% of failed 4C events, but only 16.5% of failed 5C events
- The final acceptance was approximately 78.7% for 4C and 24.7% for 5C

These numbers prove that the main deterioration appears when the fifth constraint is introduced. It is not simply a general failure of the four-momentum fit or an inability to run the ISR minimisation.

The fact that every triggered ISR fit was reported as valid means the minimiser found finite solutions satisfying the implemented constraint tolerance. However, "valid fit" does not mean "good fit": most 5C+ISR solutions still had P < 0.03.

## Test of the Fifth-Constraint Jacobian

The fifth constraint is:

**f<sub>5</sub> = M(j<sub>a</sub> + j<sub>b</sub>) − M(j<sub>c</sub> + j<sub>d</sub>) = 0**

The supervisor's first hypothesis was that the numerical derivatives of this condition might be incorrect, which would corrupt the 5C pull denominator.

The analytic and finite-difference derivatives agreed extremely closely over the event population:

- 99th-percentile relative discrepancy = 2.897 × 10<sup>−7</sup>
- 99th-percentile absolute discrepancy = 5.145 × 10<sup>−8</sup> GeV/rad

This effectively rules out a general error in the numerical fifth-row derivative as the cause of the broad 5C pull population.

A few of the most extreme events showed larger discrepancies when a fitted angle was close to a parameter boundary, because a symmetric finite difference is no longer reliable after the angle is clamped. These are genuine numerical edge cases, but they are rare and cannot explain the overall 5C behavior.

## Test of AVA<sup>T</sup> Conditioning

If **AVA<sup>T</sup>** were almost singular, its inverse could create artificially large pull values. The diagnostic compared the maximum angular pull with the matrix condition number.

The measured association was weak:

- **ρ<sub>Spearman</sub> = 0.026**

The median condition number in broad-pull events was only 1.06 times that in ordinary-pull events.

Therefore, some individual events are ill-conditioned, but matrix conditioning does not explain the broad 5C pull distribution across the full sample.

## Test of the Equal-Mass Condition

The strongest observed relationship was between broad 5C pulls and the mass difference already present after the 4C fit:

**|ΔM<sub>4C</sub>| = |M<sub>1,4C</sub> − M<sub>2,4C</sub>|**

The worst 5C-pull events frequently had very large pre-5C mass differences—for example approximately 49, 74, 86, 109, 115, 119, 137, and 143 GeV.

The 5C fit must force these two masses to become exactly equal. When the masses initially differ so strongly, the fit can satisfy the condition only by applying very large changes to the jet parameters. The estimated uncertainty of these changes may remain small, producing extremely large pulls and a very large χ<sup>2</sup>.

This is the main explanation supported by the diagnostic results. It does not yet prove why the two reconstructed masses differ so strongly; it proves that the broad pulls arise primarily when the equal-mass condition is imposed on such events.

## Additional Numerical Issue: Nested-Fit Ordering

Because 5C contains all four 4C constraints plus one additional condition, an ideally converged fit should satisfy:

**χ<sup>2</sup><sub>5C</sub> ≥ χ<sup>2</sup><sub>4C</sub>**

Nevertheless, 181 of the 45,560 events with both fits valid had:

**χ<sup>2</sup><sub>5C</sub> − χ<sup>2</sup><sub>4C</sub> < −0.001**

This is approximately 0.4% of the valid sample. Some differences were very small and compatible with limited minimisation precision, but a smaller number were extremely negative and associated with parameter boundaries or abnormal solutions.

Therefore, a real numerical pathology exists in a small subset and those events should be flagged or rejected. However, their number is far too small to explain why most 5C events have low probability.

## ISR Reconstruction Result

The 365 GeV mass comparison was:

| Reconstruction | Mean [GeV] | Width [GeV] |
|---|---|---|
| Without ISR treatment | 81.25 | 11.85 |
| With ISR treatment | 79.04 | 10.34 |
| Thesis, without ISR treatment | 82.18 | 10.15 |
| Thesis, with ISR treatment | 79.0 | 8.49 |

ISR treatment moves the mean in the same direction as the thesis and reproduces its ISR-treated mean very closely. It also narrows the distribution, showing that the implemented ISR recovery has the expected physical effect.

The remaining concern is the width: the reconstructed ISR-treated distribution is still broader than the thesis result. The fitted ISR parameter reached the photon-energy endpoint in only 422 events, and |y<sub>ISR</sub>| ≥ 4.9 occurred in only 105 events. Therefore, widespread saturation of the ISR parameter is not the main explanation.

## Overall Conclusion

The diagnostics established that:

- The fifth-row Jacobian is correct for the general event population
- Instability of **AVA<sup>T</sup>** is not the population-wide cause
- ISR treatment is operating and improves both the mass distribution and fit probabilities
- The main 5C difficulty appears when events with very different 4C dijet masses are forced to satisfy exact mass equality
- A small separate group of numerical failures exists, identified by χ<sup>2</sup><sub>5C</sub> < χ<sup>2</sup><sub>4C</sub> and parameter-boundary behavior

Consequently, the 5C result should not yet be called fully validated. The most probable underlying causes worth investigating are:

- Incorrect jet pairing in a subset of events
- Non-Gaussian jet response or missing correlations in the covariance model
- Physical radiation or reconstruction effects not completely represented by one collinear ISR photon

The present diagnostics distinguish these from the two initially suspected causes—Jacobian error and general matrix singularity—which were not supported by the data.

The important question is whether the large pre-5C mass differences are expected for the selected reconstructed sample, or whether an additional pairing/quality selection or covariance refinement should be applied before the equal-mass 5C hypothesis.
