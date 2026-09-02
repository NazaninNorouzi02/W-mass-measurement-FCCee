# Angular Resolution Photon and ISR Study

This work started because the
\[
\delta_\alpha = \frac{E_{\mathrm{reco}}-E_{\mathrm{parton}}}{E_{\mathrm{parton}}}
\]
distribution, especially for jet 4, was strongly shifted to negative values and did not look like the expected resolution shape.

The first assumption was that this could come from ISR, the event selection, the selected quarks, or the reconstructed energy being much lower than the parton energy. The diagnostic study showed that the largest effect was caused by removing all reconstructed photons before jet clustering.

## What we found

When all reconstructed photons were removed, the total reconstructed jet energy was much lower than the total parton energy. This produced a large negative shift in delta alpha, especially for jet 4.

When photons were kept, the reconstructed jet-energy sum became much closer to the parton-energy sum and the delta-alpha distribution improved significantly.

This is expected because photons from hadron decays, especially neutral mesons such as pi0 -> gamma gamma, are part of the visible energy of a hadronic jet. Removing the photon objects removes that measured energy. The parent pi0 is normally not present as another reconstructed object that can replace the removed photon energy.

The ISR/no-ISR comparison also showed that ISR was not the main reason for the large negative shift in the case where all photons were removed. Even the no-ISR sample remained strongly shifted.

The current baseline is therefore to keep normal jet photons and remove only reconstructed photons that pass the current ISR-candidate selection.

## Files used in this study

### `energyDiagnostics.h`

This header contains the helper functions used for the energy and photon diagnostic study. It includes reconstructed and MC photon-energy sums, ISR-candidate counts, reconstructed jet-energy sums, and safe access to diagnostic values.

It was used to compare:

- photons kept;
- all photons removed;
- reconstructed jet energy against parton energy;
- event-level energy closure;
- photon energy fraction;
- ISR and no-ISR samples;
- delta alpha for each jet.

### `selectISRPhotons.h`

This header selects reconstructed photons that satisfy the current ISR-like definition:

```text
E_gamma > 1 GeV
|cos(theta_gamma)| > 0.95
```

Only the selected photons are removed before reconstructed jet clustering. Other reconstructed photons remain in the event and can be clustered into the hadronic jets.

This is currently a reconstructed-level ISR approximation and still needs to be checked.

### `getElement_v2.h`

This version is used for safer access to individual jet values. Invalid or unavailable values are returned as NaN instead of using a numerical placeholder that could enter plots or later calculations.

### Energy diagnostic analysis

The energy diagnostic code compares photon-kept and photon-removed reconstruction. It produces the energy-closure, photon-fraction, reconstructed-versus-parton energy, and delta-alpha diagnostic plots.

These plots were the main evidence that removing every photon was causing the large reconstructed-energy loss.

### Photon-treatment analysis

The photon-treatment study was used to compare the different reconstruction choices and understand which photon treatment caused the change in the delta-alpha shape.

The main cases were:

- all photons removed;
- all photons kept;
- only ISR-like photons removed.

### `angular_resolution_ISR.py`

This is the multi-energy angular-resolution analysis where only the current ISR-candidate photons are removed.

It runs for:

- 160 GeV;
- 240 GeV;
- 340 GeV;
- 345 GeV;
- 350 GeV;
- 355 GeV;
- 365 GeV.

The output branch names are kept the same as in the original angular-resolution analysis so the same fitting structure can be used.

The code now uses `getElement_v2.h`. For filtered scalar values, rejected entries are stored as:

```cpp
std::numeric_limits<float>::quiet_NaN()
```

instead of:

```cpp
-999.0f
```

The processed fraction was also increased from `1e-6` to `1e-4`, so the number of processed and selected events is expected to change.

### `crystalball_ISR.py`

This code applies the Double-Sided Crystal Ball fit to the ISR-removed angular-resolution outputs for all energy levels.

The output contains only:

- PNG fit plots;
- JSON fit results.

PDF output was removed because it made the result directory too heavy.

The finite-value selection in the fit was changed from:

```python
f"({branch} == {branch})"
```

to:

```python
f"TMath::Finite({branch}) && "
```

Both can reject NaN values, but the second form makes the intention clearer and also explicitly rejects infinite values.

### `MultiE_Plots_ISR`

This directory contains the current multi-energy Crystal Ball fit plots and fit results from the ISR-removed analysis.

## Important note about ROOT entries and NaN

Replacing an unwanted value with NaN does not remove the event from the ROOT tree. The tree still has one entry for that event, but the selected branch contains NaN.

Because of this, the raw number shown as `Entries` in a ROOT tree or histogram should not automatically be used as the number of valid events after filtering.

For a cut-flow table or reported analysis statistics, the count should be made after applying the actual condition, for example:

```cpp
TMath::Finite(delta_alpha_j4)
```

together with any ISR, matching, range, or event-selection cuts.

This should be handled before reporting the final event statistics.

## Fit status 4000

Plots show `Fit status = 4000` even though the curve and chi-square look reasonable.

The current fit option contains `M`:

```python
histogram.Fit(fit_function, "SMRL")
```

In ROOT, `M` runs an additional `IMPROVE` step after the main minimization. A status of 4000 is connected to this optional improvement stage not finishing cleanly. It does not automatically mean that the fitted curve or main minimum is wrong, but the fit should not be treated as fully checked only from the plot.

A simple comparison is to run the same fit without `M`:

```python
histogram.Fit(fit_function, "SRL")
```

and compare the parameters. The fit validity, covariance status, and stability of mu and sigma should also be checked. This status was not visible in the original plots because the status was either not printed or the fitting options/output were different.

## Things I still want to check

- Reimplement the fitting with RooFit.
- Check the ISR definition and confirm whether the energy and angular cuts are suitable.
- Compare the reconstructed ISR selection with truth information if the required parent history is available.
- Make a proper cut-flow table using finite-value selections instead of raw ROOT entry numbers.
- Check the x variable again. Its distribution looks unusually narrow and the current range from -0.6 to 0.6 looks too large for the observed values.
- Confirm that `getXGen` and `getXReco` use the intended definition and subtraction order.
- Check the Crystal Ball fit stability with RooFit, different initial parameters, binning, and fit ranges.
- Add pull plots or another fit-quality check when the final fitting method is selected.
- Chack the delta_phi_j3 and delta_phi_j4 bad fitting.

The current results are encouraging because the main negative delta-alpha shift was traced to the removal of all photon energy. The remaining checks are mostly about validating the ISR definition, the x variable, event counting, and the final fitting method.
