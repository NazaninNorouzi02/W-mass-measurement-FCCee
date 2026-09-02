# Stage 3 x-resolution and jet-matching validation note

## Context

This note documents the checks performed after the suggestion to review the reco-gen jet matching in Stage 3. The main concern was that the angular and x-resolution plots looked too narrow and too clean compared with the reference thesis plots.

The variable under study is

```text
x = log(p/m)
```

This definition is kept because it is the definition used in the jet four-vector parametrization of the thesis. The goal of this note is not to change the definition of `x`, but to check whether the Stage 3 matching and filtering logic was artificially narrowing the distributions.

The Stage 3 file currently reads the Stage 2 jet-clustering output and builds `jets_gen4` and `jets_reco4` from the stored gen/reco jet four-vectors. It then performs reco-gen jet matching, computes the resolution variables, and writes diagnostic branches for the matching and for the momentum/mass components entering `x`.

---

## Original concern

The original reco-gen matching used a greedy one-to-one matching based on `DeltaR`, but it also applied a hard threshold:

```cpp
if (p.deltaR < dRThreshold) {
    matched_indices[p.reco_index] = p.gen_index;
}
```

and Stage 3 called it with

```python
greedyJetMatching(jets_reco4, jets_gen4, 0.1)
```

This means a reco-gen pair was accepted only if the matched distance satisfied

```text
DeltaR < 0.1
```

The concern was that this cut could remove the tails of the resolution distributions before the histograms were even filled. Since the resolution plots are supposed to measure the spread of the reco-gen differences, such a tight matching cut could make the result too narrow.

---

## Important distinction: DeltaR as a score vs DeltaR as a cut

In the updated logic, `DeltaR` is still used, but only as a matching score.

That means:

```text
DeltaR is used to choose the nearest gen jet for each reco jet.
DeltaR is no longer used to reject the match.
```

This is different from the original logic, where `DeltaR` was used twice:

```text
1. to rank the possible reco-gen pairs;
2. to reject matches with DeltaR larger than the threshold.
```

The second use was removed.

---

## Updated header: `greedyJetMatching_v2.h`

The updated header removes the matching threshold, adds proper phi wrapping in the `DeltaR` calculation, and adds diagnostic helper functions for checking the matched `DeltaR` values.

### Difference from the original header

The original header accepted a third argument:

```cpp
double dRThreshold = 0.1
```

and only accepted a match if

```cpp
p.deltaR < dRThreshold
```

The updated version removes that argument and always accepts the selected one-to-one greedy match. `DeltaR` is therefore still used to choose the nearest pairing, but it is not used as an acceptance cut.

A second change is that `dphi` is now wrapped to the interval `[-pi, pi]`. This avoids a fake large `DeltaR` when one jet has `phi` close to `+pi` and the other is close to `-pi`.

The added helper branches are only for diagnostics:

```text
matched_deltaR              all four matched DeltaR values per event
max_matched_deltaR          largest matched DeltaR in the event
n_matched_deltaR_gt_0p1     number of matched jets with DeltaR > 0.1
n_matched_deltaR_gt_1p0     number of matched jets with DeltaR > 1.0
```

---

## Stage 3 changes applied

The Stage 3 matching call was changed from

```python
df = df.Define("jet_match_indices", "greedyJetMatching(jets_reco4, jets_gen4, 0.1)")
```

to

```python
df = df.Define("jet_match_indices", "greedyJetMatching_v2(jets_reco4, jets_gen4)")
```

The event-level filter

```python
df = df.Filter("n_matched_jets == 4", "Require all 4 jets matched")
```

was commented out. After removing the `DeltaR` acceptance threshold, `n_matched_jets` is always 4 for this sample, so this filter no longer changes the event count.

The following diagnostic branches were added:

```python
df = df.Define("matched_deltaR", "getMatchedDeltaR(jets_reco4, jets_gen4, jet_match_indices)")
df = df.Define("max_matched_deltaR", "getMaxMatchedDeltaR(matched_deltaR)")
df = df.Define("n_matched_deltaR_gt_0p1", "countMatchedDeltaRAbove(matched_deltaR, 0.1)")
df = df.Define("n_matched_deltaR_gt_1p0", "countMatchedDeltaRAbove(matched_deltaR, 1.0)")
```

To diagnose the `x = log(p/m)` behavior, the following per-jet momentum and mass branches were also added:

```python
jet_gen_p_j1 ... jet_gen_p_j4
jet_reco_p_j1 ... jet_reco_p_j4
jet_gen_m_j1 ... jet_gen_m_j4
jet_reco_m_j1 ... jet_reco_m_j4
delta_logp_j1 ... delta_logp_j4
delta_logm_j1 ... delta_logm_j4
```

where

```text
delta_logp = log(p_reco / p_gen)
delta_logm = log(m_reco / m_gen)
```

Since

```text
delta_x = log(p_reco/m_reco) - log(p_gen/m_gen)
        = delta_logp - delta_logm
```

these branches help determine whether the shape of `delta_x` comes mainly from the momentum term, the mass term, or a cancellation between them.

---

## Fit-side fix applied

The fit code originally used

```python
bin_min = histogram.FindBin(xmin)
bin_max = histogram.FindBin(xmax)
```

For ROOT histograms, using `FindBin(xmax)` can reach the overflow bin when `xmax` is exactly the upper edge. This affected the initialization of the fit parameters and previously caused raw and filtered versions of the same variable to give different results.

The fit initialization was changed to use only the normal visible bins:

```python
bin_min = 1
bin_max = histogram.GetNbinsX()
```

After this change, raw and filtered fit results became consistent for the same visible range.

---

## Results from the tests

### 1. Original tight reco-gen matching: `DeltaR < 0.1`

With the old matching threshold, the `x = log(p/m)` distributions were already very narrow and had a sharp central peak with non-Gaussian tails. For example, for `j4` the fit gave approximately

```text
delta_x_j4:          N_fit = 37331, sigma around 0.0049
filtered_delta_x_j4: N_fit = 37331, sigma around 0.0025 before the fit-bin fix
```

The filtered/raw difference was later traced to the fit-bin initialization issue, not to the Stage 3 data itself.

### 2. Loosening the threshold to `DeltaR < 1.0`

Changing the matching call to

```python
greedyJetMatching(jets_reco4, jets_gen4, 1.0)
```

increased the number of accepted events, but did not significantly change the `x` shape. For example:

```text
delta_x_j4: N_fit = 42970, sigma = 0.00497 +/- 0.00382
```

This showed that the `DeltaR < 0.1` cut was rejecting events, but it was not the main reason for the narrow `x` shape.

### 3. Removing the reco-gen matching threshold completely

After changing to `greedyJetMatching_v2.h`, the matching no longer applies a `DeltaR` threshold. The result was:

```text
delta_x_j1: N_fit = 44865, sigma = 0.00925 +/- 0.00047
delta_x_j2: N_fit = 44660, sigma = 0.00300 +/- 0.00169
delta_x_j3: N_fit = 44554, sigma = 0.00456 +/- 0.00692
delta_x_j4: N_fit = 44377, sigma = 0.00445 +/- 0.00378
```

The corresponding filtered results were identical, confirming that the filtered branches were not changing the visible data after the fit-bin fix.

The angular resolutions remained very narrow:

```text
delta_theta_j1: sigma = 0.00059 +/- 0.00044
delta_theta_j2: sigma = 0.00067 +/- 0.00012
delta_theta_j3: sigma = 0.00047 +/- 0.00045
delta_theta_j4: sigma = 0.00050 +/- 0.00017

delta_phi_j1: sigma = 0.00092 +/- 0.00010
delta_phi_j2: sigma = 0.00081 +/- 0.00011
delta_phi_j3: sigma = 0.00060 +/- 0.00004
delta_phi_j4: sigma = 0.00055 +/- 0.00042
```

The alpha issue seen earlier for `j4` was also fixed after the fit-bin correction:

```text
delta_alpha_j4:          N_fit = 25369, sigma = 0.04538 +/- 0.00148
filtered_delta_alpha_j4: N_fit = 25369, sigma = 0.04538 +/- 0.00148
```

### 4. Removing `n_matched_jets == 4`

After removing the `DeltaR` threshold from the matching, the branch `n_matched_jets` was always 4:

```text
n_matched_jets: mean = 4.000, std dev = 0
```

Therefore, removing the filter

```python
df.Filter("n_matched_jets == 4")
```

did not change the results or the event count. This confirms that the Stage 3 matching filter is now redundant.

### 5. Stage 2 event count

The Stage 3 output has

```text
45854 entries
```

and the Stage 2 file also has

```text
45854 entries
```

Therefore, after the Stage 3 matching changes, Stage 3 is no longer reducing the event count. The current sample size and the clean jet behavior are already present in the Stage 2 output.

---

## ROOT/TBrowser commands used for diagnostics

These commands were used to inspect the matching, momentum, mass, and `x` components. Use them in TBrowser of the ROOT file.

### Matching diagnostics

```cpp
events->Draw("matched_deltaR>>h_matched_deltaR(100,0,1.2)")
events->Draw("max_matched_deltaR>>h_max_matched_deltaR(100,0,3.5)")
events->Draw("n_matched_deltaR_gt_0p1>>h_n_matched_deltaR_gt_0p1(5,-0.5,4.5)")
events->Draw("n_matched_deltaR_gt_1p0>>h_n_matched_deltaR_gt_1p0(5,-0.5,4.5)")
events->Draw("n_matched_jets>>h_n_matched_jets(5,-0.5,4.5)")
```

Observed behavior:

```text
matched_deltaR: mean = 0.029, std dev = 0.084
max_matched_deltaR: mean = 0.12, std dev = 0.39
n_matched_deltaR_gt_0p1: mean = 0.2777
n_matched_deltaR_gt_1p0: mean = 0.03585
n_matched_jets: mean = 4.000, std dev = 0
```

This confirms that most matched jets are close in `DeltaR`, even without the threshold. Some matches above 0.1 exist, which explains why removing the old threshold increased the statistics, but the shape of the resolution plots did not change significantly.

### Jet momentum and mass checks for j3

```cpp
events->Draw("jet_reco_p_j3:jet_gen_p_j3", "", "colz")
events->Draw("jet_reco_m_j3:jet_gen_m_j3", "", "colz")
events->Draw("jet_reco_p_j3 - jet_gen_p_j3>>h_delta_p_j3(100,-10,10)")
events->Draw("jet_reco_m_j3 - jet_gen_m_j3>>h_delta_m_j3(100,-10,10)")
events->Draw("delta_logp_j3>>h_delta_logp_j3(100,-0.5,0.5)")
events->Draw("delta_logm_j3>>h_delta_logm_j3(100,-0.5,0.5)")
```

Observed behavior for j3:

```text
jet_gen_p_j3:  mean = 34.32, std dev = 7.478
jet_reco_p_j3: mean = 34.16, std dev = 7.404
jet_gen_m_j3:  mean = 10.80, std dev = 5.357
jet_reco_m_j3: mean = 10.37, std dev = 5.208

jet_reco_p_j3 - jet_gen_p_j3: mean = -0.065, std dev = 1.781
jet_reco_m_j3 - jet_gen_m_j3: mean = -0.366, std dev = 2.015

delta_logp_j3: mean = -0.00229, std dev = 0.0742
delta_logm_j3: mean = -0.0311, std dev = 0.1274
```

The 2D plots show a strong diagonal structure between reco and gen momentum/mass. This explains why the difference variables are sharply peaked. Expected!

### Jet momentum and mass checks for j4

```cpp
events->Draw("jet_reco_p_j4:jet_gen_p_j4", "", "colz")
events->Draw("jet_reco_m_j4:jet_gen_m_j4", "", "colz")
events->Draw("jet_reco_p_j4 - jet_gen_p_j4>>h_delta_p_j4(100,-10,10)")
events->Draw("jet_reco_m_j4 - jet_gen_m_j4>>h_delta_m_j4(100,-10,10)")
events->Draw("delta_logp_j4>>h_delta_logp_j4(100,-0.5,0.5)")
events->Draw("delta_logm_j4>>h_delta_logm_j4(100,-0.5,0.5)")
```

Observed behavior for j4:

```text
jet_gen_p_j4:  mean = 29.40, std dev = 9.129
jet_reco_p_j4: mean = 28.90, std dev = 8.989
jet_gen_m_j4:  mean = 9.668, std dev = 5.041
jet_reco_m_j4: mean = 9.180, std dev = 4.885

jet_reco_p_j4 - jet_gen_p_j4: mean = -0.2045, std dev = 1.835
jet_reco_m_j4 - jet_gen_m_j4: mean = -0.4360, std dev = 1.905

delta_logp_j4: mean = -0.00793, std dev = 0.08398
delta_logm_j4: mean = -0.03869, std dev = 0.1326
```

Again, reco and gen are strongly correlated event-by-event. Expected!

---

## Interpretation so far

The tests show that the original `DeltaR < 0.1` matching threshold did remove events, but it is not the main reason for the narrow resolution shapes.

The strongest evidence is:

```text
1. Removing the threshold increased the number of entries.
2. The central shapes and fitted sigmas changed only mildly.
3. Removing n_matched_jets == 4 changed nothing.
4. Stage 3 now has the same number of entries as Stage 2.
5. Reco and gen jet momentum/mass are strongly correlated already at Stage 2.
```

So the current conclusion is:

```text
The narrow x and angular resolution shapes are not mainly caused by the Stage 3 matching threshold. They appear to come from the jet collections entering Stage 3, i.e. from the Stage 2 construction of gen/reco jets and the strong event-by-event similarity between those jets.
```

The next step should therefore be to audit Stage 2: how `jet_gen_*` and `jet_reco_*` are built, from which input collections, and whether the reco jet collection is truly detector-level and independent enough from the gen jet collection.

---

## Important plots

```text
01_max_matched_deltaR_per_event.png
02_matched_deltaR_all_jets.png
03_n_matched_deltaR_above_0p1_per_event.png
04_j3_reco_vs_gen_momentum.png
05_j3_reco_vs_gen_mass.png
06_j3_delta_logp.png
07_j3_delta_logm.png
08_j3_delta_momentum.png
09_j3_delta_mass.png
```

These plots show:

```text
- the old DeltaR threshold was not the dominant issue;
- the no-threshold matcher still gives mostly close matches;
- reco and gen jet p/m are strongly correlated;
- the narrow x behavior is already present in the momentum/mass components.
```
