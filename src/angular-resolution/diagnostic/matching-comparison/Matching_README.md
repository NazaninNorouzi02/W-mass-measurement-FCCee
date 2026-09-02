# Jet-Parton Matching Strategy Comparison

This analysis studies the effect of three different jet-parton matching strategies on the reconstructed jet observables in an FCC-ee $WW \to 4q$ sample at $\sqrt{s} = 160$ GeV.  
The goal is to check whether the choice of matching algorithm changes the final physics result, especially the jet energy response observable $\Delta \alpha$ and the matching efficiency.

---

## 1. Matching Strategies

Three strategies were implemented and compared:

### Strategy A: Original Non-Exclusive Matching
This is the baseline matching method.

**Core idea:**
- Each reconstructed jet is matched independently to the nearest available parton.
- The matching is **non-exclusive**, meaning that more than one jet can be associated with the same parton.

**Interpretation:**
- This is the simplest and most local method.
- It provides a reference point for comparison with more constrained algorithms.
- Because it does not enforce a global one-to-one assignment, it can in principle produce duplicate assignments in ambiguous events.

---

### Strategy B: Unified Greedy Matching
This is the exclusive direct matching strategy.

**Core idea:**
- Reconstructed jets and partons are matched with a **greedy global algorithm**.
- The algorithm chooses the best available reco-parton pair step by step, prioritizing the smallest angular distance.
- Once a jet or parton is assigned, it is removed from further consideration.

**Interpretation:**
- This enforces a **one-to-one exclusive matching**.
- It is more physically constrained than Strategy A.
- It is useful when one wants to avoid multiple jets being assigned to the same parton.

---

### Strategy C: Transitive Matching
This strategy matches through an intermediate generator-level object.

**Core idea:**
- First match **reco jet $\to$ gen jet**.
- Then match **gen jet $\to$ parton**.
- The final reco-to-parton association is obtained through this chain.

**Interpretation:**
- This method reflects the event shower history more explicitly.
- It is the most structured of the three strategies.
- It can reduce direct reco-parton ambiguity by using the gen-level jet as a bridge.

---

## 2. Code Implementation Concept

The matching strategies are implemented inside the `RDFanalysis` workflow and the comparison plot script reads the corresponding branches from the ROOT output.

A simplified conceptual form of the implementation is:


---
## 3. ROOT Output Branches

| Branch | Description |
|--------|-------------|
| n_matched_A | Number of matched partons per event (Strategy A) |
| n_matched_B | Number of matched partons per event (Strategy B) |
| n_matched_C | Number of matched partons per event (Strategy C) |
| delta_alpha_A | Jet energy response values (Strategy A) |
| delta_alpha_B | Jet energy response values (Strategy B) |
| delta_alpha_C | Jet energy response values (Strategy C) |

The main observable is the fractional energy response:

Delta_alpha = (E_reco - E_parton) / E_parton

- Delta_alpha = 0: reconstructed jet energy equals parton energy
- Delta_alpha < 0: reconstructed jet has less energy than matched parton
- Delta_alpha > 0: reconstructed jet has more energy than matched parton

---

## 4. Diagnostic Quantities

**Energy balance:** (sum E_reco_jets) / (sum E_partons)

**Minimum jet separation:** min Delta_R(reco jet i, reco jet j)

---

## 5. Running the Comparison (Plot_Matching_comparison.py)



The script saves plots to the outputs/ directory in both PNG and PDF formats.

---

## 6. Results

**Event sample:** 23,892 events, sqrt(s) = 160 GeV, WW -> 4q

### Energy Balance

| Quantity | Value |
|----------|-------|
| Mean | 0.847 |
| Std Dev | 0.081 |

### Minimum Jet Separation

| Quantity | Value |
|----------|-------|
| Mean | 1.165 |
| Std Dev | 0.422 |

### Full 4-Parton Matching Efficiency

| Strategy | Efficiency |
|----------|------------|
| A | 74.9% |
| B | 74.6% |
| C | 74.8% |

### Jet 4 Delta_alpha Results

| Strategy | Mean | RMS |
|----------|------|-----|
| A | -0.286 | 0.20 |
| B | -0.285 | 0.20 |
| C | -0.28 | 0.20 |

### Mean Delta_alpha Per Jet (All Strategies Combined)

| Jet | Mean Delta_alpha |
|-----|------------------|
| 1 | -0.20 |
| 2 | -0.12 |
| 3 | -0.19 |
| 4 | -0.29 |

### RMS Delta_alpha Per Jet (All Strategies Combined)

| Jet | RMS |
|-----|-----|
| 1 | 0.21 |
| 2 | 0.17 |
| 3 | 0.18 |
| 4 | 0.20 |

---

## 7. Conclusion

All three matching strategies produce nearly identical results.
The agreement shows that for this clean WW->4q sample at sqrt(s) = 160 GeV, the choice of jet-parton matching strategy has negligible impact on the Delta_alpha observable.
Matching-strategy dependence is not a dominant systematic uncertainty for this analysis configuration.

---

## 8. Limitations

- Strategy A: Can assign more than one jet to the same parton.
- Strategy B: Greedy matching is exclusive, but not always globally optimal in every possible event.
- Strategy C: Depends on the correctness of both matching steps: reco to gen and gen to parton.


```python
# Strategy A: non-exclusive direct matching
df = df.Define("match_A", "jetPartonMatching(jets_reco4, parton_eta, parton_phi)")

# Strategy B: unified greedy direct matching
df = df.Define("match_B", "greedyMatchToPseudoJets(jets_reco4, pseudo_jets_parton, 0.4)")

# Strategy C: transitive reco -> gen -> parton matching
df = df.Define("jet_match_indices", "greedyJetMatching(jets_reco4, jets_gen4, 0.1)")
df = df.Define("gen_to_parton_C", "greedyMatchEtaPhi(gen_eta, gen_phi, parton_eta, parton_phi, 0.4)")
df = df.Define("match_C", "transitiveMatch(jet_match_indices, gen_to_parton_C)")```









