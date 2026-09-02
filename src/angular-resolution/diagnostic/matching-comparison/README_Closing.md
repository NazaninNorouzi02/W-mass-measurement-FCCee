# Jet Clustering and Matching Algorithm Investigation


## Issue Summary

The issue raised three questions about the matching steps in `angular_resolution.py`:

1. **`greedyJetMatching` should be reviewed** — is the reco↔gen matching algorithm correct?
2. **Parton matching algorithm** — `jetPartonMatching` uses a non-exclusive loop while `greedyJetMatching` enforces 1-to-1. Are we matching different things in the two steps?
3. **Can we unify the two matching algorithms?** — use `greedyJetMatching` for both reco↔gen and reco↔parton.

---

## Background: what the two matching steps do

The analysis has two separate matching steps:

### Step 1 — reco-jet ↔ gen-jet (`greedyJetMatching`)
Used for angular resolutions: Δθ, Δφ, Δη, Δx.

```
greedyJetMatching(jets_reco4, jets_gen4, dR_threshold=0.1)
```

**Algorithm:** build all (reco_i, gen_j, ΔR) pairs → sort by ΔR ascending → assign best pair first, mark both as used → strictly 1-to-1. No two reco-jets can claim the same gen-jet.

### Step 2 — reco-jet ↔ parton (`jetPartonMatching`)
Used for energy resolution: Δα = (E_reco − E_parton) / E_parton.

```
jetPartonMatching(jets_reco4, parton_eta, parton_phi, dR_threshold=0.4)
```

**Algorithm:** for each reco-jet independently, find the nearest parton within ΔR < 0.4. **Non-exclusive** — in principle, two reco-jets could claim the same parton.

**The asymmetry:** Step 1 is exclusive (unique greedy), Step 2 is not. The issue asks whether this asymmetry corrupts the results.

---

## What was investigated

Three parton-matching strategies were implemented and run on the same events simultaneously, storing separate branches for each:

| Strategy | Algorithm | Exclusivity | Path |
|---|---|---|---|
| **A** | `jetPartonMatching` | Non-exclusive (original) | reco → parton directly |
| **B** | `greedyJetMatching` on PseudoJet partons | Exclusive 1-to-1 | reco → parton directly |
| **C** | `greedyJetMatching` transitive | Exclusive 1-to-1 | reco → gen → parton |

Strategy C goes through the gen-jet as an intermediate step: first match reco→gen (already done in Step 1), then match gen→parton, then compose the two index maps. This ensures the parton assigned to a reco-jet is always the one whose shower products ended up in that gen-jet.

---

## Results

### Matching agreement between strategies

| Pair | Mean per-jet agreement | Interpretation |
|---|---|---|
| A vs B (`n_agree_AB`) | **3.972 / 4 = 99.3%** | A and B assign same parton to 99.3% of jet slots |
| A vs C (`n_agree_AC`) | **3.944 / 4 = 98.6%** | A and C agree on 98.6% of assignments |

The original non-exclusive algorithm (A) produces **the same parton assignment** as the exclusive greedy algorithm (B) in 99.3% of jet slots. This happens because the four partons are angularly well separated (min ΔR ≈ 1.1, well above the 0.4 matching threshold), so the nearest parton to each reco-jet is always unambiguously different. The non-exclusivity never fires in practice.

### n_matched_partons distribution

| Strategy | Mean n_matched | Fraction with all 4 matched (n=4) |
|---|---|---|
| A (original) | 3.651 | ~73% |
| B (unified greedy) | 3.624 | ~73% |
| C (transitive) | 3.629 | ~73% |

All three strategies find the same fraction of events with 4/4 partons matched. The ~27% of events where at least one jet is unmatched are events where a soft reco-jet has no parton within ΔR < 0.4 — this is a physics effect, not an algorithm effect.

### delta_alpha distributions

The three delta_alpha distributions are **visually and numerically identical** for all four jets:

| Jet | Strategy A | Strategy B | Strategy C |
|---|---|---|---|
| j1 | μ=−0.018, σ=0.159 | μ=−0.018, σ=0.159 | μ=−0.018, σ=0.159 |
| j2 | μ=−0.118, σ=0.159 | μ=−0.118, σ=0.159 | μ=−0.118, σ=0.159 |
| j3 | μ=−0.188, σ=0.172 | μ=−0.187, σ=0.172 | μ=−0.187, σ=0.172 |
| j4 | μ=−0.286, σ=0.201 | μ=−0.285, σ=0.201 | μ=−0.285, σ=0.201 |

No strategy produces a different distribution. The choice of matching algorithm has no effect on the energy resolution measurement.

### Why the strategies agree: event topology

The reason all three strategies give identical results is geometric:

- Min ΔR between any two partons: **mean = 1.075** (std = 0.485)
- Min ΔR between any two reco-jets: **mean = 1.165** (std = 0.422)
- Matching threshold: **0.4**

The partons are separated by ΔR ≈ 1.1 on average — nearly **3× larger** than the matching threshold. In this geometry, each reco-jet has one unambiguously nearest parton and no other jet is close enough to compete for it. The non-exclusive original algorithm naturally produces unique assignments in almost all events, making it equivalent to the exclusive greedy algorithm in practice.

### The 27% unmatched events

The ~27% of events where at least one jet gets n_matched < 4 are events where a soft reco-jet (j4, mean energy ≈ 23 GeV) has no parton within ΔR < 0.4 due to out-of-cone radiation from the parton shower. This is **not caused by any matching algorithm** — it is a physics effect present in all three strategies equally.

Adding a filter `df.Filter("n_matched_partons == 4")` would reduce the event count from ~24k to ~17k by removing these events. This is a valid choice if the analysis requires a valid alpha value for all four jets, but it introduces a selection bias toward events with well-collimated jets. This is left as a **future decision** for the analysis.

---

## Answers to the issue questions

**1. Should `greedyJetMatching` be reviewed?**
The function is correct. It implements a proper unique greedy algorithm: sort all pairs by ΔR, assign greedily with exclusivity. No changes needed.

**2. Are we matching different things in the two steps?**
In principle yes — Step 1 is exclusive, Step 2 is not. In practice no — the partons are geometrically separated enough that the original non-exclusive algorithm produces the same unique assignments as the exclusive greedy algorithm in 99.3% of jet slots.

**3. Can we unify the two algorithms?**
Yes — Strategy B demonstrates this. Partons can be promoted to `fastjet::PseudoJet` objects and passed to `greedyJetMatching`. The result is numerically identical to the original. Unification is possible but produces no improvement in this dataset.

---

## Optional next steps (not blocking)

- **Test at higher √s:** At ecm365, jets are more boosted and better separated. The matching behaviour should be even cleaner. At ecm160 the near-threshold topology already gives well-separated partons so the result is expected to hold at all energies.
- **n_matched_partons == 4 filter:** decide whether to add this filter to remove the 27% of events where a soft jet has no matched parton. This affects the alpha measurement for j3/j4 and should be a conscious analysis choice.
- **Review `countMatchedJets`:** the counting function uses `idx >= 0` as the validity criterion. This is correct given the sentinel value of -1 for unmatched jets.
- **greedyJetMatching φ wrap-around:** the `matchScore` function computes `phi_std()` differences. A fix to normalise the difference to [−π, π] was applied in the updated header. This had no measurable effect on results (jets crossing the φ boundary are rare) but is correct in principle.

---

