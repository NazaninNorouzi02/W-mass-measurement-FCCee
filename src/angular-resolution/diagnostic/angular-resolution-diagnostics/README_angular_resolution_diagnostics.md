# Angular Resolution Diagnostic Analysis — FCC-ee WW at √s = 160 GeV
## IDEA Detector / Delphes / FCCAnalyses / Winter2023 Campaign

---

## Overview

This repository contains the Step 1 analysis script for studying jet angular and energy resolution in hadronic WW events at the FCC-ee. The physics goal is to extract the resolution parameters (σ_θ, σ_φ, σ_α, σ_x) needed to build the covariance matrix for a kinematic fit to measure the W boson mass.

The analysis evolved through multiple diagnostic iterations to understand a systematic non-Gaussian behaviour in the energy resolution variable α = (E_reco − E_parton) / E_parton for the softer jets (j3, j4). This README documents every version of the code, every diagnostic test run, all ROOT browser commands used, and all conclusions drawn.

---

## Physics Context

- **Process:** e+e− → W+W− → qqqq (fully hadronic) at √s = 160 GeV  
- **Detector:** IDEA concept, simulated with Delphes (winter2023 campaign)  
- **Generator:** Pythia 8 via FCCAnalyses  
- **Jet algorithm:** Exclusive ee-kt clustering forced to N = 4 jets  
- **Reference thesis:** W mass measurement methodology from FCC-ee thesis (papas detector simulation, CLD concept) — this work reproduces the methodology for the IDEA/Delphes simulation  

Jets are ordered by energy: j1 is hardest, j4 is softest. Each event contains two W bosons, each decaying to two quarks (partons). The four reconstructed jets are matched to four partons via a greedy ΔR-based matching algorithm.

---

## Repository Structure

```
.
├── angular_resolution_ecm160_diagnostics.py   # Main analysis script (current)
├── headers/
│   ├── greedyJetMatching.h          # Reco-to-gen jet matching (greedy ΔR)
│   ├── getDeltaTheta.h              # θ resolution: reco − gen
│   ├── getDeltaPhi.h                # φ resolution: reco − gen
│   ├── getDeltaEta.h                # η resolution: reco − gen
│   ├── getDeltaMass.h               # mass resolution
│   ├── getXGen.h                    # x = 2E_gen / √s (gen level)
│   ├── getXReco.h                   # x = 2E_reco / √s (reco level)
│   ├── getElement.h                 # safe RVec element accessor
│   ├── jetPartonMatching.h          # reco jet → parton matching by ΔR
│   ├── getDeltaAlphaParton.h        # ORIGINAL alpha header (DO NOT USE — broken sentinel)
│   ├── getDeltaAlphaPartonFixed.h   # FIXED alpha header — always pushes -999 sentinel
│   ├── filterValues.h               # RVec range filter
│   └── selectQuarks.h              # selects status=23 quarks from MC truth
└── README_angular_resolution_diagnostics.md   # This file
```

**Critical note on headers:** `getDeltaAlphaParton.h` (original) silently skips `push_back` for unmatched jets, producing a short RVec. When `getElement` reads past the end, it returns 0 or causes undefined behaviour, producing a fake spike at α = 0 and α = −1000 in plots. **Always use `getDeltaAlphaPartonFixed.h`**, which always pushes −999 as a sentinel. Filter with `delta_alpha_j4 > -998` in plots.

---


The script processes a fraction `1e-6` of the full `p8_ee_WW_ecm160` dataset (4 chunks, 4 CPUs). To run on more statistics, increase `fractions`.

---

## What the Analysis Script Does — Step by Step

### 1. Lepton/Photon Removal
Isolated electrons, muons, and photons are removed from `ReconstructedParticles` to produce `reco_clean` — the input to jet clustering.

```python
reco_clean = ReconstructedParticles - photons - electrons - muons
```

### 2. Parton Selection
`selectQuarks(Particle)` extracts the four hard-scatter quarks (status = 23) from the MC truth record. Events with exactly 4 partons are kept.

### 3. Generator-Level Jet Clustering
Final-state MC particles (status = 1) are clustered with `ee-kt` exclusive N = 4 to produce `jets_gen4`.

### 4. Reco-Level Jet Clustering
`reco_clean` particles are clustered with the same algorithm to produce `jets_reco4`.

### 5. y34 Diagnostic Variables
`y34_gen` and `y34_reco` are the Durham jet-resolution parameter at the 3 → 4 jet transition. Larger y34 means the 4th jet is more robustly separated from the 3-jet topology.

### 6. Jet Matching
`greedyJetMatching(jets_reco4, jets_gen4, 0.1)` matches each reco jet to a gen jet by minimising ΔR, with a 0.1 rad maximum. Events where all 4 jets are matched are kept.

### 7. Resolution Variables
For each matched reco–gen pair:
- `delta_theta_ji` = θ_reco − θ_gen  
- `delta_phi_ji` = φ_reco − φ_gen  
- `delta_eta_ji` = η_reco − η_gen  
- `delta_x_ji` = x_reco − x_gen, where x = 2E/√s  
- `delta_alpha_ji` = (E_reco − E_parton) / E_parton  **[main alpha — thesis definition]**  
- `delta_alpha_gen_ji` = (E_reco − E_gen) / E_gen  **[diagnostic: reco vs gen-jet]**  
- `delta_alpha_gen_vs_parton_ji` = (E_gen − E_parton) / E_parton  **[diagnostic: gen-jet vs parton]**  
- `delta_E_parton_ji` = E_reco − E_parton  **[absolute energy difference]**  

### 8. Charged/Neutral Energy Decomposition (per jet, cone ΔR < 0.4)
For each jet j1..j4:
- `E_charged_ji`, `E_neutral_ji` — energy from charged/neutral reco_clean particles  
- `n_charged_ji`, `n_neutral_ji` — multiplicities  
- `E_charged_over_gen_ji`, `E_neutral_over_gen_ji` — normalised to matched gen-jet energy  
- `E_charged_over_reco_ji`, `E_neutral_over_reco_ji` — normalised to reco-jet energy  
- `neutral_fraction_clean_ji` — E_neutral / (E_charged + E_neutral)  

### 9. Jet Geometry / Overlap Variables
- `dR_reco_j4j3`, `dR_reco_j4j1` — ΔR between j4 and its two nearest neighbours  
- `min_dR_reco_j4_others` — minimum ΔR of j4 with any of j1, j2, j3  
- `dtheta_reco_j4j3`, `dphi_reco_j4j3` — angular separations between j4 and j3  

### 10. Scalar Angle Branches
`theta_reco_ji`, `phi_reco_ji`, `eta_reco_ji`, `theta_gen_ji`, `phi_gen_ji`, `eta_gen_ji` — per-jet scalar values extracted from the vector branches for convenient ROOT plotting.

### 11. Event-Level Energy Diagnostics
- `E_reco_total_4j` = sum of all four reco-jet energies  
- `E_gen_total_4j` = sum of all four matched gen-jet energies  
- `E_ratio_total_4j` = E_reco_total / E_gen_total  

---

## All ROOT TBrowser Commands Used in the Investigation

These commands were run directly in the ROOT interactive session on the output TTree (`t` or `events`). Copy-paste them to reproduce any plot.

### Basic alpha distributions (fixed sentinel) 1-4
```cpp
t->Draw("delta_alpha_j4",     "delta_alpha_j4 > -998")
t->Draw("delta_alpha_j3",     "delta_alpha_j3 > -998")
t->Draw("delta_alpha_gen_j4", "delta_alpha_gen_j4 > -998")
t->Draw("delta_alpha_gen_j3", "delta_alpha_gen_j3 > -998")
```

### Option 1 — reco vs gen-jet energy resolution 5-8
```cpp
t->Draw("delta_alpha_gen_j1")
t->Draw("delta_alpha_gen_j2")
t->Draw("delta_alpha_gen_j3")
t->Draw("delta_alpha_gen_j4")
```

### Option 4 — gen-jet vs parton (pure MC, no detector) 9-11
```cpp
t->Draw("delta_alpha_gen_vs_parton_j4", "delta_alpha_gen_vs_parton_j4 > -998")
t->Draw("delta_alpha_gen_vs_parton_j4", "delta_alpha_gen_j4 < -0.3 && delta_alpha_gen_vs_parton_j4 > -998")
t->Draw("delta_alpha_gen_vs_parton_j4", "delta_alpha_gen_j4 > -0.1 && delta_alpha_gen_vs_parton_j4 > -998")
```

### Option 5 — absolute energy difference 12
```cpp
t->Draw("delta_E_parton_j4", "delta_E_parton_j4 > -998")
```

### Bimodal diagnostic — perfect vs imperfect reco populations 13-14
```cpp
t->Draw("delta_alpha_gen_j4", "is_perfect_reco_j4 == 1")
t->Draw("delta_alpha_gen_j4", "is_perfect_reco_j4 == 0")
```

### Particle multiplicity in j4 — good vs bad events 15-18
```cpp
t->Draw("n_reco_particles_j4", "delta_alpha_gen_j4 > -0.1")
t->Draw("n_reco_particles_j4", "delta_alpha_gen_j4 < -0.3")
t->Draw("n_reco_particles_j3")
t->Draw("n_reco_particles_j4")
```

### Angular resolution vs energy resolution (2D) 19
```cpp
t->Draw("delta_alpha_gen_j4:delta_theta_j4", "", "colz")
```

### Total 4-jet energy ratio — good vs bad events 20-21
```cpp
t->Draw("(E_reco_j1+E_reco_j2+E_reco_j3+E_reco_j4)/(E_gen_j1+E_gen_j2+E_gen_j3+E_gen_j4)", "delta_alpha_gen_j4 > -0.1")
t->Draw("(E_reco_j1+E_reco_j2+E_reco_j3+E_gen_j4)/(E_gen_j1+E_gen_j2+E_gen_j3+E_gen_j4)", "delta_alpha_gen_j4 < -0.3")
```

### J3 energy response — does j3 gain energy when j4 loses it? 22-23
```cpp
t->Draw("E_reco_j3/E_gen_j3", "delta_alpha_gen_j4 > -0.1")
t->Draw("E_reco_j3/E_gen_j3", "delta_alpha_gen_j4 < -0.3")
```

### J4 energy response — reco/gen ratio 24-25
```cpp
t->Draw("E_reco_j4/E_gen_j4", "delta_alpha_gen_j4 > -0.1")
t->Draw("E_reco_j4/E_gen_j4", "delta_alpha_gen_j4 < -0.3")
```

### Gen-jet energy — good vs bad events (is j4 harder or softer in bad events?) 26-27
```cpp
t->Draw("E_gen_j4", "delta_alpha_gen_j4 > -0.1")
t->Draw("E_gen_j4", "delta_alpha_gen_j4 < -0.3")
```

### Alpha_gen vs E_gen_j4 (2D — does energy loss scale with jet energy?) 28-29
```cpp
t->Draw("delta_alpha_gen_j4:E_gen_j4", "", "colz")
t->Draw("(E_reco_j4-E_gen_j4)/E_gen_j4:E_gen_j4", "delta_alpha_gen_j4 < -0.3", "colz")
```

### Gen-vs-parton alpha for bad events 30-31
```cpp
t->Draw("delta_alpha_gen_vs_parton_j4", "delta_alpha_gen_j4 < -0.3")
t->Draw("delta_alpha_gen_vs_parton_j4", "delta_alpha_gen_j4 < -0.3 && delta_alpha_gen_vs_parton_j4 > -998")
```

### Neutral/charged energy fractions — key diagnostic 32-36
```cpp
t->Draw("E_neutral_over_gen_j4", "delta_alpha_gen_j4 > -0.1")
t->Draw("E_neutral_over_gen_j4", "delta_alpha_gen_j4 < -0.3")
t->Draw("E_charged_over_gen_j4", "delta_alpha_gen_j4 > -0.1")
t->Draw("E_charged_over_gen_j4", "delta_alpha_gen_j4 < -0.3")
t->Draw("E_neutral_over_gen_j4:E_charged_over_gen_j4", "delta_alpha_gen_j4 < -0.3 && E_neutral_over_gen_j4 > -900", "colz")
```

### Detector geometry — is j4 in a bad detector region? 37-39
```cpp
t->Draw("theta_reco_all[3]", "delta_alpha_gen_j4 > -0.1")
t->Draw("theta_reco_all[3]", "delta_alpha_gen_j4 < -0.3")
t->Draw("delta_alpha_gen_j4:theta_reco_all[3]", "", "colz")
```

### Jet separation — is j4 close to j3? 40-43
```cpp
t->Draw("theta_reco_all[3]-theta_reco_all[2]", "delta_alpha_gen_j4 < -0.3")
t->Draw("phi_reco_all[3]-phi_reco_all[2]",     "delta_alpha_gen_j4 < -0.3")
t->Draw("min_dR_reco_j4_others", "delta_alpha_gen_j4 < -0.3")
t->Draw("min_dR_reco_j4_others", "delta_alpha_gen_j4 > -0.1")
```

### y34 — is the 4th jet marginal? 44-46
```cpp
t->Draw("y34_reco", "delta_alpha_gen_j4 < -0.3")
t->Draw("y34_gen",  "delta_alpha_gen_j4 > -0.1")
t->Draw("y34_gen",  "delta_alpha_gen_j4 < -0.3")
```

### n_reco_particles neutral/charged 2D 47
```cpp
t->Draw("n_neutral_j4:n_charged_j4", "delta_alpha_gen_j4 < -0.3", "colz")
```

---

## Complete Summary of Diagnostic Results

### What was proven step by step

| Step | Variable | Good (α_gen > −0.1) | Bad (α_gen < −0.3) | Conclusion |
|------|----------|---------------------|---------------------|------------|
| 1 | `delta_alpha_gen_vs_parton_j4` | mean ≈ −0.037, narrow | same | Gen-jet captures parton energy well. Problem is NOT parton shower. |
| 2 | `delta_alpha_gen_j4` | mean ≈ 0 by definition | mean ≈ −0.45 | Reco-jet loses ~45% of gen-jet energy. Problem IS in reco. |
| 3 | `E_reco_j4/E_gen_j4` | peak at 0.98 | ramp 0.05→0.70 | Hard upper bound at 0.70 for bad events. Not a smooth Gaussian. |
| 4 | Total 4j energy ratio | mean 0.890 | mean 0.789 | Total reco energy also low. Energy is lost, not migrated to j3. |
| 5 | `E_reco_j3/E_gen_j3` | mean 0.849 | mean 0.817 | J3 is not boosted when j4 is bad. No inter-jet migration. |
| 6 | `E_gen_j4` good vs bad | mean 29.2 GeV | mean 37.6 GeV | Bad jets are HARDER, not softer. Rules out soft-jet threshold. |
| 7 | `n_reco_particles_j4` | mean 6.25 | mean 7.26 | Bad events have MORE particles. Missing particles is not the cause. |
| 8 | `delta_alpha_gen_j4 : delta_theta_j4` | — | no correlation | Jet direction is correct even when energy is wrong. Tracking is fine. |
| 9 | `theta_reco_all[3]` good vs bad | flat 0 to π | flat 0 to π, identical | No geometric / detector-region preference. Crack hypothesis ruled out. |
| 10 | `phi/theta separation j4−j3` | — | back-to-back (Δφ ≈ ±π) | Jets are well-separated. Overlap hypothesis ruled out. |
| 11 | `y34_gen` good vs bad | mean 771 | mean 995 | Bad events have HIGHER y34 — more robust 4th jet, not a marginal one. |
| 12 | `E_neutral/charged_over_gen_j4` 2D | — | both suppressed simultaneously | Both sectors lose energy at ~same fractional level (~0.3–0.5 of gen). |
| 13 | `alpha_gen : E_gen_j4` (2D) | — | flat loss across all energies | Loss is not threshold-dependent. Same ~40% loss from 10 to 65 GeV. |

### All hypotheses tested and their status

| Hypothesis | Status | Key evidence |
|---|---|---|
| Parton shower / out-of-cone radiation as reference issue | **RULED OUT** | delta_alpha_gen_vs_parton_j4 mean = −0.037, narrow |
| Soft jet threshold / HCAL energy minimum | **RULED OUT** | Bad events have harder gen jets (37.6 vs 29.2 GeV) and more particles |
| Missing particles from detector thresholds | **RULED OUT** | Bad events have MORE reco particles (7.26 vs 6.25) |
| Jet clustering boundary / particle migration between jets | **RULED OUT** | No angular correlation; energy lost but direction preserved |
| Barrel-endcap crack / geometric detector acceptance | **RULED OUT** | theta_reco_all[3] distributions identical for good and bad events |
| Jet-jet calorimeter overlap / PF double-subtraction | **RULED OUT** | J4 and J3 are back-to-back, well-separated; j3 not boosted |
| Marginal 4th jet (low y34, barely separated topology) | **RULED OUT** | Bad events have HIGHER y34 than good events |
| Energy migration from j4 to j3 | **RULED OUT** | Total 4-jet energy is also lower in bad events; j3 not over-reconstructed |
| Energy loss scaling with jet energy (threshold behaviour) | **RULED OUT** | Alpha_gen vs E_gen_j4 shows flat ~40% loss at all energies |
| Parton matching failure | **RULED OUT** | Proven in earlier sessions; n_matched_jets == 4 filter applied |

### What has NOT been ruled out (remaining open)

The only mechanism consistent with all 13 observations simultaneously is: **a global energy scale underestimation in the Delphes particle-flow algorithm for the specific particle composition of j4 in these events, affecting both charged and neutral sectors proportionally.** The 2D plot of E_neutral_over_gen vs E_charged_over_gen for bad events shows both suppressed to 30–50% of gen simultaneously, which is the signature of either:

1. A systematic calibration bias in the Delphes IDEA card for jets of this topology at this energy scale, or  
2. A specific population of events where the WW decay topology places j4's particles in a configuration where the PF algorithm underestimates the energy of multiple particle types simultaneously.

---

## Next Tests to Run

### Priority 1 — neutral fraction vs energy response 48-50
```cpp
t->Draw("E_reco_j4/E_gen_j4:neutral_fraction_clean_j4", "delta_alpha_gen_j4 < -0.3 && neutral_fraction_clean_j4 > -998", "colz")
t->Draw("neutral_fraction_clean_j4", "delta_alpha_gen_j4 > -0.1 && neutral_fraction_clean_j4 > -998")
t->Draw("neutral_fraction_clean_j4", "delta_alpha_gen_j4 < -0.3 && neutral_fraction_clean_j4 > -998")
```
**Expected result if Delphes HCAL is the cause:** bad events will have systematically higher neutral fraction. The HCAL response in Delphes is parametric with larger downward fluctuations than the tracker, so neutral-heavy jets lose more energy.

### Priority 2 — check j1 and j2 for the same pattern 51-52
```cpp
t->Draw("E_reco_j1/E_gen_j1", "delta_alpha_gen_j1 < -0.3")
t->Draw("E_reco_j2/E_gen_j2", "delta_alpha_gen_j2 < -0.3")
```
If j1 and j2 show the same hard cutoff at ~0.70, the effect is not specific to the soft jet — it is a global event-level energy loss for a specific class of WW events.

### Priority 3 — event-level total energy vs alpha_gen of all jets 53
```cpp
t->Draw("E_ratio_total_4j:delta_alpha_gen_j4", "", "colz")
```
If the total 4-jet energy is the same regardless of which individual jet has bad alpha_gen, the problem is a different jet in different events. If the total energy is always low whenever any jet has bad alpha_gen, there is a common event-level cause.

### Priority 4 — parton energy ordering vs jet energy ordering
Check whether the "bad" j4 events are events where parton ordering and jet ordering disagree (i.e. the parton assigned to j4 position is not actually the softest parton in the event). This would indicate the `jetPartonMatching` assigns j4 to a parton that is not its natural partner.

---
## Notes on ROOT File Contents

The output ROOT file contains a TTree named `events` with the following branch groups:

- `delta_theta_j1..j4`, `delta_phi_j1..j4`, `delta_eta_j1..j4`, `delta_x_j1..j4` — angular/x resolution per jet  
- `delta_alpha_j1..j4` — main alpha (reco vs parton), uses fixed sentinel  
- `delta_alpha_gen_j1..j4` — reco vs gen-jet energy resolution  
- `delta_alpha_gen_vs_parton_j1..j4` — gen-jet vs parton (pure MC)  
- `delta_E_parton_j1..j4` — absolute energy difference reco − parton  
- `is_perfect_reco_j1..j4` — flag: 1 if |E_reco − E_gen| / E_gen < 0.001  
- `E_reco_j1..j4`, `E_gen_j1..j4` — jet energies  
- `theta/phi/eta_reco_j1..j4`, `theta/phi/eta_gen_j1..j4` — jet angles  
- `E_charged_j1..j4`, `E_neutral_j1..j4` — cone energies by particle type  
- `E_charged_over_gen_j1..j4`, `E_neutral_over_gen_j1..j4` — normalised fractions  
- `neutral_fraction_clean_j1..j4` — E_neutral / (E_charged + E_neutral) in cone  
- `n_charged_j1..j4`, `n_neutral_j1..j4`, `n_reco_particles_j1..j4` — multiplicities  
- `dR_reco_j4j3`, `dR_reco_j4j1`, `min_dR_reco_j4_others` — jet-pair geometry  
- `E_reco_total_4j`, `E_gen_total_4j`, `E_ratio_total_4j` — event-level energy  
- `y34_gen`, `y34_reco` — Durham jet-resolution parameter at 3→4 transition  
- `filtered_delta_*` — range-filtered versions of resolution variables  
- `parton_energies`, `parton_eta`, `parton_phi`, `parton_y` — truth parton kinematics  

---

*Last updated: 15 June 2026. Analysis ongoing.*
