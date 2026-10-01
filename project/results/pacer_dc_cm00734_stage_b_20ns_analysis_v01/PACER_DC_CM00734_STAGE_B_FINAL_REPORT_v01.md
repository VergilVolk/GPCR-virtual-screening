# PACER-DC CM00734 Stage B 20 ns Frozen PACER-FKG Evaluation

## 1. Purpose

Stage B is a hard-negative specificity test of the frozen PACER-FKG v02 workflow.

Compound:
- ID: `CM00734`
- experimental role: strict experimental inactive hard negative
- tested contexts:
  - `CM00734__candidate_no_probe`
  - `CM00734__candidate_probe`

This stage is not used for calibration, threshold selection, encoder tuning, bandwidth selection, RFF fitting, graph fitting, region fitting, or any other outcome-driven optimization.

## 2. Matched MD design

Frozen Stage-B production design:

- receptor: M4 / 7TRS
- 2 candidate contexts
- 3 paired replicas per context
- 20 ns per trajectory
- 50 ps per stored frame
- 400 stored frames per trajectory
- total production sampling: 120 ns
- total stored frames: 2400
- paired seeds:
  - R1 = 27101
  - R2 = 38201
  - R3 = 49301

All six trajectories completed 20.0 ns and passed production QC.

## 3. Production integrity

Local production outputs were frozen with SHA256 after transfer.

Representative trajectory SHA256 values:

### CM00734 candidate no probe
- R1 trajectory: `f60bf18d4a2c060808353b67288d74a6694a5dcff29f4bc8cfc78970f89eecf4`
- R2 trajectory: `cfbc8e1deb6c745eed00f846d555b5f75b617bb094914dbd7eed4a36348570cf`
- R3 trajectory: `29f82c14235c81c164bd6d51eb229e81183a17b6002df856297ff782a44b6df7`

### CM00734 candidate + probe
- R1 trajectory: `ac53213cb02e1988f567590e09f839e1d7d5150f6a21916aed13b154a367760e`
- R2 trajectory: `a0e6a80b66ce7552713b41d2e298c889842a679ba3cf60027f8756cc189a4651`
- R3 trajectory: `f0f533ed91750f88d2bf14327855ccf514519cef7212f7139f82b8881a503b77`

The complete manifest is stored as:

`project/results/pacer_dc_cm00734_stage_b_20ns_v01/PRODUCTION_SHA256SUMS.txt`

## 4. Frozen PACER-FKG anchor

Historical PACER-FKG v02 numerical freeze:

`b48bc74a757a3d1421acb5c4bc0544ce5ca590e3959b971dad698e32146e11bd`

The historical Phase-3 verifier returned:

- `PHASE3_VERIFIED`
- 12 authenticated Phase-1 pairs
- 600 historical trajectory blocks
- 200 calibration blocks
- fitting replica: R2 only
- excluded fitting replicas: R1 and R3
- 41 protected v01 files verified

No frozen numerical state was modified for Stage B.

## 5. Stage-B frozen application workflow

The Stage-B workflow mirrors the previously frozen 20 ns closure analysis.

### Phase 1: frozen BS256 extraction

Six Stage-B trajectories were encoded using the frozen PACER-FKG v02 encoder.

Output per trajectory:

`[400, 270, 256]`

Phase-1 status:

- `PHASE1_FROZEN`
- `PHASE1_VERIFIED`
- cache count = 6

Phase-1 freeze receipt SHA256:

`f8b5b945c81ea51dd4c5fc7d8ee0ee04114cde07043fed877b25da162ecdc8db`

### Phase 2A: frozen numerical application

Each trajectory was divided into 20 contiguous 1 ns blocks:

- 20 frames per block
- 50 ps per frame
- 20 blocks per trajectory

Frozen branches:

- `STATE_MOTION`
- `SIGNED_DRIFT`

Frozen normalization, bandwidths, RFF weights, RFF biases, graph, and region definitions were reused without refitting.

Phase-2A status:

- `PHASE2A_FROZEN`
- `PHASE2A_VERIFIED`
- outputs = 12

Phase-2A freeze receipt SHA256:

`73f1e7e5f4ef3532daf7d797a54e4d60a22cf3c65a0d16477b3871a1a7a815b8`

### Phase 2B: matched graph/region evaluation

Matched four-context panels use the same historical first-20-ns controls:

- A = apo
- P = probe_only
- C = candidate_no_probe
- CP = candidate_probe

Contrasts remain exactly:

- `Delta_PAM = CP - P`
- `Delta_AGO = C - A`
- `Delta_INT = CP - P - C + A`

Stage-B result status:

- `STAGE_B_MATCHED_20NS_COMPLETE`
- `STAGE_B_MATCHED_20NS_VERIFIED`

Final result SHA256:

`be663c7d45cc4ade1e053029e6aa09e65ca24c1e2fc22dae5398a26bf49d7d99`

## 6. Primary specificity result

The main pre-existing Stage-A signal was the `STATE_MOTION / Delta_INT / compound110_extension` representation.

Matched 20 ns values:

| Panel | R1 magnitude | R2 magnitude | R3 magnitude | R1/R3 direction cosine |
|---|---:|---:|---:|---:|
| compound110 reference | 0.140523 | 0.084188 | 0.105948 | +0.583656 |
| LY2119620 | 0.104348 | 0.085594 | 0.104993 | +0.634304 |
| CM00734 | 0.067244 | 0.068286 | 0.055084 | -0.326839 |

CM00734 magnitude relative to LY2119620:

- R1: 0.644
- R2: 0.798
- R3: 0.525

The hard negative therefore does not reproduce the positive cross-replica directional agreement seen for LY2119620 and compound110 in the principal Stage-A signal.

## 7. Region-level directional behavior

For `STATE_MOTION / Delta_INT`, CM00734 has negative R1/R3 direction cosine in every reported region.

Selected comparisons:

| Region | CM00734 | LY2119620 | compound110 |
|---|---:|---:|---:|
| compound110_extension | -0.327 | +0.634 | +0.584 |
| cooperativity_mutagenesis | -0.095 | +0.235 | -0.246 |
| intracellular_microswitches | -0.006 | +0.180 | -0.111 |
| orthosteric_activation_core | -0.385 | +0.321 | -0.074 |
| orthosteric_contact_union | -0.376 | +0.300 | -0.052 |
| pam_contact_consensus | -0.251 | +0.375 | -0.308 |
| pam_contact_union | -0.241 | +0.306 | -0.183 |

The strongest evidence for hard-negative discrimination is therefore directional reproducibility, not raw magnitude.

## 8. Important limitation: magnitude is not specific

CM00734 still generates substantial `Delta_INT` magnitudes in multiple regions.

Examples in `STATE_MOTION`:

- cooperativity_mutagenesis:
  - CM/LY = 1.148, 0.984, 1.320
- intracellular_microswitches:
  - CM/LY = 1.065, 1.177, 1.033
- orthosteric_activation_core:
  - CM/LY = 1.255, 1.089, 0.972

In `SIGNED_DRIFT`, CM00734, LY2119620, and compound110 magnitudes are frequently of the same order, while R1/R3 direction agreement is generally weak.

Therefore:

**Delta_INT magnitude alone is not a specific PAM/cooperativity marker in the current frozen representation.**

## 9. Scientific interpretation

Stage B provides **partial specificity support** for frozen PACER-FKG v02.

The experimental inactive hard negative CM00734 produces substantial representation-space Delta_INT magnitude, so the method does not support a simple magnitude-based PAM classifier.

However, CM00734 fails to reproduce the positive R1/R3 directional agreement observed for LY2119620 in the principal `STATE_MOTION / Delta_INT / compound110_extension` signal and in several additional mechanism-related regions.

The most defensible interpretation is therefore:

> Frozen PACER-FKG v02 shows partial hard-negative specificity primarily through cross-replica consistency of Delta_INT direction in the STATE_MOTION branch, rather than through Delta_INT magnitude alone.

This is a representation-space mechanistic result. It is not an efficacy classifier, potency predictor, independent-replicate significance test, or proof of pharmacological PAM activity.

## 10. Outcome

Stage-B outcome:

`PARTIAL_SPECIFICITY_SUPPORT`

This outcome was accepted without any post-hoc tuning.

No changes were made to:

- encoder
- normalization
- bandwidth
- RFF state
- graph
- region definitions
- contrast formulas
- magnitude thresholds
- classifier thresholds

The hard-negative result is retained unchanged as part of the final PACER-DC scientific audit trail.
