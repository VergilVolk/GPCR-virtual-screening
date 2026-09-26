# PACER-DC G0–G4 Freeze and H Handoff v01

Date: 2026-09-26
Branch: audit/oneprot-pacer-dc-g0
Last archived G3 commit: c820188
Training gate: CLOSED

## 1. Scope and freeze decision

G0–G4 constitute the completed short-MD audit and
representation-mechanism diagnostic phase.

No further exploratory G analyses, region reselection,
model modification or classifier training are authorized
under this frozen analysis plan.

Independent cloud MD validation is assigned to the
new H task family.

This document records the analytical freeze.
Git archival is complete only after commit and
remote synchronization have been verified.

## 2. Frozen experimental design

Target: compound110.

Four MD contexts:
- CA: compound110 + ACh
- A: ACh only
- C: compound110 only
- 0: apo

Contrasts:
- dPAM = z_CA - z_A
- dAGO = z_C - z_0
- dINT = dPAM - dAGO

Existing short-MD dataset:
- Independent replicas: R2 and R3
- Windows: W0–W4
- Duration: 1 ns per window
- Frames: 100 per window
- Saving interval: 10 ps

Adjacent windows from the same trajectory are
not independent experimental replicas.

The contrasts are operational representation-space
definitions, not established pharmacological effects.

## 3. Completed G tasks

G0:
Data provenance, checkpoint and preprocessing audit.

G1:
Five-layer representation extraction and QC.

G2:
Frozen residue mapping and regional L0/L1 analysis.

G3:
Four-context replica-bias decomposition,
cancellation analysis and structural diagnostics.

G4:
Exact FinalLayer reconstruction, four-stage
mechanism decomposition and frozen-region comparison.

The original G2 region map remains authoritative:
G2_REGION_MAP_v02.json.

## 4. G4 frozen analysis design

Four regions:
- ACh_pocket: 12 residues
- compound110_pocket: 12 residues
- distal_control: 12 residues
- global: 270 residues

Four representation stages:
- S0: original L0, 384 dimensions
- S1: token-wise LayerNorm, 384 dimensions
- S2: common conditional scaling, 384 dimensions
- S3: final L1 projection, 21 dimensions

The same checkpoint and t=0 condition are
used for all existing contexts and replicas.

The common modulation shift and final linear bias
cancel in the four-context contrasts.

Stage differences are sequential descriptive
effects, not independent causal contributions.

## 5. G4 numerical validation

All 40 L0/L1 input pairs passed reconstruction QC.

Formal output:
G4_four_stage_regions_v01.csv

Formal result rows: 80
Historical G2 contrast checks: 120

Maximum G2 matching error:
3.8510479560827093e-07

Maximum L1 reconstruction error:
2.0034512285271866e-07

An independent CSV inspection confirmed:
- 80 unique records
- 20 records per region
- 20 records per stage
- 16 records per window
- Reconstruction errors populated only for S3

All reported numerical QC checks passed.

## 6. Frozen scientific results

Five-window mean dINT cross-replica cosine:

Region                 L0         L1
ACh_pocket             0.304992   0.548038
compound110_pocket     0.128167  -0.069925
distal_control         0.129616  -0.042170
global                 0.308376   0.010757

Mean sequential stage changes:

Region                LayerNorm   Scaling   Projection
ACh_pocket            +0.007610  +0.132472  +0.102963
compound110_pocket    +0.029241  -0.150276  -0.077057
distal_control        +0.010466  -0.202948  +0.020695
global                +0.012491  -0.116105  -0.194004

The ACh pocket had positive L1 dINT cosine
in all five existing matched windows.

FinalLayer effects depended on region and window.
LayerNorm had relatively small effects on
the five-window mean cosine.

Conditional scaling and projection accounted
for larger sequential changes.

At W3, ACh pocket dINT cosine changed from
-0.293925 at S0 to +0.628047 at S3.

At W0, however, ACh pocket cosine decreased
from +0.817856 at S0 to +0.512402 at S3.

These observations establish an exact computational
description for the existing short-MD representations.

They do not establish PAM specificity or a
validated receptor conformational mechanism.

## 7. Remaining scientific limitations

Only two independent short-MD replicas were analyzed.

Five contiguous windows cannot be treated as five
independent replicates.

Direction agreement does not establish
magnitude convergence.

Replica-bias cancellation partly explains some
of the observed L1 consistency.

The previously analyzed C-alpha geometry does
not yet provide robust independent support for
the representation-space mechanism.

The contribution of insufficient MD sampling
remains unresolved.

## 8. G4 source files

The six G4 scripts are:

- g4_l1_reconstruction_smoke_v01.py
- g4_l1_matrix_qc_v01.py
- g4_stage_w0_smoke_v01.py
- g4_stage_window_smoke_v02.py
- g4_stage_region_window_v01.py
- g4_stage_summary_v01.py

The authoritative formal G4 CSV is:

project/results/pacer_dc_four_context_v01/
compound110/G4_four_stage_regions_v01.csv

The smoke-test and intermediate scripts are
retained for reproducibility.

The formal summary script and CSV are the
primary G4 computational artifacts.

## 9. H task definition

H is the independent cloud-MD validation phase.

Proposed workstreams:

H0: Predefine validation metrics and verify
    incoming trajectory provenance and QC.

H1: Assess sampling quality, temporal
    correlations and conformational states.

H2: Test independent representation stability.

H3: Reproduce cross-layer mechanism diagnostics.

H4: Test predefined structure–representation
    relationships.

H5: Integrate evidence and scientific conclusions.

The detailed H validation protocol has not
yet been frozen.

It must be finalized before inspecting new
cloud-MD representation or mechanistic results.

The frozen G model, contrasts, residue mapping
and existing scientific claims remain unchanged.

## 10. Training and archival boundaries

TRAINING_GATE = CLOSED.

Do not initiate classifier training merely
because additional cloud sampling is available.

Do not retroactively alter G region definitions,
metrics or historical result files.

G0–G4 archival and H validation must remain
distinguishable in scripts, outputs and Git history.
