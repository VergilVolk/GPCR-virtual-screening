# PACER-DC matched-context pilot MD production contract

Date frozen: 2026-09-01

## Objective

This pilot asks whether a candidate behaves like a probe-dependent M4 PAM rather than
merely occupying the allosteric pocket or activating M4 by itself. It is not a potency
prediction experiment and cannot replace a functional assay.

## Frozen systems

The run ledger is `project/config/pacer_dc_pilot_md_manifest.csv`:

- positive-control PAM: LY2119620;
- specificity control: compound-110 allosteric agonist;
- prospective hypotheses: PACER0076 and PACER0057;
- matched shared controls: ACh-only and apo M4;
- three paired seeds per condition;
- candidate + ACh and candidate without ACh for every allosteric ligand.

This is 30 initial 100 ns trajectories (3.0 microseconds total). A condition is extended
to 500 ns only when its endpoint confidence interval or state effective-sample-size fails
the frozen quality gate; extensions are not chosen by whether the result looks favorable.

## Common setup

1. Use the same 7TRS receptor, G-protein, membrane composition, protonation policy,
   water/ion model, box geometry, force-field family and equilibration schedule in every
   context.
2. Use paired velocity seeds across `candidate_probe`, `candidate_no_probe`, `probe_only`
   and `apo` to reduce nuisance variance while retaining three independent seed groups.
3. Transfer crystallographic control poses by receptor alignment; prospective candidates
   use the already frozen 7TRS allosteric-pocket pose proposal. Pose source is recorded in
   the manifest.
4. Do not delete unstable replicas. Setup failures, ligand escape and topology failures
   remain in the ledger and are reported as outcomes/QC failures.

## Frame-level extraction, replica-level inference

Every trajectory produces time series for:

- frozen internal-distance `coupling_coordinate`;
- sticky-soft prototype posterior and transitions (interpretation only);
- ACh pose RMSD and native-contact persistence when ACh is present;
- candidate pocket-contact persistence, heavy-atom RMSD and ligand strain proxy;
- receptor/G-protein stability and membrane/box QC.

Frames estimate a trajectory endpoint; they are never treated as independent samples.
PACER-DC bootstrap and comparisons operate on the three replica endpoints.

## Frozen candidate evidence vector

The scoring program is `project/scripts/pacer_dc_score.py`. It reports:

1. `CoupledShift = candidate+ACh coupling - ACh-only coupling`;
2. `OrthostericStabilization = ACh-only RMSD - candidate+ACh RMSD`;
3. `IntrinsicActivationRisk = candidate-without-ACh coupling - apo coupling`;
4. `BindingCompatibility` in both ligand contexts;
5. replica-level uncertainty and missing-context audit.

Higher values are desirable for the first two and binding compatibility; lower intrinsic
activation risk is desirable. Complete candidates are compared by Pareto dominance, not
an optimized weighted sum. The program never outputs a PAM probability.

The frozen coordinate artifact is
`project/config/pacer_dc_coupling_coordinate_v01.json`. Raw trajectories are converted
to frame traces and replica-level evidence with
`project/scripts/extract_pacer_dc_trajectory_endpoints.py`; the scorer accepts one or
multiple resulting evidence CSV files. A public MK-97 run-1 smoke test reproduced the
previous coupling-coordinate mean (`0.2508562`) from 500 frames while independently
extracting ACh pose RMSD and allosteric-pocket contact coverage.

## Baselines and falsification

The pilot has failed to add value if any of these hold:

- compound-110 is not separated from LY2119620 by the no-ACh intrinsic-risk endpoint;
- PACER-DC does not outperform static coupling, single-structure docking, ensemble docking,
  PCA/tICA and the non-triplet continuous coordinate on control discrimination;
- conclusions depend on one seed or disappear after initial-transient removal;
- candidate ranking changes under reasonable block sizes or one-replica deletion;
- a prospective candidate leaves the pocket in most replicas.

Passing this pilot supports prioritization for assay, not the statement that a generated
molecule is an effective PAM. Confirmation requires at minimum M4 ACh EC20 modulation,
candidate-alone intrinsic agonism, ACh EC80/probe-dependence and M2 counter-screening.

## Current host readiness

The Windows host has OpenMM CPU/OpenCL, PDBFixer, MDAnalysis and the trajectory extraction
stack. It does not currently have OpenFF Toolkit or AmberTools/tleap, so reproducible
candidate-ligand parameterization is not yet available. The exact environment is frozen in
`project/environment_pacer_dc_md.yml`; `project/scripts/audit_pacer_dc_md_environment.py`
must report `production_ready=true` before any ledger row changes to `running`.
