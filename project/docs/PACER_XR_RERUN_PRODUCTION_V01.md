# PACER-XR rerun production contract

This implementation uses only the byte-locked Module-2 rerun top200 from
`origin/codex/drugclip-pacer-handoff` (commit `b1f5787b95138a5f56df84878e1458b793a819dd`).
The roster SHA256 is `41c5a69e6c54825d0cbe9149a2f1fe9dbe5ea4ae2de137c161450607c3fa28eb`.
Canonical molecular identity is calculated with RDKit. Original `gen_id`, `smiles`,
`score`, `rankpct`, and `final_rank` fields remain present in assembled outputs.
The normalized fields are `candidate_id`, `canonical_smiles`, and `rerun_final_rank`.

## Execution modes

The master script requires one of `--preflight`, `--acceptance`, or `--run`.
`--preflight` verifies the full roster, anchor structures, asset/tool versions,
and frozen reference distributions without docking. `--acceptance` runs only
rerun ranks 1, 2, and 3, then stops after six-channel raw assembly. It cannot
invoke fusion or create a natural/final shortlist. `--run` is the full production
mode and requires separate user authorization.

The master stages are Glide, Vina, six-channel assembly, reference verification,
PACER-XR, cascade, full200, natural top20, and final experimental shortlist20.
References are also checked during preflight so an unrecoverable reference
dependency cannot waste a complete docking campaign.

## Frozen assets and software

`project/config/pacer_xr_assets_v01.json` locks the 11 real, successful Glide
grid archives, prepared/source receptors, frozen receptor PDBQT files, PMF,
smoke manifests, and executable identities. Normal production never regenerates
grids or minimizes receptors. The lock retains the ensemble smoke's ACH bond
and formal-charge repair provenance. Centers remain unchanged; Glide's inner
boxes preserve the original centroid search regions and outer boxes accommodate
ligand extent (7TRS 22/44 Å; ensemble 30/60 Å).

Expected software: Schrödinger 2023-1 build 128, Glide 98128, Vina 1.2.7,
Molscrub 0.1.1, RDKit 2025.3.5, and Meeko 0.7.1. The validated local environment
is `C:\Users\barne\anaconda3\envs\chrm4_vs\python.exe`.

## Preparation and raw scores

Glide uses Molscrub pH7 with at most 16 states, the existing enumeration order,
RDKit ETKDGv3 seed 42+state index and MMFF (300 iterations), and LigPrep
`-i 0 -nt -ns -bff 14`. Every returned variant is retained and assigned an
explicit parent/state/variant ID. Loss of an entire source state during LigPrep
fails the preparation stage. Rank 1 reuses its byte-exact validated smoke ligand,
which was made with this same preparation chain. Each candidate's prepared input
is docked against all 11 grids in SP with five poses per ligand and five
postdocking poses. Reduction is the lowest finite `r_i_glide_gscore`, then
prepared variant ID, then poseviewer entry order. Emodel is never a raw channel.
Returned poses need not cover every prepared variant; returned/accepted variant
counts are recorded, and a candidate/receptor with no finite pose fails.

Vina PDB preserves the actual `dock_candidate_portfolio.py` contract: parent
canonical SMILES, RDKit/Meeko preparation, seed42, exhaustiveness4, one mode,
one CPU, and the local 7TRS 22 Å box. The historical PDB implementation does
not enumerate pH7 states; this is explicitly distinguished from the ensemble.
Vina ensemble reuses Molscrub pH7 state enumeration, RDKit/Meeko preparation,
seed42, exhaustiveness8, nine modes, energy range3 kcal/mol,
one CPU, and the frozen 30 Å boxes. Each accepted state must finish successfully.
The rank-1 Vina affinity is parsed for each state; the lowest finite affinity
across states is the candidate/cluster score, with state ID resolving ties.

For either engine, `BE_i = cluster_score + PMF_i`. BEmin is the minimum;
BEavg is the arithmetic mean of exactly ten clusters. No weighting or imputation
is permitted. The four historical overlaps are rerun because a complete
protocol/state-preparation fingerprint for reuse was not established.

## References and fusion

`project/config/pacer_xr_references_v01.json` locks six exact published score
tables from `tylerdt1/gpcr-am-ensemble-docking` commit
`44798c841ee77230b1f89fe41074e41b458c7070`, the commit already frozen in
`THOMPSON_MIAO_2026_EXACT_REPRODUCTION.md`. Recovery validates upstream Git blob
identity and records SHA256. Tables live under
`project/results/pacer_xr_rerun_v01/references/{PDB,Ensemble}`. Production checks
hashes and expected finite-array sizes before loading. Candidate scores never
serve as reference distributions.

Fusion directly imports the existing `empirical_rank`: for sorted reference
energies, percentile = `1 - (less + 0.5*equal)/N`. Higher is better. PACER_XR
is the equal mean of all six ranks. Stage1 uses the inclusive Glide-BEmin raw
threshold `quantile(reference, 0.01, method="higher")` and orders by Glide-BEmin
percentile descending. Stage2 orders remaining candidates by PACER_XR descending.
Stable ties preserve frozen rerun order. Both energy engines supply all200
PDB evaluations and all200×10 ensemble evaluations before production fusion.

Natural top20 is cascade ranks1..20. The final experimental shortlist inserts
missing structural anchors (PACER0010→PACERGEN01755, PACER0073→PACERGEN00094,
PACER0027→PACERGEN00123) and removes the worst cascade-ranked non-anchor after
each insertion. It retains exactly20 unique canonical structures without
changing the full200 ranking. Selection reason, legacy ID, and anchor-required
fields distinguish natural selections from forced anchors.

## Resume and failures

Each preparation/docking job has an input/code/tool fingerprint, atomic current
status, immutable numbered attempts, output SHA256s, and a JSONL ledger.
Only matching PASS records with intact outputs resume. Failed or changed jobs
receive new attempts; old inputs and logs remain available. Missing states,
poses, clusters, molecules, or channels prevent channel/fusion output. The
assembler verifies the complete expected job-identity set and all job hashes.
The master resumes valid completed engine stages and verifies assembly/fusion
outputs through its own stage ledger. Outputs from acceptance are isolated under
`acceptance_rank1_3`; they are not promoted into full200 scientific outputs.

## Outputs and tests

Each engine writes ligand-state manifests, detailed pose/state and cluster
tables, three raw channels, completion summaries, and channel provenance. The
master's production outputs are `candidate_six_channel_raw_scores.csv`,
`pacer_xr_rerun_full200_ranking.csv`, `pacer_xr_natural_top20.csv`,
`pacer_xr_final_shortlist20.csv`, and `fusion_audit.json` under
`project/results/pacer_xr_rerun_v01`.

Contract tests: `project/tests/test_pacer_xr_production_v01.py`. They use
synthetic energies in temporary fixtures and cannot produce a real shortlist.
The separate rank1–3 acceptance exercises real installed docking tools.
Existing Vina analysis, candidate scoring algorithms, and historical smoke
outputs remain unchanged. No commit/push is part of this implementation task.
