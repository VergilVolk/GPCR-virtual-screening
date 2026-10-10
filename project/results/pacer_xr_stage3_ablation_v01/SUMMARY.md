# PACER-XR Stage3 retrospective ablation (read-only)

- Frozen full200 SHA256: `6fff7e34d94466ff288b52281b1c3a64fef4ea6cd9da84ca45bc703429a04559`
- **Frozen Cascade reproduced exactly before all ablations.**
- Cascade Top20: PACERGEN02432, PACERGEN00123, PACERGEN00462, PACERGEN02015, PACERGEN02420, PACERGEN01803, PACERGEN02355, PACERGEN02164, PACERGEN02464, PACERGEN01495, PACERGEN01580, PACERGEN01593, PACERGEN01682, PACERGEN01954, PACERGEN01494, PACERGEN01409, PACERGEN01876, PACERGEN01579, PACERGEN02247, PACERGEN02496
- Pure six-channel XR Top3: PACERGEN01495, PACERGEN01580, PACERGEN01593
- This is **rank stability / structural evidence**, not functional PAM ground-truth performance.
- Anchor-forced final_shortlist20.csv was NOT used.

## Scenario comparisons versus frozen Cascade

| Scenario | Spearman | Kendall | Top10 shared | Top20 shared | Mean absolute rank movement |
|---|---:|---:|---:|---:|---:|
| cascade_v01 | 1.000 | 1.000 | 10/10 | 20/20 | 0.0 |
| pure_six_mean | 0.971 | 0.957 | 3/10 | 16/20 | 4.0 |
| glide_only | 0.777 | 0.583 | 2/10 | 11/20 | 29.5 |
| vina_only | 0.788 | 0.608 | 2/10 | 11/20 | 26.9 |
| pdb_only | 0.786 | 0.618 | 2/10 | 13/20 | 25.8 |
| ensemble_only | 0.896 | 0.732 | 5/10 | 11/20 | 19.3 |
| bemin_only | 0.846 | 0.664 | 4/10 | 9/20 | 23.6 |
| beavg_only | 0.881 | 0.707 | 3/10 | 11/20 | 20.3 |
| loo_Glide_PDB | 0.923 | 0.778 | 4/10 | 14/20 | 16.3 |
| loo_Glide_BEmin | 0.927 | 0.815 | 2/10 | 16/20 | 12.9 |
| loo_Glide_BEavg | 0.961 | 0.900 | 3/10 | 16/20 | 7.2 |
| loo_Vina_PDB | 0.940 | 0.818 | 3/10 | 12/20 | 13.1 |
| loo_Vina_BEmin | 0.959 | 0.866 | 3/10 | 16/20 | 9.5 |
| loo_Vina_BEavg | 0.963 | 0.890 | 3/10 | 16/20 | 8.0 |
| cascade_gate_0.5pct | 0.986 | 0.978 | 5/10 | 18/20 | 2.1 |
| cascade_gate_2pct | 0.976 | 0.948 | 9/10 | 11/20 | 4.9 |
| cascade_gate_5pct | 0.898 | 0.820 | 9/10 | 11/20 | 15.4 |

## Selected molecule channel-dependence diagnostics

| Candidate | Cascade | Pure XR | Min channel rank | Glide-Vina mean gap | Ensemble-PDB gap | Diagnostic flags |
|---|---:|---:|---:|---:|---:|---|
| PACERGEN02432 | 1 | 10 | 0.541 | -0.060 | +0.357 | channel_range_ge_0.40;PDB_ensemble_mean_gap_ge_0.25 |
| PACERGEN00123 | 2 | 100 | 0.255 | +0.112 | +0.624 | channel_below_0.50;channel_range_ge_0.40;PDB_ensemble_mean_gap_ge_0.25;LOO_pure_rank_shift_ge_20 |
| PACERGEN00462 | 3 | 104 | 0.231 | +0.077 | +0.684 | channel_below_0.50;channel_range_ge_0.40;PDB_ensemble_mean_gap_ge_0.25;LOO_pure_rank_shift_ge_20 |
| PACERGEN02015 | 4 | 122 | 0.200 | +0.217 | +0.564 | channel_below_0.50;channel_range_ge_0.40;Glide_Vina_mean_gap_ge_0.20;PDB_ensemble_mean_gap_ge_0.25;LOO_pure_rank_shift_ge_20 |
| PACERGEN02420 | 5 | 11 | 0.412 | +0.149 | +0.361 | channel_below_0.50;channel_range_ge_0.40;PDB_ensemble_mean_gap_ge_0.25 |
| PACERGEN01803 | 6 | 7 | 0.542 | +0.138 | +0.260 | channel_range_ge_0.40;PDB_ensemble_mean_gap_ge_0.25 |
| PACERGEN02355 | 7 | 66 | 0.339 | +0.071 | +0.577 | channel_below_0.50;channel_range_ge_0.40;PDB_ensemble_mean_gap_ge_0.25 |
| PACERGEN02164 | 8 | 15 | 0.452 | -0.094 | +0.407 | channel_below_0.50;channel_range_ge_0.40;PDB_ensemble_mean_gap_ge_0.25 |
| PACERGEN02464 | 9 | 14 | 0.529 | -0.039 | +0.383 | channel_range_ge_0.40;PDB_ensemble_mean_gap_ge_0.25 |
| PACERGEN01495 | 10 | 1 | 0.712 | -0.068 | +0.192 | none |
| PACERGEN01580 | 11 | 2 | 0.799 | -0.044 | +0.144 | none |
| PACERGEN01593 | 12 | 3 | 0.726 | +0.034 | +0.210 | none |
| PACERGEN01682 | 13 | 4 | 0.617 | +0.101 | +0.206 | none |
| PACERGEN01954 | 14 | 5 | 0.712 | +0.009 | +0.202 | none |
| PACERGEN01494 | 15 | 6 | 0.699 | -0.111 | +0.143 | none |
| PACERGEN01409 | 16 | 8 | 0.528 | +0.134 | +0.241 | channel_range_ge_0.40 |
| PACERGEN01876 | 17 | 9 | 0.738 | +0.019 | +0.181 | none |
| PACERGEN01579 | 18 | 12 | 0.601 | +0.014 | +0.325 | PDB_ensemble_mean_gap_ge_0.25 |
| PACERGEN02247 | 19 | 13 | 0.490 | +0.095 | +0.338 | channel_below_0.50;channel_range_ge_0.40;PDB_ensemble_mean_gap_ge_0.25 |
| PACERGEN02496 | 20 | 16 | 0.366 | +0.189 | +0.347 | channel_below_0.50;channel_range_ge_0.40;PDB_ensemble_mean_gap_ge_0.25;LOO_pure_rank_shift_ge_20 |

## Notes
- PDB and ensemble have separately calibrated reference score distributions; compare channel percentiles rather than raw energies.
- BEmin represents one PMF-corrected winning ensemble; BEavg represents a ten-conformation statistic without a unique pose.
- Contact residues are geometric heavy-atom proximity at <=4.0 A only; they do not establish hydrogen bonds or modulation mechanism.
- Ensemble receptor residue numbering may differ from PDB; do not assume residue-number correspondence without mapping.
- Changing gate/weights AFTER seeing ranking is exploratory: no experimentally labeled enrichment or prospective validation was performed.
- `pose_availability.csv` records missing pose files; no fabricated contact rows.
- Glide contact analysis requires the separate Schrodinger helper, which reads existing MAEGZ without docking.

## Output tables
- `selected_candidate_diagnostics.csv`: per-candidate six-channel diagnostics and leave-one-out rank shifts
- `six_channel_spearman.csv`: six-channel correlation matrix
- `within200_six_channel_ranks.csv`: all six channel ranks within the current 200 molecules (different from official reference percentiles)
- `ablation_summary.csv`: all ranking/threshold scenarios
- `ablation_all_ranks.csv`: 200 candidates under every scenario
- `pose_availability.csv`, `vina_pose_contacts.csv`, `glide_pose_requests.json`: local pose evidence
- `glide_pose_contacts.csv` and `glide_pose_status.csv`: available after optional Schrodinger helper
