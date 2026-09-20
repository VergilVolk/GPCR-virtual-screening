# PACER-DC evidence table contract

Each row is one observed endpoint from one independent structure or trajectory replica.
Trajectory frames must never be entered as independent replicas.

Required columns:

| column | allowed meaning |
|---|---|
| `candidate_id` | stable molecule identifier |
| `context` | `candidate_probe`, `candidate_no_probe`, `probe_only`, or `apo` |
| `replicate_id` | independent MD replica or independent experimental structure |
| `metric` | `coupling_coordinate`, `orthosteric_pose_rmsd_A`, or `binding_compatibility` |
| `value` | finite numeric observation |
| `evidence_level` | e.g. `trajectory`, `static_structure`, `docking`; retained for audit |
| `source_id` | dataset, PDB, or simulation-run identifier |

Endpoint signs are frozen:

- `CoupledShift = mean(candidate_probe coupling) - mean(probe_only coupling)`;
- `OrthostericStabilization = mean(probe_only RMSD) - mean(candidate_probe RMSD)`;
- `IntrinsicActivationRisk = mean(candidate_no_probe coupling) - mean(apo coupling)`;
- larger `BindingCompatibility` must always mean better compatibility.

Only candidates with all four contexts and required metrics are Pareto eligible. The
program never emits a PAM probability and never treats a complete vector as experimental
confirmation.

Identical `probe_only` and `apo` simulations may use `candidate_id=SHARED_CONTROL`.
Candidate-specific control rows override a shared row with the same context and metric.
