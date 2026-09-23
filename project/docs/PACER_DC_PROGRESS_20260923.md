# PACER-DC Progress Report
Date: 2026-09-23
Branch: audit/oneprot-pacer-dc-g0

## 1. Project objective

Evaluate reproducibility of frozen OneProt-MD trajectory
embeddings across four matched CHRM4 simulation contexts.

d_PAM = z_CA - z_A
d_AGO = z_C - z_0

These are descriptive differential features, not validated
biological activity predictions.

## 2. Completed four-context technical tests

LY2119620:
- Known PAM positive control.
- Replica 01, window 000 completed.
- Initial embedding-level differential calculated.

compound110:
- Allosteric agonist specificity control.
- Replica 01, window 000 completed.
- Four-context re-audit with extractor 498a2d2: PASS.
- Independent topology provenance verified.
- Checkpoint loading audit: PASS.
- All six old/new vectors exactly identical.

Initial descriptive differential results:

                  LY2119620    compound110
||d_PAM||2          0.630557       0.544822
||d_AGO||2          0.448729       0.744751
cosine             -0.564953       0.385785

These observations do not establish functional classification.

## 3. Apo replica 2

- Production: 5 ns complete.
- Recovery: binary checkpoint.
- 500 frames, 10 ps spacing.
- 227685 atoms.
- All trajectory coordinates finite.
- Mean receptor CA RMSD: 1.223 Angstrom.
- Mean temperature: 300.31 K.
- Preliminary structural and numerical QC: PASS.

## 4. Apo replica 3

- Production: 5 ns complete.
- Recovery from verified 50000-step binary checkpoint.
- Original DCD had 14 frames through 140 ps.
- Original XML/checkpoint state was at step 50000.
- Original files preserved and backed up.
- Recovery performed in an independent output directory.
- Recovery DCD: 490 frames.
- Final merged DCD: 500 frames, 10-5000 ps.
- Final merged CSV: 500 records, steps 5000-2500000.
- All merged coordinates match their source frames exactly.
- Maximum source-versus-merged coordinate difference: 0.
- Recovery-boundary adjacent CA RMSD: 0.566 Angstrom.
- Mean receptor CA RMSD: 1.181 Angstrom.
- Mean recovery-stage temperature: 300.32 K.
- Preliminary structural and numerical QC: PASS.

Canonical merged trajectory:
project/results/pacer_dc_recovery_merged_v01/apo/replica_03/trajectory_corrected.dcd

The recovery-stage thermodynamic statistics exclude
the original first 10 records.

## 5. Current limitations

- Five-nanosecond pilots do not prove conformational convergence.
- Replica 2 and 3 receptor RMSDs were computed relative
  to their respective first frames.
- The apo replica 3 recovery boundary passed CA RMSD
  continuity checks, not exhaustive PBC or per-atom checks.
- Remaining compound110 contexts require additional replicas.
- Four-molecule functional comparison is incomplete.
- LY2119620 re-audit compatibility remains to be confirmed.
- Training gate remains CLOSED.

## 6. Next execution steps

1. Freeze and commit small apo replica 2/3 QC manifests.
2. Confirm Windows compound110 re-audit assets.
3. Inventory remaining compound110 production start states.
4. Complete matched four-context, three-replica 5 ns pilots.
5. Audit cross-replica and cross-window differential stability.
6. Extend to additional predefined functional controls.

Large DCD files, checkpoints, membrane systems,
and temporary recovery files must not be committed.

Windows is the sole Git authority.
