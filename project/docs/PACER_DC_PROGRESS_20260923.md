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

## 2026-09-23 最新里程碑：compound110 replicas 2 and 3

### Git milestone
Latest pushed commit: 8d51112. Branch: audit/oneprot-pacer-dc-g0. Two candidate_probe QC manifests committed and pushed.

### Production and QC
compound110 + ACh replicas 2 and 3: both completed 5 ns CUDA production, 500 frames each, 227742 atoms, 10 ps interval. Trajectory integrity, preliminary thermodynamic QC, receptor CA RMSD and centroid PBC audits PASS. Mean receptor CA RMSD: replica 2 = 1.102387 A; replica 3 = 1.200738 A. Mean temperature: 300.243 K and 300.277 K.

### Cross-replica observations
Mean compound110 RMSD: 4.379 A and 3.836 A. Final initial-pocket contacts: 26 and 28. Final core-centroid distances: 4.656 A and 4.485 A. Both replicas retain core contacts and lose contacts with GLY65, SER232, ILE234 and PRO235. THR148, VAL149 and PRO150 contact gains occur in replica 2 but not replica 3. A common final binding pose or biological mechanism is NOT established.

### Four-context progress
compound110 + ACh: replica 1 = 1 ns; replicas 2 and 3 = 5 ns + preliminary QC.
compound110 only: replica 1 = 1 ns; replicas 2 and 3 not started.
ACh only: replica 1 = 1 ns; replicas 2 and 3 not started.
apo: replica 1 = 1 ns; replicas 2 and 3 = 5 ns + preliminary QC.

### Next actions
1. Start probe_only replica 2 and 3 production after checking existing files.
2. Complete compound110 candidate_no_probe replica 2 and 3.
3. Audit replica 1 production status and complete matched 5 ns windows.
4. Perform matched four-context atom14 and frozen OneProt-MD analysis.
5. Expand the functional control matrix after QC.

### Limitations
Independent replicas share a common starting structure. Five nanoseconds does not establish conformational convergence. Simulation E-chain residue indices are not original CHRM4 sequence indices. Large DCD files, checkpoints and membrane systems remain local, not on GitHub. TRAINING_GATE = CLOSED.

## 2026-09-25 Progress Update

### Completed milestones
- compound110 R2 and R3: all four contexts completed 5 ns production, preliminary QC and provenance audit.
- R3 candidate-only QC archived in commit a189942.
- Extracted and verified 40 raw trajectory windows: 2 replicas x 5 windows x 4 contexts.
- All 90 transferred raw, topology and manifest files passed WSL/Windows SHA256 checks.
- R3 apo extraction uses the corrected recovery trajectory with verified provenance.
- Completed 40 frozen OneProt-MD embeddings (1024 dimensions) and 10 matched dPAM vectors.
- Added reusable Windows embedding batch runner and an apo trajectory override to the extractor.

### dPAM stability results
- Ten complete differential units and 25 differential pair comparisons.
- Matched R2/R3 window cosine similarities: 0.502745, 0.712310, -0.353644, 0.262047, 0.085430.
- R2 mean dPAM L2: 0.223579; R3 mean dPAM L2: 0.259751.
- Mean dPAM cross-replica cosine: 0.339248; L2 distance: 0.279383.
- Cancellation ratios: R2 0.353899; R3 0.399608.
- Window-level and cross-replica directional stability have not been established.

### Current decision
- TRAINING_GATE = CLOSED. No classifier training or fine-tuning.
- Descriptive results do not establish PAM classification or biological mechanism.
- Pause further computation pending the afternoon project meeting.
- Next scientific steps will be decided after discussion.
