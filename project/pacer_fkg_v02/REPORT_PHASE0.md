# PACER-FKG v02 long-MD - Phase 0 report

**Gate: PASS_WITH_DOCUMENTED_CAVEAT**

Phase 0 established isolation and audited provenance, completion, trajectory metadata, PBC behavior, and the frozen encoder input mapping. No ViSNet inference or downstream PACER-FKG computation was performed.

## 1. v01 protection coverage

The pre-existing manifest was preserved byte-for-byte. It contains 41 files across the seven protected historical result roots: G1 diagnostics (6), G2A PBC (5), G2A (3), G2B robustness (5), G2B (5), G2C (5), and R2/R3 (12). End-of-phase recomputation: `V01_RESULT_HASHES_UNCHANGED = TRUE`.

## 2. Baseline tracked-tree protection

`BASELINE_TRACKED_TREE_SHA256.json` records SHA256 and size for all 630 Git-tracked files at `0a4b4f2`, branch `experiment/pacer-fkg-v02-longmd-v01`, tag `encoder-candidate-ah-frozen-20260929`. End-of-phase recomputation: `BASELINE_TRACKED_FILES_UNCHANGED = TRUE`.

## 3. Server/local checksums

The SHA256 of `SERVER_SHA256SUMS.txt` matches its companion file. All 60 server-recorded production files match locally; no mismatch occurred. The inventory contains 75 backup files total, including validation and rsync evidence not recorded by the server manifest.

## 4-5. Exact trajectory inventory and metadata

| System | Rep | Atoms | Residues | Receptor residues | Frames | First ps | Last ps | dt ps | PBC |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| apo | 01 | 227685 | 21864 | 270 | 1000 | 50.000000 | 50000.000149 | 50.000000 | yes |
| apo | 02 | 227685 | 21864 | 270 | 1000 | 50.000000 | 50000.000149 | 50.000000 | yes |
| apo | 03 | 227685 | 21864 | 270 | 1000 | 50.000000 | 50000.000149 | 50.000000 | yes |
| probe_only | 01 | 227711 | 21865 | 270 | 1000 | 50.000000 | 50000.000149 | 50.000000 | yes |
| probe_only | 02 | 227711 | 21865 | 270 | 1000 | 50.000000 | 50000.000149 | 50.000000 | yes |
| probe_only | 03 | 227711 | 21865 | 270 | 1000 | 50.000000 | 50000.000149 | 50.000000 | yes |
| compound110__candidate_no_probe | 01 | 227716 | 21865 | 270 | 1000 | 50.000000 | 50000.000149 | 50.000000 | yes |
| compound110__candidate_no_probe | 02 | 227716 | 21865 | 270 | 1000 | 50.000000 | 50000.000149 | 50.000000 | yes |
| compound110__candidate_no_probe | 03 | 227716 | 21865 | 270 | 1000 | 50.000000 | 50000.000149 | 50.000000 | yes |
| compound110__candidate_probe | 01 | 227742 | 21866 | 270 | 1000 | 50.000000 | 50000.000149 | 50.000000 | yes |
| compound110__candidate_probe | 02 | 227742 | 21866 | 270 | 1000 | 50.000000 | 50000.000149 | 50.000000 | yes |
| compound110__candidate_probe | 03 | 227742 | 21866 | 270 | 1000 | 50.000000 | 50000.000149 | 50.000000 | yes |

Each trajectory is consistent with 50 ns from the DCD header and local completion ledger. Five deterministic full-coordinate frames per trajectory are finite; all frames were scanned at receptor-CA level for PBC discontinuities.

## 6. Completion/QC

All 12 `progress.json` files independently record `complete`, target/completed 50.0 ns. Each trajectory and `state.csv` has 1,000 stored records. Paired seeds are locally recorded and verified: R1=27101, R2=38201, R3=49301.

## 7. PBC/preprocessing route

Established short-data route: DCD -> per-system minimized PDB -> `chainID E and protein` -> raw stored receptor coordinates -> `pacer_dc_atom14_v0.1` -> expected heavy-atom point cloud -> ViSNet. The route does not make molecules whole, center, align, unwrap, or wrap. DCD box vectors are retained but no coordinate transformation is applied. Across 12,000 frames, the audit found zero adjacent-CA PBC splits, zero severe adjacent-CA breaks, and zero temporal CA wrapping events, so the established route is usable unchanged.

## 8. Atom/residue compatibility

All systems select 270 receptor residues in identical order and match the authenticated short-data sequence. Expected N/CA/C/O and side-chain atom14 heavy atoms are present; Gly uses N/CA/C/O only; HSE/HSD/HSP -> HIS and CHARMM ILE CD -> CD1 rules are preserved; terminal OXT/hydrogens are ignored. All 12 trajectories satisfy the frozen input contract without repair.

## 9. Path guard

`python -m unittest discover -s project/pacer_fkg_v02 -p 'test_*.py' -v`: 4/4 tests passed. The guard uses resolved paths plus component-aware containment checks and rejects every protected historical root and all paths outside the three allowed v02 roots.

## 10. Gate decision

**PASS_WITH_DOCUMENTED_CAVEAT**

Documented non-blocking caveats:

- The long-MD backup contains no topology. Read-only, checksummed per-system sanitized PDBs from C:\projects\PACER_DC_MD_delivery_20260927 were used; their atom counts match every DCD and their receptor sequence/mapping matches the authenticated short-data atom14 CSVs.
- All state.csv timestamps are offset from DCD header times by about 3.999860 ps at the final frame (DCD: 50.0-50000.0 ps; state.csv: 54.0-50004.0 ps). Completion is independently recorded as 50.0 ns in each progress.json; no timing value was inferred from filenames.
- MDAnalysis warns that some non-receptor PDB records lack resid information and defaults those records to resid 1. The chain-E protein selection remains 270 residues, matches the short-data sequence exactly, and has complete expected atom14 heavy atoms.

## 11. Git diff/status summary

Only new v02 audit/guard/report files are present; no baseline tracked file is modified. The provenance result root is intentionally ignored by the repository's existing rules. Current concise status:

```text
## experiment/pacer-fkg-v02-longmd-v01
?? project/pacer_fkg_v02/
!! project/results/pacer_fkg_v02_longmd_v01/
```

Phase 0 stops here. No encoder inference, BS256 cache, v02 normalization, bandwidth fit, RFF draw, PACER-FKG, or biological analysis was performed.
