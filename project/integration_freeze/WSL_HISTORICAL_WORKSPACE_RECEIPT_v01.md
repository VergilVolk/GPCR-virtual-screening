# WSL historical execution workspace - read-only receipt v01

Recorded: 2026-10-02 (integration preflight, step A)
Status: **PROTECTED AS FOUND. Nothing in this workspace was reverted, overwritten, moved or deleted.**

## 1. Why this receipt exists

Correction 1 of the integration preflight reclassifies
`/root/work/GPCR-virtual-screening` from "the MD workspace to be cleaned up" to
**historical execution workspace**. Its dirty code state, its Git state and its
unique result assets are evidence, not cruft. The final production MD will use a
clean copy re-synchronised from the code-freeze worktree onto a native Linux
filesystem; it will **not** reuse this tree.

## 2. Identity

| Field | Value |
|---|---|
| Windows path | `\\\\wsl.localhost\\Ubuntu-24.04\\root\\work\\GPCR-virtual-screening` |
| Distro | `Ubuntu-24.04` (WSL 2.7.12.0, kernel 6.18.33.2-microsoft-standard-WSL2) |
| Branch | `main` |
| HEAD | `5ef9904cc3ca688fac7e4521473d54f3cb09c345` |
| HEAD date / subject | Sun Sep 20 17:51:46 2026 +0800 - "Use molecule-level validation for PACER-DC early stopping" |
| Remote | `origin https://github.com/VergilVolk/GPCR-virtual-screening.git` |
| Relationship to the Windows repo | **Same repository.** `5ef9904` is an ancestor of every current branch, so this is an old checkout, not a fork. |
| Stash | empty |
| Working-tree size | 24,266,102,850 bytes (~22.6 GiB) |
| `project/results` | **514 files, 24,250,250,806 bytes** |
| Tracked files with a diff | **327** |
| Untracked entries | **6** |

## 3. Git state

- `git status --porcelain=v1` -> 333 lines (327 ` M`, 6 `??`), captured verbatim in
  `WSL_GIT_STATUS_PORCELAIN_v01.txt`.
- `git diff --numstat` -> captured in `WSL_GIT_DIFF_NUMSTAT_v01.txt`.

### 3.1 Line-ending noise vs real content

Of the 327 modified tracked files, **322 differ only by line endings / file mode**
(the checkout is mode 0777 with CRLF, the Git blobs are LF).

**Real content differences - 5 files** (`git diff --ignore-cr-at-eol`):

| File | +/- |
|---|---|
| `project/environment_pacer_dc_md.yml` | +15 |
| `project/scripts/build_pacer_dc_membrane_reference.py` | +64/-43 |
| `project/scripts/build_pacer_dc_openmm_reference_systems.py` | +49/-43 |
| `project/scripts/pacer_dc_score.py` | +2/-2 |
| `project/tests/test_pacer_dc_score.py` | +69 |

The full patch is preserved in `WSL_DIRTY_CODE_DIFF_v01.patch`, and the tracked-blob
vs working-tree blob hashes are in `WSL_DIRTY_FILES_BLOB_HASHES_v01.tsv`.

**Note on the two MD builders:** their working-tree copies have the ligand
`definitions` block removed (leaving an empty dict), i.e. they are abandoned
in-progress edits that exist in no commit. They are **preserved here as found**
and are **not** the code that produced the historical membrane references. The
committed versions are canonical.

### 3.2 Untracked entries (also preserved)

```
?? project/docs/PACER_DC_PROGRESS_20260923.md
?? project/pacer_dc_training/audit_pacer_dc_training_gate.py
?? project/pacer_dc_training/extract_pacer_dc_four_context_embeddings.py
?? project/pacer_dc_training/extract_pacer_dc_four_context_embeddings.py.pre_498a2d2.bak
?? project/pacer_dc_training/extract_pacer_dc_four_context_embeddings.py.pre_trajectory_override.bak
?? project/scripts/analyze_pacer_dc_trajectory_qc.py
```

## 4. Unique scientific assets that exist ONLY here

| Directory | Files | Bytes | Note |
|---|---:|---:|---|
| `pacer_dc_production_v01` | 70 | 13,338,540,946 | 6 systems x 3 replicas, **5 ns pilot** generation (the 50 ns set is on Windows in `PACER_DC_MD_backup`) |
| `pacer_dc_recovery_merged_v01` | 4 | 2,732,366,151 | apo recovery lineage |
| `pacer_dc_restraint_release_v01` | 84 | 1,786,517,339 | **the only complete 6-system set**; the Windows main worktree holds CM00734 only |
| `pacer_dc_recovery_v01` | 5 | 1,410,508,104 | apo recovery |
| `pacer_dc_four_context_v01` | 131 | 1,255,229,035 | LY2119620 / apo_only / compound110 + `training_gate.json` |
| `pacer_dc_short_equilibration_v01` | 54 | 1,018,462,037 | |
| `pacer_dc_membrane_reference_v01` | 22 | 991,556,263 | **the only complete 6-system set**, plus `shared_OPM_POPC_water.pdb` |
| `pacer_dc_membrane_smoke_v01` | 54 | 735,208,789 | |
| `pacer_dc_membrane_reference_compound110_test` | 8 | 343,903,645 | |
| `pacer_dc_membrane_reference_compound110_formalbase_test` | 8 | 323,708,956 | |
| `pacer_dc_openmm_reference_qc_v01` | 20 | 144,453,196 | |
| `pacer_dc_recovery_backup_v01` | 5 | 109,862,448 | |
| `pacer_dc_openmm_reference_qc_compound110_test` | 6 | 48,210,912 | |
| `pacer_dc_unified_qc_v01` | 14 | 6,049,106 | |
| `pacer_dc_reference_complexes_v01` | 7 | 4,175,647 | |
| `pacer_dc_qc_v01` | 14 | 70,391 | |
| `pacer_dc_common_protein_v01` | 2 | 1,417,069 | `7TRS_common_protein_pH74.pdb` + audit |
| `pacer_dc_phase1_audit_v2` | 2 | 5,777 | |
| `pacer_dc_score_smoke` | 2 | 1,637 | |
| `pacer_dc_ligand_parameterization_audit.json` | 1 | 4,096 | |
| `pacer_dc_md_environment_audit.json` | 1 | 4,096 | `production_ready: true`, `blocking_reason: "none"` |

Also under `/root/work` (outside the repo):

| File | Bytes | SHA256 |
|---|---:|---|
| `gpcr_dcd_20260921_0703.tar.gz` | 2,557,996,597 | `6d162039ca87a2f550b0754b07348e1fc887ff35cc0326ad0bd7e5ba718dc436` |
| `gpcr_migration_20260921_0703.tar.gz` | 731,661,248 | `b23deeb95f4cb22d24534f7f89095b015ac587b8e6a4e1009dc902812f43c326` |
| `pacer-dc-md_export.yml` | 11,897 | `2b54378df67371f2e15a6b0c1e7947bc38c018682982147879c121ddc97aada0` |
| `pacer-dc-md_explicit.txt` | 28,399 | `7480e06fbf689255986e3f68e7a788291a981b7d559f821312e05b99159185a6` |
| `env_backup_environment_pacer_dc_md.yml` | 1,127 | `0167dae634d31d14f00c26c3e64fec9da0e9cb41736315bd5e6be20bac8e840b` |
| `env_backup_pacer-dc-md_pre_cuda_fix.yml` | 11,897 | `3c85a8dff2304756d300e6862a07587acb3b7c9564a7a82a0d2f09e9ec637510` |
| `pacer_dc_three_contexts_progress.txt` | 422 | `167778936ae3d8ab872be48ffe3667e1d8169e0aaa3eaf4baf35a17eee61bb71` |

### 4.1 Production completeness caveat inside `pacer_dc_production_v01`

| System | R1 | R2 | R3 |
|---|---|---|---|
| `apo` | 5.0 ns complete | 5.0 ns complete | **`status: running`, 0.1 ns - INCOMPLETE** |
| `probe_only` | 1.0 ns | 5.0 ns | 5.0 ns |
| `compound110__candidate_no_probe` | 1.0 ns | 5.0 ns | 5.0 ns |
| `compound110__candidate_probe` | 1.0 ns | 5.0 ns | 5.0 ns |
| `LY2119620__candidate_no_probe` | 1.0 ns | - | - |
| `LY2119620__candidate_probe` | 1.0 ns | - | - |

## 5. Environment captured alongside

`pacer-dc-md` (verified independently in the same preflight):

```
python 3.11.16 | OpenMM 8.6.1 | platforms: Reference, CPU, CUDA
openmmforcefields 0.14.1 | openff.toolkit 0.16.5 | openff.interchange 0.3.29
pdbfixer present | ParmEd 4.3.1 | MDAnalysis 2.10.0 | mdtraj 1.11.1
RDKit 2026.03.1 | numpy 1.26.4 | scipy 1.17.1 | pandas 2.3.3
tleap / sander present ; GROMACS and pmemd absent
nvidia-smi: NVIDIA GeForce RTX 4070 Laptop GPU, 610.60, 8188 MiB, 0 MiB used
```

Other WSL environments: `base` (/root/miniforge3), `pacer_dc_prep`.
**`drugclip-cpu` does not exist yet.**

## 6. Integrity evidence

- `WSL_RESULTS_SHA256_v01.txt` - SHA256 of **all 514 files** under `project/results`.
- `WSL_RESULTS_FILES_v01.tsv` - path, bytes, mtime for the same 514 files.
- `WSL_TOPFILES_v01.tsv` - repo-root files.
- `WSL_WORK_EXTRAS_SHA256_v01.txt` - the `/root/work` extras above.
- `WSL_RESULTS_TOTAL_BYTES_v01.txt` - `24250250806` and `514`.

## 7. Handling rules for future phases

1. **Do not revert, clean, or overwrite anything in this workspace.** It is a
   historical evidence record.
2. Production MD must run from **a clean copy of the code-freeze worktree**
   re-synchronised onto a native Linux filesystem (per the frozen
   `setup_wsl_pacer_dc_md.sh` contract, which explicitly rejects 9p/drvfs for
   trajectory I/O).
3. If any asset listed in section 4 is needed downstream, copy it out; never
   operate in place.
