# PACER-DC four-context window contract (v0.1)

Status: **draft, frozen for the smoke run described at the bottom.**
Date: 2026-09-21
Owner file: `project/pacer_dc_training/extract_pacer_dc_four_context_embeddings.py`

## 1. Analysis unit

The unit of analysis is one **`candidate × replica × window`**, instantiated across
**four matched contexts**:

| context | system directory | meaning |
|---|---|---|
| `candidate_probe` | `<candidate>__candidate_probe` | M4 + candidate + ACh |
| `candidate_no_probe` | `<candidate>__candidate_no_probe` | M4 + candidate |
| `probe_only` | `probe_only` | M4 + ACh (shared control) |
| `apo` | `apo` | M4 (shared control) |

`probe_only` and `apo` are **shared** control trajectories: every candidate reuses
the *same* control trajectory and the *same* window indices, which is what makes
the differential a matched comparison.

## 2. Window definition

* `window_frames = 100` (the PACER-DC main analysis window; the 39-frame G0 smoke
  is historical validation only and **is not** this protocol).
* `start_frame = window_index * window_frames`, `end_frame = start_frame + 99`
  (inclusive). Windows are **non-overlapping**.
* A window is *materialised* only if the context's trajectory has at least
  `end_frame + 1` frames.
* All four contexts of a unit share `candidate_id`, `replica_id`, `window_id`,
  `start_frame`, `end_frame`. A context may not use a different window, a
  different frame count, or different frame sampling — even if a "better" window
  exists elsewhere in its trajectory.

## 3. Per-row metadata (long-form table)

`four_context_long.csv`, one row per (unit, context):

```
candidate_id, replica_id, window_id, start_frame, end_frame, frame_count,
context, system, trajectory_path, topology_path, trajectory_frames_total,
frame_spacing_ps, receptor_selection, preprocessing_version,
checkpoint_sha256, embedding_path, embedding_dim, finite, norm, status
```

`status ∈ {ok, missing_trajectory, insufficient_frames}`.

## 4. Differential table

`four_context_differential.csv`, one row per unit:

```
candidate_id, replica_id, window_id, start_frame, end_frame,
d_PAM_path, d_AGO_path, d_PAM_norm, d_AGO_norm, cosine_dPAM_dAGO,
n_contexts_present, completeness, missing_contexts
```

* `d_PAM = z_CA - z_A`, `d_AGO = z_C - z_0`.
* `completeness ∈ {complete, incomplete}`. **If any context is missing or short,
  the row is written with `completeness=incomplete`, `d_*_path` empty, and no
  differential is computed.** Missing contexts are named in `missing_contexts`.

## 5. Shared preprocessing (invariant across contexts)

Every context goes through exactly the same path:

1. topology = `project/results/pacer_dc_membrane_reference_v01/<system>/minimized.pdb`
   trajectory = `project/results/pacer_dc_production_v01/<system>/replica_<NN>/trajectory.dcd`
2. receptor selection = `chainID E and protein` — chain E is the M4 receptor
   (270 residues) in every PACER-DC membrane system; chains A–D are the G protein
   and nanobody, F is POPC, G/H water, I ions.
3. CHARMM → atom14 mapping identical to the G0 converter: `HSE/HSD/HSP → HIS`,
   `ILE CD → CD1`, slot names from mdgen `residue_constants`.
4. coordinates in ångström; atom14 shape `(100, L, 14, 3)` float32.
5. sequence for the encoder comes from the *same* topology selection, so the
   four contexts share one `seqres` and one residue mask.
6. `build_batch` is imported verbatim from the G0 audit package
   (`oneprot_g0_audit/scripts/run_g0_embedding.py`) — the same code that produced
   the G0 evidence, not a re-implementation.
7. one checkpoint for all contexts; its sha256 is recorded in every row.

`receptor_selection` and `preprocessing_version` are written into each row so a
reviewer can detect any row that deviated.

## 6. Baselines (Phase 4; currently blocked)

PCA / tICA / VAMP must use **the same atom14 window** as the OneProt embedding —
same frames, same residues and same window boundaries. They must also use one
frozen basis fitted on training systems only. Fitting a separate basis per
context makes the resulting coordinates incomparable.

The current extractor's historical per-context implementation is therefore
disabled by a hard error. It may not be used as a reported baseline until a
train-only shared-basis implementation and transform-only evaluation path are
added.

**Leakage rule.** Any comparative claim requires a frozen training manifest: fit
the basis on training trajectories and only transform validation/test units.
With the current molecule count no comparative claim is permitted (see §7).

## 7. Data gate

`TRAINING_GATE` is decided from the molecule-level inventory, not from the number
of windows:

* independent candidates, independent chemotypes,
* complete replicas per molecule, complete four-context windows per molecule,
* functional PAM labels, binder-but-not-PAM / agonist-like controls,
* molecules splittable into molecule-level, chemotype-aware, leak-free
  train/val/test.

Windows and frames are **not** independent samples. If any requirement is unmet,
the extractor writes `TRAINING_GATE = CLOSED` and no training may start.
