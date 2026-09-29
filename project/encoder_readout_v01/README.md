# Standalone encoder readout Experiment A v01

This implementation reads cached Geom2Vec NPZs and the existing frozen region map.
It does not import the encoder, train anything, calculate biological contrasts,
run kernels, or alter PACER-FKG methods or results. External inputs are read-only.

## Commands (from the isolated worktree root)

Synthetic-only tests (scratch files are confined to the test directory and cleaned):

```powershell
python -B -m unittest discover -s project/encoder_readout_v01/tests -v
```

Read-only real-data preflight:

```powershell
python -B project/encoder_readout_v01/experiment_a.py --input-root C:/projects/pacer_encoder_data/R2R3_full_v01 --preflight-only
```

Separate real-data decomposition and validation command:

```powershell
python -B project/encoder_readout_v01/experiment_a.py --input-root C:/projects/pacer_encoder_data/R2R3_full_v01 --output-root project/results/encoder_readout_A_v01/run_001
```

An output run must be new and beneath `project/results/encoder_readout_A_v01/` in
this worktree. Existing outputs are never overwritten. All 40 inputs must pass
preflight before any output is created. A failure after output creation leaves
an incomplete directory without a PASS report; use a fresh run name after
resolving the failure. Do not treat partial output as validated.

## Input authentication and metadata

The authoritative NPZ manifest is
`project/results/pacer_dc_fkg_R2R3_v01/PACER_FKG_R2R3_AUDIT.json`, field
`provenance.input_files`. Its SHA256 and the archived batch-audit SHA256 are
pinned in the standalone module. The graph hash must match the manifest.
Exactly the 40 declared replica/window/context paths are required. Each NPZ is
hashed and parsed from the same in-memory bytes, both in preflight and immediately
before decomposition. Missing, hash-mismatched and nonfinite inputs fail closed.

Required data are float32 `residue_features (100,270,128)`, float32
`global_features (100,128)`, scalar Unicode `sequence` identical to the frozen
manifest and graph, and int64 `frame_ids (100,)` exactly 0 through 99. These IDs
are local to each input window, not absolute simulation times. No absolute time
or DCD offset is invented. Optional `residue_index` must be int64 0 through 269.
All numeric NPZ fields must be finite. Archived global features must equal the
historical FP32 residue mean. These strict checks describe this frozen archive,
not an interchangeable generic embedding format.

## Arithmetic and output contract

Convert H to a C-contiguous FP64 working array. Compute `G = mean(H, axis=1,
dtype=float64)` and `R = H - G[:,None,:]`. There is no robust scaling, L2
normalization, fitted parameter, clipping or final FP32 cast. Preserve all
frames, channels and residue indices. The output names and shapes per file are:

| Field | Shape / meaning |
|---|---|
| `global_mean` | FP64 `(100,128)` |
| `local_residual` | FP64 `(100,270,128)` |
| `region_local_mean` | FP64 `(100,9,128)`, mean residual over each region |
| `region_original_mean` | FP64 `(100,9,128)`, independent original-H region mean |
| `frame_ids`, `sequence` | Original arrays, unchanged |
| `residue_index`, `residue_sites` | Canonical embedding indices and frozen graph site identities |
| `region_names` | Original frozen JSON region order |
| `region_indices`, `region_offsets` | Flat indices plus offsets; no object arrays |
| source path/hash, graph hash, version | Per-output provenance |

The logical per-region two-branch readout is `[region_local_mean[:,r,:], G]`,
shape `(100,256)`. G is stored only once. R/G preserve original residue features
up to FP64 rounding; pooled regions alone do not. This cannot restore atom-level
information already lost in extraction. Original source files remain the exact
FP32 archive; they are never replaced by reconstructed values.

The frozen graph is read directly from
`project/results/pacer_dc_geom2vec_pilot_v01/M4_MULTISTRUCTURE_GRAPH_v01.json`.
Its nine region sizes are **12,16,8,14,13,13,9,12,12** in JSON order:
**109 memberships, 69 unique residues**. The earlier design's 119 was an
arithmetic error. Regions overlap; membership count is not independent evidence.

## Validation and provenance

Validate `R+G` against H, residual zero mean, and `region_local_mean+G` against
each original regional mean. A conservative numerical bound is
`128 * eps(float64) * max(1,max(abs(H)))`; this is an arithmetic acceptance bound,
not a statistical/biological gate. Record maximum absolute error, RMSE and
relative L2 error (null when the reference norm is zero), separately by file and
region. Recompute all branches independently and compare array bytes, then
reload every written NPZ and compare all arrays byte-for-byte. NPZ container
hashes are provenance, not a guarantee of identical ZIP bytes across reruns.

`validation_report.json` is written only after all files pass. It includes input
and output hashes, schema, the frozen graph nodes and region memberships, full
archived extraction audit, numerical metrics, runtime versions, source-code
hashes, Git revision/status, paths and timestamps. Historical checkpoint metadata
is explicitly attributed to the archived audit: this run does not verify weight
coverage or installed-source equivalence. Frame/time metadata absent from the
archive remains unknown.

Generated NPZs are already excluded by the repository's `*.npz` ignore rule.
Retain the small JSON report for review; do not force-add large datasets.
No attention, encoder inference, kernel computation, biological model selection,
or long-timescale MD is involved. Frames/windows are not independent biological
samples; the validation reports numerical properties only.
