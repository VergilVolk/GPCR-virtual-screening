# ENCODER_OPTIMIZATION_DESIGN_v01 — Read-only architecture audit

## 1. Executive summary

**The smallest supported experiment is a representation-only numerical precision comparison using existing per-frame residue embeddings, with pretrained weights, pooling, kernels, and downstream statistics unchanged.** It requires no retraining or re-embedding, but the required R2/R3 embedding files are absent from this worktree.

Workspace verification:

| Item | Verified value |
|---|---|
| Working directory | `C:\projects\GPCR-virtual-screening-encoder` |
| Branch | `experiment/pacer-encoder-opt-v01` |
| Initial working-tree status | Clean |
| Final working-tree status | Clean; no diff |
| Actions performed | Read-only source/document inspection, file inventory, checkpoint hashing |
| Actions not performed | File writes, training, extraction, MD analysis, experiments, original-worktree or cloud access |

The implemented Geom2Vec path is:

**Receptor coordinates → atom14 → heavy-atom point cloud → frozen ViSNet outputs → scalar channels concatenated with vector norms → atom mean per residue → per-frame residue features.**

The downstream paths differ materially:

- **Historical linear pilot:** time mean, then region mean.
- **Original PACER-FKG:** joint-region RBF statistics and graph propagation of nonnegative node magnitudes. Its magnitude gate is withdrawn.
- **G2-B:** shared preprocessing across replicas, but a separate shared kernel for each region; signed kernel-mean contrasts.
- **G2-C:** one shared residue-descriptor kernel; signed node contrasts, optional graph propagation, then region averaging.
- **U2 erratum v02:** shared-node, bag-of-residue contrast norms and directions; no replacement gate.

The G2-B, G2-C, robustness, and U2 erratum reports are present. Their historical input-verification claims could be inspected but **could not be independently repeated**, because the underlying 40 R2/R3 NPZ files are missing.

## 2. Verified current architecture

### Coordinate and encoder interface

The production preparation script selects `chainID E and protein`, retaining receptor coordinates and excluding direct ligand, solvent, lipid, ion, and other protein-chain inputs. Ligand conditions therefore enter Geom2Vec through receptor geometry, not explicit ligand features.

Preparation uses 100-frame windows by default, CHARMM histidine-name normalization, and `ILE CD → CD1`. The atom14 representation is float32, with coordinates documented in ångströms. See [preparation defaults and mappings](/C:/projects/GPCR-virtual-screening-encoder/project/pacer_dc_training/extract_pacer_dc_four_context_embeddings.py:50) and [window contract](/C:/projects/GPCR-virtual-screening-encoder/project/pacer_dc_training/FOUR_CONTEXT_WINDOW_CONTRACT.md:64).

The local extraction wrapper imports `geom2vec.create_model_from_checkpoint`; it does not implement ViSNet itself. The environment pins Geom2Vec to commit:

`371d642ec1061664f16e49fcac702d07fc8d0b51`

Other numerical dependencies mostly use minimum versions or remain unpinned. See [environment definition](/C:/projects/GPCR-virtual-screening-encoder/project/environment_pacer_dc_geom2vec.yml:5).

### Checkpoint

Present file:

`project/tools/geom2vec_checkpoints/visnet_l6_h64_rbf64_r75.pth`

Its independently computed SHA256 matches the extraction audit:

`b8f1ef9b591c57f7687566bd60d3664800956ac280cf2ae8e64f54196cd8d417`

The filename suggests six layers, 64 hidden channels, 64 radial basis functions, and a 7.5-radius configuration. **These are naming implications, not independently verified instantiated model settings.** The package source and loader implementation are absent from this worktree; the checkpoint was hashed, not deserialized.

The wrapper has no intermediate-layer selection option. It consumes the returned `x, v, _` from one model forward pass. See [encoder invocation](/C:/projects/GPCR-virtual-screening-encoder/project/pacer_dc_training/extract_geom2vec_atom14.py:138).

### Representation readout

For each atom, the wrapper concatenates scalar features and vector-channel Euclidean norms, then averages over atoms belonging to each residue:

\[
h_r=\frac{1}{|A_r|}\sum_{a\in A_r}[x_a,\|v_a\|_{\mathrm{xyz}}].
\]

This is **mean of vector norms**, not norm of mean vectors. See [invariant residue readout](/C:/projects/GPCR-virtual-screening-encoder/project/pacer_dc_training/extract_geom2vec_atom14.py:128).

The archived extraction audit records:

- 40 files: two replicas × five windows × four contexts;
- 100 frames per file;
- `(100, 270, 128)` residue features;
- CUDA, batch size 1, stride 1;
- compiled neighborhood backend, with fallback disabled;
- approximately 289 seconds total extraction time.

These are archived execution records, not measurements reproduced here. See [batch audit](/C:/projects/GPCR-virtual-screening-encoder/project/results/pacer_dc_geom2vec_R2R3_full_v01/batch_audit.json:1).

### Historical OneProt branch

OneProt remains a separate implementation, not an additional stage inside Geom2Vec:

- `build_batch` constructs relative rigid-frame offsets and torsions, producing `(1,T,L,21)`.
- The extractor instantiates `TrajectoryEncoder`, loads `network.md.*` weights, checks missing/unexpected keys, and produces a 1024-dimensional window embedding.
- Separate diagnostic hooks target L0 `(1,100,270,384)`, L1 `(1,100,270,21)`, pooling L2 `(1,21)`, projection L3 `(1,1024)`, and normalization L4 `(1,1024)`.

References: [OneProt preprocessing](/C:/projects/GPCR-virtual-screening-encoder/project/pacer_dc_training/oneprot_g0_audit/scripts/run_g0_embedding.py:46), [loading and extraction](/C:/projects/GPCR-virtual-screening-encoder/project/pacer_dc_training/extract_pacer_dc_four_context_embeddings.py:331), [diagnostic hooks and expected shapes](/C:/projects/GPCR-virtual-screening-encoder/project/pacer_dc_training/diagnose_oneprot_five_layers.py:101).

These hooks establish a historical OneProt layer-audit mechanism. They do not establish an equivalent Geom2Vec mechanism.

## 3. End-to-end data flow with file/line references

Notation: `T` frames, `L` residues, `N` selected heavy atoms, `K` region residues, `D` residue channels, `M` RFF channels.

### 3.1 Trajectory → raw receptor window

**Implementation:** `stage_traj`.

- Input: topology plus production DCD.
- Default window: `start = window_id × 100`, inclusive end `start + 99`; explicit start override exists.
- Output: `coords (T, selected_atoms, 3)` plus atom names, residue names, original residue IDs, and atom-to-residue indices.
- Frame order follows the trajectory slice.
- Residue order follows the selected topology residues.
- No coordinate alignment, centering, unwrapping, or reimaging is performed here.
- The sanitized topology copy removes CONECT records for parser compatibility.
- Frame spacing metadata is assigned the constant `10 ps`, rather than derived from trajectory timestamps.

**Information omitted:** non-receptor coordinates, unit-cell data in the raw NPZ, and direct ligand identity/geometry.

References: [metadata construction](/C:/projects/GPCR-virtual-screening-encoder/project/pacer_dc_training/extract_pacer_dc_four_context_embeddings.py:152), [trajectory extraction](/C:/projects/GPCR-virtual-screening-encoder/project/pacer_dc_training/extract_pacer_dc_four_context_embeddings.py:176), [default start calculation](/C:/projects/GPCR-virtual-screening-encoder/project/pacer_dc_training/extract_pacer_dc_four_context_embeddings.py:499).

### 3.2 Raw window → atom14

**Implementation:** `atom14_index_map`, `stage_atom14`.

- Input: raw coordinates and topology-derived atom/residue identities.
- Output: float32 `(T,L,14,3)` NPY and sequence CSV.
- Residue order is preserved through ascending `res_index`.
- Slot order comes from MDGen residue constants.
- Absent slots are zero-filled.
- The validity mask and original residue IDs are not included in the atom14 NPY/sequence pair.

**Potential information loss:** actual atom presence becomes indistinguishable from padding without upstream metadata. Sequence and coordinates alone do not preserve full chain/residue identity.

Reference: [atom14 mapping and storage](/C:/projects/GPCR-virtual-screening-encoder/project/pacer_dc_training/extract_pacer_dc_four_context_embeddings.py:241).

The public-trajectory G0 converter is a separate historical entry point, using PSF/DCD and segments A+B. It should not be confused with the production chain-E selector. See [G0 converter](/C:/projects/GPCR-virtual-screening-encoder/project/pacer_dc_training/oneprot_g0_audit/scripts/convert_dcd_to_atom14.py:62).

### 3.3 Atom14 → encoder inputs

**Implementation:** `flatten_atom14`.

- Input: `(T,L,14,3)` and sequence.
- Output: `xyz (T,N,3)`, atomic numbers `(N,)`, residue index `(N,)`.
- Atoms are enumerated residue-first, then in the hardcoded atom14 slot order.
- Only chemically expected heavy-atom slots are selected.
- Checks: rank/shape, sequence length, supported amino acids, finite coordinates.
- No coordinate normalization or unit conversion occurs.

Archived pilot metadata records `L=270`, `N=2139`. See [flattening implementation](/C:/projects/GPCR-virtual-screening-encoder/project/pacer_dc_training/extract_geom2vec_atom14.py:106) and [pilot input dimensions](/C:/projects/GPCR-virtual-screening-encoder/project/results/pacer_dc_geom2vec_pilot_v01/apo_stride20.audit.json:8).

### 3.4 Frame selection and frozen forward pass

**Implementation:** `embed`; extraction entry points.

- Frame selection: `arange(0,T,stride)`, starting at local frame zero.
- Single-file default stride: 10; batch default: 1; historical pilot used 20.
- Batched positions become `(B×N,3)`, with graph membership distinguishing frames.
- Returned scalar features are reshaped to `(B,N,Hx)`.
- Vector features are reshaped to `(B,N,3,Hv)`.
- Model runs in evaluation and inference mode.
- Output frames are appended in original selected-frame order.

The wrapper itself contains no temporal mixing. Exact internal encoder behavior remains dependent on the missing package implementation.

Reference: [forward pass and frame selection](/C:/projects/GPCR-virtual-screening-encoder/project/pacer_dc_training/extract_geom2vec_atom14.py:138).

### 3.5 Atomic outputs → residue embeddings → storage

**Implementation:** `invariant_residue_features`; NPZ writers.

- Atom invariant descriptor: `(N,Hx+Hv)`.
- Residue feature: `(L,Hx+Hv)`.
- Stored trajectory: `(Tselected,L,D)`, archived `D=128`.
- Global feature: residue mean `(Tselected,D)`.
- All stored embeddings are float32.
- No L2 normalization is applied by the extractor.
- No atomic outputs or intermediate layers are retained.

The single-file NPZ contains `residue_index`; the batch NPZ does not. Both contain sequence and local `frame_ids`.

References: [readout](/C:/projects/GPCR-virtual-screening-encoder/project/pacer_dc_training/extract_geom2vec_atom14.py:128), [single-file storage](/C:/projects/GPCR-virtual-screening-encoder/project/pacer_dc_training/extract_geom2vec_atom14.py:178), [batch storage and audit](/C:/projects/GPCR-virtual-screening-encoder/project/pacer_dc_training/run_geom2vec_four_context_batch.py:49).

**Irreversible losses:** vector orientation, atomic heterogeneity within residues, intermediate activations. Global pooling is also lossy, but the NPZ retains residue features separately.

### 3.6 Normalization and dimensionality reduction

Three distinct behaviors exist:

| Path | Transformation | Calibration |
|---|---|---|
| Linear Geom2Vec pilot | Raw feature means | No channel scaling in this script |
| Original FKG | Optional per-frame residue-global subtraction; channel median/MAD scaling | Refit on each loaded four-context unit |
| G2-B/G2-C | Per-frame residue-global subtraction; `(h−median)/max(1.4826×MAD,1e−6)` | G2-B fits R2 only; G2-C reuses its parameters |

G2-B’s `unit_normalize` is misleadingly named: it subtracts the residue mean; it does **not** produce unit-length vectors.

References: [historical linear aggregation](/C:/projects/GPCR-virtual-screening-encoder/project/pacer_dc_training/audit_geom2vec_cross_replica.py:29), [original FKG preprocessing](/C:/projects/GPCR-virtual-screening-encoder/project/pacer_dc_training/pacer_factorial_kernel_graph.py:34), [G2-B normalization](/C:/projects/GPCR-virtual-screening-encoder/project/pacer_dc_training/analyze_pacer_fkg_g2b.py:38).

G2-B calibration uses eight evenly spaced frames from each of five R2 windows and four contexts: 160 frame-level samples, pooled over residues for channel statistics. R3 receives the frozen transform. See [calibration](/C:/projects/GPCR-virtual-screening-encoder/project/pacer_dc_training/analyze_pacer_fkg_g2b.py:154).

**Loss:** residue-global subtraction removes any channel shift shared across all residues in a frame. Whether that removes nuisance or meaningful response is unresolved.

No PCA/tICA/VAMP stage was found in the traced Geom2Vec→FKG path. Historical CA-coordinate baseline code exists, but its CLI is explicitly disabled because it fits separate context-specific bases. RFF is a kernel feature approximation, not PCA compression. See [disabled baseline path](/C:/projects/GPCR-virtual-screening-encoder/project/pacer_dc_training/extract_pacer_dc_four_context_embeddings.py:489).

### 3.7 Residue/region interface

There is no single universal region representation.

**Linear pilot:** `(T,L,128) → (L,128) → (128,)` through time and region means. Both means discard heterogeneity. It uses `G2_REGION_MAP_v02.json`.

**G2-B:** select ordered region indices, then concatenate:

`(100,270,128) → (100,K,128) → (100,K×128) → (100,512)`.

Region flattening preserves residue identity and simultaneous within-frame configurations. Each region gets its own R2 bandwidth and seeded RFF map, shared across contexts/windows/replicas. It is not one globally shared region space. See [G2-B region interface](/C:/projects/GPCR-virtual-screening-encoder/project/pacer_dc_training/analyze_pacer_fkg_g2b.py:187).

**G2-C:** apply one shared 128-input-dimensional node kernel:

`(100,270,128) → (100,270,256)`.

Regional averaging occurs after signed node contrasts and optional graph propagation. This loses within-region residue identity in the final vector. See [shared-node transform](/C:/projects/GPCR-virtual-screening-encoder/project/pacer_dc_training/analyze_pacer_fkg_g2c.py:132).

The FKG map is `M4_MULTISTRUCTURE_GRAPH_v01.json`, with nine regions and 270 nodes. Its builder maps structure residue IDs to embedding indices; region entries are sorted by structure residue ID. The stored map records ligand-contact cutoff 4.5 Å and graph cutoff 8 Å.

References: [mapping and ordered entries](/C:/projects/GPCR-virtual-screening-encoder/project/pacer_dc_training/build_m4_multistructure_graph.py:34), [graph artifact](/C:/projects/GPCR-virtual-screening-encoder/project/results/pacer_dc_geom2vec_pilot_v01/M4_MULTISTRUCTURE_GRAPH_v01.json:1).

### 3.8 Temporal aggregation, signed contrasts, and statistics

G2-B averages each 100-frame RFF sequence into five contiguous 20-frame blocks. Per-context storage becomes `(5 windows,5 blocks,512)`. G2-C retains nodes: `(5,5,270,256)`.

The contrast definitions are:

- ΔPAM: `candidate_probe − probe_only`;
- ΔAGO: `candidate_no_probe − apo`;
- ΔINT: `candidate_probe − probe_only − candidate_no_probe + apo`.

These combine context distribution summaries, not physically paired conformations.

G2-B compares signed 512-dimensional contrast vectors between replicas. G2-C propagates signed `(270,256)` vectors with:

\[
V^{(k+1)}=(1-\alpha)V^{(0)}+\alpha PV^{(k)},
\]

using default `α=0.65`, 20 steps, then averages region nodes. References: [G2-B aggregation](/C:/projects/GPCR-virtual-screening-encoder/project/pacer_dc_training/analyze_pacer_fkg_g2b.py:202), [G2-C propagation](/C:/projects/GPCR-virtual-screening-encoder/project/pacer_dc_training/analyze_pacer_fkg_g2c.py:38), [G2-C contrasts](/C:/projects/GPCR-virtual-screening-encoder/project/pacer_dc_training/analyze_pacer_fkg_g2c.py:158).

The pooled point estimates discard temporal order. Block-based sensitivity calculations still depend on block membership; arbitrary frame shuffling can therefore change sensitivity results even when pooled point estimates remain unchanged.

U2 v02 uses the G2-C shared-node representation without restoring the old gate. Its archived 27-item G2-C crosscheck passed, with maximum absolute difference approximately `1.85e−8`. See [erratum audit](/C:/projects/GPCR-virtual-screening-encoder/project/results/pacer_dc_fkg_U2_erratum_v02/FKG_U2_ERRATUM_AUDIT_v02.json:1).

## 4. Encoder versus post-processing versus downstream boundaries

| Boundary | Current implementation | Permissible optimization boundary |
|---|---|---|
| **A. Encoder architecture/pretrained weights** | External Geom2Vec model, local checkpoint | Frozen-layer readout experiments are possible; architecture/weight changes would need separate justification. |
| **B. Extraction/preprocessing** | Atom14 mapping, atom selection, radius backend, stride, inference | Validation, provenance, precision and extraction consistency are representation work. |
| **C. Post-processing/normalization** | Vector norms; global residualization; robust channel scaling | Separate, versioned representation experiments are permissible. Frozen G2-B/G2-C implementations and calibrations remain intact. |
| **D. Residue/region pooling** | Atom mean; historical region mean; G2-B concatenation; G2-C node mean | Alternative readouts may be studied separately. Region membership and frozen downstream estimands must not be silently changed. |
| **E. Temporal aggregation** | Independent frame encoding; block/window means | Representation-side temporal descriptors can be investigated separately; altering frozen block statistics is outside scope. |
| **F. Replica reproducibility** | R2 calibration, R3 transform; descriptive comparisons | A validation boundary, not a target to maximize. |
| **G. Kernel/RFF configuration** | Region-specific 512D G2-B; shared-node 256D G2-C | Hold fixed. Existing approximation audits inform interpretation, not a seed/dimension search. |
| **H. PACER-FKG statistics** | Contrasts, graph propagation, norms/cosines, sensitivity summaries | Frozen. No gate restoration, replacement, threshold tuning, or statistical-definition changes. |

Code location alone does not determine scientific ownership: channel normalization is representation processing even though it currently resides inside analysis scripts. A new experiment should use a separate adapter and must avoid applying normalization twice.

## 5. Confirmed issues and unresolved questions

### Confirmed implementation limitations

1. **Missing-atom validation is incomplete.** Atom14 preparation zero-fills missing slots; Geom2Vec flattening selects expected atoms from sequence without reading an actual presence mask. An unexpectedly absent heavy atom could enter inference at zero coordinates. This is a confirmed validation gap, **not confirmed corruption of historical data**.

2. **Finite-output checks do not fail extraction.** The batch writer records `finite`, but assigns `status="ok"` unconditionally. Consequently, `files_ok` alone does not establish numerical validity. Later G2-B/G2-C preflights explicitly reject nonfinite embeddings. See [batch status assignment](/C:/projects/GPCR-virtual-screening-encoder/project/pacer_dc_training/run_geom2vec_four_context_batch.py:67).

3. **Embedding files are not self-contained provenance records.** Batch NPZs omit original residue IDs, chain IDs, absolute trajectory frames, timestamps, atom masks, checkpoint identity, and preprocessing version. Some information exists upstream or in the batch audit, but reconstructing the full chain requires joining artifacts.

4. **Extraction provenance is incomplete for exact reproduction.** The batch audit records input hashes, checkpoint hash, device, stride, batch size and backend flag, but not complete package versions, hardware details, script hashes, or output-file hashes. Later FKG manifests add NPZ hashes.

5. **Representation losses are real; their harm is unproven.** Vector norms remove orientation; atom means remove within-residue variation; global residualization removes frame-wide shifts; temporal means remove order. None alone proves poor scientific performance.

6. **Hardcoded analysis shapes constrain compatibility.** G2-B/G2-C expect 100 frames, 270 residues, 128 channels, two named replicas and five windows. Changing layer dimension or frame count is not a drop-in change to frozen analyses.

### Unresolved implementation questions

- Exact ViSNet layers, internal normalization, neighbor cap, cutoff behavior, checkpoint-loading strictness, and accessible intermediate activations remain unverified because package source is missing.
- The fallback claims mathematical equivalence, but implements nearest-neighbor truncation with `topk`, uses `<=r`, and ignores `flow`. Equivalence to the compiled backend under neighbor saturation, distance ties, and boundary cases is untested here. See [fallback implementation](/C:/projects/GPCR-virtual-screening-encoder/project/pacer_dc_training/extract_geom2vec_atom14.py:50).
- Float32 global subtraction may lose precision when residuals are small relative to channel offsets. No measured instability was established.
- A `1e−6` MAD floor prevents division by zero but does not establish that nearly constant channels are harmless.
- The fixed 10-ps metadata assumption needs reconciliation with actual source trajectory metadata before time-based representation work.

### Historical evidence and its limits

The archived G2-A audit reports 40 geometry-screen passes, matching masks, and zero mapping/missing-slot errors. The separate PBC report is explicitly a screen, not a universal reconstruction proof. These reduce concern about historical corruption but do not replace missing-input verification. See [G2-A audit](/C:/projects/GPCR-virtual-screening-encoder/project/results/pacer_dc_fkg_G2A_v02/G2A_AUDIT.json:1) and [PBC audit](/C:/projects/GPCR-virtual-screening-encoder/project/results/pacer_dc_fkg_G2A_PBC_v01/G2A_PBC_AUDIT.json:1).

G2-B documents selected repeated ΔPAM directions and limitations of two replicas. The supplementary exact-kernel comparison preserves some historical signs but explicitly distinguishes subsample estimates from full-data estimates. G2-C reports mixed graph effects. These are background observations, not configuration-selection objectives.

References: [G2-B report](/C:/projects/GPCR-virtual-screening-encoder/project/results/pacer_dc_fkg_G2B_v01/G2B_DIRECTION_REPORT.md:37), [robustness report](/C:/projects/GPCR-virtual-screening-encoder/project/results/pacer_dc_fkg_G2B_robustness_v01/G2B_ROBUSTNESS_REPORT.md:30), [G2-C report](/C:/projects/GPCR-virtual-screening-encoder/project/results/pacer_dc_fkg_G2C_v01/G2C_GRAPH_ABLATION_REPORT.md:37).

### Explicitly missing artifacts

- All 40 R2/R3 `.geom2vec.npz` files under `project/results/pacer_dc_geom2vec_R2R3_full_v01/`; only `batch_audit.json` is present.
- Pilot Geom2Vec NPZ files; pilot audit JSONs remain.
- R2/R3 atom14 and raw coordinate payloads needed for re-extraction.
- `project/tools/oneprot-embeddings/`, including MDGen constants, OneProt implementation and checkpoint.
- A repository-local Geom2Vec/ViSNet package source checkout.
- `project/results/pacer_dc_G2A_source_v01/`.
- Production trajectory and membrane-reference directories referenced by preprocessing.

The R1/W0 atom14 files, historical OneProt embeddings, `G2_REGION_MAP_v02.json`, frozen FKG graph, Geom2Vec checkpoint, and requested historical reports are present.

## 6. Candidate optimization comparison table

These are hypotheses and possible future experiments, not changes made by this audit. “Reuse” assumes the missing original embedding payloads become available.

| Candidate and concrete location | Evidence / suspected issue | Retraining | Reuse / re-embedding | Approximate cost | Risks and independent validation |
|---|---|---|---|---|---|
| **Precision-controlled normalization** — G2-B `unit_normalize`/application; G2-C `load_contrasts` | Float32 residue-global subtraction precedes higher-precision scaling. Cancellation is possible, not demonstrated. | No | Reuse residue NPZs; no re-embedding | Low CPU cost; linear in stored features | Cannot recover precision lost during inference. Compare FP32/FP64 arithmetic against a reference with frozen parameters; evaluate numerical errors, not Δ-axis improvement. |
| **Channel/block normalization diagnostic** — adapter after NPZ loading, before downstream interface | Scalar and vector-norm channels have different constructions; current median/MAD may already compensate. No demonstrated imbalance defect. | No | Reuse; no re-embedding | Low CPU cost | New scaling changes geometry and interacts with kernel calibration. Initially assess channel dispersion/conditioning only on a predefined split; do not substitute into frozen analyses. |
| **Intermediate ViSNet readout** — `embed` around model call, line 149 | Only returned outputs are retained. Earlier features might preserve different information; no superiority evidence. | No for frozen hooks | Existing final embeddings insufficient; re-embedding required | One inference pass plus activation/storage overhead | Layer scales and semantics may differ. First inspect pinned package source; predeclare one layer comparison and assess invariance and independent structural fidelity. |
| **Alternative atom-to-residue readout** — `invariant_residue_features` | Equal atom means erase atom-level heterogeneity and mix backbone/side-chain contributions. | No | Existing residue means insufficient; re-embedding required | Approximately one extraction pass; modest readout overhead | CA-only, backbone-only or split pooling changes residue chemistry weighting. Test a fixed readout with rigid-motion and atom-order consistency plus held-out structural probes. |
| **Region pooling diagnostic** — separate adapter alongside G2-B region selection / G2-C `region_vectors` | Joint concatenation and node averaging preserve different information. Opposing residue changes can cancel under averaging. | No for fixed pooling | Reuse; no re-embedding | Low–moderate CPU cost | Changes the representation’s meaning; cannot be treated as the same frozen estimand. Keep region membership fixed; inspect residue-level retention and predetermined reconstruction diagnostics. |
| **Temporal representation diagnostic** — adapter before existing block means | Current pooled estimates do not encode transitions. | No for fixed lag descriptors | Reuse ordered full-frame NPZs; no re-embedding | Linear feature pass; lag concatenation roughly doubles channels | Short trajectories, uncertain timestamps, boundary handling and new dimensions. Use one predefined lag and chronological holdout; evaluate order sensitivity separately from frozen FKG. |
| **Frame-selection sensitivity** — stored `frame_ids`; extraction stride | Defaults differ between entry points; deterministic offset-zero sampling may be sensitive to sampling phase. | No | Full stride-1 NPZs support subsampling; discarded frames require re-embedding | Low CPU cost when reusable | Selecting favorable windows induces bias. Predeclare offsets/equal frame counts, retain all contexts, and report all results without choosing a winner by contrast cosine. |
| **Backend/batch numerical consistency** — fallback and `embed` | CPU fallback and compiled graph equivalence is asserted but unverified here. | No | Existing embeddings are reference only; small re-embedding required | A few frames; much less than full extraction | GPU reduction order and neighbor truncation can differ. Compare edge sets where accessible, feature errors, rigid-motion invariance and batch-size consistency on fixed inputs. |
| **Replica robustness protocol** — calibration boundary and immutable representation manifest | R2 is both calibration and evaluated replica; only R3 is preprocessing-held-out. | No | Reuse for descriptive checks; independent confirmation needs independent short-trajectory data | Low analysis cost; new data availability dominates | R2/R3 are already inspected. Do not optimize their ΔINT/ΔPAM outcomes; lock choices before a new independent short-replica evaluation. |
| **Metadata and validity contract** — atom14 mapper, NPZ writers, batch audit | Missing masks, incomplete identity metadata, and unconditional “ok” are verified gaps. | No | Existing features reusable if provenance reconstructs uniquely; ambiguity may require re-extraction | Low I/O/hash cost; no inference for validation alone | Never invent identity mappings. Independently reconcile sequence, structure IDs, frame mapping, hashes and checkpoint; reject unresolved mismatches. |

Historical cost context: the archived CUDA run processed 4,000 frames in approximately five minutes; the CPU pilot reports about 3–5 seconds/frame. These are hardware-specific historical observations, not forecasts for this workstation.

## 7. Proposed minimal validation strategy

1. **Restore auditability before evaluation.** Obtain the exact archived embedding payloads through a separately authorized transfer. Verify them against the frozen NPZ manifest. Do not regenerate substitutes and label them original.

2. **Freeze a representation experiment contract.** Specify input hashes, frame IDs, residue order, arithmetic precision, calibration parameters, and expected output shape. Keep checkpoint, region map, RFF settings and downstream files fixed.

3. **Use label-independent engineering criteria first.** Check finite values, deterministic replay, metadata consistency, feature error against a numerical reference, channel degeneracy and preservation of residue/frame identity. For any future re-embedding, add rigid-motion and batch/backend consistency checks.

4. **Keep calibration separate from evaluation.** Reusing frozen R2 parameters allows an engineering replay. R3 is not a fresh scientific holdout after these historical analyses. Any later scientific comparison needs a locked representation and independent short-trajectory replica data.

5. **Respect temporal dependence.** Evaluate contiguous windows and fixed subsampling patterns without treating frames or windows as independent biological replicates. Report numerical sensitivity separately from trajectory sampling variability.

6. **Report unsuccessful and negligible changes.** A precision experiment may show no meaningful difference. That is a valid outcome; it must not trigger an open-ended search for favorable biological contrasts.

Any later downstream evaluation must preserve the frozen definitions and clearly identify new representation inputs. It must not overwrite archived reports or introduce a new qualification gate.

## 8. Scientific constraints and prohibited changes

- The old PACER-FKG `kernel_u2` magnitude gate remains withdrawn.
- `kernel_u2` was a squared-kernel-norm estimator; the problem was not simply “forgetting a square root.”
- Region-minus-stable comparisons from independently calibrated kernel spaces cannot be reinstated as functional qualification.
- G2-B region kernels and G2-C shared-node kernels represent different objects; their cosine magnitudes are not interchangeable.
- U2 v02 explicitly defines **no new gate**. See [erratum conclusions](/C:/projects/GPCR-virtual-screening-encoder/project/results/pacer_dc_fkg_U2_erratum_v02/FKG_U2_ERRATUM_REPORT_v02.md:3).
- Do not select layers, normalization, pooling, frames, seeds or other settings to improve ΔINT cosine, maximize ΔPAM separation, or recover historical findings.
- Do not use long-timescale MD for development or hyperparameter selection.
- Do not modify frozen statistics, graph definitions, gates, archived results or analysis implementations.
- Shared kernels, signed representations and replica direction consistency remain methodological considerations, not evidence of efficacy.
- Positive direction cosine is not pharmacological sign; repeated selected-region signals are not proof of PAM function.
- Historical documentation proposing new controls or transition kernels is not authorization to change frozen downstream methods in this task.

## 9. Recommended next experimental step

After the exact archived NPZs are available, perform **one paired normalization-precision experiment**:

- Baseline: existing float32 residue-global subtraction.
- Variant: float64 subtraction and scaling using the **same frozen median/MAD parameters**, followed by the same final float32 representation format.
- Keep frame/residue order, channels, checkpoint, pooling and all downstream settings unchanged.
- Compare representation-level numerical error, finite values and repeatability only; do not select the variant using ΔINT or ΔPAM outcomes.
- Store future results separately, with immutable provenance.

**Recommendation: start with this numerical normalization comparison. It is reversible, requires neither retraining nor re-embedding, and tests a concrete uncertainty without changing frozen PACER-FKG statistics.**
