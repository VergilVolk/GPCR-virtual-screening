# OneProt-MD G0 reproduction — audit package

Self-contained audit material for the **G0** deliverable of
`START_HERE_ONEPROT_PACER_DC.md`: reproduce the OneProt-MD trajectory encoder and
show that it produces an embedding from a real M4 trajectory.

This package exists because the work was performed inside
`project/tools/oneprot-embeddings/`, a **nested third-party Git checkout that the
main repository ignores wholesale** (`.gitignore`: `project/tools/oneprot-embeddings/`).
That tree is 5.9 GB (3.6 GB checkpoint + public trajectory + caches), so none of
it is reviewable from GitHub. Here we copy out only the small, load-bearing
pieces: the six scripts that did the work, the 4.6 MB atom14 input, and four
outputs — plus a verification script written for this audit that re-derives the
embedding from the package's own inputs.

---

## 1. Purpose

Make the G0 claim checkable by a reviewer who has only this Git branch:

| Question a reviewer should be able to answer | Where the answer lives |
|---|---|
| Does the official checkpoint contain the MD branch? | `../../../runs/oneprot_checkpoint_audit.json` |
| Does it load into the encoder without missing weights? | `../../../runs/oneprot_md_branch_runtime_audit.json` |
| Was the input a real MD trajectory, or synthetic data? | `input/`, `scripts/convert_dcd_to_atom14.py` |
| Did a real forward pass actually run? | `scripts/run_g0_embedding.py`, `outputs/run_log.txt` |
| Is the embedding finite, non-zero, reproducible? | `outputs/embedding_smoke.npy`, `outputs/forward_reverification.txt` |

Not in scope: embedding *quality*, PAM prediction, G1 four-context trajectories.

---

## 2. Upstream OneProt provenance

| item | value |
|---|---|
| repository | `https://github.com/oneprot-models/oneprot-embeddings.git` |
| commit | `53fa9c0f9150aa44a9258db6948277a4602036ab` ("Update README.md", Alina Bazarova, 2026-07-09) |
| local commits ahead of upstream | 0 |
| submodule | `external/mdgen` @ `1b82cc13f18ea6a4a2f680d5c086bd07715b89c5` (`https://github.com/bjing2016/mdgen.git`) |
| checkpoint | Zenodo record 20997998, member `checkpoints/Pocket_Text_ST_SG_MD/epoch_012_01100-v1.ckpt`, 3,614,862,456 B, sha256 `a79bce2e…1997`, zip CRC32 `2f07819f` |

**No upstream OneProt source files were modified.**
**No mdgen source files were modified.**

Verified for this package:

```
git -C project/tools/oneprot-embeddings diff --stat                          -> empty
git -C project/tools/oneprot-embeddings rev-list --count origin/main..HEAD   -> 0
git -C project/tools/oneprot-embeddings/external/mdgen diff --stat           -> empty
```

The only additions inside the nested checkout are untracked: `artifacts/` (this
work) and `**/__pycache__/` byte-code caches. Details, including the two
hard-coded cluster paths that remain in upstream source and why they are inert
here, are in `upstream/oneprot_provenance.txt`.

---

## 3. What was done

Six steps, in order. Each row states the producing script, the input, the output
and the verdict.

1. **Checkpoint retrieval and integrity.**
   `scripts/fetch_md_ckpt.py` pulls only the needed member out of the 14.73 GB
   Zenodo `checkpoints.zip` using HTTP Range requests (resumable, keeps a
   `*.deflate.part` while in flight); `scripts/verify_md_ckpt.py` then checks
   size, recomputes the member CRC32 and prints local md5/sha256.
   → `outputs/checkpoint_sha256.txt`. **Result: size and CRC32 match; the
   zip-level md5 published by Zenodo is for the whole archive and is explicitly
   not used as a per-member check.**

2. **Checkpoint MD-branch inspection.**
   `project/pacer_dc_training/inspect_oneprot_md_checkpoint.py` (already tracked,
   upstream of this package) counts the MD tensors inside the checkpoint.
   → `../../../runs/oneprot_checkpoint_audit.json`.
   **Result: 1442 top-level keys, 125 `network.md.transformer.*` tensors
   (~61 MB) → `embedded_md_transformer_candidate: true`, so the separate
   `forward_sim.ckpt` is not required.**

3. **Runtime weight-load audit.**
   `scripts/audit_md_branch_runtime.py` instantiates
   `TrajectoryEncoder(pretrained=False)` — **no `model_path`, no
   `forward_sim.ckpt`** — and loads the checkpoint's `network.md.*` weights
   after stripping the prefix. It does not trust `strict=False`; it re-derives
   matched / absent / shape-mismatched sets per key and cross-checks them against
   the `load_state_dict` return value.
   → `../../../runs/oneprot_md_branch_runtime_audit.json` (and the sibling
   `.log`, which ends in `RUNTIME_AUDIT_PASSED`).
   **Result: configuration A (the one the smoke actually uses) supplies 131/131
   encoder keys with 0 missing and 0 shape mismatches under `transformer.*`; the
   single tolerated extra is the 4-byte `norm.1.log_logit_scale` buffer, which
   only exists when `use_logit_scale=True` and is never read by
   `TrajectoryEncoder.forward`. Configuration B (counterfactual
   `use_logit_scale=True`) loads with zero missing and zero unexpected keys.
   `"all_checks_passed": true` in both.**

4. **Real M4 trajectory → atom14.**
   `scripts/convert_dcd_to_atom14.py` reads the public M4 trajectory
   (PSF + DCD) with MDAnalysis and maps it onto the atom14 representation MDGen
   consumes.
   → `input/m4xan_i1.npy` (4,586,528 B) + `input/m4xan_test.csv`.
   **Result: `(100, 273, 14, 3)` float32, ångström; 273 residues; the CSV's
   seqres is 273 characters and matches the array.**

5. **OneProt-MD forward smoke.**
   `scripts/run_g0_embedding.py` builds the model inputs from step 4, loads the
   checkpoint, and calls the encoder twice.
   → `outputs/embedding_smoke.npy` (1024 float32) and `outputs/run_log.txt`.
   **Result: output `(1, 1024)`, finite, non-zero, bit-identical between the two
   calls, L2 norm exactly 1.0.**

6. **Independent re-derivation (added for this audit).**
   `scripts/reverify_forward.py` re-runs step 5 in a fresh process using **this
   package's** `input/` files and the original `build_batch` preprocessing, then
   compares against **this package's** `outputs/embedding_smoke.npy`.
   → `outputs/forward_reverification.txt`.
   **Result: `REVERIFY_PASSED` — `identical_to_package: True`,
   `max_abs_diff_vs_saved: 0.0`, and identical sha256 of the raw array bytes
   (`4f484f49d353af01d831eb7525b185863f79769e138b9e29b20cf7363e036e25`).**

---

## 4. Real trajectory provenance

* **Source:** Zenodo record **8136971** (DOI `10.5281/zenodo.8136971`),
  dataset "M4 allosteric xanomeline stability", **replica 1**.
* **Files used:** `M4_xan_ortho_and_allo.psf` (990,472 B) and
  `M4_xan_ortho_and_allo_1.dcd` (~10.1 MB). **Not committed** — size, and the
  package is a code/evidence audit, not a trajectory archive.
* **Selection:** `segid A B and not resname ACE XAN` → 273 protein residues
  (the ACE cap and the XAN ligands are dropped).
* **CHARMM → atom14 fixes applied by the converter:** `HSE/HSD/HSP → HIS`,
  and isoleucine's `CD → CD1`. Atom14 slot names come from mdgen's
  `residue_constants`.
* **Output:** `input/m4xan_i1.npy`, shape `(100, 273, 14, 3)`, float32,
  ångström, 100 frames (`NUM_FRAMES = 100`, matching the encoder's
  `num_frames=100`; `convert_dcd_to_atom14.py:32`).
* **Why ~43 % of the array is zero:** atom14 reserves 14 slots per residue;
  residues that lack a given atom (e.g. GLY has no CB, and hydrogens absent from
  the trajectory) are zero-padded by construction. Verified: nonzero fraction
  0.567, and **no all-zero frames** — i.e. the padding is per-atom, not a
  placeholder trajectory.

This is a genuine MD trajectory, not synthetic or random data: the array is
derived frame-by-frame from a DCD file with per-residue atom index mapping, and
the values span −35.2 to +34.4 Å with a mean of 0.013 Å.

---

## 5. Forward path

Confirmed by reading the source, not inferred from filenames.

```
input/m4xan_i1.npy                     (100, 273, 14, 3) float32, Å
  └─ scripts/run_g0_embedding.py::build_batch()          lines 45-95
       mirrors MDDataset.collate_fn + NewMDGenWrapper.prep_batch in plain torch
       (the upstream dataset classes import pytorch_lightning, which this
        project never uses)
       → atom14_to_frames → offsets (7) + torsions (14) → latents
  └─ latents [1, 100, 273, 21]  +  model_kwargs {start_frames, end_frames,
                                                mask, aatype, x_cond, x_cond_mask}
  └─ src/models/components/md_encoder.py::TrajectoryEncoder.forward()   lines 182-216
       self.transformer(x=latents, t=0, mask, start_frames, end_frames,
                        x_cond, x_cond_mask, aatype)
       == mdgen LatentMDGenModel  ← the OneProt MD branch
  └─ hidden states                                                     lines 218-251
       4-D output → mean over time, then masked pooling over residues
       3-D output → masked pooling over residues
  └─ self.proj(pooled) → self.norm(projected)                          lines 256-257
       the final L2 normalisation is why the norm is exactly 1.0
  └─ outputs/embedding_smoke.npy        (1024,) float32
```

* **Weights actually used:** the encoder is built with `pretrained=False`, which
  initialises a fresh `LatentMDGenModel(MDGenArgs(), latent_dim=21)`; the
  checkpoint's `network.md.*` tensors are then loaded on top (prefix stripped).
  Because the runtime audit shows *every* encoder key was supplied, the
  transformer, the projection head and the normalisation all come from the
  checkpoint — **no layer in the embedding path is left randomly initialised**.
* **Not bypassed:** there is no precomputed-embedding lookup and no early exit;
  the transformer runs and the output flows through projection and
  normalisation.
* **Device:** CPU (the `oneprot-md` environment has a CPU-only torch 2.4.1 build).

**Verdict — real-trajectory OneProt-MD forward: YES.**

Two honest caveats, so the verdict is not read as more than it is:

1. The **preprocessing** is a faithful but *re-implemented* copy of OneProt's
   dataset collation (`build_batch`), because the upstream classes require
   `pytorch_lightning`. The model call itself is upstream code, unmodified.
2. The `t=0` argument is the flow-matching time; this is an
   **encoder-style read-out** of the MD branch at t=0, which is what the
   PACER-DC adapter consumes. It is not an MDGen sampling/generation run.

---

## 6. Outputs

| file | bytes | what it is |
|---|---|---|
| `outputs/checkpoint_sha256.txt` | 769 | checkpoint provenance: source URL, archive member, size, sha256, local md5, zip CRC32 |
| `outputs/embedding_smoke.npy` | 4,224 | the G0 embedding: `(1024,)` float32, L2-normalised, finite, all 1024 entries non-zero |
| `outputs/oneprot_environment_audit.json` | 909 | environment/provenance record: upstream commit, mdgen submodule commit, checkpoint sha256, torch/MDAnalysis versions, input shape, output shape, finite/non-zero/reproducible flags |
| `outputs/run_log.txt` | 1,264 | stdout of step 5, ending in `G0_EMBEDDING_OK`. **Encoding: UTF-16LE** (verbatim copy of the file produced at the time) — read it with `iconv -f UTF-16 -t UTF-8 outputs/run_log.txt` |
| `outputs/forward_reverification.txt` | 1,345 | UTF-8 report of the independent re-derivation (§3 step 6), ending in `REVERIFY_PASSED` |

Companion evidence lives outside this directory, already committed on this
branch (commit `2e23271`):

* `runs/oneprot_checkpoint_audit.json` — step 2 output
* `runs/oneprot_md_branch_runtime_audit.json` — step 3 output
* `runs/oneprot_md_branch_runtime_audit.log` — step 3 stdout, ends in `RUNTIME_AUDIT_PASSED`

They are referenced rather than duplicated so that the two copies cannot drift.

---

## 7. Claim boundary

This package **does** support:

* the official OneProt checkpoint contains a complete MD transformer branch;
* the encoder can be instantiated **without** the cluster-only `forward_sim.ckpt`
  and loaded from that checkpoint with zero missing and zero shape-mismatched
  keys under `network.md.transformer.*`;
* a real, public M4 trajectory can be converted into the encoder's input format
  by the shipped script;
* the MD branch executes a genuine forward pass on that trajectory and returns a
  finite, non-zero, L2-normalised 1024-d embedding;
* that embedding is deterministic within a process and **bit-identical when
  re-derived in a fresh process from this package's own inputs**.

This package **does not** support, and must not be cited as:

* any statement about embedding *quality* or about PAM functional prediction —
  no biological endpoint is involved anywhere in this chain;
* any G1 claim: the four-context production trajectories are **not** produced by
  this work, and no four-context candidate exists;
* any claim that a candidate molecule is a PAM — there is no wet-lab evidence in
  this repository at all;
* any claim that OneProt-MD is the PACER-DC backbone. It is a candidate encoder
  that still has to beat tICA/VAMP and the manual differential endpoints on
  held-out chemotypes;
* generality beyond the single public trajectory and the single 100-frame window
  used here (one replica, one system).

---

## 8. Large files intentionally omitted

| omitted | size | why |
|---|---|---|
| `epoch_012_01100-v1.ckpt` | 3.37 GB | model weights; excluded by `.gitignore` (`*.ckpt`) and far beyond any code-audit need — re-fetchable via `scripts/fetch_md_ckpt.py` |
| `epoch_012_01100-v1.ckpt.deflate.part` | 2.30 GB | leftover partial download from the range fetch |
| raw M4 PSF + DCD | ~11 MB | source trajectory; the derived 4.6 MB atom14 array is shipped instead |
| nested `artifacts/` tree | 5.9 GB | the whole nested working tree (checkpoint, trajectory, caches) |
| `project/tools/oneprot-embeddings/**` | — | third-party checkout; ignored by the main repository by design |
| WSL production trajectories / checkpoints / state XMLs | 4.90 GB | G1 runtime state, lives outside Git entirely |
| migration tarballs (`*.tar.gz`) | 3.08 GB | transfer artefacts; excluded by `.gitignore` |

None of these is required to review the logic or the evidence above.

---

## 9. Reproducing this package

Prerequisites (deliberately not vendored):

1. `project/tools/oneprot-embeddings` at `53fa9c0` with `external/mdgen` at
   `1b82cc1` — `project/scripts/setup_wsl_pacer_dc_md.sh` and
   `upstream/oneprot_provenance.txt` document how it is obtained;
2. the checkpoint at
   `project/tools/oneprot-embeddings/artifacts/Pocket_Text_ST_SG_MD/epoch_012_01100-v1.ckpt`
   — fetch/inflate with `scripts/fetch_md_ckpt.py`, then check with
   `scripts/verify_md_ckpt.py`;
3. a CPU torch environment (validated with python 3.11.16, torch 2.4.1,
   numpy 2.4.6; `project/environment_pacer_dc_train_cpu.yml`).

Then, from the repository root:

```bash
# weights are not needed for this step
python project/pacer_dc_training/oneprot_g0_audit/scripts/audit_md_branch_runtime.py \
  --output runs/oneprot_md_branch_runtime_audit.json

# the forward itself (overwrites outputs/embedding_smoke.npy - copy it first if you
# want a diff)
python project/pacer_dc_training/oneprot_g0_audit/scripts/run_g0_embedding.py

# and the check that needs neither the raw DCD nor a second copy of anything
python project/pacer_dc_training/oneprot_g0_audit/scripts/reverify_forward.py
```

`reverify_forward.py` exits non-zero unless it reproduces
`outputs/embedding_smoke.npy` bit-for-bit.

---

## 10. Package layout

```
project/pacer_dc_training/oneprot_g0_audit/
├── README.md                        this file
├── scripts/
│   ├── fetch_md_ckpt.py             range-extract the checkpoint member from Zenodo
│   ├── verify_md_ckpt.py            size + CRC32 + hashes of the extracted member
│   ├── load_md_branch.py            first successful load (observational)
│   ├── audit_md_branch_runtime.py   assertive load audit -> runs/*.json
│   ├── convert_dcd_to_atom14.py     PSF+DCD -> atom14 npy (step 4)
│   ├── run_g0_embedding.py          the forward smoke (step 5)
│   └── reverify_forward.py          added for this audit: bit-exact re-derivation
├── input/
│   ├── m4xan_i1.npy                 (100, 273, 14, 3) float32
│   └── m4xan_test.csv               name,seqres (273 aa)
├── outputs/
│   ├── checkpoint_sha256.txt
│   ├── embedding_smoke.npy          (1024,) float32
│   ├── oneprot_environment_audit.json
│   ├── run_log.txt                  UTF-16LE stdout of step 5
│   └── forward_reverification.txt   UTF-8 stdout of step 6
└── upstream/
    └── oneprot_provenance.txt       remote/commit/submodule + no-modification record
```

The six `scripts/` files other than `reverify_forward.py` are byte-identical
copies of the files that produced the evidence; the originals remain in the
nested checkout, which this package does not modify or delete from.
