# Frozen DrugCLIP apply-only runbook (v01)

**Every command below is inference-only. Training, fine-tuning, adapter fitting,
PocketRealign and any parameter update are forbidden in this path.**

## 0. Frozen identities (must match exactly)

| Item | Identity |
|---|---|
| DrugCLIP source | `bowen-gao/DrugClip` @ `7a3a3fa33673f8668c811790f2e4681c98af44ef` |
| Uni-Core source | `dptech-corp/Uni-Core` @ `44f6386f4dcd7137fc1e5d5e768117d635d64a26` |
| Official checkpoint | `checkpoint_best.pt`, 1,183,713,459 bytes, sha256 `dc2c76d0f02f9bb079a613f09d538dcda1bf9075f2952d91dc1bea55571f667e` |
| Embedding / representation dims | 128 / 512 |
| M4 pockets (frozen, chain R) | `7TRQ_M4_allosteric`, `7TRP_M4_allosteric`, `7TRS_M4_allosteric` |
| Source PDBs | `project/data/pdb/{7TRQ,7TRP,7TRS}.pdb` - hashes pinned in `model_a.EXPECTED_PDB_SHA256` |
| Molecule conformer policy | 1 conformer, seed `20260922` |
| Model B deployed weights | `project/results/drugclip_science2026/family_aug_v01/science2026_13target_ep80_famaug.seed{20260925,20260926,20260927}.projection.pt` |
| Dual-model rule | `old_top50 INTERSECT new_top25` (**no score fusion**) |

## 1. What has to be restored before Model A can run

| # | Asset | Target | Verification |
|---|---|---|---|
| 1 | DrugCLIP checkout @ pinned commit | `project/tools/DrugCLIP` | `git rev-parse HEAD` |
| 2 | Uni-Core checkout @ pinned commit | `project/tools/Uni-Core` | `git rev-parse HEAD` |
| 3 | Official checkpoint | `project/tools/DrugCLIP/artifacts/checkpoint_best.pt` | size **and** SHA256 |
| 4 | `drugclip-cpu` WSL env + micromamba | `bash project/scripts/setup_drugclip_cpu_wsl.sh` | `micromamba run -n drugclip-cpu python -V` |

The molecule and pocket LMDBs are **not** restored - they are rebuilt
deterministically from frozen inputs. The three PDBs are already present and
hash-verified.

## 2. Model A - official frozen apply-only

```bash
# WSL, from the code-freeze worktree
bash project/scripts/setup_drugclip_cpu_wsl.sh        # pins commits + audited CPU device-placement patch
bash project/scripts/run_drugclip_m4_cpu_baseline.sh  # PDB hash gate -> build LMDBs -> extract -> evaluate
```

`prepare_drugclip_cpu_checkout.py` asserts the checkout HEAD equals the pinned
commit and applies **exactly four** `.cuda()` -> `device=` substitutions in
`unimol/models/drugclip.py`, each required to occur exactly once, recording an
audit. Scope is device placement only; architecture and tensors are untouched.

`extract_drugclip_embeddings_cpu.py` runs `model.cpu().float().eval()` under
`torch.inference_mode()`, loads with `strict=False`, and serialises
`missing_checkpoint_keys` / `unexpected_checkpoint_keys` into the archive.
**Both must be empty.**

### Candidate scoring (new molecules)

```bash
micromamba run -n drugclip-cpu python project/scripts/build_drugclip_inputs.py molecules \
  --csv <candidates.csv> --output out/molecules.lmdb \
  --id-column <id col> --smiles-column <smiles col> --conformers 1 --seed 20260922

micromamba run -n drugclip-cpu python project/scripts/extract_drugclip_embeddings_cpu.py \
  --drugclip project/tools/DrugCLIP --unicore project/tools/Uni-Core \
  --checkpoint project/tools/DrugCLIP/artifacts/checkpoint_best.pt \
  --molecules out/molecules.lmdb \
  --pockets out/7TRQ_pocket.lmdb out/7TRP_pocket.lmdb out/7TRS_pocket.lmdb \
  --output out/embeddings.npz --batch-size 16
```

### Validate the archive (fail closed)

```bash
python -m project.drugclip_freeze.apply validate --archive out/embeddings.npz
```

Checks: required keys, unique molecule ids, shapes `(n,128)`/`(n,512)`,
`scores == (n_pockets, n_molecules)`, finiteness, non-empty M4 pocket set, and
**empty missing/unexpected checkpoint keys** when the audit JSON is present.

Report all three pockets plus `state_mean` and `state_max`. **Never select the
best pocket by labels.**

## 3. Model B - frozen projection-only

```bash
python -m project.drugclip_freeze.apply model-b \
  --archive out/embeddings.npz --output out/modelb_scores.csv
```

It consumes Model A's frozen 512-D `molecule_representations` and
`pocket_representations`, applies the three deployed MLP heads
(`linear1: 512x512` -> ReLU -> `linear2: 128x512` -> L2-normalise), and takes
the **official max over the declared M4 pockets**, then averages the three seeds.

- The backbone is **not** re-run.
- No gradient is computed anywhere.
- No training module (`finetune_*`, `train_*`) is imported: the forward pass is
  re-expressed functionally in `projection.py` and every state dict is validated
  key-by-key and shape-by-shape before use.

## 4. Dual-model decision

```bash
python -m project.drugclip_freeze.apply decide --pairs pairs.csv --output out/dual_model.json
python -m project.drugclip_freeze.apply verify-shortlist \
  --shortlist project/results/project_wide_integration_benchmark_v02/candidate_dual_model_shortlist_v02.csv
```

`pairs.csv` columns: `candidate_id, old_rank, new_rank`.

- Model A (official weights) is the **M4 decision line**.
- Model B (family-augmented) is a **second opinion**.
- Decision = `old_rank <= 50 AND new_rank <= 25`.
- **Score fusion is forbidden.** `assert_no_fusion` guards the decision inputs;
  descriptive fusion columns may exist on benchmark artifacts but may never enter
  the decision. Rank fusion was evaluated upstream and rejected
  (ROC +0.002 but EF1% 1.49 vs 1.75).

## 5. Explicitly forbidden

| Script family | Why |
|---|---|
| `screen_m4_candidates_drugclip.py` | fits a `PocketRealign` adapter (300 epochs x 3 seeds) before ranking |
| `finetune_drugclip_*` (20 scripts) | trains projections |
| `train_drugclip_*`, `train_muscarinic_*` | adapter fitting |
| `fit_gpcr_drugclip_ce_checkpoint.py` | fits a checkpoint |

`evaluate_drugclip_*` / `summarize_drugclip_*` / `compare_drugclip_*` are
read-only benchmarks over frozen artifacts and are permitted for audit only.

## 6. Claim boundary

DrugCLIP is a **binding-compatibility** representation. Its cosine is not a PAM
probability, not efficacy, and not potency. The dual-model intersection is a
prioritization device. Every downstream statement must stay inside
"prospective representation-space mechanistic evidence".
