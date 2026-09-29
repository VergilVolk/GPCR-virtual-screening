#!/usr/bin/env python3
"""Strict leave-one-target-out audit of DrugCLIP internal projection LoRA.

The held target, its labels, and every held-target scaffold are absent from the
optimization objective.  Retrieval is normalized over seen pockets only.  The
experiment tests target transfer, not potency or PAM efficacy.
"""
from __future__ import annotations

import argparse, copy, json
from pathlib import Path
import numpy as np
import pandas as pd
import torch
from torch.nn import functional as F
from sklearn.linear_model import LogisticRegression

from finetune_drugclip_gpcr_screening import TARGETS, FastProjectionLoRA, balanced_weights, screening_metrics, fp_matrix
from finetune_drugclip_muscarinic_triplet import Proj


class Model(torch.nn.Module):
    def __init__(self, initial, rank, seed):
        super().__init__()
        self.mol = FastProjectionLoRA(initial["mol_project"], rank, seed)
        self.pocket = FastProjectionLoRA(initial["pocket_project"], rank, seed + 17)


def fit(mol, pocket, pair_mol, pair_target, labels, seen, initial, seed, args, random_labels=False):
    model = Model(initial, args.rank, seed)
    with torch.no_grad():
        hm = model.mol.hidden(mol); hp = model.pocket.hidden(pocket)
        z0m = model.mol.base_from_hidden(hm); z0p = model.pocket.base_from_hidden(hp)
    labels = np.asarray(labels, dtype=int).copy()
    if random_labels:
        rng = np.random.default_rng(seed)
        for target in seen:
            keep = np.flatnonzero(pair_target == target)
            labels[keep] = rng.permutation(labels[keep])
    pm = torch.as_tensor(np.asarray(pair_mol, dtype=np.int64).copy())
    pt = torch.as_tensor(np.asarray(pair_target, dtype=np.int64).copy())
    y = torch.as_tensor(labels, dtype=torch.float32)
    weights = torch.as_tensor(balanced_weights(pair_target, labels))
    train_mol = torch.unique(pm); seen_t = torch.as_tensor(seen, dtype=torch.long)
    remap = {target: j for j, target in enumerate(seen)}
    retrieval = {}
    for m, t, label in zip(pair_mol, pair_target, labels):
        if label:
            retrieval.setdefault(int(m), set()).add(int(t))
    retrieval = [(m, remap[next(iter(ts))]) for m, ts in retrieval.items()
                 if len(ts) == 1 and next(iter(ts)) in remap]
    rm = torch.as_tensor([row[0] for row in retrieval], dtype=torch.long)
    rt = torch.as_tensor([row[1] for row in retrieval], dtype=torch.long)
    optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    for _ in range(args.epochs):
        optimizer.zero_grad()
        zm = model.mol.from_hidden(hm); zp = model.pocket.from_hidden(hp)
        scores = zm @ zp.T
        bce = (F.binary_cross_entropy_with_logits(scores[pm, pt] / args.temperature, y,
                                                  reduction="none") * weights).mean()
        retrieval_loss = F.cross_entropy(scores[rm][:, seen_t] / args.temperature, rt)
        preserve = ((1 - (zm[train_mol] * z0m[train_mol]).sum(1)).mean()
                    + (1 - (zp[seen_t] * z0p[seen_t]).sum(1)).mean())
        loss = bce + args.retrieval_weight * retrieval_loss + args.preserve_weight * preserve
        loss.backward(); torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0); optimizer.step()
    return model


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--representations", type=Path, required=True)
    p.add_argument("--projection", type=Path, required=True)
    p.add_argument("--pairs", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--seeds", default="20260925,20260926,20260927")
    p.add_argument("--rank", type=int, default=8)
    p.add_argument("--epochs", type=int, default=60)
    p.add_argument("--lr", type=float, default=5e-3)
    p.add_argument("--weight-decay", type=float, default=1e-3)
    p.add_argument("--temperature", type=float, default=.07)
    p.add_argument("--retrieval-weight", type=float, default=.25)
    p.add_argument("--preserve-weight", type=float, default=.2)
    p.add_argument("--bootstrap", type=int, default=500)
    a = p.parse_args(); torch.set_num_threads(max(1, min(8, torch.get_num_threads())))
    arc = np.load(a.representations, allow_pickle=False); initial = torch.load(a.projection, map_location="cpu")
    pairs0 = pd.read_csv(a.pairs); mi = {str(v): i for i, v in enumerate(arc["molecule_ids"])}
    pairs = pairs0[pairs0.canonical_smiles.astype(str).isin(mi)].copy().reset_index(drop=True)
    pair_mol = np.asarray([mi[str(v)] for v in pairs.canonical_smiles], dtype=int)
    pair_target = np.asarray([TARGETS.index(v) for v in pairs.target], dtype=int)
    labels = pairs.label.to_numpy(int); scaffolds = pairs.murcko_scaffold.astype(str).to_numpy()
    ids = list(map(str, arc["pocket_ids"])); porder = [ids.index(f"{target}_cluster0") for target in TARGETS]
    mol = torch.as_tensor(arc["molecule_representations"].astype(np.float32))
    pocket = torch.as_tensor(arc["pocket_representations"][porder].astype(np.float32))
    frozen_m = Proj(copy.deepcopy(initial["mol_project"])).eval(); frozen_p = Proj(copy.deepcopy(initial["pocket_project"])).eval()
    with torch.inference_mode(): official_matrix = (frozen_m(mol) @ frozen_p(pocket).T).numpy()
    official = official_matrix[pair_mol, pair_target]
    seeds = [int(v) for v in a.seeds.split(",")]
    tuned_seeds = [np.full(len(pairs), np.nan) for _ in seeds]
    random_seeds = [np.full(len(pairs), np.nan) for _ in seeds]
    ecfp = np.full(len(pairs), np.nan); fingerprints = fp_matrix(list(map(str, arc["molecule_ids"])))
    audit = {}
    for held in range(len(TARGETS)):
        test = pair_target == held; held_scaffolds = set(scaffolds[test])
        train = (pair_target != held) & ~np.asarray([s in held_scaffolds for s in scaffolds])
        seen = [j for j in range(len(TARGETS)) if j != held]
        audit[TARGETS[held]] = {"train_pairs": int(train.sum()), "test_pairs": int(test.sum()),
                                "purged_rows": int(((pair_target != held) & ~train).sum()),
                                "held_target_in_loss": False, "scaffold_overlap": 0}
        classifier = LogisticRegression(C=1.0, class_weight="balanced", max_iter=3000,
                                        solver="liblinear", random_state=seeds[0])
        classifier.fit(fingerprints[pair_mol[train]], labels[train])
        ecfp[test] = classifier.predict_proba(fingerprints[pair_mol[test]])[:, 1]
        for seed_index, seed in enumerate(seeds):
            for random_labels, output in [(False, tuned_seeds[seed_index]), (True, random_seeds[seed_index])]:
                model = fit(mol, pocket, pair_mol[train], pair_target[train], labels[train], seen,
                            initial, seed + held * 1009, a, random_labels)
                with torch.inference_mode(): matrix = model.mol.from_hidden(model.mol.hidden(mol)) @ model.pocket.from_hidden(model.pocket.hidden(pocket)).T
                output[test] = matrix.numpy()[pair_mol[test], pair_target[test]]
    tuned = np.mean(tuned_seeds, axis=0); random = np.mean(random_seeds, axis=0)
    metric_names = ["roc_auc", "pr_auc", "bedroc_alpha20", "ef1pct", "ef5pct"]
    groups = np.unique(scaffolds); group_idx = {g: np.flatnonzero(scaffolds == g) for g in groups}
    rng = np.random.default_rng(20260926); boot = {name: [] for name in metric_names}
    for _ in range(a.bootstrap):
        idx = np.concatenate([group_idx[g] for g in rng.choice(groups, len(groups), replace=True)])
        try:
            sample = pairs.iloc[idx].reset_index(drop=True)
            left = screening_metrics(sample, tuned[idx])["macro"]
            right = screening_metrics(sample, official[idx])["macro"]
        except ValueError:
            continue
        for name in metric_names: boot[name].append(left[name] - right[name])
    ci = {name: list(map(float, np.quantile(values, [0.025, 0.5, 0.975]))) for name, values in boot.items()}
    report = {"protocol": "strict leave-one-entire-GPCR-target-out; held-target scaffolds purged; held pocket absent from all losses",
              "metrics": {"official_drugclip": screening_metrics(pairs, official),
                          "pooled_ecfp4_logistic_transfer": screening_metrics(pairs, ecfp),
                          "three_target_bce_retrieval_transfer": screening_metrics(pairs, tuned),
                          "three_target_random_label_control": screening_metrics(pairs, random)},
              "tuned_minus_official_scaffold_bootstrap_95ci": ci,
              "per_seed_tuned_metrics": [screening_metrics(pairs, value)["macro"] for value in tuned_seeds],
              "training_audit": audit,
              "claim_boundary": "Retrospective four-target transfer audit; no broad GPCR, affinity, efficacy, or PAM claim.",
              "hyperparameters": {k: getattr(a, k) for k in ["rank", "epochs", "lr", "temperature", "retrieval_weight", "preserve_weight", "seeds", "bootstrap"]}}
    a.output.parent.mkdir(parents=True, exist_ok=True); a.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    out = pairs[["pair_id", "target", "label", "murcko_scaffold"]].copy(); out["official"] = official; out["ecfp"] = ecfp; out["tuned"] = tuned; out["random"] = random
    out.to_csv(a.output.with_suffix(".predictions.csv"), index=False); print(json.dumps(report, indent=2))


if __name__ == "__main__": main()
