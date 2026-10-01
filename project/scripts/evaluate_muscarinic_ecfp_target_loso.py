#!/usr/bin/env python3
"""ECFP-to-frozen-pocket baseline under strict muscarinic target LOSO."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from torch.nn import functional as F

from evaluate_muscarinic_ecfp_pocket_baseline import ECFPPocket, fps
from finetune_drugclip_muscarinic_triplet import clustered_bootstrap


TARGETS = ["M1", "M2", "M3", "M4", "M5"]


def fit_model(x, triplets, pocket, mi, pi, seed, args, random_orientation=False):
    torch.manual_seed(seed)
    model = ECFPPocket(seed)
    opt = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    mol = torch.as_tensor([mi[str(v)] for v in triplets.canonical_smiles], dtype=torch.long)
    pos = torch.as_tensor([pi[str(v)] for v in triplets.positive_subtype], dtype=torch.long)
    neg = torch.as_tensor([pi[str(v)] for v in triplets.negative_subtype], dtype=torch.long)
    if random_orientation:
        swap = torch.rand(len(triplets)) < 0.5
        old = pos.clone(); pos[swap] = neg[swap]; neg[swap] = old[swap]
    for _ in range(args.epochs):
        opt.zero_grad()
        score = model(x) @ pocket.T
        loss = F.relu(args.margin - score[mol, pos] + score[mol, neg]).mean()
        loss.backward(); torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0); opt.step()
    return model


def pair_scores(model, x, pocket, frame, mi, pi):
    with torch.inference_mode():
        matrix = (model(x) @ pocket.T).numpy()
    return np.asarray([matrix[mi[str(r.canonical_smiles)], pi[str(r.positive_subtype)]] -
                       matrix[mi[str(r.canonical_smiles)], pi[str(r.negative_subtype)]]
                       for r in frame.itertuples()])


def metrics(frame, score, held):
    preferred = frame.positive_subtype.to_numpy() == held
    disfavored = frame.negative_subtype.to_numpy() == held
    a = float((score[preferred] > 0).mean()); b = float((score[disfavored] > 0).mean())
    return {"accuracy": float((score > 0).mean()), "held_preferred_accuracy": a,
            "held_disfavored_accuracy": b, "direction_balanced_accuracy": (a + b) / 2}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--representations", type=Path, required=True)
    p.add_argument("--benchmark-dir", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--seeds", default="20260924,20260925,20260926")
    p.add_argument("--epochs", type=int, default=300)
    p.add_argument("--lr", type=float, default=1e-3)
    p.add_argument("--weight-decay", type=float, default=1e-2)
    p.add_argument("--margin", type=float, default=0.1)
    p.add_argument("--fingerprint-radius", type=int, default=2)
    p.add_argument("--skip-random-control", action="store_true")
    args = p.parse_args(); torch.set_num_threads(max(1, min(8, torch.get_num_threads())))

    arc = np.load(args.representations, allow_pickle=False)
    ids = list(map(str, arc["molecule_ids"])); mi = {v: i for i, v in enumerate(ids)}
    pids = [str(v).split("_")[0] for v in arc["pocket_ids"]]; pi = {v: i for i, v in enumerate(pids)}
    x = fps(ids, radius=args.fingerprint_radius); pocket = F.normalize(torch.as_tensor(arc["pocket_embeddings"].astype(np.float32)), dim=-1)
    reports, rows = {}, []
    for held in TARGETS:
        train_raw = pd.read_csv(args.benchmark_dir / f"{held}_train.csv")
        test_raw = pd.read_csv(args.benchmark_dir / f"{held}_test.csv")
        train = train_raw[train_raw.canonical_smiles.astype(str).isin(mi)].reset_index(drop=True)
        test = test_raw[test_raw.canonical_smiles.astype(str).isin(mi)].reset_index(drop=True)
        targeted, controls = [], []
        for seed in map(int, args.seeds.split(",")):
            model = fit_model(x, train, pocket, mi, pi, seed, args, False)
            ts = pair_scores(model, x, pocket, test, mi, pi)
            if args.skip_random_control:
                cs = np.full(len(test), np.nan)
            else:
                control = fit_model(x, train, pocket, mi, pi, seed, args, True)
                cs = pair_scores(control, x, pocket, test, mi, pi)
                controls.append(cs)
            targeted.append(ts)
            for i, row in test.iterrows():
                rows.append({"held_subtype": held, "canonical_smiles": row.canonical_smiles,
                             "positive_subtype": row.positive_subtype, "negative_subtype": row.negative_subtype,
                             "seed": seed, "random_control_delta": float(cs[i]), "targeted_delta": float(ts[i])})
        target = np.mean(targeted, axis=0)
        control = np.mean(controls, axis=0) if controls else None
        reports[held] = {
            "train_pairs": int(len(train)), "test_pairs": int(len(test)),
            "test_molecules": int(test.canonical_smiles.nunique()),
            "targeted": metrics(test, target, held),
            "random_control": metrics(test, control, held) if control is not None else None,
            "targeted_minus_random_molecule_cluster_95ci": clustered_bootstrap(test, control, target, n=1000) if control is not None else None,
            "per_seed_targeted_accuracy": [float((v > 0).mean()) for v in targeted],
        }
    names = ["targeted"] + ([] if args.skip_random_control else ["random_control"])
    macro = {name: {key: float(np.mean([reports[h][name][key] for h in TARGETS]))
                    for key in ["accuracy", "direction_balanced_accuracy"]}
             for name in names}
    report = {"method": f"ECFP radius-{args.fingerprint_radius} linear projection into frozen DrugCLIP pocket space",
              "split": "strict leave-one-muscarinic-target-out", "folds": reports,
              "macro_average": macro,
              "hyperparameters": {"epochs": args.epochs, "lr": args.lr, "margin": args.margin,
                                      "seeds": args.seeds},
              "claim_boundary": "Muscarinic subtype activity transfer; not PAM efficacy or binding affinity."}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    pd.DataFrame(rows).to_csv(args.output.with_suffix(".predictions.csv"), index=False)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
