#!/usr/bin/env python3
"""One-shot external benchmark evaluation of source-trained CGDA checkpoints."""
from __future__ import annotations

import argparse
import json
import pickle
from pathlib import Path

import numpy as np
import torch

from evaluate_drugclip_cgda_loso import CGDA
from evaluate_drugclip_official_litpcba_embeddings import metrics


METRIC_NAMES = ["roc_auc", "pr_auc", "bedroc_alpha80_5", "ef0.005", "ef0.01", "ef0.02", "ef0.05"]


def load_contexts(path: Path):
    archive = np.load(path, allow_pickle=False)
    targets = archive["targets"].astype(str).tolist()
    p, po = archive["pocket_embeddings"], archive["pocket_offsets"]
    r, ro = archive["reference_embeddings"], archive["reference_offsets"]
    return {target: {"pockets": p[po[i] : po[i + 1]], "references": r[ro[i] : ro[i + 1]]}
            for i, target in enumerate(targets)}


def load_model(path: Path):
    payload = torch.load(path, map_location="cpu", weights_only=False)
    architecture = payload["architecture"]
    model = CGDA(architecture["dimension"], architecture["rank"], architecture["experts"], 0)
    model.load_state_dict(payload["model_state"])
    model.eval()
    return model, payload


def bootstrap(per_target, iterations, seed):
    targets = sorted(per_target)
    delta = np.asarray([[per_target[t]["cgda"][name] - per_target[t]["reference"][name]
                         for name in METRIC_NAMES] for t in targets])
    rng = np.random.default_rng(seed)
    samples = np.empty((iterations, len(METRIC_NAMES)))
    for index in range(iterations):
        samples[index] = delta[rng.integers(0, len(targets), len(targets))].mean(axis=0)
    return {name: {"mean_delta": float(delta[:, i].mean()),
                   "ci95": [float(v) for v in np.quantile(samples[:, i], [0.025, 0.975])],
                   "targets_improved": int((delta[:, i] > 0).sum()),
                   "targets_total": len(targets)} for i, name in enumerate(METRIC_NAMES)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--contexts", type=Path, required=True)
    parser.add_argument("--checkpoints", type=Path, nargs="+", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--candidate-pickle", default="drugclip_emb/mols.lmdb.pkl")
    parser.add_argument("--bootstrap", type=int, default=20000)
    parser.add_argument("--seed", type=int, default=20260926)
    args = parser.parse_args()

    contexts = load_contexts(args.contexts)
    loaded = [load_model(path) for path in args.checkpoints]
    models = [item[0] for item in loaded]
    training_targets = [item[1]["training_targets"] for item in loaded]
    per_target = {}
    score_dir = args.output.with_suffix("")
    score_dir.mkdir(parents=True, exist_ok=True)
    for index, target in enumerate(sorted(contexts)):
        path = args.root / target / args.candidate_pickle
        molecules, _, labels = pickle.load(path.open("rb"))
        molecules = np.asarray(molecules, dtype=np.float32)
        labels = np.asarray(labels, dtype=np.int64)
        context = contexts[target]
        pocket_score = (context["pockets"] @ molecules.T).max(axis=0)
        reference_score = (context["references"] @ molecules.T).max(axis=0)
        seed_scores = []
        with torch.no_grad():
            pockets = torch.as_tensor(context["pockets"], dtype=torch.float32)
            references = torch.as_tensor(context["references"], dtype=torch.float32)
            molecule_tensor = torch.as_tensor(molecules, dtype=torch.float32)
            for model in models:
                adapted, _ = model.adapted_references(pockets, references)
                seed_scores.append((molecule_tensor @ adapted.T).max(dim=1).values.numpy())
        cgda_score = np.mean(seed_scores, axis=0)
        per_target[target] = {"n": int(len(labels)), "positives": int(labels.sum()),
                              "pocket": metrics(labels, pocket_score),
                              "reference": metrics(labels, reference_score),
                              "cgda": metrics(labels, cgda_score)}
        np.savez_compressed(score_dir / f"{target}.npz", labels=labels, pocket=pocket_score,
                            reference=reference_score, cgda=cgda_score)
        print(f"completed {index + 1}/{len(contexts)} {target}", flush=True)

    macro = {method: {name: float(np.mean([per_target[t][method][name] for t in per_target]))
                      for name in METRIC_NAMES} for method in ("pocket", "reference", "cgda")}
    report = {
        "protocol": "Frozen source-trained CGDA applied once to an external benchmark",
        "source_training_targets": training_targets,
        "external_targets": len(per_target),
        "macro": macro,
        "paired_cgda_vs_reference": bootstrap(per_target, args.bootstrap, args.seed),
        "per_target": per_target,
        "parameters": {key: [str(x) for x in value] if isinstance(value, list) else
                       str(value) if isinstance(value, Path) else value for key, value in vars(args).items()},
        "claim_boundary": "No external labels used for fitting or selection; benchmark-specific decoy bias remains possible.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps({"macro": macro, "paired": report["paired_cgda_vs_reference"]}, indent=2))


if __name__ == "__main__":
    main()
