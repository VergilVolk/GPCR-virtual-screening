#!/usr/bin/env python3
"""GaMD population-aware ensemble DrugCLIP fine-tuning for GPCR screening.

Ten published cluster representatives per target are embedded by DrugCLIP.
The proposed aggregator is a population-prior log-sum-exp over conformations:

  S(m,t) = tau_t * log sum_c pi_tc * exp(s(m,t,c) / tau_t)

It continuously interpolates between max-like and population-average scoring.
Tau is learned only on each outer fold's training data.  Encoders stay frozen;
rank-r residuals update the two final DrugCLIP projections.
"""
from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from torch import nn
from torch.nn import functional as F

from finetune_drugclip_gpcr_screening import (
    TARGETS, FastProjectionLoRA, balanced_weights, screening_metrics,
)
from finetune_drugclip_muscarinic_triplet import Proj


class EnsembleAdapter(nn.Module):
    def __init__(self, initial, rank, seed, initial_tau):
        super().__init__()
        self.molecule = FastProjectionLoRA(initial["mol_project"], rank, seed)
        self.pocket = FastProjectionLoRA(initial["pocket_project"], rank, seed + 17)
        raw = np.log(np.expm1(initial_tau))
        self.raw_tau = nn.Parameter(torch.full((len(TARGETS),), float(raw)))
        # Per-target corrections to the GaMD population prior.  Zero means the
        # published populations are used exactly; only 4 x 10 parameters are
        # added and their deviation is explicitly regularised.
        self.population_delta = nn.Parameter(torch.zeros(len(TARGETS), 10))
        # Evidence gate: dynamic evidence is an optional residual over the
        # cluster-0 score, rather than an unconditional replacement.
        self.raw_gate = nn.Parameter(torch.full((len(TARGETS),), -2.1972246))  # sigmoid = 0.1

    def tau(self):
        return F.softplus(self.raw_tau) + 1e-4

    def gate(self):
        return torch.sigmoid(self.raw_gate)


def corrected_populations(populations, delta):
    return torch.softmax(populations.clamp_min(1e-8).log() + delta, dim=1)


def aggregate(cluster_scores, populations, pair_targets, mode, tau=None, population_delta=None):
    # cluster_scores: pair x cluster; populations: target x cluster
    if mode == "cluster0":
        return cluster_scores[:, 0]
    if mode == "population_mean":
        return (cluster_scores * populations[pair_targets]).sum(1)
    if mode == "max":
        return cluster_scores.max(1).values
    if mode in {"population_lse", "adaptive_population_lse"}:
        pair_tau = tau[pair_targets]
        effective = (populations if population_delta is None
                     else corrected_populations(populations, population_delta))
        log_pi = effective[pair_targets].clamp_min(1e-8).log()
        return pair_tau * torch.logsumexp(log_pi + cluster_scores / pair_tau[:, None], dim=1)
    raise ValueError(mode)


def all_cluster_scores(zm, zp, pair_mol, pair_target):
    selected_molecule = zm[pair_mol]
    selected_pockets = zp.view(len(TARGETS), 10, -1)[pair_target]
    return torch.einsum("pd,pcd->pc", selected_molecule, selected_pockets)


def fit(mode, train_pair_mol, train_pair_target, train_labels, molecule_rep, pocket_rep,
        populations, initial, seed, args, random_labels=False):
    model = EnsembleAdapter(initial, args.rank, seed, args.initial_tau)
    with torch.no_grad():
        hm = model.molecule.hidden(molecule_rep); hp = model.pocket.hidden(pocket_rep)
        z0m = model.molecule.base_from_hidden(hm); z0p = model.pocket.base_from_hidden(hp)
    labels = train_labels.copy()
    if random_labels:
        rng = np.random.default_rng(seed)
        for target in range(len(TARGETS)):
            keep = np.flatnonzero(train_pair_target == target)
            labels[keep] = rng.permutation(labels[keep])
    pair_mol = torch.as_tensor(train_pair_mol, dtype=torch.long)
    pair_target = torch.as_tensor(train_pair_target, dtype=torch.long)
    # Keep the outer test fold completely unseen by the train-time
    # preservation term; otherwise the evaluation would be transductive.
    preserve_mol = torch.as_tensor(np.unique(train_pair_mol), dtype=torch.long)
    preserve_target = torch.as_tensor(
        np.concatenate([np.arange(target * 10, target * 10 + 10)
                        for target in np.unique(train_pair_target)]),
        dtype=torch.long,
    )
    y = torch.as_tensor(labels, dtype=torch.float32)
    weights = torch.as_tensor(balanced_weights(train_pair_target, labels))
    optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    for _ in range(args.epochs):
        optimizer.zero_grad()
        zm = model.molecule.from_hidden(hm); zp = model.pocket.from_hidden(hp)
        clusters = all_cluster_scores(zm, zp, pair_mol, pair_target)
        adaptive_mode = mode in {"adaptive_population_lse", "gated_adaptive_lse"}
        delta = model.population_delta if adaptive_mode else None
        aggregate_mode = "adaptive_population_lse" if adaptive_mode else mode
        ensemble_scores = aggregate(
            clusters, populations, pair_target, aggregate_mode, model.tau(), delta
        )
        if mode == "gated_adaptive_lse":
            gate = model.gate()[pair_target]
            scores = clusters[:, 0] + gate * (ensemble_scores - clusters[:, 0])
        else:
            scores = ensemble_scores
        scores = scores / args.temperature
        bce = (F.binary_cross_entropy_with_logits(scores, y, reduction="none") * weights).mean()
        preserve = ((1 - (zm[preserve_mol] * z0m[preserve_mol]).sum(1)).mean()
                    + (1 - (zp[preserve_target] * z0p[preserve_target]).sum(1)).mean())
        if adaptive_mode:
            q = corrected_populations(populations, model.population_delta)
            population_kl = (q * (q.clamp_min(1e-8).log()
                                  - populations.clamp_min(1e-8).log())).sum(1).mean()
        else:
            population_kl = torch.zeros((), dtype=bce.dtype)
        if mode == "gated_adaptive_lse":
            positive_counts = torch.stack([(y[pair_target == target] > 0.5).sum()
                                           for target in range(len(TARGETS))]).float().clamp_min(1)
            scarcity = torch.sqrt(positive_counts.max() / positive_counts)
            gate_penalty = (scarcity * model.gate().square()).mean()
        else:
            gate_penalty = torch.zeros((), dtype=bce.dtype)
        loss = (bce + args.preserve_weight * preserve
                + args.population_kl_weight * population_kl
                + args.gate_weight * gate_penalty)
        loss.backward(); torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0); optimizer.step()
    return model


def predict(model, mode, molecule_rep, pocket_rep, pair_mol, pair_target, populations):
    with torch.inference_mode():
        zm = model.molecule.from_hidden(model.molecule.hidden(molecule_rep))
        zp = model.pocket.from_hidden(model.pocket.hidden(pocket_rep))
        pair_mol_t = torch.as_tensor(pair_mol, dtype=torch.long)
        pair_target_t = torch.as_tensor(pair_target, dtype=torch.long)
        clusters = all_cluster_scores(zm, zp, pair_mol_t, pair_target_t)
        adaptive_mode = mode in {"adaptive_population_lse", "gated_adaptive_lse"}
        delta = model.population_delta if adaptive_mode else None
        aggregate_mode = "adaptive_population_lse" if adaptive_mode else mode
        ensemble_scores = aggregate(
            clusters, populations, pair_target_t, aggregate_mode, model.tau(), delta
        )
        if mode == "gated_adaptive_lse":
            gate = model.gate()[pair_target_t]
            return (clusters[:, 0] + gate * (ensemble_scores - clusters[:, 0])).numpy()
        return ensemble_scores.numpy()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--representations", type=Path, required=True)
    parser.add_argument("--projection", type=Path, required=True)
    parser.add_argument("--pairs", type=Path, required=True)
    parser.add_argument("--pocket-metadata", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seeds", default="20260925")
    parser.add_argument("--rank", type=int, default=4)
    parser.add_argument("--epochs", type=int, default=100)
    parser.add_argument("--lr", type=float, default=5e-3)
    parser.add_argument("--weight-decay", type=float, default=1e-3)
    parser.add_argument("--temperature", type=float, default=0.07)
    parser.add_argument("--initial-tau", type=float, default=0.10)
    parser.add_argument("--preserve-weight", type=float, default=0.2)
    parser.add_argument("--population-kl-weight", type=float, default=0.05)
    parser.add_argument("--gate-weight", type=float, default=0.02)
    args = parser.parse_args(); torch.set_num_threads(max(1, min(8, torch.get_num_threads())))

    archive = np.load(args.representations, allow_pickle=False)
    initial = torch.load(args.projection, map_location="cpu")
    pairs_raw = pd.read_csv(args.pairs)
    mol_index = {value: i for i, value in enumerate(map(str, archive["molecule_ids"]))}
    pairs = pairs_raw[pairs_raw.canonical_smiles.astype(str).isin(mol_index)].copy().reset_index(drop=True)
    pair_mol = np.asarray([mol_index[value] for value in pairs.canonical_smiles.astype(str)], dtype=int)
    pair_target = np.asarray([TARGETS.index(value) for value in pairs.target], dtype=int)
    labels = pairs.label.to_numpy(int); folds = pairs.scaffold_fold.to_numpy(int)
    molecule_rep = torch.as_tensor(archive["molecule_representations"].astype(np.float32))
    metadata = pd.read_csv(args.pocket_metadata).sort_values(["target", "cluster"])
    expected_ids = [f"{target}_cluster{cluster}" for target in TARGETS for cluster in range(10)]
    pocket_ids = list(map(str, archive["pocket_ids"]))
    if set(pocket_ids) != set(expected_ids) or len(pocket_ids) != len(expected_ids):
        raise ValueError(f"Pocket ID set mismatch: {pocket_ids[:5]} ...")
    pocket_order = [pocket_ids.index(value) for value in expected_ids]
    pocket_rep = torch.as_tensor(
        archive["pocket_representations"][pocket_order].astype(np.float32)
    )
    metadata = metadata.set_index("pocket_id").loc[expected_ids].reset_index()
    populations = torch.as_tensor(metadata.population.to_numpy(np.float32).reshape(len(TARGETS), 10))
    populations = populations / populations.sum(1, keepdim=True)

    frozen_mol = Proj(copy.deepcopy(initial["mol_project"])).eval()
    frozen_pocket = Proj(copy.deepcopy(initial["pocket_project"])).eval()
    with torch.inference_mode():
        zm0 = frozen_mol(molecule_rep); zp0 = frozen_pocket(pocket_rep)
        clusters0 = all_cluster_scores(zm0, zp0, torch.as_tensor(pair_mol), torch.as_tensor(pair_target))
        fixed_tau = torch.full((len(TARGETS),), args.initial_tau)
        official = {
            "official_cluster0": aggregate(clusters0, populations, torch.as_tensor(pair_target), "cluster0").numpy(),
            "official_population_mean": aggregate(clusters0, populations, torch.as_tensor(pair_target), "population_mean").numpy(),
            "official_max": aggregate(clusters0, populations, torch.as_tensor(pair_target), "max").numpy(),
            "official_population_lse": aggregate(clusters0, populations, torch.as_tensor(pair_target), "population_lse", fixed_tau).numpy(),
        }

    modes = ["cluster0", "population_lse", "adaptive_population_lse", "gated_adaptive_lse"]
    all_predictions = {mode: [] for mode in modes}; random_predictions = []
    tau_records, output_rows = [], []
    for seed in map(int, args.seeds.split(",")):
        predictions = {mode: np.full(len(pairs), np.nan) for mode in modes}
        random = np.full(len(pairs), np.nan)
        for fold in sorted(np.unique(folds)):
            train = folds != fold; test = folds == fold
            for mode in modes:
                model = fit(mode, pair_mol[train], pair_target[train], labels[train], molecule_rep, pocket_rep,
                            populations, initial, seed + 101 * int(fold), args, False)
                predictions[mode][test] = predict(model, mode, molecule_rep, pocket_rep,
                                                   pair_mol[test], pair_target[test], populations)
                tau_records.append({"seed": seed, "fold": int(fold), "mode": mode,
                                    **{f"tau_{target}": float(model.tau().detach()[i])
                                       for i, target in enumerate(TARGETS)},
                                    **{f"population_kl_{target}": float((
                                        corrected_populations(populations, model.population_delta)[i]
                                        * (corrected_populations(populations, model.population_delta)[i].clamp_min(1e-8).log()
                                           - populations[i].clamp_min(1e-8).log())).sum().detach())
                                       for i, target in enumerate(TARGETS)},
                                    **{f"gate_{target}": float(model.gate().detach()[i])
                                       for i, target in enumerate(TARGETS)}})
            random_model = fit("gated_adaptive_lse", pair_mol[train], pair_target[train], labels[train],
                               molecule_rep, pocket_rep, populations, initial,
                               seed + 101 * int(fold), args, True)
            random[test] = predict(random_model, "gated_adaptive_lse", molecule_rep, pocket_rep,
                                   pair_mol[test], pair_target[test], populations)
        for mode in modes: all_predictions[mode].append(predictions[mode])
        random_predictions.append(random)
        for i, row in pairs.iterrows():
            record = {"pair_id": row.pair_id, "seed": seed, "target": row.target,
                      "label": int(row.label), "scaffold_fold": int(row.scaffold_fold),
                      "random_population_lse": float(random[i])}
            record.update({name: float(values[i]) for name, values in official.items()})
            record.update({f"finetuned_{mode}": float(predictions[mode][i]) for mode in modes})
            output_rows.append(record)
    ensemble = {mode: np.mean(values, axis=0) for mode, values in all_predictions.items()}
    random = np.mean(random_predictions, axis=0)
    metrics = {name: screening_metrics(pairs, values) for name, values in official.items()}
    metrics.update({f"finetuned_{mode}": screening_metrics(pairs, values) for mode, values in ensemble.items()})
    metrics["random_population_lse"] = screening_metrics(pairs, random)
    report = {
        "method": "GaMD population-prior log-sum-exp ensemble DrugCLIP with dual-projection rank-4 LoRA",
        "protocol": "Five-fold global Murcko-scaffold OOF; tau learned only on outer-training folds.",
        "pairs_evaluated": int(len(pairs)), "input_pairs_excluded": int(len(pairs_raw) - len(pairs)),
        "metrics": metrics, "tau_by_fold": tau_records,
        "hyperparameters": {"rank": args.rank, "epochs": args.epochs, "lr": args.lr,
                            "temperature": args.temperature, "initial_tau": args.initial_tau,
                            "preserve_weight": args.preserve_weight,
                            "population_kl_weight": args.population_kl_weight,
                            "gate_weight": args.gate_weight,
                            "seeds": args.seeds},
        "claim_boundary": "Retrospective ensemble active/decoy screening; snapshots are structural samples, not a PAM efficacy model.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    pd.DataFrame(output_rows).to_csv(args.output.with_suffix(".predictions.csv"), index=False)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
