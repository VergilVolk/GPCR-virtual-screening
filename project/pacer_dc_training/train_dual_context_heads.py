#!/usr/bin/env python
"""Train PACER-DC PAM-synergy and intrinsic-agonism heads on frozen embeddings.

Input rows are four matched contexts for one candidate/replica/window. Windows
are augmentation only: losses and metrics are aggregated to unique molecules.
The script refuses candidate or chemotype leakage across train/val/test splits.
"""
from __future__ import annotations

import argparse
import json
import random
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from sklearn.metrics import average_precision_score, balanced_accuracy_score, roc_auc_score
from torch import nn
from torch.nn import functional as F

CONTEXTS = {"candidate_probe", "candidate_no_probe", "probe_only", "apo"}


class DualHead(nn.Module):
    def __init__(self, n_features: int, hidden: int, projection: int):
        super().__init__()
        self.pam_encoder = nn.Sequential(nn.Linear(n_features, hidden), nn.ReLU(), nn.LayerNorm(hidden))
        self.ago_encoder = nn.Sequential(nn.Linear(n_features, hidden), nn.ReLU(), nn.LayerNorm(hidden))
        self.pam_projection = nn.Linear(hidden, projection)
        self.ago_projection = nn.Linear(hidden, projection)
        self.pam_classifier = nn.Linear(hidden, 1)
        self.ago_classifier = nn.Linear(hidden, 1)

    def forward(self, d_pam: torch.Tensor, d_ago: torch.Tensor):
        hp = self.pam_encoder(d_pam)
        ha = self.ago_encoder(d_ago)
        return {
            "pam_embedding": F.normalize(self.pam_projection(hp), dim=1),
            "ago_embedding": F.normalize(self.ago_projection(ha), dim=1),
            "pam_logit": self.pam_classifier(hp).squeeze(1),
            "ago_logit": self.ago_classifier(ha).squeeze(1),
        }


def read_config(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def validate_and_build(frame: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    required = {"candidate_id", "chemotype", "replicate_id", "window_id", "context",
                "pam_label", "agonism_label", "split"}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")
    feature_cols = sorted(c for c in frame.columns if c.startswith("f_"))
    if not feature_cols:
        raise ValueError("No embedding columns found; expected f_000, f_001, ...")
    if not set(frame.context).issubset(CONTEXTS):
        raise ValueError(f"Unknown contexts: {sorted(set(frame.context) - CONTEXTS)}")
    if not set(frame.split).issubset({"train", "val", "test"}):
        raise ValueError("split must contain only train, val, or test")
    for column in ["candidate_id", "chemotype"]:
        leakage = frame.groupby(column).split.nunique()
        leakage = leakage[leakage > 1]
        if len(leakage):
            raise ValueError(f"{column} leakage across splits: {list(leakage.index[:10])}")

    keys = ["candidate_id", "chemotype", "replicate_id", "window_id", "split"]
    rows = []
    for key, part in frame.groupby(keys, sort=False):
        seen = set(part.context)
        if seen != CONTEXTS or part.context.duplicated().any():
            raise ValueError(f"Incomplete/duplicate four-context unit {key}: {sorted(seen)}")
        vectors = {r.context: r[feature_cols].to_numpy(np.float32) for _, r in part.iterrows()}
        pam_labels = set(part.pam_label.dropna().astype(int))
        ago_labels = set(part.agonism_label.dropna().astype(int))
        if len(pam_labels) > 1 or len(ago_labels) > 1:
            raise ValueError(f"Conflicting labels within {key}")
        rows.append({
            **dict(zip(keys, key)),
            "pam_label": next(iter(pam_labels)) if pam_labels else np.nan,
            "agonism_label": next(iter(ago_labels)) if ago_labels else np.nan,
            "d_pam": vectors["candidate_probe"] - vectors["probe_only"],
            "d_ago": vectors["candidate_no_probe"] - vectors["apo"],
        })
    units = pd.DataFrame(rows)
    if not {"train", "val", "test"}.issubset(set(units.split)):
        raise ValueError("All three splits train/val/test are required")
    return units, feature_cols


def candidate_means(values: torch.Tensor, ids: list[str]) -> tuple[torch.Tensor, list[str]]:
    ordered = list(dict.fromkeys(ids))
    means = [values[[i for i, cid in enumerate(ids) if cid == candidate]].mean(0) for candidate in ordered]
    return torch.stack(means), ordered


def candidate_labels(units: pd.DataFrame, ids: list[str], column: str) -> torch.Tensor:
    values = []
    for candidate in ids:
        labels = units.loc[units.candidate_id.eq(candidate), column].dropna().astype(int).unique()
        if len(labels) != 1:
            raise ValueError(f"Candidate {candidate} needs one consistent {column}")
        values.append(int(labels[0]))
    return torch.tensor(values, dtype=torch.float32)


def cross_chemotype_triplet(embedding: torch.Tensor, labels: torch.Tensor,
                            candidates: list[str], chemotype: dict[str, str], margin: float) -> torch.Tensor:
    terms = []
    for a in range(len(candidates)):
        if labels[a].item() != 1:
            continue
        positives = [p for p in range(len(candidates)) if p != a and labels[p].item() == 1
                     and chemotype[candidates[p]] != chemotype[candidates[a]]]
        negatives = [n for n in range(len(candidates)) if labels[n].item() == 0]
        for p in positives:
            for n in negatives:
                dap = 1 - torch.dot(embedding[a], embedding[p])
                dan = 1 - torch.dot(embedding[a], embedding[n])
                terms.append(F.relu(margin + dap - dan))
    if not terms:
        return embedding.sum() * 0.0
    return torch.stack(terms).mean()


def metrics(y: np.ndarray, p: np.ndarray) -> dict:
    out = {"n": int(len(y))}
    if len(np.unique(y)) == 2:
        out.update({
            "roc_auc": float(roc_auc_score(y, p)),
            "pr_auc": float(average_precision_score(y, p)),
            "balanced_accuracy_0p5": float(balanced_accuracy_score(y, p >= 0.5)),
        })
    else:
        out["not_evaluable_reason"] = "only one class"
    return out


def validation_bce(model: DualHead, units: pd.DataFrame, device: torch.device) -> float:
    """Candidate-level validation objective; windows never get independent weight."""
    model.eval()
    x_pam = torch.tensor(np.stack(units.d_pam), device=device)
    x_ago = torch.tensor(np.stack(units.d_ago), device=device)
    with torch.no_grad():
        out = model(x_pam, x_ago)
        pam_logit, ids = candidate_means(out["pam_logit"].unsqueeze(1), units.candidate_id.tolist())
        ago_logit, ids2 = candidate_means(out["ago_logit"].unsqueeze(1), units.candidate_id.tolist())
        if ids != ids2:
            raise RuntimeError("Validation candidate aggregation mismatch")
        yp = candidate_labels(units, ids, "pam_label").to(device)
        ya = candidate_labels(units, ids, "agonism_label").to(device)
        loss = F.binary_cross_entropy_with_logits(pam_logit.squeeze(1), yp)
        loss = loss + F.binary_cross_entropy_with_logits(ago_logit.squeeze(1), ya)
    return float(loss.cpu())


def evaluate(model: DualHead, units: pd.DataFrame, device: torch.device) -> tuple[dict, pd.DataFrame]:
    model.eval()
    x_pam = torch.tensor(np.stack(units.d_pam), device=device)
    x_ago = torch.tensor(np.stack(units.d_ago), device=device)
    with torch.no_grad():
        out = model(x_pam, x_ago)
    rows = []
    for candidate in units.candidate_id.unique():
        idx = np.flatnonzero(units.candidate_id.to_numpy() == candidate)
        rows.append({
            "candidate_id": candidate,
            "chemotype": units.iloc[idx[0]].chemotype,
            "split": units.iloc[idx[0]].split,
            "pam_label": int(units.iloc[idx[0]].pam_label),
            "agonism_label": int(units.iloc[idx[0]].agonism_label),
            "pam_probability": float(torch.sigmoid(out["pam_logit"][idx]).mean().cpu()),
            "agonism_probability": float(torch.sigmoid(out["ago_logit"][idx]).mean().cpu()),
            "n_windows": int(len(idx)),
        })
    pred = pd.DataFrame(rows)
    summary = {}
    for split in ["val", "test"]:
        part = pred[pred.split.eq(split)]
        summary[split] = {
            "pam": metrics(part.pam_label.to_numpy(), part.pam_probability.to_numpy()),
            "agonism": metrics(part.agonism_label.to_numpy(), part.agonism_probability.to_numpy()),
        }
    return summary, pred


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--features", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    cfg = read_config(args.config)
    seed = int(cfg["seed"])
    random.seed(seed); np.random.seed(seed); torch.manual_seed(seed)
    torch.set_num_threads(max(1, min(8, torch.get_num_threads())))
    device = torch.device("cpu")

    units, feature_cols = validate_and_build(pd.read_csv(args.features))
    train = units[units.split.eq("train")].reset_index(drop=True)
    validation = units[units.split.eq("val")].reset_index(drop=True)
    if train.candidate_id.nunique() < 6:
        raise ValueError("Training gate failed: fewer than 6 independent training molecules")
    for label in ["pam_label", "agonism_label"]:
        counts = train.drop_duplicates("candidate_id")[label].value_counts()
        if set(counts.index) != {0, 1} or counts.min() < 2:
            raise ValueError(f"Training gate failed for {label}: need >=2 independent molecules per class")

    x_pam = torch.tensor(np.stack(train.d_pam), device=device)
    x_ago = torch.tensor(np.stack(train.d_ago), device=device)
    model = DualHead(len(feature_cols), int(cfg["hidden_dim"]), int(cfg["projection_dim"])).to(device)
    optimizer = torch.optim.AdamW(model.parameters(), lr=float(cfg["learning_rate"]),
                                  weight_decay=float(cfg["weight_decay"]))
    chemotype = dict(zip(train.candidate_id, train.chemotype))
    best_state, best_val_loss, stale = None, float("inf"), 0
    history = []
    for epoch in range(int(cfg["epochs"])):
        model.train(); optimizer.zero_grad()
        out = model(x_pam, x_ago)
        pam_emb, ids = candidate_means(out["pam_embedding"], train.candidate_id.tolist())
        ago_emb, ids2 = candidate_means(out["ago_embedding"], train.candidate_id.tolist())
        pam_logit, _ = candidate_means(out["pam_logit"].unsqueeze(1), train.candidate_id.tolist())
        ago_logit, _ = candidate_means(out["ago_logit"].unsqueeze(1), train.candidate_id.tolist())
        if ids != ids2:
            raise RuntimeError("Candidate aggregation mismatch")
        yp = candidate_labels(train, ids, "pam_label").to(device)
        ya = candidate_labels(train, ids, "agonism_label").to(device)
        bce = F.binary_cross_entropy_with_logits(pam_logit.squeeze(1), yp)
        bce = bce + F.binary_cross_entropy_with_logits(ago_logit.squeeze(1), ya)
        tp = cross_chemotype_triplet(pam_emb, yp, ids, chemotype, float(cfg["triplet_margin"]))
        ta = cross_chemotype_triplet(ago_emb, ya, ids, chemotype, float(cfg["triplet_margin"]))
        loss = bce + float(cfg["pam_triplet_weight"]) * tp + float(cfg["ago_triplet_weight"]) * ta
        loss.backward(); optimizer.step()
        value = float(loss.detach())
        val_loss = validation_bce(model, validation, device)
        history.append({"epoch": epoch, "loss": value, "bce": float(bce.detach()),
                        "pam_triplet": float(tp.detach()), "ago_triplet": float(ta.detach()),
                        "validation_bce": val_loss})
        if val_loss < best_val_loss - 1e-6:
            best_val_loss = val_loss
            best_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}
            stale = 0
        else:
            stale += 1
        if stale >= int(cfg["early_stopping_patience"]):
            break
    if best_state is None:
        raise RuntimeError("No model state produced")
    model.load_state_dict(best_state)
    summary, predictions = evaluate(model, units, device)

    args.output.mkdir(parents=True, exist_ok=True)
    torch.save({"model_state": best_state, "config": cfg, "feature_columns": feature_cols},
               args.output / "model.pt")
    pd.DataFrame(history).to_csv(args.output / "training_history.csv", index=False)
    predictions.to_csv(args.output / "candidate_predictions.csv", index=False)
    audit = {
        "model_id": cfg["model_id"], "device": "cpu", "seed": seed,
        "n_embedding_features": len(feature_cols),
        "n_four_context_units": int(len(units)),
        "candidate_counts": units.drop_duplicates("candidate_id").split.value_counts().to_dict(),
        "epochs_completed": len(history), "best_validation_bce": best_val_loss,
        "metrics": summary,
        "candidate_level_metrics": True,
        "windows_are_augmentation_not_independent_samples": True,
        "claim_boundary": cfg["claim_boundary"],
    }
    (args.output / "audit.json").write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(audit, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
