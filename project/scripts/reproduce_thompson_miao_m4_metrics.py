# -*- coding: utf-8 -*-
"""Exact M4 metric reproduction using the authors' released score tables.

Outputs the source-code metrics alongside Table 4 values transcribed from the
2026-08-18 bioRxiv PDF, and flags differences above rounding tolerance.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd


P = Path(__file__).resolve().parents[1]
REPO = P / "tools" / "gpcr-am-ensemble-docking"
SCORES = REPO / "docking_scores" / "M4R"
OUT = P / "results" / "m4_official_exact_reproduction"
OUT.mkdir(parents=True, exist_ok=True)


def module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    obj = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj


AUC = module(REPO / "scripts" / "auc_calculator.py", "official_auc")
ENRICH = module(REPO / "scripts" / "enrichment_PDB.py", "official_enrich")


FILES = {
    "Glide_PDB": (SCORES / "PDB" / "M4R_PDB_Glide_scores.csv", "glide_gscore"),
    "Glide_BEmin": (SCORES / "Ensemble" / "M4R_Ensemble_Glide_BEmin_ranked.csv", "BE_min"),
    "Glide_BEavg": (SCORES / "Ensemble" / "M4R_Ensemble_Glide_BEavg_ranked.csv", "BE_avg"),
    "Vina_PDB": (SCORES / "PDB" / "M4R_PDB_Vina_scores.csv", "vina_score"),
    "Vina_BEmin": (SCORES / "Ensemble" / "M4R_Ensemble_Vina_BEmin_ranked.csv", "BE_min"),
    "Vina_BEavg": (SCORES / "Ensemble" / "M4R_Ensemble_Vina_BEavg_ranked.csv", "BE_avg"),
}


PAPER = {
    "Glide_PDB": dict(EF05=11.499, EFp05=10.374, EF1=8.646, EFp1=10.251, AUC=67.40, logAUC=11.61),
    "Glide_BEmin": dict(EF05=20.370, EFp05=21.768, EF1=13.692, EFp1=19.187, AUC=70.04, logAUC=16.10),
    "Glide_BEavg": dict(EF05=6.128, EFp05=5.808, EF1=5.097, EFp1=5.548, AUC=71.83, logAUC=11.63),
    "Vina_PDB": dict(EF05=.778, EFp05=.589, EF1=1.512, EFp1=1.236, AUC=72.53, logAUC=10.62),
    "Vina_BEmin": dict(EF05=.261, EFp05=.199, EF1=.784, EFp1=.587, AUC=73.38, logAUC=10.18),
    "Vina_BEavg": dict(EF05=.959, EFp05=.766, EF1=1.482, EFp1=1.136, AUC=69.58, logAUC=8.97),
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    rows = []
    for method, (path, score_col) in FILES.items():
        d = pd.read_csv(path)
        y = d.ligand_id.astype(str).str.startswith("ASD").astype(int).to_numpy()
        raw = d[score_col].to_numpy(float)
        pct_d, pct_a = AUC.roc_curve(y, -raw)
        auc = AUC.calc_auc(pct_d, pct_a) * 100
        logauc = AUC.calc_logauc(pct_d, pct_a) * 100
        sorted_d = ENRICH.load_and_sort(str(path), "ligand_id", score_col, True)
        n, hits = ENRICH.get_compound_counts(sorted_d, "ligand_id", "ASD")
        e = ENRICH.compute_ef(sorted_d, "ligand_id", "ASD", n, hits, [.5, 1.0])
        ep = ENRICH.compute_ef_prime(sorted_d, "ligand_id", "ASD", n, hits, [.5, 1.0])
        got = dict(
            EF05=e[.5]["EF"], EFp05=ep[.5]["EF_prime"],
            EF1=e[1.0]["EF"], EFp1=ep[1.0]["EF_prime"], AUC=auc, logAUC=logauc,
        )
        row = {"method": method, "n": n, "actives": hits, "sha256": sha(path)}
        for metric, value in got.items():
            row[f"code_{metric}"] = value
            row[f"paper_{metric}"] = PAPER[method][metric]
            row[f"delta_{metric}"] = value - PAPER[method][metric]
        rows.append(row)

    out = pd.DataFrame(rows)
    out.to_csv(OUT / "m4_table4_exact_reproduction.csv", index=False)
    discrepancies = []
    for _, row in out.iterrows():
        for metric in ("EF05", "EFp05", "EF1", "EFp1", "AUC", "logAUC"):
            tolerance = .05 if metric in {"AUC", "logAUC"} else .02
            if abs(row[f"delta_{metric}"]) > tolerance:
                discrepancies.append({
                    "method": row.method, "metric": metric,
                    "code": float(row[f"code_{metric}"]),
                    "paper": float(row[f"paper_{metric}"]),
                    "delta": float(row[f"delta_{metric}"]),
                })
    audit = {
        "paper": "Thompson and Miao, bioRxiv, posted 2026-08-18",
        "doi": "10.64898/2026.08.12.744492",
        "repository_commit": "44798c841ee77230b1f89fe41074e41b458c7070",
        "official_analysis_code_executed": True,
        "n_methods": len(out),
        "discrepancies_above_rounding_tolerance": discrepancies,
        "reproduction_scope": (
            "Exact re-analysis of released docking scores. Raw Glide docking and 3x500 ns GaMD production "
            "were not regenerated locally because they require a commercial license and GPU/HPC trajectories."
        ),
    }
    (OUT / "reproduction_audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(out[["method", "n", "actives", "code_EF05", "code_EF1", "code_AUC", "code_logAUC"]].to_string(index=False))
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
