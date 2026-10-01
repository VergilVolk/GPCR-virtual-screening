#!/usr/bin/env python3
"""Build the project-wide baseline registry with protocol and execution status."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[2])
    args = p.parse_args()
    root = args.repo_root
    out = root / "project/results/project_wide_integration_benchmark_v02"
    binding = json.loads((out / "binding_recomputed_v02.json").read_text())
    master = json.loads((root / "project/results/project_wide_integration_benchmark_v01/project_wide_benchmark_master_v01.json").read_text())
    rows = []

    def add(module, task, method, family, status, comparable, metric="", value="", role="baseline", source=""):
        rows.append(dict(module=module, task=task, method=method, family=family, status=status,
                         same_protocol_comparable=comparable, primary_metric=metric, primary_value=value,
                         pipeline_role=role, source=source))

    for protocol in ("13T", "20T"):
        for method, family, role in [
            ("official_drugclip", "3D pocket-molecule contrastive retrieval", "binding baseline"),
            ("ecfp4_logistic", "2D fingerprint linear model", "strong ligand-only baseline"),
            ("ep80_score_ensemble", "DrugCLIP dual projection fine-tune", "binding model"),
            ("famaug_score_ensemble", "family-augmented DrugCLIP fine-tune", "13T primary / 20T ablation"),
        ]:
            node = binding["methods"][f"{protocol}/{method}"]["macro"]
            add("binding", f"{protocol} strict target-LOSO", method, family, "RECOMPUTED_RAW", True,
                "macro ROC-AUC", round(node["roc_auc"], 6), role,
                "binding_recomputed_v02.json")
    add("binding", "13T strict target-LOSO", "random-label control", "negative control", "COMMITTED_RESULT", True,
        "macro ROC-AUC", round(master["binding"]["m13_summary"]["random"]["macro"]["roc_auc"], 6), "leakage control", "v01 master")
    add("binding", "LIT-PCBA 15-target LOTO", "CGDA", "context-gated low-rank DrugCLIP adapter", "COMMITTED_RESULT", True,
        "macro ROC-AUC", 0.5788, "external transfer baseline", "DRUGCLIP_CGDA_RESULT.md")
    add("binding", "DUD-E GPCR external", "CGDA", "context-gated low-rank DrugCLIP adapter", "COMMITTED_RESULT", True,
        "macro ROC-AUC", 0.7437, "external stress test", "DRUGCLIP_CGDA_RESULT.md")

    qsar = master["qsar_static"]["baselines"]["pam_vs_inactive"]
    for split in ("scaffold", "source", "series"):
        for method in ("TanimotoKNN", "RandomForest", "ExtraTrees"):
            key = f"{split}|{method}"
            add("2D/QSAR", f"PAM vs inactive {split}-holdout", method, "ECFP classical ML", "COMMITTED_RESULT", True,
                "ROC-AUC", round(qsar[key]["aggregate_oof"]["ROC_AUC"], 6), "generalization baseline", "pacer_baselines_v1")

    structure = master["qsar_static"]["structure_loso_aggregate"]
    for method, node in structure.items():
        add("static structure", "430 compounds leave-one-source-out", method, "docking/IFP/2D fusion", "COMMITTED_RESULT", True,
            "macro Spearman", round(node["macro_Spearman"], 6), "static rescoring baseline", "pacer_structure_loso_v01")
    add("static structure", "external M4 enrichment", "PACER-FS centered", "series-centered static fusion", "COMMITTED_RESULT", True,
        "Spearman delta", round(float(master["qsar_static"]["pacer_fs"]["centered_vs_absolute_spearman"]["estimate"]), 6),
        "static-line best", "pacer_fs_final")

    for row in master["qsar_static"].get("representation_bakeoff", []):
        add("MD representation", "public M4 PAM replica holdout", row["method"], "linear slow-coordinate model", "COMMITTED_RESULT", True,
            "mean lagged R2", round(float(row["mean_lagged_R2"]), 6), "dynamic baseline", "pacer_dynamic_representation_bakeoff_v01")
    for method, value, status, role in [
        ("OneProt-MD final embedding", 0.25, "COMMITTED_NEGATIVE", "failed pretrained baseline"),
        ("Geom2Vec C0-M128 + common kernel", 0.467, "COMMITTED_RESULT", "historical representation baseline"),
        ("C1-BS256 + PACER-FKG v02 compound110", 0.583656, "COMMITTED_RESULT", "functional direction model"),
        ("C1-BS256 + PACER-FKG v02 LY2119620", 0.634304, "COMMITTED_RESULT", "positive-control generalization"),
        ("C1-BS256 + PACER-FKG v02 CM00734", -0.326839, "COMMITTED_RESULT", "hard-negative challenge"),
        ("PACER-MCV physical CV", "", "SOFTWARE_SMOKE_ONLY", "interpretable no-encoder baseline"),
    ]:
        add("functional MD", "four-context DeltaINT direction", method, "trajectory representation", status, method != "PACER-MCV physical CV",
            "R1/R3 direction cosine" if value != "" else "", value, role, "PACER-DC Stage A/B + MCV audit")

    # Relevant published methods that are legitimate future baselines but have not
    # been run under our exact inputs.  They are never plotted as achieved scores.
    for module, method, family, reason, source in [
        ("static structure", "AutoDock Vina", "empirical docking", "RUN_AS_VinaOnly", "https://vina.scripps.edu/"),
        ("static structure", "GNINA 1.3", "CNN docking/rescoring", "NOT_RUN_NO_MATCHED_BINARY_GPU", "https://github.com/gnina/gnina"),
        ("static structure", "RTMScore", "residue-atom graph transformer", "NOT_RUN_LEGACY_CUDA_STACK", "https://doi.org/10.1021/acs.jcim.2c01226"),
        ("static structure", "DeepRLI", "multi-objective graph transformer", "PILOT_ONLY_NOT_COMMON_PROTOCOL", "https://doi.org/10.1039/D4DD00403E"),
        ("static structure", "EquiScore", "equivariant heterogeneous GNN", "NOT_RUN_NO_COMMON_PREDICTIONS", "https://doi.org/10.1038/s42256-024-00849-z"),
        ("MD representation", "VAMPnet", "deep kinetic model", "NOT_RUN_DATA_GATE", "https://doi.org/10.1038/s41467-017-02388-1"),
        ("MD representation", "SPIB", "predictive information bottleneck", "NOT_RUN_DATA_GATE", "https://doi.org/10.1063/5.0038198"),
        ("MD representation", "Geom2Vec pretrained", "geometric MD featurizer", "RUN_SUPERSEDED_BY_C1_BS256", "https://github.com/dinner-group/geom2vec"),
    ]:
        add(module, "future/common-protocol baseline", method, family, reason, False, "", "", "candidate baseline", source)

    frame = pd.DataFrame(rows)
    frame.to_csv(out / "benchmark_method_registry_v02.csv", index=False)
    summary = {
        "rows": len(frame), "modules": frame.module.value_counts().to_dict(),
        "status_counts": frame.status.value_counts().to_dict(),
        "rule": "Only same_protocol_comparable=true rows may enter head-to-head winner claims.",
    }
    (out / "benchmark_method_registry_v02.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
