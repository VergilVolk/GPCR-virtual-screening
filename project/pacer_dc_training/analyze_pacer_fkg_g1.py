#!/usr/bin/env python3
"""PACER-FKG G1: read-only diagnostics over *frozen* R2/R3 JSON results.

No embedding recomputation, graph reconstruction, training, parameter tuning, or
redefinition of previously reported gates. This program writes only to a new
output directory and refuses to overwrite it.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics as stats
from pathlib import Path

AXES = ("synergy_interaction", "intrinsic_agonism", "conditional_pam_effect")
CONTROLS = ("stable_core_control", "distal_control")
EXPECTED_REGIONS = 9
EXPECTED_REPLICAS = (2, 3)
EXPECTED_WINDOWS = (0, 1, 2, 3, 4)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def median(values):
    return float(stats.median(values))


def rankdata(values):
    """Average ranks for ties; uses no third-party dependencies."""
    ordered = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)
    start = 0
    while start < len(values):
        end = start + 1
        while end < len(values) and values[ordered[end]] == values[ordered[start]]:
            end += 1
        rank = (start + 1 + end) / 2.0
        for pos in range(start, end):
            ranks[ordered[pos]] = rank
        start = end
    return ranks


def spearman(x, y):
    if len(x) != len(y) or len(x) < 2:
        return None
    rx, ry = rankdata(x), rankdata(y)
    mx, my = stats.mean(rx), stats.mean(ry)
    dx = [a - mx for a in rx]
    dy = [b - my for b in ry]
    denom = math.sqrt(sum(z * z for z in dx) * sum(z * z for z in dy))
    return sum(a * b for a, b in zip(dx, dy)) / denom if denom else None


def approx(a, b, tol=1e-9):
    return a is not None and b is not None and math.isclose(a, b, rel_tol=tol, abs_tol=tol)


def write_csv(path, rows, fields):
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def validate_and_load(input_root: Path, audit_path: Path):
    if not input_root.is_dir() or not audit_path.is_file():
        raise FileNotFoundError("Frozen input directory or aggregate audit JSON not found")
    d = json.loads(audit_path.read_text(encoding="utf-8"))
    if tuple(d["replicas"]) != EXPECTED_REPLICAS or tuple(d["windows"]) != EXPECTED_WINDOWS:
        raise ValueError("G1 only accepts frozen R2/R3 W0-W4 audit")
    if d["provenance"]["n_inputs"] != 40 or len(d["provenance"]["input_files"]) != 40:
        raise ValueError("Frozen source manifest does not contain 40 Embedding inputs")
    if len(d["window_results"]) != 90 or len(d["cross_replica_summary"]) != 27:
        raise ValueError("Unexpected aggregate result shape")
    if d["provenance"].get("preprocess") != "global_residual_channel_robust":
        raise ValueError("Unexpected frozen preprocessing")
    if d["provenance"].get("bandwidth") != "per-call median heuristic":
        raise ValueError("Unexpected frozen kernel bandwidth specification")
    graph_hash = d["provenance"].get("graph_sha256", "")
    if len(graph_hash) != 64:
        raise ValueError("Missing graph SHA256")
    windows = {}
    window_hashes = []
    for r in EXPECTED_REPLICAS:
        for w in EXPECTED_WINDOWS:
            p = input_root / f"R{r}_W{w:03d}.fkg.json"
            if not p.is_file():
                raise FileNotFoundError(p)
            content = json.loads(p.read_text(encoding="utf-8"))
            if len(content["region_scores"]) != EXPECTED_REGIONS or len(content["node_scores"]) != 270:
                raise ValueError(f"Unexpected region or residue count: {p}")
            region_names = [x["region"] for x in content["region_scores"]]
            if len(set(region_names)) != EXPECTED_REGIONS or not set(CONTROLS).issubset(region_names):
                raise ValueError(f"Incomplete/duplicate region names: {p}")
            if [n["embedding_index"] for n in content["node_scores"]] != list(range(270)):
                raise ValueError(f"Wrong node ordering: {p}")
            windows[(r, w)] = content
            window_hashes.append({"replica": r, "window": w, "path": p.name, "sha256": sha256(p)})
    aggregates = {(int(row["replica"]), int(row["window"]), row["region"]): row
                  for row in d["window_results"]}
    if len(aggregates) != 90:
        raise ValueError("Duplicate aggregate keys")
    # Cross-check all source window-level region results with frozen aggregate.
    for (r, w), content in windows.items():
        for src in content["region_scores"]:
            dest = aggregates.get((r, w, src["region"]))
            if dest is None:
                raise ValueError(f"Missing aggregate row for {r}/{w}/{src['region']}")
            for axis in AXES:
                for suffix in ("kernel_u2", "kernel_biased2", "bandwidth", "linear", "graph_diffused_mean"):
                    key = f"{axis}_{suffix}"
                    if not approx(src[key], dest[key]):
                        raise ValueError(f"Aggregate/window mismatch: R{r} W{w} {src['region']} {key}")
    return d, windows, aggregates, window_hashes


def compute(audit, aggregates):
    all_regions = sorted({region for _, _, region in aggregates})
    summary_index = {(s["axis"], s["region"]): s for s in audit["cross_replica_summary"]}
    if len(summary_index) != len(AXES) * len(all_regions):
        raise ValueError("Incomplete/duplicate axis/region summary")
    window_rows, summary_rows, notes = [], [], []
    for axis in AXES:
        for region in all_regions:
            summary = summary_index[(axis, region)]
            per_replica = {}
            for r in EXPECTED_REPLICAS:
                series = [aggregates[(r, w, region)] for w in EXPECTED_WINDOWS]
                stable = [aggregates[(r, w, "stable_core_control")] for w in EXPECTED_WINDOWS]
                u = [float(x[f"{axis}_kernel_u2"]) for x in series]
                b = [float(x[f"{axis}_bandwidth"]) for x in series]
                gs = [float(x[f"{axis}_graph_diffused_mean"]) for x in series]
                du = [float(x[f"{axis}_kernel_u2"] - s[f"{axis}_kernel_u2"]) for x, s in zip(series, stable)]
                dg = [float(x[f"{axis}_graph_diffused_mean"] - s[f"{axis}_graph_diffused_mean"]) for x, s in zip(series, stable)]
                raw = [float(x[f"{axis}_kernel_biased2"]) for x in series]
                cv = stats.pstdev(b) / stats.mean(b) if stats.mean(b) else None
                per_replica[r] = dict(du=du, dg=dg, bandwidth=b, graph=gs, u2=u)
                frozen = summary["replicas"][str(r)]
                if not approx(median(du), frozen["median_u2_minus_stable"]):
                    raise ValueError(f"Stored/recomputed u2 median differs: {axis}/{region}/R{r}")
                if not approx(median(dg), frozen["median_graph_diffused_minus_stable"]):
                    raise ValueError(f"Stored/recomputed graph median differs: {axis}/{region}/R{r}")
                if sum(v > 0 for v in du) != frozen["positive_windows"]:
                    raise ValueError("Stored/recomputed positive-window counts differ")
                for i, w in enumerate(EXPECTED_WINDOWS):
                    window_rows.append({
                        "axis": axis, "region": region, "replica": r, "window": w,
                        "n_residues": series[i]["n_residues"], "kernel_u2": u[i],
                        "kernel_biased2": raw[i], "bandwidth": b[i],
                        "kernel_u2_minus_stable": du[i],
                        "graph_diffused_mean": gs[i],
                        "graph_diffused_minus_stable": dg[i],
                        "linear_norm_descriptive": series[i][f"{axis}_linear"],
                    })
                summary_rows.append({
                    "axis": axis, "region": region, "replica": r,
                    "n_residues": series[0]["n_residues"],
                    "median_u2_minus_stable": median(du),
                    "positive_windows": sum(v > 0 for v in du),
                    "median_graph_minus_stable": median(dg),
                    "graph_positive_windows": sum(v > 0 for v in dg),
                    "median_bandwidth": median(b),
                    "min_bandwidth": min(b), "max_bandwidth": max(b),
                    "bandwidth_cv": cv,
                    "median_kernel_u2": median(u),
                    "gate_in_frozen_report": summary["preregistered_specificity_gate"],
                })
            rho_u = spearman(per_replica[2]["du"], per_replica[3]["du"])
            rho_g = spearman(per_replica[2]["dg"], per_replica[3]["dg"])
            if (rho_u is None) != (summary["matched_window_spearman"] is None) or (rho_u is not None and not approx(rho_u, summary["matched_window_spearman"])):
                raise ValueError(f"Matched-window rho mismatch: {axis}/{region}")
            if (rho_g is None) != (summary["matched_window_graph_diffused_spearman_descriptive"] is None) or (rho_g is not None and not approx(rho_g, summary["matched_window_graph_diffused_spearman_descriptive"])):
                raise ValueError(f"Graph rho mismatch: {axis}/{region}")
            notes.append({
                "axis": axis, "region": region, "rho_u2": rho_u, "rho_graph": rho_g,
                "rho_graph_minus_u2": None if rho_g is None or rho_u is None else rho_g - rho_u,
                "gate_in_frozen_report": summary["preregistered_specificity_gate"],
                "median_bandwidth_r2": median(per_replica[2]["bandwidth"]),
                "median_bandwidth_r3": median(per_replica[3]["bandwidth"]),
                "bandwidth_median_ratio_r3_r2": median(per_replica[3]["bandwidth"]) / median(per_replica[2]["bandwidth"]),
                "r2_median_u2_minus_stable": median(per_replica[2]["du"]),
                "r3_median_u2_minus_stable": median(per_replica[3]["du"]),
                "r2_median_graph_minus_stable": median(per_replica[2]["dg"]),
                "r3_median_graph_minus_stable": median(per_replica[3]["dg"]),
            })
    return window_rows, summary_rows, notes


def build_report(audit, windows, window_hashes, summary_rows, comparison):
    counts = {}
    for axis in AXES:
        active = [r for r in comparison if r["axis"] == axis and r["region"] not in CONTROLS]
        counts[axis] = {
            "frozen_gate_passing_regions": [r["region"] for r in active if r["gate_in_frozen_report"]],
            "graph_rho_increased": sum(r["rho_graph_minus_u2"] is not None and r["rho_graph_minus_u2"] > 1e-9 for r in active),
            "graph_rho_decreased": sum(r["rho_graph_minus_u2"] is not None and r["rho_graph_minus_u2"] < -1e-9 for r in active),
            "graph_rho_unchanged": sum(r["rho_graph_minus_u2"] is not None and abs(r["rho_graph_minus_u2"]) <= 1e-9 for r in active),
        }
    distal = [r for r in comparison if r["region"] == "distal_control"]
    return {
        "method": "PACER-FKG G1 read-only frozen-result diagnostic v01",
        "claim_boundary": "Descriptive two-replica method diagnostics only; no PAM pharmacology, confidence intervals, directional RKHS agreement, or independent graph advantage established.",
        "input_audit_sha256": audit["_input_sha256"],
        "input_graph_sha256": audit["provenance"]["graph_sha256"],
        "input_window_sha256": window_hashes,
        "input_manifest_sha256": [item["sha256"] for item in audit["provenance"]["input_files"]],
        "replicas": list(EXPECTED_REPLICAS), "windows": list(EXPECTED_WINDOWS),
        "axes": counts,
        "distal_control": distal,
        "comparison": comparison,
        "interpretation_constraints": [
            "Five contiguous windows are not five independent replicates.",
            "Per-region/per-window and per-axis median-heuristic bandwidths limit cross-region raw kernel magnitude comparisons.",
            "Squared RKHS contrast norm is not the direction or pharmacological sign of a response.",
            "Graph-diffused residue kernel-root means and joint-region kernel U2 are distinct estimands; their numerical magnitudes are not directly comparable.",
            "No inference on graph causality from pre/post graph Spearman alone; no shuffled-graph control here.",
            "The frozen gate excludes distal control by definition, so passing it does not certify distal specificity.",
            "The prior-registration status of the frozen numerical gate thresholds has not been independently verified.",
        ],
    }


def make_markdown(report):
    lines = ["# PACER-FKG G1 · 冻结结果只读专项审计", "",
             "## 范围与数据完整性", "",
             "R2/R3 × W0–W4 × 四上下文；输入为已归档的 10 份窗口 JSON 及综合报告。",
             f"综合报告 SHA256：`{report['input_audit_sha256']}`；图 SHA256：`{report['input_graph_sha256']}`。", "",
             "没有重新提取 Embedding、训练编码器、重新计算核统计或更改原有门禁。", "",
             "## 原门禁及图传播描述性结果", "",
             "| 轴 | 原门禁通过区域 | 图传播后窗口 rho 上升/下降/持平（七功能区域） |", "|---|---|---|" ]
    for axis in AXES:
        rec = report["axes"][axis]
        lines.append(f"| {axis} | {', '.join(rec['frozen_gate_passing_regions']) or '无'} | {rec['graph_rho_increased']}/{rec['graph_rho_decreased']}/{rec['graph_rho_unchanged']} |")
    lines.extend(["", "## 分区核带宽与跨 replica 描述性一致性", "",
                  "详细数据见 `G1_axis_region_comparison.csv`、`G1_replica_region_summary.csv` 和 `G1_window_diagnostics.csv`。", "",
                  "特别注意：当前区域核 U² 和先计算残基级核距离再图传播的区域均值不是同一估计量。", "",
                  "## 远端对照", ""])
    for rec in report["distal_control"]:
        lines.append(f"- {rec['axis']}：R2 中位差 {rec['r2_median_u2_minus_stable']:.6f}，R3 {rec['r3_median_u2_minus_stable']:.6f}；原 rho={rec['rho_u2']:.3f}，传播后 rho={rec['rho_graph']:.3f}。")
    lines.extend(["", "## 解释限制", ""])
    lines.extend("- " + s for s in report["interpretation_constraints"])
    lines.extend(["", "## 后续工作（不属于本次 G1 结果）", "",
                  "先核验 atom14/PBC 的独立 QC，再在新冻结方案下开展共同核空间方向比较、时间块重采样，以及无图/置乱图对照。", ""])
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", type=Path, required=True, help="Frozen R2/R3 FKG result directory")
    parser.add_argument("--output-root", type=Path, required=True, help="New G1 result directory; must not exist")
    parser.add_argument("--preflight-only", action="store_true")
    args = parser.parse_args()
    inp, out = args.input_root.resolve(), args.output_root.resolve()
    if out.exists() or out == inp or inp in out.parents:
        raise FileExistsError("G1 output must be a NEW, separate directory (no overwrite or nesting)")
    audit_path = inp / "PACER_FKG_R2R3_AUDIT.json"
    audit, windows, aggregate, hashes = validate_and_load(inp, audit_path)
    audit["_input_sha256"] = sha256(audit_path)
    window_rows, summary_rows, comparison = compute(audit, aggregate)
    report = build_report(audit, windows, hashes, summary_rows, comparison)
    if args.preflight_only:
        print(json.dumps({"G1_PREFLIGHT": "PASS", "windows_checked": len(windows), "aggregate_regions": len(aggregate), "axes": len(AXES), "output_untouched": str(out)}, indent=2))
        return
    out.mkdir(parents=True, exist_ok=False)
    write_csv(out / "G1_window_diagnostics.csv", window_rows, list(window_rows[0]))
    write_csv(out / "G1_replica_region_summary.csv", summary_rows, list(summary_rows[0]))
    write_csv(out / "G1_axis_region_comparison.csv", comparison, list(comparison[0]))
    (out / "G1_DIAGNOSTIC_AUDIT.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    (out / "G1_DIAGNOSTIC_REPORT.md").write_text(make_markdown(report), encoding="utf-8")
    print(json.dumps({"G1_STATUS": "PASS", "windows_checked": len(windows), "window_rows": len(window_rows), "replica_summary_rows": len(summary_rows), "axis_region_comparisons": len(comparison), "output": str(out)}, indent=2))


if __name__ == "__main__":
    main()
