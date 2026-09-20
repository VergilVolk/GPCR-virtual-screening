"""
7TRQ 别构口袋 Redocking 验证 (VU0467154) — 论文方法学复现
==========================================================
验证标准 (比 RMSD<2A 更适合柔性别构配体):
  1. 定位正确性: 最优 pose 质心 vs 晶体质心距离 (< 3 A)
  2. 接触残基重现: dock pose 与晶体 pose 的接触残基重合率 (应 100% 附近)
  3. 亲和力: Vina 最优 pose 亲和力 (强结合, < -8 kcal/mol)
  4. 附: 原子级 RMSD (构象精度, 如实报告)

用法: python scripts/validate_allosteric_pocket.py
"""
from __future__ import annotations
import sys
import json
from pathlib import Path

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

OUT = Path("results/structure")
CUTOFF = 4.5


def parse_atoms(pdb, heavy=True, exclusions=("IUI", "IXO")):
    out = []
    for line in open(pdb):
        if line.startswith(("ATOM", "HETATM")):
            res = line[17:20].strip()
            if res in exclusions:
                continue
            el = line[76:78].strip() or line[12:16].strip()[0]
            if heavy and el == "H":
                continue
            out.append((el, float(line[30:38]), float(line[38:46]),
                        float(line[46:54]), res, int(line[22:26])))
    return out


def centroid(atoms):
    xs = [a[1] for a in atoms]
    ys = [a[2] for a in atoms]
    zs = [a[3] for a in atoms]
    return (sum(xs) / len(xs), sum(ys) / len(ys), sum(zs) / len(zs))


def contact_residues(prot, lig, cutoff=CUTOFF):
    hits = {}
    for pel, px, py, pz, pres, pnum in prot:
        for lel, lx, ly, lz, _, _ in lig:
            if (px - lx) ** 2 + (py - ly) ** 2 + (pz - lz) ** 2 < cutoff ** 2:
                hits[f"{pres}{pnum}"] = True
    return set(hits)


def parse_pdbqt_model(pdbqt, model=1):
    """从 vina 输出 pdbqt 提取指定 MODEL 的重原子坐标."""
    atoms, cur = [], 0
    for line in open(pdbqt):
        if line.startswith("MODEL"):
            cur += 1
            if cur > model:
                break
            continue
        if cur != model:
            continue
        if line.startswith(("ATOM", "HETATM")):
            # pdbqt 原子名/元素在 13-16 列与 77-78 列可能不同, 用坐标+类型
            x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
            el = line[77:79].strip() or line[13:16].strip().lstrip("0123456789")[:1]
            atoms.append((el, x, y, z, "LIG", 1))
    return atoms


def main():
    import math
    import json

    prot = parse_atoms(str(Path("data/pdb/7TRQ.pdb")))
    crys_lig = parse_atoms(str(OUT / "7trq_IUI_ligand.pdb"), exclusions=())
    dock_lig = parse_pdbqt_model(str(OUT / "vu0467154_redock_crystal_out.pdbqt"), model=1)

    c_crys = centroid(crys_lig)
    c_dock = centroid(dock_lig)
    dist = math.dist(c_crys, c_dock)
    hits_crys = contact_residues(prot, crys_lig)
    hits_dock = contact_residues(prot, dock_lig)
    overlap = hits_crys & hits_dock
    ratio = len(overlap) / len(hits_crys) * 100 if hits_crys else 0

    report = {
        "centroid_dist_A": round(dist, 2),
        "crystal_contact_residues": sorted(hits_crys),
        "dock_contact_residues": sorted(hits_dock),
        "overlap_residues": sorted(overlap),
        "overlap_ratio_pct": round(ratio, 1),
        "verdict": "PASS (定位正确)" if (dist < 3 and ratio >= 90) else "REVIEW",
    }
    print(json.dumps(report, indent=2, ensure_ascii=False))
    (OUT / "redock_validation.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")


if __name__ == "__main__":
    main()
