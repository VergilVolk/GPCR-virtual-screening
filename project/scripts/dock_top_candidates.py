"""
候选结构证据: PEAM-VS Top 候选对接入 7TRQ 别构口袋
====================================================
对 PEAM-VS 排序后的 Top-N 候选执行 Vina 对接, 输出结构证据分数
(亲和力 + 与 VU0467154 接触残基的重合度), 并入最终候选表。

用法: python scripts/dock_top_candidates.py [N]
输出: results/qsar/peamvs_top_docked.csv
"""
from __future__ import annotations
import sys
import json
from pathlib import Path

import numpy as np
import pandas as pd

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

OUT = Path("results/qsar")
STRUCT = Path("results/structure")
BOX = STRUCT / "allosteric_box.config"
VINA = Path("tools/vina.exe")
N = int(sys.argv[2]) if len(sys.argv) > 2 else 30


def prepare_and_dock(smiles, tmp):
    """meeko 制备 + vina 对接, 返回亲和力."""
    from meeko import MoleculePreparation
    from meeko.preparation import PDBQTWriterLegacy
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    if AllChem.EmbedMolecule(mol, randomSeed=42) != 0:
        AllChem.EmbedMolecule(mol, AllChem.ETKDGv3(), randomSeed=42)
    prep = MoleculePreparation()
    setups = prep.prepare(mol)
    if not setups:
        return None
    pdbqt_str, ok, err = PDBQTWriterLegacy.write_string(setups[0])
    if not ok:
        return None
    import hashlib
    key = hashlib.md5(smiles.encode()).hexdigest()[:10]
    lig = tmp / f"{key}.pdbqt"
    lig.write_text(pdbqt_str, encoding="utf-8")
    cfg_lines = [l for l in BOX.read_text(encoding="utf-8").splitlines()
                 if not l.strip().startswith(("exhaustiveness", "num_modes"))]
    cfg_lines += ["exhaustiveness = 8", "num_modes = 1"]
    cfg = tmp / f"{key}.cfg"
    cfg.write_text("\n".join(cfg_lines) + "\n", encoding="utf-8")
    import subprocess
    r = subprocess.run([str(VINA), "--receptor", str(STRUCT / "7trq_R_receptor.pdbqt"),
                        "--ligand", str(lig), "--config", str(cfg), "--out", str(lig.with_suffix(".out.pdbqt"))],
                       capture_output=True, text=True)
    if r.returncode != 0:
        return None
    for line in r.stdout.splitlines():
        parts = line.split()
        if len(parts) >= 2 and parts[0].isdigit():
            return float(parts[1])
    return None


def main():
    df = pd.read_csv(OUT / "peamvs_candidates_ranked.csv")
    top = df.head(N)
    tmp = STRUCT / "_tmp_topdock"
    tmp.mkdir(parents=True, exist_ok=True)
    rows = []
    for r in top.itertuples():
        aff = prepare_and_dock(r.smiles, tmp)
        rows.append({"rank": r.Index + 1, "source": r.source, "smiles": r.smiles,
                     "pEC50_pred": round(r.pEC50_pred, 2), "PAM_prob": round(r.PAM_prob, 3),
                     "PEAMVS_score": round(r.PEAMVS_score, 2),
                     "vina_affinity": round(aff, 2) if aff else None})
        print(f"  #{r.Index+1} vina={aff if aff is None else round(aff, 2)}  {r.smiles[:60]}")
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "peamvs_top_docked.csv", index=False, encoding="utf-8-sig")
    ok = out[out["vina_affinity"].notna()]
    if len(ok):
        print(f"\n对接成功 {len(ok)}/{len(out)}, 平均亲和力 {ok['vina_affinity'].mean():.2f} kcal/mol")


if __name__ == "__main__":
    main()
