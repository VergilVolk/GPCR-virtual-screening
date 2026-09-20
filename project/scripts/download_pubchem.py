"""
下载 PubChem 类药 SMILES 语料 (用于 SMILES LSTM 语言先验)
========================================================
网络说明: 本机 Python urllib/requests 的 TLS 指纹被中间网络重置, 但 curl 可通,
故所有请求走 curl subprocess。

策略: 扫 CID 区间 (1~MAX_CID), 每批 BATCH 个用 PUG-REST POST 取 SMILES+MW,
本地 RDKit 规范化 + 类药过滤, 存成一行一个 SMILES。

用法: python scripts/download_pubchem.py
"""
from __future__ import annotations
import subprocess
import json
import sys
import time
from pathlib import Path
from rdkit import Chem
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

OUT = Path("data/generated")
OUT.mkdir(parents=True, exist_ok=True)
CORPUS = OUT / "corpus.smi"

MAX_CID = 100000
BATCH = 400
MW_RANGE = (150, 600)
TARGET = 40000        # 收满即停
ALLOWED = {"C", "N", "O", "S", "P", "F", "Cl", "Br", "I"}

URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/property/CanonicalSMILES,MolecularWeight/JSON"


def fetch_batch(cids, retries=3):
    data = "cid=" + ",".join(str(c) for c in cids)
    for _ in range(retries):
        try:
            out = subprocess.run(
                ["curl", "-s", "-m", "40", "-X", "POST", URL, "-d", data],
                capture_output=True, text=True, timeout=50,
            )
            d = json.loads(out.stdout)
            return d.get("PropertyTable", {}).get("Properties", [])
        except Exception:
            time.sleep(0.5)
    return []


def is_druglike(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return False
    if Chem.GetFormalCharge(mol) != 0:
        return False
    atoms = {a.GetSymbol() for a in mol.GetAtoms()}
    if not atoms <= ALLOWED or "C" not in atoms:
        return False
    return True


def main():
    seen = set()
    collected = 0
    n_req = 0
    t0 = time.time()
    for start in range(1, MAX_CID + 1, BATCH):
        cids = list(range(start, min(start + BATCH, MAX_CID + 1)))
        props = fetch_batch(cids)
        n_req += 1
        for p in props:
            smi = p.get("CanonicalSMILES") or p.get("ConnectivitySMILES")
            mw = p.get("MolecularWeight")
            if not smi or not mw:
                continue
            try:
                mw = float(mw)
            except ValueError:
                continue
            if not (MW_RANGE[0] <= mw <= MW_RANGE[1]):
                continue
            if not is_druglike(smi):
                continue
            mol = Chem.MolFromSmiles(smi)
            canon = Chem.MolToSmiles(mol)
            if canon in seen:
                continue
            seen.add(canon)
            collected += 1

        if n_req % 20 == 0:
            print(f"  CID ~{start} | 请求 {n_req} | 已收 {collected} | "
                  f"耗时 {time.time()-t0:.0f}s", flush=True)

        if collected >= TARGET:
            print(f"  达到目标 {TARGET}，停止")
            break
        time.sleep(0.1)

    # 写盘
    with open(CORPUS, "w", encoding="utf-8") as fh:
        fh.write("\n".join(sorted(seen)) + "\n")
    print(f"完成: {len(seen)} 个类药 SMILES -> {CORPUS} (请求 {n_req}, {time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()
