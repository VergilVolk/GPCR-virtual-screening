"""
别构-正构耦合机制分析 (7TRQ 双配体结构) — 方向 1
==================================================
科学问题: PAM 的效力(增强正构激动剂响应)在结构上如何体现?
7TRQ 是 M4 + iperoxo(正构) + VU0467154(别构) + Gi1 的双配体复合物,
允许我们量化"别构占据对正构口袋"的影响。

四步:
  1. 双口袋定义: iperoxo 接触残基(正构) / VU0467154 接触残基(别构, <4.5A)
  2. 变构通路: 蛋白残基-残基接触图(原子<4.5A)中, 别构口袋残基 -> 正构口袋残基
     的最短路径(桥接残基) = 变构传导通路
  3. 双配体能量耦合: vina score_only iperoxo 晶体 pose,
     受体含/不含 VU0467154 原子 -> 别构占据对正构结合能的刚性近似
  4. 通路残基药效团: 候选/基准分子对接 pose 对通路残基的接触数 ->
     作为"效力感知结构证据"重打分, 对比原始 Vina 的富集

用法: python scripts/allosteric_coupling_analysis.py
输出: results/structure/coupling_analysis.json + 通路残基
"""
from __future__ import annotations
import sys
import json
import subprocess
from collections import deque
from pathlib import Path

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

PDB = Path("data/pdb/7TRQ.pdb")
OUT = Path("results/structure")
CUTOFF = 4.5  # 接触距离 (A)
CACUT = 8.0   # 残基-残基接触用 Cα 距离阈值


def parse_atoms(pdb_path, exclusions=()):
    out = []
    for line in open(pdb_path):
        if line.startswith(("ATOM", "HETATM")):
            res = line[17:20].strip()
            if res in exclusions:
                continue
            el = line[76:78].strip() or line[12:16].strip()[0]
            if el == "H":
                continue
            out.append((el, float(line[30:38]), float(line[38:46]), float(line[46:54]),
                        res, int(line[22:26]), line[12:16].strip()))
    return out


def contact_residues(prot, lig, cutoff=CUTOFF):
    """配体原子 <cutoff A 的受体残基集合."""
    hits = {}
    for pel, px, py, pz, pres, pnum, pname in prot:
        for lel, lx, ly, lz, _, _, _ in lig:
            if (px - lx) ** 2 + (py - ly) ** 2 + (pz - lz) ** 2 < cutoff ** 2:
                hits[(pnum, pres)] = True
    return set(hits)


def ca_centroid(prot, res):
    """残基 Cα 坐标 (取 'CA' 原子)."""
    for el, x, y, z, rname, rnum, aname in prot:
        if rnum == res and aname == "CA":
            return (x, y, z)
    return None


def residue_contact_graph(prot, residues):
    """残基间接触图: Cα 距离 < CACUT 的残基对连边."""
    nodes = list(residues)
    coords = {r: ca_centroid(prot, r) for r in nodes}
    graph = {r: [] for r in nodes}
    for i in range(len(nodes)):
        for j in range(i + 1, len(nodes)):
            a, b = nodes[i], nodes[j]
            ca, cb = coords.get(a), coords.get(b)
            if ca is None or cb is None:
                continue
            d2 = (ca[0] - cb[0]) ** 2 + (ca[1] - cb[1]) ** 2 + (ca[2] - cb[2]) ** 2
            if d2 < CACUT ** 2:
                graph[a].append(b)
                graph[b].append(a)
    return graph


def shortest_paths(graph, sources, targets):
    """从任一起点到任一点的最短路径 (BFS), 返回所有 (path, 长度)."""
    best_paths = []
    best_len = None
    for s in sources:
        dist = {s: 0}
        prev = {s: None}
        q = deque([s])
        while q:
            u = q.popleft()
            if u in targets:
                # 重建路径
                path = []
                cur = u
                while cur is not None:
                    path.append(cur)
                    cur = prev[cur]
                path.reverse()
                if best_len is None or len(path) < best_len:
                    best_len = len(path)
                    best_paths = [path]
                elif len(path) == best_len:
                    best_paths.append(path)
                break  # BFS 首个到达即最短
            for v in graph.get(u, []):
                if v not in dist:
                    dist[v] = dist[u] + 1
                    prev[v] = u
                    q.append(v)
    return best_paths, best_len


def atom_contact_graph(prot, residues, cutoff=CUTOFF):
    """残基间原子接触图: 任意原子对 < cutoff A 即连边 (比 Cα 更细)."""
    nodes = list(residues)
    atoms_by_res = {}
    for el, x, y, z, rname, rnum, aname in prot:
        if rnum in residues:
            atoms_by_res.setdefault(rnum, []).append((x, y, z))
    graph = {r: [] for r in nodes}
    for i in range(len(nodes)):
        for j in range(i + 1, len(nodes)):
            a, b = nodes[i], nodes[j]
            hit = False
            for ax, ay, az in atoms_by_res.get(a, ()):
                for bx, by, bz in atoms_by_res.get(b, ()):
                    if (ax - bx) ** 2 + (ay - by) ** 2 + (az - bz) ** 2 < cutoff ** 2:
                        hit = True
                        break
                if hit:
                    break
            if hit:
                graph[a].append(b)
                graph[b].append(a)
    return graph


def all_shortest_paths(graph, sources, targets, max_len=5):
    """别构->正构 的所有最短路径 (BFS, 限长)."""
    found = []
    best = None
    for s in sources:
        dist = {s: 0}
        prev = {s: []}
        q = deque([s])
        while q:
            u = q.popleft()
            if best is not None and dist[u] >= best:
                continue
            if u in targets:
                # 回溯所有等长路径 (简化: 记录一条)
                path = []
                cur = u
                while cur is not None:
                    path.append(cur)
                    cur = prev[cur][0] if prev[cur] else None
                path.reverse()
                if best is None or len(path) < best:
                    best = len(path)
                    found = [path]
                elif len(path) == best:
                    found.append(path)
                continue
            for v in graph.get(u, []):
                if v not in dist:
                    dist[v] = dist[u] + 1
                    prev[v] = [u]
                    q.append(v)
                elif dist[v] == dist[u] + 1:
                    prev[v].append(u)
    return [p for p in found if len(p) <= max_len], best


def main():
    prot = parse_atoms(PDB, exclusions=("IUI", "IXO"))
    iui = parse_atoms(PDB, exclusions=())  # 全部, 下面过滤
    iui = [a for a in iui if a[4] == "IUI"]
    ixo = [a for a in parse_atoms(PDB, exclusions=()) if a[4] == "IXO"]

    # 1) 双口袋
    allosteric = contact_residues(prot, iui)   # 别构口袋 (VU0467154)
    orthosteric = contact_residues(prot, ixo)  # 正构口袋 (iperoxo)
    print(f"别构口袋残基 ({len(allosteric)}): {sorted(allosteric)}")
    print(f"正构口袋残基 ({len(orthosteric)}): {sorted(orthosteric)}")

    # 2) 变构通路 (原子级接触图)
    all_res = allosteric | orthosteric
    graph = atom_contact_graph(prot, all_res)
    paths, plen = all_shortest_paths(graph, list(allosteric), list(orthosteric))
    print(f"\n变构通路 (别构->正构, 原子级接触图, 最短 {plen} 个残基):")
    bridge = set()
    for p in paths:
        print("  " + " -> ".join(f"{r[1]}{r[0]}" for r in p))
        for r in p:
            if r not in allosteric and r not in orthosteric:
                bridge.add(r)
    # 共享残基 (同时在两个口袋)
    shared = allosteric & orthosteric
    print(f"共享残基(同时接触两配体): {sorted(shared) if shared else '无'}")
    print(f"桥接残基(传导中间体): {sorted(bridge) if bridge else '无'}")

    # 两配体最短原子距离 (验证不直接接触)
    lig_all = parse_atoms(PDB, exclusions=())
    iui_all = [a for a in lig_all if a[4] == "IUI"]
    ixo_all = [a for a in lig_all if a[4] == "IXO"]
    mind = min((a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2 + (a[3] - b[3]) ** 2
               for a in iui_all for b in ixo_all) ** 0.5
    print(f"两配体最短原子距离: {mind:.2f} A ({'不直接接触' if mind > 5 else '直接接触'})")

    # 3) 双配体能量耦合 (vina score_only): 受体含/不含 VU0467154
    vina = Path("tools/vina.exe")
    # iperoxo pdbqt: 从 7TRQ IXO 片段 (带 CONECT, 序号重映射) obabel 转换
    ixo_serials = sorted(int(line[6:11]) for line in open(PDB)
                         if line.startswith("HETATM") and line[17:20].strip() == "IXO")
    idx_map = {s: i + 1 for i, s in enumerate(ixo_serials)}
    ixo_lines, conect = [], []
    with open(PDB) as f:
        for line in f:
            if line.startswith("HETATM") and line[17:20].strip() == "IXO":
                ixo_lines.append(line[:6] + "{:5d}".format(idx_map[int(line[6:11])]) + line[11:])
            elif line.startswith("CONECT"):
                nums = [int(line[i:i + 5]) for i in range(6, len(line.rstrip()), 5)]
                if all(n in idx_map for n in nums):
                    conect.append("CONECT" + "".join("{:5d}".format(idx_map[n]) for n in nums) + "\n")
    ixo_pdb = OUT / "7trq_IXO_ligand.pdb"
    ixo_pdb.write_text("".join(ixo_lines) + "".join(conect), encoding="utf-8")
    lig_pdbqt = OUT / "iperoxo_crystal.pdbqt"
    subprocess.run(["obabel", str(ixo_pdb), "-O", str(lig_pdbqt)], capture_output=True, text=True)
    print(f"iperoxo pdbqt: {lig_pdbqt.exists()}")

    box = OUT / "orthosteric_box.config"
    cx = sum(a[1] for a in ixo) / len(ixo)
    cy = sum(a[2] for a in ixo) / len(ixo)
    cz = sum(a[3] for a in ixo) / len(ixo)
    box.write_text(f"center_x = {cx:.3f}\ncenter_y = {cy:.3f}\ncenter_z = {cz:.3f}\n"
                   f"size_x = 18.0\nsize_y = 18.0\nsize_z = 18.0\n", encoding="utf-8")

    def score_only(receptor_pdbqt):
        r = subprocess.run([str(vina), "--receptor", str(receptor_pdbqt),
                            "--ligand", str(lig_pdbqt), "--config", str(box),
                            "--score_only"], capture_output=True, text=True, timeout=300)
        for line in r.stdout.splitlines():
            if "Estimated Free Energy of Binding" in line:
                return float(line.split(":")[1].split("(")[0].strip())
        return None

    e_no_pam = score_only(OUT / "7trq_R_receptor.pdbqt")  # 别构空态
    print(f"iperoxo 正构结合能 (无 VU0467154 占据): {e_no_pam} kcal/mol")

    # 受体含 VU0467154: 蛋白 pdbqt + IUI pdbqt 原子行拼接
    prot_pdbqt = OUT / "7trq_R_receptor.pdbqt"
    iui_pdbqt = OUT / "vu0467154_crystal.pdbqt"
    prot_lines = [l for l in prot_pdbqt.read_text(encoding="utf-8").splitlines()
                  if l.startswith(("ATOM", "HETATM"))]
    iui_lines = [l for l in iui_pdbqt.read_text(encoding="utf-8").splitlines()
                 if l.startswith(("ATOM", "HETATM"))]
    with_pam = OUT / "7trq_R_with_VU0467154.pdbqt"
    with_pam.write_text("\n".join(prot_lines + iui_lines) + "\n", encoding="utf-8")
    e_with_pam = score_only(with_pam)
    print(f"iperoxo 正构结合能 (含 VU0467154 占据): {e_with_pam} kcal/mol")
    if e_no_pam is not None and e_with_pam is not None:
        delta = e_with_pam - e_no_pam
        print(f"别构占据对正构结合能的影响: {delta:+.2f} kcal/mol "
              f"({'增强' if delta < 0 else '削弱' if delta > 0 else '无变化'})")

    # 4) 通路残基药效团: 输出供重打分脚本用
    report = {
        "orthosteric_pocket": sorted([f"{r[1]}{r[0]}" for r in orthosteric]),
        "allosteric_pocket": sorted([f"{r[1]}{r[0]}" for r in allosteric]),
        "coupling_paths": [" -> ".join(f"{r[1]}{r[0]}" for r in p) for p in paths],
        "bridge_residues": sorted([f"{r[1]}{r[0]}" for r in bridge]),
        "orthosteric_centroid": [round(cx, 2), round(cy, 2), round(cz, 2)],
    }
    (OUT / "coupling_analysis.json").write_text(json.dumps(report, indent=2, ensure_ascii=False),
                                                encoding="utf-8")
    print(f"\n已保存 -> {OUT/'coupling_analysis.json'}")


if __name__ == "__main__":
    main()
