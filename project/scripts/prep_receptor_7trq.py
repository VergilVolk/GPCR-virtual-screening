"""
7TRQ 受体制备 + M4 别构口袋定义 + Redocking 验证 (纯 CPU)
==========================================================
PDB 7TRQ: M4–iperoxo(IXO)–VU0467154(IUI)–Gi1 复合物 (2.5 A cryo-EM)

步骤:
  1. 提取 chain R (M4 受体) + IUI (VU0467154 PAM) -> receptor.pdb / vu0467154.pdb
  2. 接触残基分析: IUI 各原子 < 4.5 A 的受体残基 -> 别构口袋残基表 (论文证据)
  3. 盒子定义: IUI 质心 +- 11 A (Vina 搜索空间)
  4. 受体制备: openbabel 加氢 + Gasteiger 电荷 -> receptor.pdbqt
  5. 配体制备: meeko (SMILES -> PDBQT)
  6. Redocking VU0467154: vina 对接 -> 与晶体 pose RMSD (应 < 2 A)
  7. 输出: results/structure/ 下的口袋报告 + config

用法: python scripts/prep_receptor_7trq.py
"""
from __future__ import annotations
import sys
import subprocess
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

PDB = Path("data/pdb/7TRQ.pdb")
OUT = Path("results/structure")
OUT.mkdir(parents=True, exist_ok=True)

# VU0467154 SMILES (从 7TRQ IUI 配体经 obabel 提取并校验):
# thieno[2,3-d]pyrimidine 核心 + 3-氨基 + 2-甲基 + 5-酰胺-NHCH2-4-(trifluoromethylsulfonyl)苯基
VU_SMILES = "Cc1c(N)c2sc(C(=O)NCc3ccc(cc3)S(=O)(=O)C(F)(F)F)c(C)c2nc1"


def extract_chains():
    """提取 chain R 蛋白原子 + IUI 配体原子 (重写序号 1-N + 有效 CONECT)."""
    prot_atoms, iui_atoms = [], []
    iui_serials = []
    with open(PDB) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                ch = line[21]
                if line.startswith("ATOM") and ch == "R":
                    prot_atoms.append(line)
                elif line.startswith("HETATM"):
                    name = line[17:20].strip()
                    if name == "IUI":
                        iui_serials.append(int(line[6:11]))
                        iui_atoms.append(line)
    # IUI 序号重映射为 1..N (保证 CONECT 有效, RDKit/obabel 才能正确建键)
    idx_map = {s: i + 1 for i, s in enumerate(iui_serials)}
    iui_fixed = []
    for line in iui_atoms:
        iui_fixed.append(line[:6] + "{:5d}".format(idx_map[int(line[6:11])]) + line[11:])
    # 收集 IUI 相关 CONECT 并重映射
    iui_set = set(iui_serials)
    conect_lines = []
    with open(PDB) as f:
        for line in f:
            if line.startswith("CONECT"):
                nums = [int(line[i:i + 5]) for i in range(6, len(line.rstrip()), 5)]
                if all(n in iui_set for n in nums):
                    remap = [idx_map[n] for n in nums]
                    conect_lines.append("CONECT" + "".join("{:5d}".format(n) for n in remap) + "\n")
    (OUT / "7trq_R_receptor.pdb").write_text("".join(prot_atoms), encoding="utf-8")
    (OUT / "7trq_IUI_ligand.pdb").write_text("".join(iui_fixed) + "".join(conect_lines),
                                             encoding="utf-8")
    print(f"chain R ATOM: {len(prot_atoms)} | IUI HETATM: {len(iui_fixed)} (CONECT {len(conect_lines)})")
    return prot_atoms, iui_fixed


def parse_xyz(atoms, only_heavy=True):
    """从 PDB 行提取 (element, x, y, z, resname, resnum)."""
    out = []
    for line in atoms:
        el = line[76:78].strip() or line[12:16].strip()[0]
        if only_heavy and el == "H":
            continue
        out.append((el, float(line[30:38]), float(line[38:46]), float(line[46:54]),
                    line[17:20].strip(), int(line[22:26])))
    return out


def contact_residues(prot, lig, cutoff=4.5):
    """IUI 各原子 <cutoff A 的受体残基集合 (残基号+名称)."""
    contacts = {}
    for pel, px, py, pz, pres, pnum in prot:
        for lel, lx, ly, lz, _, _ in lig:
            d2 = (px - lx) ** 2 + (py - ly) ** 2 + (pz - lz) ** 2
            if d2 < cutoff ** 2:
                contacts.setdefault((pnum, pres), []).append(round(d2 ** 0.5, 2))
    return contacts


def centroid(lig):
    xs = [a[1] for a in lig]
    ys = [a[2] for a in lig]
    zs = [a[3] for a in lig]
    return sum(xs) / len(xs), sum(ys) / len(ys), sum(zs) / len(zs)


def main():
    prot, iui = extract_chains()
    prot_h = parse_xyz(prot)
    iui_h = parse_xyz(iui)
    cx, cy, cz = centroid(iui_h)
    print(f"\nVU0467154 质心 (盒子中心): ({cx:.2f}, {cy:.2f}, {cz:.2f})")

    contacts = contact_residues(prot_h, iui_h)
    print(f"\n别构口袋接触残基 (<4.5 A, N={len(contacts)}):")
    res_lines = []
    for (num, name), ds in sorted(contacts.items()):
        line = f"  {name}{num}: min {min(ds):.2f} A, {len(ds)} 原子接触"
        print(line)
        res_lines.append(f"{name}{num}\t{min(ds):.2f}\t{len(ds)}")
    (OUT / "allosteric_pocket_residues.tsv").write_text(
        "residue\tmin_dist_A\tn_contacts\n" + "\n".join(res_lines), encoding="utf-8")

    # 盒子 config (vina)
    box = 22.0
    cfg = (f"center_x = {cx:.3f}\ncenter_y = {cy:.3f}\ncenter_z = {cz:.3f}\n"
           f"size_x = {box}\nsize_y = {box}\nsize_z = {box}\nexhaustiveness = 16\nnum_modes = 9\n")
    (OUT / "allosteric_box.config").write_text(cfg, encoding="utf-8")
    print(f"\n盒子 config 已写 -> {OUT/'allosteric_box.config'}")

    # 受体制备: Meeko 保留 PDB 残基身份并赋 Gasteiger 电荷。
    # OpenBabel 在本 Windows 环境会把部分蛋白残基错误改名，禁止用于正式结果。
    receptor_pdbqt = OUT / "7trq_R_meeko.pdbqt"
    mk = subprocess.run(
        ["mk_prepare_receptor.exe", "--read_pdb", str(OUT / "7trq_R_receptor.pdb"),
         "-o", str(OUT / "7trq_R_meeko"), "-p", "--charge_model", "gasteiger", "-a"],
        capture_output=True, text=True)
    print("\nMeeko receptor:", mk.returncode, mk.stdout.strip()[-300:] or "ok")
    if mk.returncode != 0 or not receptor_pdbqt.exists():
        raise RuntimeError(mk.stderr[-1000:])

    # 配体制备 (meeko): 需显式 H + 3D 构象; v0.7+ 用 PDBQTWriterLegacy
    from meeko import MoleculePreparation
    from meeko.preparation import PDBQTWriterLegacy
    mol = Chem.MolFromSmiles(VU_SMILES)
    if mol is None:
        print("VU SMILES 无效")
        return
    mol = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol, randomSeed=42) != 0:
        AllChem.EmbedMolecule(mol, AllChem.ETKDGv3(), randomSeed=42)
    prep = MoleculePreparation()
    setups = prep.prepare(mol)
    with open(OUT / "vu0467154_ligand.pdbqt", "w") as fh:
        for setup in setups:
            pdbqt_string, is_ok, err_msg = PDBQTWriterLegacy.write_string(setup)
            if not is_ok:
                raise RuntimeError(f"PDBQTWriterLegacy: {err_msg}")
            fh.write(pdbqt_string)
    print("meeko ligand pdbqt 已生成")

    # redocking
    vina = Path("tools/vina.exe")
    dock_out = OUT / "vu0467154_redock_out.pdbqt"
    r = subprocess.run(
        [str(vina), "--receptor", str(receptor_pdbqt),
         "--ligand", str(OUT / "vu0467154_ligand.pdbqt"),
         "--config", str(OUT / "allosteric_box.config"),
         "--out", str(dock_out)],
        capture_output=True, text=True)
    print("\nvina:", r.returncode)
    print(r.stdout[-800:] if r.returncode == 0 else r.stderr[-800:])

    if dock_out.exists():
        # RMSD: 晶体 pose (IUI pdb) vs 最优 redock pose (openbabel OBAlign, 自动处理对称性)
        dock_sdf = OUT / "vu0467154_redock_out.sdf"
        subprocess.run(["obabel", str(dock_out), "-O", str(dock_sdf), "-d"], capture_output=True, text=True)
        from openbabel import openbabel as ob
        ref_ob = ob.OBMol()
        dock_ob = ob.OBMol()
        conv = ob.OBConversion()
        conv.SetInAndOutFormats("pdb", "sdf")
        conv.ReadFile(ref_ob, str(OUT / "7trq_IUI_ligand.pdb"))
        conv.SetInAndOutFormats("sdf", "pdb")
        conv.ReadFile(dock_ob, str(dock_sdf))
        print(f"OBAlign: ref {ref_ob.NumAtoms()} atoms, dock {dock_ob.NumAtoms()} atoms")
        align = ob.OBAlign(False, False)
        align.SetRefMol(ref_ob)
        align.SetTargetMol(dock_ob)
        align.Align()
        rmsd = align.GetRMSD()
        print(f"\nRedocking RMSD (VU0467154): {rmsd:.2f} A  {'✅ <2A' if rmsd < 2 else '❌'}")


if __name__ == "__main__":
    main()
