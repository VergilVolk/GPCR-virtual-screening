"""
方向1: CHRM4 PAM 类似物片段生成 (结构驱动, 纯 CPU)
==============================================
从冻结数据 (CHRM4_PAM_Modeling_Handoff_v1.0) 中 pEC50 最高的 CHRM4 PAM 种子出发,
BRICS 拆分 → 受控片段重组 → 药效团+类药性过滤, 产出 PAM 类似物候选库。

设计要点:
  - 结构驱动, 不依赖大规模活性数据 (适配 PAM 功能活性数据稀缺)
  - PAM 药效团约束来自 M4 别构口袋 (7TRQ: VU0467154 结合位点, 胞外 vestibule)
    的接触模式: 芳环簇 (π-π/疏水 vs 口袋芳香残基) + 极性簇 (H键 vs 极性残基)
  - 重组用「递归生长」而非 RDKit BRICSBuild: BRICSBuild 对含
    [3*]O[3*]/[5*]N[5*] 这类双连接点同标签 linker 会无限聚合导致组合爆炸,
    这里用「片段数预算 + 结果上限」严格限界。

用法: python scripts/fragment_generation.py
"""
from __future__ import annotations
import sys
from collections import defaultdict
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import BRICS, Descriptors, rdMolDescriptors
from rdkit.Chem import Crippen, Lipinski
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

OUT = Path("results/generated")
ACTIVES = Path("data/generated/actives.smi")

# 已知 CHRM4 (M4) PAM 阳性种子: 优先读 data/generated/actives.smi (12 个高活性,
# 取自冻结数据 modeling_potency_exact_calcium.csv pEC50 最高的一批),
# 否则用内置默认 (thienopyrimidinone 系 VU 类 PAM).
_DEFAULT_SEEDS = [
    ("Cc1nnc2sc(C(=O)NC3CN(c4ccncc4Cl)C3)c(N)c2c1C", "M4_PAM_01", "thienopyrimidinone"),
    ("Cc1c(Cl)c2nnc(C)n2c2sc(C(=O)NC3Cc4ccccc4C3)c(N)c12", "M4_PAM_02", "thienopyrimidinone"),
    ("Cc1c(Cl)nnc2sc(C(=O)NC3Cc4ccccc4C3)c(N)c12", "M4_PAM_03", "thienopyridazinone"),
    ("COc1cc(N2CC(NC(=O)c3sc4nnc(Cl)c(C)c4c3N)C2)c(Cl)cn1", "M4_PAM_04", "thienopyridazinone"),
    ("Cc1nnc2sc(C(=O)NC3CN(c4cc(F)nc(F)c4)C3)c(N)c2c1C", "M4_PAM_05", "thienopyrimidinone"),
    ("Cc1c(Cl)c2nncn2c2sc(C(=O)NC3Cc4ccccc4C3)c(N)c12", "M4_PAM_06", "thienopyrimidinone"),
]


def load_seeds():
    """读 data/generated/actives.smi (SMILES 每行一个), 带名字回退为内置默认."""
    smis = []
    if ACTIVES.exists():
        smis = [l.strip() for l in ACTIVES.read_text(encoding="utf-8").splitlines()
                if l.strip()]
    if smis:
        return [(s, f"M4_PAM_{i:02d}", "dataset_top_active") for i, s in enumerate(smis, 1)]
    return _DEFAULT_SEEDS

# 药物样性边界
MW_RANGE = (250, 550)
LOGP_MAX = 5.5
HBD_MAX = 5
HBA_MAX = 10
ROTB_MAX = 10
TPSA_RANGE = (40, 140)

# PAM 药效团 (简化自 M4 别构口袋 VU0467154 接触模式): ≥2 芳环 + ≥2 HBA + ≥1 HBD
MIN_AROM_RINGS = 2
MIN_HBA = 2
MIN_HBD = 1

MAX_FRAGS = 5          # 单分子最大片段数 (自由重组在 ≥8 时组合爆炸; 5 已覆盖种子 MW 范围)
MAX_GENERATE = 5000    # 去重后分子总数上限
MAX_ATTEMPTS = 200000  # 原始 join 尝试次数上限 (安全阀, 防止新种子接入后卡死)

# BRICS 键规则: 从 reactionDefs 提取「标签对 -> 键型」的兼容表
COMPAT = defaultdict(set)   # label(str) -> 可兼容的 label 集合
BOND_TYPE = {}              # (label_a, label_b) -> bond char ('-' / '=')
for _defn in BRICS.reactionDefs:
    for _a, _b, _bt in _defn:
        COMPAT[_a].add(_b)
        COMPAT[_b].add(_a)
        BOND_TYPE[(_a, _b)] = _bt
        BOND_TYPE[(_b, _a)] = _bt


def attachments(mol):
    """返回 (原子idx, 连接点标签) 列表; 标签=假原子同位素字符串 (如 [16*] -> '16')."""
    return [(a.GetIdx(), str(a.GetIsotope())) for a in mol.GetAtoms() if a.GetAtomicNum() == 0]


def join(mol_a, mol_b, da, db, bond_char):
    """把 mol_a[da] 与 mol_b[db] 两个连接点对接 (锚原子按 BRICS 键型成键, 删假原子)."""
    aa = mol_a.GetAtomWithIdx(da).GetNeighbors()[0].GetIdx()
    bb = mol_b.GetAtomWithIdx(db).GetNeighbors()[0].GetIdx()
    combo = Chem.RWMol(Chem.CombineMols(mol_a, mol_b))
    off = mol_a.GetNumAtoms()
    bt = Chem.BondType.DOUBLE if bond_char == "=" else Chem.BondType.SINGLE
    combo.AddBond(aa, bb + off, bt)
    for idx in sorted([db + off, da], reverse=True):
        combo.RemoveAtom(idx)
    try:
        mol = combo.GetMol()
        Chem.SanitizeMol(mol)
        return mol
    except Exception:
        return None


def grow(mol, frag_pool, budget, results, attempts):
    """递归生长: 每次把最左连接点用兼容片段补上, 无连接点即完整分子入库."""
    dummies = attachments(mol)
    if not dummies:
        try:
            smi = Chem.MolToSmiles(mol)
            if smi not in results and len(results) < MAX_GENERATE:
                results.add(smi)
        except Exception:
            pass
        return
    if budget <= 0 or len(results) >= MAX_GENERATE or attempts[0] >= MAX_ATTEMPTS:
        return
    da, label = dummies[0]
    for frag, db, bond_char in frag_pool.get(label, ()):
        if len(results) >= MAX_GENERATE or attempts[0] >= MAX_ATTEMPTS:
            return
        attempts[0] += 1
        new = join(mol, frag, da, db, bond_char)
        if new is not None:
            grow(new, frag_pool, budget - 1, results, attempts)


def _guided_rebuild(cur, frag_list, target, attempts):
    """guided 重建: 只从该种子自身的 BRICS 片段集中选兼容片段装配,
    检验 join 的键型/原子处理是否正确 (能否精确重建种子)."""
    dummies = attachments(cur)
    if not dummies:
        try:
            return Chem.MolToSmiles(cur) == target
        except Exception:
            return False
    if attempts[0] >= MAX_ATTEMPTS:
        return False
    da, label = dummies[0]
    for i, frag in enumerate(frag_list):
        for db, lb in attachments(frag):
            if lb not in COMPAT[label] or (label, lb) not in BOND_TYPE:
                continue
            attempts[0] += 1
            nxt = join(cur, frag, da, db, BOND_TYPE[(label, lb)])
            if nxt is not None and _guided_rebuild(
                    nxt, frag_list[:i] + frag_list[i + 1:], target, attempts):
                return True
    return False


def validate_rebuild(seeds):
    """用 CHRM4 种子做确定性重建校验 (替代旧 MRGPRX1 ML382 校验):
    对每个种子, 从它自己的 BRICS 片段集出发 guided 装配, 必须精确重建出种子."""
    ok_all = True
    for smi, name, chemotype in seeds:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        target = Chem.MolToSmiles(mol)
        frags = [Chem.MolFromSmiles(f) for f in BRICS.BRICSDecompose(mol)]
        frags = [f for f in frags if f is not None]
        rebuilt = False
        for i, start in enumerate(frags):
            attempts = [0]
            if _guided_rebuild(start, frags[:i] + frags[i + 1:], target, attempts):
                rebuilt = True
                break
        ok_all = ok_all and rebuilt
        print(f"  {name} ({chemotype}): 重建 {'✓' if rebuilt else '✗'}")
    return ok_all


def pharmacophore_ok(mol):
    rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    hba = Lipinski.NumHAcceptors(mol)
    hbd = Lipinski.NumHDonors(mol)
    return rings >= MIN_AROM_RINGS and hba >= MIN_HBA and hbd >= MIN_HBD


def druglike_ok(mol):
    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    rot = rdMolDescriptors.CalcNumRotatableBonds(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    return (MW_RANGE[0] <= mw <= MW_RANGE[1]
            and logp <= LOGP_MAX
            and hbd <= HBD_MAX and hba <= HBA_MAX
            and rot <= ROTB_MAX
            and TPSA_RANGE[0] <= tpsa <= TPSA_RANGE[1])


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    seeds = load_seeds()
    print("=" * 60)
    print(f"步骤1: BRICS 拆分种子 (N={len(seeds)})")
    frag_mols = []
    for smi, name, chemotype in seeds:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            print(f"  {name}: 无效 SMILES, 跳过")
            continue
        fs = list(BRICS.BRICSDecompose(mol))
        print(f"  {name} ({chemotype}): {len(fs)} 个片段")
        for f in fs:
            m = Chem.MolFromSmiles(f)
            if m is not None:
                frag_mols.append(m)

    # 去重 + 建立「标签 -> 兼容片段(及其连接点+键型)」索引
    seen_frag = set()
    unique_frags = []
    for m in frag_mols:
        cs = Chem.MolToSmiles(m)
        if cs in seen_frag:
            continue
        seen_frag.add(cs)
        unique_frags.append(m)

    frag_pool = defaultdict(list)
    for m in unique_frags:
        for db, label_b in attachments(m):
            for label_a in COMPAT[label_b]:
                frag_pool[label_a].append((m, db, BOND_TYPE[(label_a, label_b)]))
    print(f"  唯一片段: {len(unique_frags)}, 兼容池大小: "
          f"{ {k: len(v) for k, v in sorted(frag_pool.items())} }")

    print("\n步骤2: 受控重组 (每分子≤%d片段, 总上限%d)" % (MAX_FRAGS, MAX_GENERATE))
    results = set()
    attempts = [0]
    for seed in unique_frags:
        grow(seed, frag_pool, MAX_FRAGS - 1, results, attempts)
    print(f"  生成 {len(results)} 个完整分子 (join 尝试 {attempts[0]} 次)")

    # join 正确性校验 (CHRM4 种子确定性重建)
    ok = validate_rebuild(seeds)
    print(f"  join 正确性校验: CHRM4 种子重建 {'✓ 全部通过' if ok else '✗ 存在失败'}")

    print("\n步骤3: 药效团 + 类药性过滤")
    passed = []
    for smi in results:
        mol = Chem.MolFromSmiles(smi)
        if mol is None or not druglike_ok(mol) or not pharmacophore_ok(mol):
            continue
        passed.append({
            "smiles": smi,
            "mw": Descriptors.MolWt(mol),
            "logp": Crippen.MolLogP(mol),
            "aromatic_rings": rdMolDescriptors.CalcNumAromaticRings(mol),
            "hba": Lipinski.NumHAcceptors(mol),
            "hbd": Lipinski.NumHDonors(mol),
            "tpsa": rdMolDescriptors.CalcTPSA(mol),
        })
    passed.sort(key=lambda r: (abs(r["logp"] - 3.0), abs(r["mw"] - 350)))
    print(f"  通过 {len(passed)} 个 (药效团+类药性)")

    import csv
    csv_path = OUT / "generated_pam_analogs.csv"
    with open(csv_path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["smiles", "mw", "logp", "aromatic_rings", "hba", "hbd", "tpsa"])
        w.writeheader()
        w.writerows(passed)

    print("\n步骤4: 输出候选 (Top 15)")
    print(f"{'#':>3} {'SMILES':<55} {'MW':>6} {'logP':>5} {'环':>2} {'HBA':>3} {'HBD':>2}")
    print("-" * 82)
    for i, r in enumerate(passed[:15], 1):
        smi = r["smiles"] if len(r["smiles"]) <= 54 else r["smiles"][:51] + "..."
        print(f"{i:>3} {smi:<55} {r['mw']:>6.1f} {r['logp']:>5.1f} {r['aromatic_rings']:>2} "
              f"{r['hba']:>3} {r['hbd']:>2}")

    print(f"\n已保存全部 {len(passed)} 个候选到 {csv_path}")


if __name__ == "__main__":
    main()
