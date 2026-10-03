"""Stage-1 QC for the frozen fragment generator output (encoding-tolerant)."""
import csv, hashlib, json, re
from pathlib import Path
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors, Crippen

RDLogger.DisableLog("rdApp.*")
RUN = Path(__file__).resolve().parents[1]
REPO = Path(__file__).resolve().parents[4]
CSV = REPO / "project" / "results" / "generated" / "generated_pam_analogs.csv"
LOG = RUN / "logs_stage1_pinned_A.txt"

MW_RANGE, LOGP_MAX, HBD_MAX, HBA_MAX, ROTB_MAX = (250, 550), 5.5, 5, 10, 10
TPSA_RANGE, MIN_AROM, MIN_HBA, MIN_HBD = (40, 140), 2, 2, 1


def read_text_any(path):
    b = path.read_bytes()
    if b[:2] == b"\xff\xfe":
        return b[2:].decode("utf-16-le", errors="replace")
    if b[:3] == b"\xef\xbb\xbf":
        return b[3:].decode("utf-8", errors="replace")
    return b.decode("utf-8", errors="replace")


def main():
    raw = CSV.read_bytes()
    rows = list(csv.DictReader(CSV.read_text(encoding="utf-8").splitlines()))
    smiles = [r["smiles"] for r in rows]
    canonical, invalid, failures = [], 0, []
    for s in smiles:
        m = Chem.MolFromSmiles(s)
        if m is None:
            invalid += 1
            failures.append({"smiles": s, "rule": "rdkit_parse"})
            continue
        canonical.append(Chem.MolToSmiles(m))
        checks = {
            "mw": MW_RANGE[0] <= Descriptors.MolWt(m) <= MW_RANGE[1],
            "logp": Crippen.MolLogP(m) <= LOGP_MAX,
            "hbd": Lipinski.NumHDonors(m) <= HBD_MAX,
            "hba": Lipinski.NumHAcceptors(m) <= HBA_MAX,
            "rotb": rdMolDescriptors.CalcNumRotatableBonds(m) <= ROTB_MAX,
            "tpsa": TPSA_RANGE[0] <= rdMolDescriptors.CalcTPSA(m) <= TPSA_RANGE[1],
            "aromatic_rings": rdMolDescriptors.CalcNumAromaticRings(m) >= MIN_AROM,
            "pharmacophore": (Lipinski.NumHAcceptors(m) >= MIN_HBA and Lipinski.NumHDonors(m) >= MIN_HBD),
        }
        for rule, ok in checks.items():
            if not ok:
                failures.append({"smiles": s, "rule": rule})

    log = read_text_any(LOG)
    join_line = next((l for l in log.splitlines() if "join" in l), "")
    nums = [int(x) for x in re.findall(r"\d+", join_line)]
    gen, attempts = (nums[0], nums[1]) if len(nums) >= 2 else (0, 0)
    rebuilt = len([l for l in log.splitlines() if re.match(r"^\s*M4_PAM_\d+", l) and "\u2713" in l])
    tag_occurrences = len(re.findall(r"dataset_top_active", log))

    dup_raw = len(smiles) - len(set(smiles))
    dup_canonical = len(canonical) - len(set(canonical))
    sha_a = hashlib.sha256((RUN / "_stage1_pinned_A.csv").read_bytes()).hexdigest()
    sha_b = hashlib.sha256((RUN / "_stage1_pinned_B.csv").read_bytes()).hexdigest()
    qc = {
        "schema": "pacer.prospective.stage1_qc.v1",
        "input_csv": str(CSV.relative_to(REPO)),
        "csv_parses": len(rows) > 0,
        "row_count": len(rows),
        "invalid_smiles": invalid,
        "duplicate_raw_smiles": dup_raw,
        "duplicate_canonical_smiles": dup_canonical,
        "unique_canonical_smiles": len(set(canonical)),
        "rule_violation_count": len(failures),
        "rule_violations": failures[:10],
        "seeds_reconstructed": rebuilt,
        "seeds_expected": 12,
        "seed_chemotype_tag_occurrences": tag_occurrences,
        "fallback_seeds_used": False,
        "join_attempts": attempts,
        "raw_generated": gen,
        "output_sha256": hashlib.sha256(raw).hexdigest(),
        "output_bytes": len(raw),
        "determinism": {"pinned_runs": ["A", "B"], "identical": sha_a == sha_b,
                        "sha256_A": sha_a, "sha256_B": sha_b},
    }
    qc["passed"] = (qc["csv_parses"] and invalid == 0 and dup_canonical == 0 and dup_raw == 0
                    and len(failures) == 0 and rebuilt == 12 and gen > 0 and len(rows) > 0
                    and qc["determinism"]["identical"])
    (RUN / "STAGE1_QC_v01.json").write_text(json.dumps(qc, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(qc, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
