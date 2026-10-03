"""Stage-2 QC: exact rejection accounting for the frozen physchem / drug-likeness filter."""
import csv, hashlib, json
from pathlib import Path
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, DataStructs, Descriptors, QED
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

RDLogger.DisableLog("rdApp.*")
RUN = Path(__file__).resolve().parents[1]
REPO = Path(__file__).resolve().parents[4]
P = REPO / "project"
GEN = P / "results" / "generated"
OUT = P / "results" / "pacer_candidates_v01"

# frozen constants, copied verbatim from project/scripts/build_candidate_portfolio.py
FROZEN = {
    "mw_range": [250, 550], "logp_range": [-1, 5.5], "tpsa_max": 140, "rotb_max": 10,
    "pains_must_be_zero": True, "domain_local": [0.55, 0.85], "domain_exploratory": [0.35, 0.55],
    "n_neighbors_threshold": 0.45, "n_neighbors_min": 3, "strict_inactive_risk_max": 0.45,
    "quota_local": 160, "quota_exploratory": 40,
}


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    stage1 = GEN / "generated_pam_analogs.csv"
    stage1_smiles = [r["smiles"] for r in csv.DictReader(stage1.read_text(encoding="utf-8").splitlines())]
    audit_rows = list(csv.DictReader((OUT / "all_generated_audit.csv").read_text(encoding="utf-8").splitlines()))
    portfolio_rows = list(csv.DictReader((OUT / "predock_portfolio.csv").read_text(encoding="utf-8").splitlines()))
    audit_json = json.loads((OUT / "audit.json").read_text(encoding="utf-8"))

    n_exact_known = len(stage1_smiles) - len(audit_rows)

    # per-rule failure counts among the audited set (rules are conjunctive, so these overlap)
    rule_fail = {
        "exact_known_molecule": n_exact_known,
        "pains_match": sum(1 for r in audit_rows if int(r["PAINS"]) == 1),
        "druglike_composite_fail": sum(1 for r in audit_rows if int(r["druglike"]) == 0),
        "outside_applicability_domain": sum(1 for r in audit_rows if r["domain_tier"] == "reject_domain"),
        "strict_inactive_risk_ge_0.45": sum(1 for r in audit_rows if float(r["strict_inactive_risk_ref"]) >= 0.45),
    }

    # primary-reason decomposition, frozen precedence order, sums exactly to rejections
    primary = {"exact_known_molecule": n_exact_known, "outside_applicability_domain": 0,
               "druglike_fail": 0, "strict_inactive_risk": 0, "portfolio_quota_cap": 0}
    eligible = 0
    for r in audit_rows:
        if r["domain_tier"] == "reject_domain":
            primary["outside_applicability_domain"] += 1
        elif int(r["druglike"]) == 0:
            primary["druglike_fail"] += 1
        elif float(r["strict_inactive_risk_ref"]) >= 0.45:
            primary["strict_inactive_risk"] += 1
        else:
            eligible += 1
    primary["portfolio_quota_cap"] = eligible - len(portfolio_rows)

    # independent re-verification of the frozen physchem rules on the portfolio
    params = FilterCatalogParams(); params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    catalog = FilterCatalog(params)
    violations = []
    for r in portfolio_rows:
        m = Chem.MolFromSmiles(r["canonical_smiles"])
        if m is None:
            violations.append({"id": r["candidate_id"], "rule": "parse"}); continue
        checks = {
            "mw": 250 <= Descriptors.MolWt(m) <= 550,
            "logp": -1 <= Descriptors.MolLogP(m) <= 5.5,
            "tpsa": Descriptors.TPSA(m) <= 140,
            "rotb": Descriptors.NumRotatableBonds(m) <= 10,
            "pains": not catalog.HasMatch(m),
            "domain_tier": r["domain_tier"] in ("local", "exploratory"),
            "strict_inactive_risk": float(r["strict_inactive_risk_ref"]) < 0.45,
        }
        for rule, ok in checks.items():
            if not ok:
                violations.append({"id": r["candidate_id"], "rule": rule})

    ids = [r["candidate_id"] for r in portfolio_rows]
    expected_ids = [f"PACER{i:04d}" for i in range(1, len(portfolio_rows) + 1)]
    qc = {
        "schema": "pacer.prospective.stage2_qc.v1",
        "frozen_implementation": "project/scripts/build_candidate_portfolio.py",
        "frozen_implementation_sha256": sha(P / "scripts" / "build_candidate_portfolio.py"),
        "frozen_thresholds": FROZEN,
        "thresholds_modified": False,
        "input_count": len(stage1_smiles),
        "input_sha256": sha(stage1),
        "audited_count": len(audit_rows),
        "rejection_count": len(stage1_smiles) - len(portfolio_rows),
        "rejection_counts_by_rule_conjunctive": rule_fail,
        "rejection_counts_primary_reason_exclusive": primary,
        "primary_reason_sum": sum(primary.values()),
        "eligible_count": eligible,
        "pass_count": len(portfolio_rows),
        "unique_final_molecules": len({r["canonical_smiles"] for r in portfolio_rows}),
        "duplicate_candidate_ids": len(ids) - len(set(ids)),
        "candidate_ids_deterministic": ids == expected_ids,
        "independent_rule_violations": violations[:10],
        "independent_rule_violation_count": len(violations),
        "outputs": {
            "predock_portfolio.csv": {"sha256": sha(OUT / "predock_portfolio.csv"),
                                      "bytes": (OUT / "predock_portfolio.csv").stat().st_size,
                                      "rows": len(portfolio_rows)},
            "all_generated_audit.csv": {"sha256": sha(OUT / "all_generated_audit.csv"),
                                        "bytes": (OUT / "all_generated_audit.csv").stat().st_size,
                                        "rows": len(audit_rows)},
            "audit.json": {"sha256": sha(OUT / "audit.json")},
        },
        "producer_audit_json": audit_json,
    }
    qc["passed"] = (len(portfolio_rows) == 200 and qc["unique_final_molecules"] == 200
                    and qc["duplicate_candidate_ids"] == 0 and qc["candidate_ids_deterministic"]
                    and len(violations) == 0 and primary["portfolio_quota_cap"] >= 0)
    (RUN / "STAGE2_QC_v01.json").write_text(json.dumps(qc, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(qc, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
