from pathlib import Path
import pandas as pd
import numpy as np

ROOT = Path(r"C:\projects\GPCR-virtual-screening")

ACC = ROOT / "project/results/pacer_xr_rerun_v01/acceptance_rank1_3"
R45 = ROOT / "project/results/pacer_xr_emergency_rank4_5_v01"
ANC = ROOT / "project/results/pacer_xr_emergency_anchors3_v01"
REF = ROOT / "project/results/pacer_xr_rerun_v01/references"
OUT = ROOT / "project/results/pacer_xr_emergency_top5_plus_anchors_v01"
OUT.mkdir(parents=True, exist_ok=True)

def load_pair(root):
    g = pd.read_csv(root / "glide/glide_channels.csv")
    v = pd.read_csv(root / "vina/vina_channels.csv")

    keys = ["candidate_id", "canonical_smiles", "rerun_final_rank"]

    return g.merge(
        v,
        on=keys,
        validate="one_to_one"
    )

# 1–3 + 4–5 + anchors
parts = [
    load_pair(ACC),
    load_pair(R45),
    load_pair(ANC),
]

df = pd.concat(parts, ignore_index=True)

# Exact expected pool
expected_ranks = {1, 2, 3, 4, 5, 17, 60, 114}

assert len(df) == 8, df
assert df["candidate_id"].nunique() == 8, df
assert set(df["rerun_final_rank"].astype(int)) == expected_ranks, df

refs = {
    "Glide_PDB": (
        REF / "PDB/M4R_PDB_Glide_scores.csv",
        "glide_gscore",
        "Glide_PDB_raw",
    ),
    "Glide_BEmin": (
        REF / "Ensemble/M4R_Ensemble_Glide_BEmin_ranked.csv",
        "BE_min",
        "Glide_BEmin_raw",
    ),
    "Glide_BEavg": (
        REF / "Ensemble/M4R_Ensemble_Glide_BEavg_ranked.csv",
        "BE_avg",
        "Glide_BEavg_raw",
    ),
    "Vina_PDB": (
        REF / "PDB/M4R_PDB_Vina_scores.csv",
        "vina_score",
        "Vina_PDB_raw",
    ),
    "Vina_BEmin": (
        REF / "Ensemble/M4R_Ensemble_Vina_BEmin_ranked.csv",
        "BE_min",
        "Vina_BEmin_raw",
    ),
    "Vina_BEavg": (
        REF / "Ensemble/M4R_Ensemble_Vina_BEavg_ranked.csv",
        "BE_avg",
        "Vina_BEavg_raw",
    ),
}

# Frozen semantics:
# lower raw energy = better
# mid-rank tie handling
# 1.0 = best end of frozen reference distribution
def empirical_rankpct(score, reference):
    a = np.asarray(reference, dtype=float)
    a = a[np.isfinite(a)]

    n = len(a)
    better = np.sum(a < score)
    equal = np.sum(a == score)

    midrank = 1.0 + better + 0.5 * equal

    return float(
        1.0 - (midrank - 1.0) / n
    )

pct_cols = []

for channel, (path, ref_col, raw_col) in refs.items():

    ref_frame = pd.read_csv(path)

    if ref_col not in ref_frame.columns:
        raise RuntimeError(
            f"{channel}: expected {ref_col}; "
            f"found {list(ref_frame.columns)}"
        )

    reference = pd.to_numeric(
        ref_frame[ref_col],
        errors="coerce"
    ).dropna().to_numpy()

    outcol = channel + "_rankpct"

    df[outcol] = [
        empirical_rankpct(float(x), reference)
        for x in df[raw_col]
    ]

    pct_cols.append(outcol)

# Six-channel PACER-XR
df["PACER_XR"] = df[pct_cols].mean(axis=1)

# rank within this emergency comparison set only
df = df.sort_values(
    ["PACER_XR", "rerun_final_rank"],
    ascending=[False, True]
).reset_index(drop=True)

df["PACER_XR_comparison_rank"] = np.arange(1, len(df) + 1)

# DrugCLIP natural top5 indicator
df["group"] = np.where(
    df["rerun_final_rank"] <= 5,
    "DrugCLIP_top5",
    "historical_anchor"
)

# Legacy labels — do NOT assign PACER0027 until identity conflict is resolved
legacy = {
    "PACERGEN01755": "PACER0010",
    "PACERGEN00094": "PACER0073",
    "PACERGEN00123": "PACER0027_MAPPING_UNRESOLVED",
}

df["legacy_id"] = df["candidate_id"].map(legacy).fillna("")

out_csv = OUT / "pacer_xr_top5_plus_anchors_ranking.csv"
df.to_csv(out_csv, index=False)

show = [
    "PACER_XR_comparison_rank",
    "candidate_id",
    "legacy_id",
    "group",
    "rerun_final_rank",
    "Glide_PDB_rankpct",
    "Glide_BEmin_rankpct",
    "Glide_BEavg_rankpct",
    "Vina_PDB_rankpct",
    "Vina_BEmin_rankpct",
    "Vina_BEavg_rankpct",
    "PACER_XR",
]

print(df[show].to_string(index=False))

print()
print("WROTE:", out_csv)