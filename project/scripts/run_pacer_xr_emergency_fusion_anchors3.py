from pathlib import Path
import pandas as pd
import numpy as np

ROOT = Path(r"C:\projects\GPCR-virtual-screening")

E = ROOT / "project/results/pacer_xr_emergency_anchors3_v01"
R = ROOT / "project/results/pacer_xr_rerun_v01/references"

g = pd.read_csv(E / "glide/glide_channels.csv")
v = pd.read_csv(E / "vina/vina_channels.csv")

keys = ["candidate_id", "canonical_smiles", "rerun_final_rank"]
df = g.merge(v, on=keys, validate="one_to_one")

refs = {
    "Glide_PDB": (
        R / "PDB/M4R_PDB_Glide_scores.csv",
        "glide_gscore",
        "Glide_PDB_raw",
    ),
    "Glide_BEmin": (
        R / "Ensemble/M4R_Ensemble_Glide_BEmin_ranked.csv",
        "BE_min",
        "Glide_BEmin_raw",
    ),
    "Glide_BEavg": (
        R / "Ensemble/M4R_Ensemble_Glide_BEavg_ranked.csv",
        "BE_avg",
        "Glide_BEavg_raw",
    ),
    "Vina_PDB": (
        R / "PDB/M4R_PDB_Vina_scores.csv",
        "vina_score",
        "Vina_PDB_raw",
    ),
    "Vina_BEmin": (
        R / "Ensemble/M4R_Ensemble_Vina_BEmin_ranked.csv",
        "BE_min",
        "Vina_BEmin_raw",
    ),
    "Vina_BEavg": (
        R / "Ensemble/M4R_Ensemble_Vina_BEavg_ranked.csv",
        "BE_avg",
        "Vina_BEavg_raw",
    ),
}

# Frozen PACER-XR convention:
# lower raw score = better;
# empirical mid-rank against frozen reference;
# 1.0 = best end of the reference distribution.
def empirical_rankpct(score, reference):
    a = np.asarray(reference, dtype=float)
    a = a[np.isfinite(a)]
    n = len(a)

    better = np.sum(a < score)
    equal = np.sum(a == score)

    # candidate empirical mid-rank
    rank = 1.0 + better + 0.5 * equal

    # normalize such that the best end is 1 and the worst end tends to 0
    pct = 1.0 - (rank - 1.0) / n

    return float(pct)

pct_cols = []

for channel, (path, ref_col, raw_col) in refs.items():
    r = pd.read_csv(path)

    if ref_col not in r.columns:
        raise RuntimeError(
            f"{channel}: expected reference column {ref_col!r}; "
            f"found {list(r.columns)}"
        )

    ref = pd.to_numeric(r[ref_col], errors="coerce").dropna().to_numpy()

    outcol = channel + "_rankpct"
    df[outcol] = [
        empirical_rankpct(float(x), ref)
        for x in df[raw_col]
    ]
    pct_cols.append(outcol)

df["PACER_XR"] = df[pct_cols].mean(axis=1)

df = df.sort_values(
    ["PACER_XR", "rerun_final_rank"],
    ascending=[False, True]
).reset_index(drop=True)

df["PACER_XR_anchor_rank"] = np.arange(1, len(df) + 1)

out = E / "anchor_pacer_xr_ranking.csv"
df.to_csv(out, index=False)

show = [
    "candidate_id",
    "rerun_final_rank",
    *pct_cols,
    "PACER_XR",
    "PACER_XR_anchor_rank",
]

print(df[show].to_string(index=False))
print()
print("WROTE:", out)