"""Create the frozen Function-Space v3 evidence figure."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


PROJECT = Path(__file__).resolve().parents[1]
OUT = PROJECT / "results" / "pacer_function_space_v3_figure"
ORDER = ["A", "B", "C", "D"]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    metrics = pd.read_csv(PROJECT / "results" / "pacer_external_us20260055116_v01" / "metrics.csv")
    ext = pd.read_csv(PROJECT / "results" / "pacer_external_us20260055116_v01" /
                      "external_unique_active_moieties.csv")
    cliffs = pd.read_csv(PROJECT / "results" / "us20260055116_functional_divergence_v01" /
                         "local_functional_cliffs.csv")
    context = pd.read_csv(PROJECT / "results" / "pacer_context_kernel_v01" / "metrics.csv")

    sns.set_theme(style="whitegrid", context="talk")
    fig = plt.figure(figsize=(16, 12), constrained_layout=True)
    grid = fig.add_gridspec(2, 2)

    ax = fig.add_subplot(grid[0, 0])
    seven = metrics[metrics.protocol == "seven_anchor"].sort_values("ordinal_concordance")
    colors = ["#C44E52" if "DeltaSAR-7" in x else "#4C72B0" for x in seven.method]
    ax.barh(seven.method.str.replace("+7AnchorOffset", "+offset", regex=False)
            .str.replace("PACER-FS-", "", regex=False),
            seven.ordinal_concordance, color=colors)
    ax.axvline(.5, color="black", linestyle="--", linewidth=1)
    ax.set_xlim(.2, .65); ax.set_xlabel("Ordinal concordance")
    ax.set_title("A  Fifth-campaign seven-anchor stress")

    ax = fig.add_subplot(grid[0, 1])
    table = pd.crosstab(ext.human_m4_perk, ext.rat_m4_perk).reindex(
        index=ORDER, columns=ORDER, fill_value=0)
    sns.heatmap(table, annot=True, fmt="d", cmap="Blues", cbar=False, ax=ax)
    ax.set_xlabel("Rat M4 pERK bin"); ax.set_ylabel("Human M4 pERK bin")
    ax.set_title("B  Species shift: 21/22 non-ties favor human")

    ax = fig.add_subplot(grid[1, 0])
    shift = cliffs.max_abs_endpoint_bin_delta.value_counts().sort_index()
    ax.bar(shift.index.astype(str), shift.values, color="#55A868")
    ax.set_xlabel("Maximum endpoint shift within local pair (bins)")
    ax.set_ylabel("Local molecular pairs")
    ax.set_title("C  Local functional cliffs (23 pairs >=2 bins)")

    ax = fig.add_subplot(grid[1, 1])
    macro = context[context.scope == "macro"].sort_values("ordinal_concordance")
    labels = macro.method.str.replace("-KRR", "", regex=False)
    ax.barh(labels, macro.ordinal_concordance, color=["#8172B3", "#CCB974", "#64B5CD"])
    ax.axvline(.5, color="black", linestyle="--", linewidth=1)
    ax.set_xlim(.45, .78); ax.set_xlabel("Macro ordinal concordance")
    ax.set_title("D  Context model did not beat pooled baseline")

    fig.suptitle("PACER Function-Space v3: PAM identity is conditional, not a docking score",
                 fontsize=22, fontweight="bold")
    fig.savefig(OUT / "PACER_FUNCTION_SPACE_V3_EVIDENCE.png", dpi=220, bbox_inches="tight")
    fig.savefig(OUT / "PACER_FUNCTION_SPACE_V3_EVIDENCE.svg", bbox_inches="tight")
    plt.close(fig)
    print(OUT / "PACER_FUNCTION_SPACE_V3_EVIDENCE.png")


if __name__ == "__main__":
    main()
