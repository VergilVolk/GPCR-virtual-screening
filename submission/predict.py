from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd
from module1_generate_filter import evaluate_library
from module3_ensemble_dock import fuse_six_channels
from module4_fkg_eval import export_evidence
from screen import run_screen

ROOT = Path(__file__).resolve().parent

def run_pipeline(library: Path, molecule_repr: Path, pocket_repr: Path,
                 docking: Path, output_dir: Path, output: Path,
                 stage4_summary: Path | None = None,
                 stage4_id_map: Path | None = None) -> pd.DataFrame:
    output_dir.mkdir(parents=True, exist_ok=True)
    chemistry = evaluate_library(library, output_dir / "module1_chemistry.csv")
    run_screen(molecule_repr, pocket_repr, output_dir / "module2_binding.csv")
    binding = pd.read_csv(output_dir / "module2_binding.csv").rename(columns={"molecule_id": "candidate_id"})
    structure = fuse_six_channels(docking, output_dir / "module3_structure.csv")
    result = structure.merge(chemistry, on="candidate_id", how="left", validate="one_to_one")
    result = result.merge(binding, on="candidate_id", how="left", validate="one_to_one")
    result["stage4_candidate_id"] = result["candidate_id"]
    if stage4_id_map:
        mapping = pd.read_csv(stage4_id_map)
        result = result.drop(columns="stage4_candidate_id").merge(mapping, on="candidate_id", how="left")
    if stage4_summary:
        function = export_evidence(stage4_summary, output_dir / "module4_function.csv")
        result = result.merge(function, left_on="stage4_candidate_id", right_on="candidate_id",
                              how="left", suffixes=("", "_stage4"))
    result["four_context_status"] = result.get("four_context_status", pd.Series(index=result.index, dtype=object)).fillna("not_run")
    result["functional_interpretation"] = result.get("functional_interpretation", pd.Series(index=result.index, dtype=object)).fillna("not evaluated by four-context MD")
    result = result.sort_values("pacer_xr_score", ascending=False).reset_index(drop=True)
    result.insert(0, "final_rank", range(1, len(result) + 1))
    output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output, index=False)
    return result

def main() -> None:
    p = argparse.ArgumentParser(description="Run the four-module M4 inference workflow")
    p.add_argument("--library", type=Path, required=True)
    p.add_argument("--molecule-representations", type=Path, required=True)
    p.add_argument("--pocket-representation", type=Path, required=True)
    p.add_argument("--docking-six-channel", type=Path, required=True)
    p.add_argument("--stage4-summary", type=Path)
    p.add_argument("--stage4-id-map", type=Path)
    p.add_argument("--output-dir", type=Path, default=ROOT / "run")
    p.add_argument("--output", type=Path, default=ROOT / "run/results.csv")
    args = p.parse_args()
    result = run_pipeline(args.library, args.molecule_representations,
                          args.pocket_representation, args.docking_six_channel,
                          args.output_dir, args.output, args.stage4_summary,
                          args.stage4_id_map)
    print(result[["final_rank", "candidate_id", "binding_rank", "pacer_xr_score", "four_context_status"]].to_string(index=False))

if __name__ == "__main__":
    main()
