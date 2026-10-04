from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd

def main() -> None:
    p = argparse.ArgumentParser(description="Build a local HTML review panel")
    p.add_argument("--results", type=Path, required=True)
    p.add_argument("--output", type=Path, default=Path("run/dashboard.html"))
    args = p.parse_args()
    table = pd.read_csv(args.results)
    columns = [c for c in ["final_rank", "candidate_id", "binding_rank", "pacer_xr_score",
                           "four_context_status", "functional_interpretation"] if c in table]
    html = """<!doctype html><meta charset='utf-8'><title>M4 screening panel</title>
<style>body{font:15px Arial;margin:36px;color:#222}table{border-collapse:collapse;width:100%}
th,td{padding:9px;border-bottom:1px solid #ddd;text-align:left}th{background:#eee}</style>
<h1>M4 virtual screening panel</h1>""" + table[columns].to_html(index=False, border=0)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(html, encoding="utf-8")
    print(args.output)

if __name__ == "__main__":
    main()
