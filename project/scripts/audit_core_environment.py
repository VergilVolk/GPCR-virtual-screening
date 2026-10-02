#!/usr/bin/env python3
"""Audit the canonical CPython 3.9 PACER-M4 core environment."""
from __future__ import annotations

import argparse
import importlib
import json
import py_compile
import sys
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any


PROJECT = Path(__file__).resolve().parents[1]
ROOT = PROJECT.parent
EXPECTED = {
    "numpy": ("numpy", "1.26.4"),
    "pandas": ("pandas", "2.3.1"),
    "scipy": ("scipy", "1.13.1"),
    "scikit-learn": ("sklearn", "1.6.1"),
    "rdkit": ("rdkit", "2025.3.5"),
    "torch": ("torch", "2.8.0"),
    "lightgbm": ("lightgbm", "4.6.0"),
    "meeko": ("meeko", "0.7.1"),
    "molscrub": ("scrubber", "0.1.1"),
    "gemmi": ("gemmi", "0.7.3"),
    "matplotlib": ("matplotlib", "3.9.4"),
    "reportlab": ("reportlab", "4.2.5"),
}
CORE_SCRIPTS = [
    "project/scripts/fragment_generation.py",
    "project/scripts/score_candidates.py",
    "project/scripts/select_pareto_candidates.py",
    "project/scripts/export_candidate_sdf.py",
    "project/scripts/dock_candidate_portfolio.py",
    "project/scripts/dock_potency_coupling_features.py",
    "project/scripts/dock_m4_gamd_ensemble.py",
    "project/scripts/score_pacer200_drugclip2023_m4_loto.py",
    "project/scripts/pacer_drugclip_router.py",
    "project/scripts/build_pacer_m4_algorithm_report_v01.py",
]


def audit() -> dict[str, Any]:
    packages: dict[str, Any] = {}
    for distribution, (module, expected_version) in EXPECTED.items():
        try:
            installed = version(distribution)
            importlib.import_module(module)
            import_ok = True
            error = None
        except PackageNotFoundError:
            installed, import_ok, error = None, False, "distribution_not_installed"
        except Exception as exc:  # retain exact import failure for handoff diagnosis
            installed, import_ok, error = None, False, f"{type(exc).__name__}: {exc}"
        packages[distribution] = {
            "module": module,
            "expected": expected_version,
            "installed": installed,
            "version_match": (
                installed is not None and installed.split("+", 1)[0] == expected_version
            ),
            "import_ok": import_ok,
            "error": error,
        }

    scripts: dict[str, Any] = {}
    for relative in CORE_SCRIPTS:
        path = ROOT / relative
        try:
            py_compile.compile(str(path), doraise=True)
            scripts[relative] = {"exists": True, "syntax_ok": True, "error": None}
        except Exception as exc:
            scripts[relative] = {
                "exists": path.exists(),
                "syntax_ok": False,
                "error": f"{type(exc).__name__}: {exc}",
            }

    python_match = sys.version_info[:2] == (3, 9)
    packages_ok = all(
        item["version_match"] and item["import_ok"] for item in packages.values()
    )
    scripts_ok = all(item["syntax_ok"] for item in scripts.values())
    return {
        "protocol": "pacer_m4_core_environment_v01",
        "python": sys.version,
        "expected_python": "3.9.x (validated environment: 3.9.23)",
        "python_match": python_match,
        "packages": packages,
        "scripts": scripts,
        "core_environment_valid": python_match and packages_ok and scripts_ok,
        "external_assets": {
            "vina_binary_in_requirements": False,
            "model_assets_in_requirements": False,
            "api_dependencies_in_requirements": False,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument(
        "--report-only", action="store_true", help="Return zero even when the environment differs"
    )
    args = parser.parse_args()
    report = audit()
    text = json.dumps(report, indent=2, ensure_ascii=False)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text + "\n", encoding="utf-8")
    print(text)
    if not args.report_only and not report["core_environment_valid"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
