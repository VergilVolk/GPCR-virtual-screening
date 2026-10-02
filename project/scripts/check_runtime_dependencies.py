#!/usr/bin/env python
"""Import and executable-presence checks only. Never import pipeline entrypoints.

Use --profile core|drugclip|md|fkg with that stage's own Python interpreter.
Optional --source NAME=PATH supplies existing DrugCLIP/Uni-Core/Geom2Vec roots.
No assets are opened, no predictions/simulations run, and no packages installed.
"""
from __future__ import annotations

import argparse
import contextlib
import importlib
import importlib.metadata as metadata
import io
import json
import os
from pathlib import Path
import platform
import shutil
import sys

ROOT = Path(__file__).resolve().parents[2]
CATALOG = ROOT / "project/config/runtime_dependency_map_v01.json"


def check_import(name, spec):
    log = io.StringIO()
    try:
        with contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
            module = importlib.import_module(spec.get("smoke_import", name))
        try:
            version = metadata.version(spec["distribution"])
        except metadata.PackageNotFoundError:
            version = getattr(module, "__version__", None)
        result = {"status": "PASS", "version": version,
                  "distribution": spec["distribution"],
                  "module_path": getattr(module, "__file__", None)}
    except Exception as exc:
        result = {"status": "FAIL", "version": None,
                  "error": f"{type(exc).__name__}: {exc}"}
    if log.getvalue():
        result["import_messages"] = log.getvalue()[-2000:]
    result["role"] = spec.get("role", "runtime")
    return result


def tool_presence(name, spec):
    candidates = []
    is_windows = os.name == "nt"
    if name == "vina" and is_windows:
        candidates.append(ROOT / "project/tools/vina_1.2.7.exe")
    bindir = Path(sys.prefix) / ("Scripts" if is_windows else "bin")
    for command in spec["commands"]:
        if not is_windows and command.endswith(".exe"):
            continue
        found = shutil.which(command)
        if found:
            candidates.append(Path(found))
        candidates.append(bindir / command)
        if is_windows:
            candidates.append(bindir / (command + ".exe"))
    found = next((p for p in candidates if p.is_file()
                  and (is_windows or os.access(p, os.X_OK))), None)
    return {"status": "PASS" if found else "FAIL", "path": str(found) if found else None,
            "version": spec.get("recorded_version"),
            "version_evidence": "recorded baseline; executable not invoked",
            "role": spec.get("role", "runtime")}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", choices=["all", "core", "drugclip", "md", "fkg"], default="all")
    parser.add_argument("--source", action="append", default=[], metavar="NAME=PATH")
    parser.add_argument("--output", type=Path, help="Optional dependency-check JSON output")
    args = parser.parse_args()
    catalog = json.loads(CATALOG.read_text(encoding="utf-8"))
    for item in args.source:
        name, separator, path = item.partition("=")
        if not separator or name not in {"drugclip", "unicore", "geom2vec"}:
            parser.error("--source must be drugclip=PATH, unicore=PATH or geom2vec=PATH")
        root = Path(path).resolve()
        sys.path.insert(0, str(root / "src" if (root / "src").is_dir() else root))

    def selected(spec):
        return args.profile == "all" or args.profile in spec["profiles"]

    deps = {name: check_import(name, spec) for name, spec in catalog["imports"].items()
            if selected(spec)}
    tools = {name: tool_presence(name, spec) for name, spec in catalog["external_tools"].items()
             if selected(spec)}
    required = [r for r in list(deps.values()) + list(tools.values()) if r["role"] != "optional"]
    passed = sum(r["status"] == "PASS" for r in required)
    expected = catalog["python_versions"]
    python_ok = args.profile == "all" or platform.python_version() == expected[args.profile]
    overall = "PASS" if passed == len(required) and python_ok else "PARTIAL" if passed else "FAIL"
    result = {"schema": "pacer.runtime_dependencies.v1", "profile": args.profile,
              "python": platform.python_version(), "python_executable": sys.executable,
              "python_version_check": {"status": "PASS" if python_ok else "FAIL",
                                       "expected": expected.get(args.profile, "separate stage runtimes")},
              "python_dependencies": deps, "external_tools": tools, "overall": overall,
              "limitation": "Imports/presence only; no clean solve, native ABI, GPU platform or asset verification."}
    output = json.dumps(result, indent=2, ensure_ascii=False)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(output + "\n", encoding="utf-8")
    print(output)
    return 0 if overall == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
