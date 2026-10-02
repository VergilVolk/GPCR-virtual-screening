#!/usr/bin/env python
"""Static dependency coverage for explicit production roots and import closure.

Reads Python ASTs and declarations. Never imports or executes scientific code.
The report explicitly separates mapping completeness from deployment readiness.
"""
from __future__ import annotations

import argparse
import ast
from collections import defaultdict
import hashlib
import json
from pathlib import Path
import pkgutil
import re
import subprocess
import sys
import sysconfig

ROOT = Path(__file__).resolve().parents[2]
MAP = ROOT / "project/config/runtime_dependency_map_v01.json"


def normalized(value):
    return re.sub(r"[-_.]+", "-", value).lower()


def imports(tree):
    """Capture nested imports, literal dynamic imports, and fallback context."""
    rows = []

    class Visitor(ast.NodeVisitor):
        optional = False
        type_only = False

        def visit_Import(self, node):
            for alias in node.names:
                rows.append((alias.name, [], 0, node.lineno, self.optional, self.type_only))

        def visit_ImportFrom(self, node):
            rows.append((node.module or "", [a.name for a in node.names], node.level,
                         node.lineno, self.optional, self.type_only))

        def visit_Try(self, node):
            old = self.optional
            caught = any(isinstance(h.type, ast.Name) and h.type.id in
                         {"ImportError", "ModuleNotFoundError"} for h in node.handlers)
            self.optional = old or caught
            for child in node.body:
                self.visit(child)
            self.optional = old
            for child in node.handlers + node.orelse + node.finalbody:
                self.visit(child)

        def visit_If(self, node):
            old = self.type_only
            test = ast.unparse(node.test)
            self.type_only = old or test.endswith("TYPE_CHECKING")
            for child in node.body:
                self.visit(child)
            self.type_only = old
            for child in node.orelse:
                self.visit(child)

        def visit_Call(self, node):
            function = ast.unparse(node.func)
            if function in {"__import__", "importlib.import_module"} and node.args:
                if isinstance(node.args[0], ast.Constant) and isinstance(node.args[0].value, str):
                    rows.append((node.args[0].value, [], 0, node.lineno, self.optional, self.type_only))
            self.generic_visit(node)

    Visitor().visit(tree)
    return rows


def local_file(module, parent):
    parts = module.split(".") if module else []
    for base in [ROOT, parent, ROOT / "project/scripts", ROOT / "project/pacer_fkg_v02"]:
        p = base.joinpath(*parts)
        for candidate in [p.with_suffix(".py"), p / "__init__.py"]:
            if candidate.is_file():
                return candidate.resolve()
    return None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "project/results/dependency_audit_v01.json")
    parser.add_argument("--geom2vec-source", type=Path, help="Existing source root for vendor import scan")
    args = parser.parse_args()
    catalog = json.loads(MAP.read_text(encoding="utf-8"))
    stdlib = set(getattr(sys, "stdlib_module_names", ())) | set(sys.builtin_module_names)
    if not getattr(sys, "stdlib_module_names", None):
        stdlib.update(x.name for x in pkgutil.iter_modules([sysconfig.get_path("stdlib")]))
    stdlib.add("__future__")
    todo = [(ROOT / name, info["stage"], name in catalog["optional_entrypoints"], "entrypoint")
            for name, info in catalog["entrypoints"].items()]
    vendor_files = catalog["vendor_scan_files"]
    if args.geom2vec_source:
        vendor_files = [str(args.geom2vec_source / "src/geom2vec" / Path(x.replace('\\', '/')).relative_to('C:/projects/geom2vec-source/src/geom2vec')) for x in vendor_files]
    todo += [(Path(x), "pacer_fkg", False, "geom2vec_source_import_surface") for x in vendor_files]
    seen = set()
    scanned, missing, syntax_errors = [], [], []
    third, unmapped, standard, internal, type_only = [], [], [], [], []
    while todo:
        path, stage, optional_root, reason = todo.pop(0)
        key = (str(path.resolve()), stage, optional_root)
        if key in seen:
            continue
        seen.add(key)
        label = path.relative_to(ROOT).as_posix() if path.is_relative_to(ROOT) else str(path)
        if not path.is_file():
            missing.append({"path": label, "stage": stage, "reason": reason})
            continue
        data = path.read_bytes()
        try:
            tree = ast.parse(data.decode("utf-8-sig"), filename=label)
        except (SyntaxError, UnicodeError) as exc:
            syntax_errors.append({"path": label, "error": str(exc)})
            continue
        scanned.append({"path": label, "stage": stage, "optional_entrypoint": optional_root,
                        "reason": reason, "sha256": hashlib.sha256(data).hexdigest()})
        for module, aliases, level, line, fallback, types in imports(tree):
            record = {"import": module, "file": label, "line": line,
                      "stage": stage, "optional": optional_root or fallback, "fallback": fallback}
            if types:
                type_only.append(record)
                continue
            base = module.split(".")[0]
            local = None
            if level:
                p = path.parent
                for _ in range(level - 1):
                    p = p.parent
                local = local_file(module, p)
                # Relative imports are internal to their project or vendor package.
                internal.append(record)
                if path.is_relative_to(ROOT) and local:
                    todo.append((local, stage, optional_root, "relative_import_closure"))
                continue
            if base in stdlib:
                standard.append(record)
                continue
            if path.is_relative_to(ROOT):
                local = local_file(module, path.parent)
            if local or base == "project":
                namespace_exists = ROOT.joinpath(*module.split(".")).is_dir()
                if base == "project" and not local and not namespace_exists:
                    unmapped.append(dict(record, kind="missing_internal_module"))
                    continue
                internal.append(record)
                if local:
                    todo.append((local, stage, optional_root, "internal_import_closure"))
                for alias in aliases:
                    child = local_file(module + "." + alias, path.parent)
                    if child:
                        todo.append((child, stage, optional_root, "from_import_closure"))
                continue
            spec = catalog["imports"].get(module) or catalog["imports"].get(base)
            if not spec:
                unmapped.append(record)
                continue
            record.update(distribution=spec["distribution"], declarations=spec["declarations"],
                          kind=spec.get("kind", "python_package"))
            third.append(record)
        for helper in catalog["explicit_subprocess_helpers"].get(label, []):
            todo.append((ROOT / helper, stage, optional_root, "subprocess_script_closure"))

    used = defaultdict(list)
    for record in third:
        used[record["distribution"]].append(record)
    mapped = []
    declaration_errors = []
    for distribution, records in sorted(used.items()):
        spec = next(v for v in catalog["imports"].values() if v["distribution"] == distribution)
        declarations = spec["declarations"]
        if spec.get("kind") != "local_source":
            present = False
            for filename in declarations:
                path = ROOT / filename
                if path.is_file():
                    for line in path.read_text(encoding="utf-8").splitlines():
                        match = re.match(r"\s*(?:-\s+)?([A-Za-z][A-Za-z0-9_.-]*)\s*(?:[=<>!~]|$)", line)
                        declared_names = spec.get("declaration_names", [distribution])
                        if match and normalized(match[1]) in {normalized(x) for x in declared_names}:
                            present = True
            if not present:
                declaration_errors.append(distribution)
        mapped.append({"distribution": distribution, "imports": sorted({r["import"] for r in records}),
                       "usage_modules": sorted({r["file"] for r in records}),
                       "required": any(not r["optional"] for r in records),
                       "declarations": declarations, "kind": spec.get("kind", "python_package"),
                       "version_evidence": spec["version_evidence"],
                       "metadata_evidence": spec["metadata_evidence"]})
    optional_metadata = [{"distribution": v["distribution"], "required": False,
                          "kind": "metadata_only", "declarations": v["declarations"],
                          "version_evidence": v["version_evidence"],
                          "metadata_evidence": v["metadata_evidence"]}
                         for v in catalog["imports"].values() if v.get("role") == "optional"
                         and v["distribution"] not in used]
    scanned_paths = {x["path"] for x in scanned}
    excluded = [{"path": p.relative_to(ROOT).as_posix(),
                 "reason": "not in production roots/import closure; historical, experimental, training, audit or test entrypoint"}
                for directory in [ROOT / "project/scripts", ROOT / "project/pacer_fkg_v02"]
                for p in sorted(directory.glob("*.py")) if p.relative_to(ROOT).as_posix() not in scanned_paths
                and p.name not in {"check_runtime_dependencies.py", "audit_runtime_dependencies.py"}]
    missing_required = [x for x in missing if x["reason"] != "geom2vec_source_import_surface"]
    unmapped_required = [r for r in unmapped if not r["optional"]]
    head = subprocess.check_output(["git", "-C", str(ROOT), "rev-parse", "HEAD"], text=True).strip()
    smoke_reports = {}
    for profile in ["core", "drugclip", "md", "fkg", "overview"]:
        receipt = ROOT / f"project/results/dependency_smoke_{profile}_v01.json"
        if receipt.is_file():
            smoke = json.loads(receipt.read_text(encoding="utf-8"))
            smoke_reports[profile] = {"path": receipt.relative_to(ROOT).as_posix(),
                                      "python": smoke["python"], "overall": smoke["overall"],
                                      "failed_imports": [k for k,v in smoke["python_dependencies"].items() if v["status"] == "FAIL"]}
    report = {"schema": "pacer.dependency_audit.v1", "git_head": head,
              "scope_notes": catalog["scope_notes"], "scanned_files": sorted(scanned, key=lambda r:(r['path'],r['stage'])),
              "scanned_file_count": len(scanned_paths), "scan_context_count": len(scanned),
              "third_party_imports": third, "mapped_dependencies": mapped,
              "mapped_import_count": len({r['import'] for r in third}),
              "unmapped_imports": unmapped, "unmapped_required_imports": unmapped_required,
              "unmapped_required_import_count": len(unmapped_required),
              "optional_dependencies": [d for d in mapped if not d["required"]] + optional_metadata,
              "optional_imports": [r for r in third if r["optional"]],
              "external_tools": catalog["external_tools"], "local_source_dependencies": catalog["local_sources"],
              "standard_library_imports": sorted({r['import'] for r in standard}),
              "internal_module_imports": internal, "type_only_imports": type_only,
              "excluded_entrypoints": excluded, "missing_files": missing,
              "syntax_errors": syntax_errors, "declaration_errors": declaration_errors,
              "coverage": "PASS" if not (unmapped_required or missing_required or syntax_errors or declaration_errors) else "FAIL",
              "deployment_readiness": "PARTIAL",
              "runtime_evidence": "project/results/dependency_runtime_evidence_v01.json",
              "scientific_assets": "project/results/dependency_scientific_assets_v01.json",
              "smoke_reports": smoke_reports,
              "python_dependency_counts": {
                  "direct_required_ast_distributions": sum(d["required"] and d["kind"] != "local_source" for d in mapped),
                  "optional_ast_distributions": sum(not d["required"] for d in mapped),
                  "optional_metadata_distributions": len(optional_metadata),
                  "optional_total_distributions": sum(not d["required"] for d in mapped) + len(optional_metadata),
                  "core_direct_or_necessary_requirements": 11},
              "validation_limit": "Mapping coverage is not proof of a clean environment solve or scientific reproducibility."}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({k:report[k] for k in ['scanned_file_count','mapped_import_count','unmapped_required_import_count','missing_files','syntax_errors','declaration_errors','coverage']},indent=2))
    return 0 if report["coverage"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
