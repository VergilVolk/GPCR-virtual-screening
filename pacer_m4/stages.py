"""Frozen registry of PACER-M4 workflow entry points.

The registry does not merge scientifically distinct scores.  It provides one
discoverable interface to the existing audited scripts while preserving each
stage's own inputs, outputs and claim boundary.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from importlib.util import find_spec
from pathlib import Path


@dataclass(frozen=True)
class StageSpec:
    stage_id: str
    script: str
    category: str
    description: str
    claim_boundary: str
    hardware: str = "CPU"
    optional_assets: tuple[str, ...] = ()
    external_requirement: str | None = None
    python_modules: tuple[str, ...] = ()

    def to_dict(self, root: Path) -> dict[str, object]:
        result = asdict(self)
        result["script_exists"] = (root / self.script).is_file()
        result["optional_asset_status"] = {
            asset: (root / asset).exists() for asset in self.optional_assets
        }
        result["python_module_status"] = {
            module: find_spec(module) is not None for module in self.python_modules
        }
        result["code_ready"] = bool(result["script_exists"])
        python_ready = bool(result["script_exists"]) and all(
            result["python_module_status"].values()
        )
        result["python_ready"] = python_ready
        if not python_ready:
            result["execution_readiness"] = "missing_code_or_python_module"
        elif self.external_requirement:
            result["execution_readiness"] = "requires_external_assets_or_software"
        else:
            result["execution_readiness"] = "ready"
        return result


def _stage(
    stage_id: str,
    script: str,
    category: str,
    description: str,
    claim_boundary: str,
    *,
    hardware: str = "CPU",
    optional_assets: tuple[str, ...] = (),
    external_requirement: str | None = None,
    python_modules: tuple[str, ...] = (),
) -> StageSpec:
    return StageSpec(
        stage_id=stage_id,
        script=script,
        category=category,
        description=description,
        claim_boundary=claim_boundary,
        hardware=hardware,
        optional_assets=optional_assets,
        external_requirement=external_requirement,
        python_modules=python_modules,
    )


STAGES = {
    spec.stage_id: spec
    for spec in [
        _stage(
            "validate-release", "project/scripts/validate_pacer_m4_release.py", "repository",
            "Validate the reviewer-facing source release and optional artifact bundle.",
            "Repository completeness check; not a scientific benchmark.",
        ),
        _stage(
            "drugclip-score-pacer200", "project/scripts/score_pacer200_drugclip2023_m4_loto.py", "binding",
            "Recompute the packaged three-seed M4-held-out DrugCLIP adapter scores.",
            "Binding retrieval only; not PAM function or potency.",
            optional_assets=("project/artifacts/drugclip2023_m4_loto_pacer200_v01",),
            python_modules=("numpy", "pandas", "torch"),
        ),
        _stage(
            "drugclip-route", "project/scripts/pacer_drugclip_router.py", "binding",
            "Apply the frozen M4-safe 2023/2026 DrugCLIP routing rule.",
            "Binding retrieval only; not PAM function or potency.",
            python_modules=("numpy", "pandas"),
        ),
        _stage(
            "drugclip-benchmark", "project/scripts/run_drugclip_generation_bakeoff.py", "benchmark",
            "Run the common-data DrugCLIP generation and fine-tuning benchmark.",
            "Retrospective target-level benchmark; not a prospective PAM claim.",
            python_modules=("numpy", "pandas", "rdkit"),
        ),
        _stage(
            "dock-single", "project/scripts/dock_candidate_portfolio.py", "docking",
            "Dock candidates against the single 7TRS ACh-state structure.",
            "Pose proposal and pocket compatibility only.",
            external_requirement="AutoDock Vina and receptor assets",
            python_modules=("numpy", "pandas", "rdkit"),
        ),
        _stage(
            "dock-three-state", "project/scripts/dock_potency_coupling_features.py", "docking",
            "Extract candidate features across 7TRQ, 7TRP and 7TRS receptor states.",
            "Static state features do not establish functional PAM efficacy.",
            external_requirement="AutoDock Vina and three prepared receptor states",
            python_modules=("numpy", "pandas", "rdkit"),
        ),
        _stage(
            "dock-gamd-ensemble", "project/scripts/run_m4_gamd_complete.py", "docking",
            "Run the ten-conformation M4 GaMD ensemble docking workflow.",
            "Ensemble binding compatibility only; not functional potentiation.",
            external_requirement="AutoDock Vina and prepared GaMD representatives",
            python_modules=("numpy", "pandas", "rdkit"),
        ),
        _stage(
            "dock-xr-six-channel", "project/scripts/apply_pacer_xr_candidates.py", "docking",
            "Apply six-channel Glide/Vina PDB, BEmin and BEavg rank fusion and cascade.",
            "Broad allosteric-modulator ranking; not functional PAM efficacy.",
            external_requirement="All six matched scores; Glide channels require a legal Schrodinger installation",
            python_modules=("numpy", "pandas"),
        ),
        _stage(
            "pacer-fs", "project/scripts/run_pacer_fs_complete.py", "potency",
            "Run series-centred potency ranking and three-anchor evaluation.",
            "Few-shot series calibration; not a zero-shot PAM classifier.",
            python_modules=("numpy", "pandas", "sklearn"),
        ),
        _stage(
            "pareto-select", "project/scripts/select_pareto_candidates.py", "selection",
            "Select a diverse non-dominated candidate portfolio without a composite efficacy score.",
            "Computational hypotheses requiring functional assay confirmation.",
            python_modules=("numpy", "pandas", "rdkit"),
        ),
        _stage(
            "md-build-complexes", "project/scripts/build_pacer_dc_reference_complexes.py", "four_context_md",
            "Build matched A, P, C and CP reference complexes.",
            "System preparation only.", external_requirement="Prepared receptor and ligand parameters",
            python_modules=("openmm",),
        ),
        _stage(
            "md-build-qc", "project/scripts/build_pacer_dc_openmm_reference_systems.py", "four_context_md",
            "Build and validate non-solvated OpenMM reference systems.",
            "System quality control only.", external_requirement="OpenMM environment",
            python_modules=("openmm",),
        ),
        _stage(
            "md-build-membrane", "project/scripts/build_pacer_dc_membrane_reference.py", "four_context_md",
            "Build matched membrane systems for selected contexts.",
            "System construction only.", hardware="CPU or GPU", external_requirement="OpenMM and membrane assets",
            python_modules=("openmm",),
        ),
        _stage(
            "md-equilibrate", "project/scripts/run_pacer_dc_short_equilibration.py", "four_context_md",
            "Run restrained short equilibration for matched systems.",
            "Equilibration pilot; cannot support a PAM claim.", hardware="GPU recommended",
            external_requirement="OpenMM and completed membrane systems",
            python_modules=("openmm",),
        ),
        _stage(
            "md-release-restraints", "project/scripts/run_pacer_dc_restraint_release.py", "four_context_md",
            "Release restraints in the frozen 500 to 100 to 10 to 0 schedule.",
            "Equilibration only; cannot support a PAM claim.", hardware="GPU recommended",
            external_requirement="OpenMM and stable equilibration checkpoints",
            python_modules=("openmm",),
        ),
        _stage(
            "md-production", "project/scripts/run_pacer_dc_production_md.py", "four_context_md",
            "Run matched four-context production trajectories.",
            "Trajectories are inputs to functional analysis, not functional labels.", hardware="GPU recommended",
            external_requirement="OpenMM and released production systems",
            python_modules=("openmm",),
        ),
        _stage(
            "md-extract", "project/scripts/extract_pacer_dc_production_endpoints.py", "four_context_md",
            "Extract quality-control and dynamic endpoints from production trajectories.",
            "Endpoint extraction alone does not classify PAM function.",
            external_requirement="MDAnalysis and completed trajectories",
            python_modules=("MDAnalysis", "pandas"),
        ),
        _stage(
            "dynamic-score", "project/scripts/pacer_dc_score.py", "dynamic_analysis",
            "Score compatible multi-replica dynamic evidence tables.",
            "Evidence summary only; PAM probability remains disabled without validation.",
            python_modules=("numpy", "pandas"),
        ),
        _stage(
            "benchmark-project", "project/scripts/run_project_benchmark_v02.py", "benchmark",
            "Assemble module-level benchmark tables and figures.",
            "Modules retain distinct tasks; no synthetic end-to-end AUC is created.",
            python_modules=("numpy", "pandas", "rdkit"),
        ),
    ]
}
