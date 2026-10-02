"""Local FastAPI service for the registered PACER-M4 stages.

Execution is disabled by default. Set ``PACER_M4_API_ALLOW_EXECUTION=1`` only
on a trusted local machine; this service is not designed for public exposure.
"""
from __future__ import annotations

import os
import threading
import uuid
from pathlib import Path
from typing import Any

from fastapi import BackgroundTasks, FastAPI, HTTPException
from pydantic import BaseModel, Field

from . import __version__
from .runner import repository_root, run_stage, stage_capabilities
from .stages import STAGES


class StageRunRequest(BaseModel):
    args: list[str] = Field(default_factory=list)
    execute: bool = False


class JobStore:
    def __init__(self) -> None:
        self._lock = threading.Lock()
        self._jobs: dict[str, dict[str, Any]] = {}

    def put(self, job_id: str, value: dict[str, Any]) -> None:
        with self._lock:
            self._jobs[job_id] = value

    def get(self, job_id: str) -> dict[str, Any] | None:
        with self._lock:
            value = self._jobs.get(job_id)
            return dict(value) if value else None


ROOT = repository_root()
JOBS = JobStore()
app = FastAPI(
    title="PACER-M4 local service",
    version=__version__,
    description="Audited access to registered PACER-M4 workflow stages.",
)


def _execution_allowed() -> bool:
    return os.environ.get("PACER_M4_API_ALLOW_EXECUTION", "0") == "1"


def _execute_job(job_id: str, stage_id: str, arguments: list[str]) -> None:
    run_dir = ROOT / "runs" / "pacer_m4_api"
    JOBS.put(job_id, {"job_id": job_id, "stage_id": stage_id, "status": "running"})
    try:
        receipt = run_stage(
            stage_id,
            arguments,
            root=ROOT,
            receipt_path=run_dir / f"{job_id}.json",
            log_path=run_dir / f"{job_id}.log",
        )
        JOBS.put(job_id, {"job_id": job_id, **receipt})
    except Exception as error:  # the job ledger must retain unexpected failures
        JOBS.put(
            job_id,
            {
                "job_id": job_id,
                "stage_id": stage_id,
                "status": "failed",
                "error_type": type(error).__name__,
                "error": str(error),
            },
        )


@app.get("/health")
def health() -> dict[str, Any]:
    return {
        "status": "ok",
        "version": __version__,
        "repository_root": str(ROOT),
        "execution_allowed": _execution_allowed(),
    }


@app.get("/v1/capabilities")
def capabilities() -> dict[str, Any]:
    return stage_capabilities(ROOT)


@app.get("/v1/stages/{stage_id}")
def stage_detail(stage_id: str) -> dict[str, Any]:
    if stage_id not in STAGES:
        raise HTTPException(status_code=404, detail="unknown stage")
    return STAGES[stage_id].to_dict(ROOT)


@app.post("/v1/stages/{stage_id}/run")
def submit_stage(
    stage_id: str,
    request: StageRunRequest,
    background_tasks: BackgroundTasks,
) -> dict[str, Any]:
    if stage_id not in STAGES:
        raise HTTPException(status_code=404, detail="unknown stage")
    if not request.execute:
        return run_stage(stage_id, request.args, root=ROOT, dry_run=True)
    if not _execution_allowed():
        raise HTTPException(
            status_code=403,
            detail="execution is disabled; set PACER_M4_API_ALLOW_EXECUTION=1 on a trusted local host",
        )
    job_id = uuid.uuid4().hex
    JOBS.put(job_id, {"job_id": job_id, "stage_id": stage_id, "status": "queued"})
    background_tasks.add_task(_execute_job, job_id, stage_id, request.args)
    return JOBS.get(job_id) or {}


@app.get("/v1/jobs/{job_id}")
def job_status(job_id: str) -> dict[str, Any]:
    job = JOBS.get(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail="unknown job")
    return job

