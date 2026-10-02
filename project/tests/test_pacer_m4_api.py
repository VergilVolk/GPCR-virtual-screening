from __future__ import annotations

from fastapi.testclient import TestClient

from pacer_m4.api import app


client = TestClient(app)


def test_health_is_read_only_by_default(monkeypatch) -> None:
    monkeypatch.delenv("PACER_M4_API_ALLOW_EXECUTION", raising=False)
    response = client.get("/health")
    assert response.status_code == 200
    assert response.json()["execution_allowed"] is False


def test_capabilities_expose_registered_stages() -> None:
    response = client.get("/v1/capabilities")
    assert response.status_code == 200
    stage_ids = {item["stage_id"] for item in response.json()["stages"]}
    assert {"drugclip-route", "dock-xr-six-channel", "md-production", "dynamic-score"} <= stage_ids


def test_stage_dry_run_returns_exact_command_without_execution() -> None:
    response = client.post(
        "/v1/stages/drugclip-route/run",
        json={"args": ["--input", "scores.csv", "--output", "routed.csv"]},
    )
    assert response.status_code == 200
    payload = response.json()
    assert payload["status"] == "dry_run"
    assert payload["command"][-4:] == ["--input", "scores.csv", "--output", "routed.csv"]


def test_execution_requires_explicit_local_opt_in(monkeypatch) -> None:
    monkeypatch.delenv("PACER_M4_API_ALLOW_EXECUTION", raising=False)
    response = client.post(
        "/v1/stages/validate-release/run",
        json={"args": [], "execute": True},
    )
    assert response.status_code == 403


def test_unknown_stage_returns_404() -> None:
    response = client.get("/v1/stages/not-a-stage")
    assert response.status_code == 404
