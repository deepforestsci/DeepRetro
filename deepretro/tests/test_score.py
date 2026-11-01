"""Offline score contracts, including optional model failures."""

from __future__ import annotations

import json
from types import SimpleNamespace
from typing import Any

import pytest
from rdkit import rdBase

from deepretro import score


def test_canonicalize_smiles_preserves_identity_and_rejects_invalid() -> None:
    assert score.canonicalize_smiles(" C(C)O ") == "CCO"
    assert score.canonicalize_smiles("not-a-molecule") is None
    assert score.canonicalize_smiles("") is None


def test_canonicalization_preserves_rdkit_log_settings() -> None:
    before = rdBase.LogStatus()
    try:
        rdBase.DisableLog("rdApp.warning")
        expected = rdBase.LogStatus()
        assert score.canonicalize_smiles("invalid") is None
        assert rdBase.LogStatus() == expected
    finally:
        for line in before.splitlines():
            channel, state = line.split(":")
            (rdBase.EnableLog if state == "enabled" else rdBase.DisableLog)(channel)


@pytest.mark.parametrize("value", [[], (), float("nan"), float("inf"), "bad", None])
def test_invalid_optional_model_output_becomes_warning(
    value: Any, monkeypatch: pytest.MonkeyPatch
) -> None:
    scorer = SimpleNamespace(get_score_from_smi=lambda smiles: value)
    monkeypatch.setattr(score, "_get_scscorer", lambda *args: (scorer, None))
    result = score.score_molecule("CCO")
    assert result["valid"] is True
    assert result["sc_score"] is None
    assert any("SCScore scoring failed" in warning for warning in result["warnings"])
    json.dumps(result, allow_nan=False)


@pytest.mark.parametrize("value", [2.5, ("CCO", 2.5), [2.5]])
def test_sc_score_accepts_supported_model_outputs(
    value: Any, monkeypatch: pytest.MonkeyPatch
) -> None:
    scorer = SimpleNamespace(get_score_from_smi=lambda smiles: value)
    monkeypatch.setattr(score, "_get_scscorer", lambda *args: (scorer, None))
    assert score.sc_score("CCO", score.SCScoreConfig(fp_len=1024)) == 2.5


def test_sa_score_uses_rdkit_contrib() -> None:
    result = score.sa_score("CCO")
    assert result is not None
    assert 1 <= result <= 10
    assert score.sa_score("invalid") is None


def test_step_and_pathway_summaries_use_available_scores(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    values = {"CCC": (4.0, 3.0), "CC": (2.0, 1.0), "O": (None, None)}
    monkeypatch.setattr(
        score, "_sa_score_with_error", lambda smi: (values[smi][0], None)
    )
    monkeypatch.setattr(
        score, "_sc_score_with_error", lambda smi, config: (values[smi][1], None)
    )
    step = score.score_step("CCC", "CC.O")
    assert step["mean_reactant_sa"] == 2.0
    assert step["delta_sa_mean"] == 2.0
    assert step["delta_sc_mean"] == 2.0
    assert step["sa_simplifies"] is True
    assert step["warnings"]  # One missing reactant metric must remain visible.

    result = score.score_pathway(
        [
            {"product_smiles": "CCC", "reactants_smiles": ["CC"]},
            {"products": [{"smiles": "CC"}], "reactants": [{"smiles": "CCC"}]},
            {"product": "CCC"},
        ]
    )
    summary = result["route_summary"]
    assert summary["n_steps"] == 3
    assert summary["mean_delta_sa"] == 0.0
    assert summary["fraction_steps_sa_simplifying"] == 0.5
    assert summary["complexity_increasing_steps"] == [1]
    assert result["errors"] == ["Malformed step 2: missing product or reactant SMILES"]


def test_empty_score_reports_preserve_schema() -> None:
    result = score.empty_pathway_scores("unavailable")
    assert result["errors"] == ["unavailable"]
    assert result["step_scores"] == []
    assert result["route_summary"]["n_steps"] == 0
    assert result["route_summary"]["mean_product_sa"] is None
    assert score.score_pathway({})["errors"] == ["Malformed pathway: missing steps"]
