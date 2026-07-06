from __future__ import annotations

from typing import Any

from deepretro import score


class FakeChem:
    @staticmethod
    def MolFromSmiles(smiles: str) -> object:
        return object()


class ExplodingSAScorer:
    @staticmethod
    def calculateScore(mol: object) -> float:
        raise ValueError("boom")


def test_score_molecule_surfaces_sa_failure_detail(monkeypatch) -> None:
    monkeypatch.setattr(score, "canonicalize_smiles", lambda smiles: "CCO")
    monkeypatch.setattr(score, "_get_chem", lambda: FakeChem())
    monkeypatch.setattr(score, "_get_sascorer", lambda: ExplodingSAScorer())
    monkeypatch.setattr(score, "_sc_score_with_error", lambda smiles, config: (None, None))

    result = score.score_molecule("CCO")

    assert "SA_Score failed: boom" in result["warnings"]
