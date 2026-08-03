"""Unit tests for :mod:`deepretro.score`.

Covers the public scoring API: SMILES canonicalisation, single-molecule and
step scoring, whole-pathway scoring, and the empty-payload helper. Optional
metrics (SA_Score / SCScore) degrade to ``None`` when their backends are
unavailable, so the tests assert the None-or-float contract and only make a
numeric assertion when the backend is present.
"""

import pytest
from deepretro import score

# ---------------------------------------------------------------------------
# canonicalize_smiles
# ---------------------------------------------------------------------------


def test_canonicalize_valid_is_idempotent():
    """Canonicalisation returns a stable canonical form."""
    once = score.canonicalize_smiles("C1=CC=CC=C1")  # benzene, Kekulé form
    assert once is not None
    assert score.canonicalize_smiles(once) == once


@pytest.mark.parametrize("bad", ["", "   ", "not_a_smiles!!!", None])
def test_canonicalize_invalid_returns_none(bad):
    """Empty, whitespace, invalid, or None inputs canonicalise to None."""
    assert score.canonicalize_smiles(bad) is None


# ---------------------------------------------------------------------------
# sa_score / sc_score contracts
# ---------------------------------------------------------------------------


def test_sa_score_none_or_float_for_valid():
    """SA score is a float when the backend is present, else None (never raises)."""
    value = score.sa_score("c1ccccc1")
    if value is None:
        pytest.skip("SA_Score backend unavailable in this environment")
    assert isinstance(value, float)
    assert value > 0


def test_sa_score_invalid_smiles_is_none():
    """An invalid SMILES yields no SA score."""
    assert score.sa_score("not_a_smiles!!!") is None


def test_sc_score_returns_none_or_float():
    """SC score honours the None-or-float contract without raising."""
    value = score.sc_score("c1ccccc1")
    assert value is None or isinstance(value, float)


# ---------------------------------------------------------------------------
# score_molecule
# ---------------------------------------------------------------------------

_MOLECULE_KEYS = {
    "input_smiles",
    "canonical_smiles",
    "valid",
    "sa_score",
    "sc_score",
    "warnings",
    "errors",
}


def test_score_molecule_valid():
    """A valid molecule is marked valid with a canonical form and full schema."""
    result = score.score_molecule("c1ccccc1")
    assert set(result) == _MOLECULE_KEYS
    assert result["valid"] is True
    assert result["canonical_smiles"] is not None
    assert result["errors"] == []


def test_score_molecule_invalid():
    """An invalid molecule is flagged with an error and no canonical form."""
    result = score.score_molecule("not_a_smiles!!!")
    assert result["valid"] is False
    assert result["canonical_smiles"] is None
    assert "Invalid SMILES" in result["errors"]


# ---------------------------------------------------------------------------
# score_step
# ---------------------------------------------------------------------------


def test_score_step_valid_shape():
    """A valid step scores its product and reactants and reports delta metrics."""
    result = score.score_step("CCO", "CC=O")
    assert result["errors"] == []
    assert result["product"]["valid"] is True
    assert len(result["reactants"]) == 1
    assert isinstance(result["sa_simplifies"], bool)
    for key in ("delta_sa_mean", "delta_sa_max", "delta_sc_mean", "delta_sc_max"):
        assert key in result


def test_score_step_accepts_list_reactants():
    """Reactants may be provided as a list rather than a dot-joined string."""
    result = score.score_step("CCO", ["CC=O", "O"])
    assert len(result["reactants"]) == 2


def test_score_step_missing_reactants_is_malformed():
    """A step with no reactant SMILES is reported as malformed."""
    result = score.score_step("CCO", "")
    assert result["errors"] == ["Malformed step: missing reactant SMILES"]


# ---------------------------------------------------------------------------
# score_pathway
# ---------------------------------------------------------------------------

_PATHWAY_KEYS = {"step_scores", "route_summary", "warnings", "errors"}


def test_score_pathway_viewer_style():
    """A viewer-style ``{'steps': [...]}`` pathway scores each step and summarises."""
    pathway = {"steps": [{"product": "CCO", "reactants": "CC=O"}]}
    result = score.score_pathway(pathway)
    assert set(result) == _PATHWAY_KEYS
    assert len(result["step_scores"]) == 1
    assert result["route_summary"]["n_steps"] == 1


def test_score_pathway_list_of_steps():
    """A bare list of step dicts is accepted as a pathway."""
    result = score.score_pathway([{"product": "CCO", "reactants": "CC=O"}])
    assert len(result["step_scores"]) == 1


def test_score_pathway_malformed_input():
    """A non-dict/non-list pathway is reported as malformed with no steps."""
    result = score.score_pathway("not a pathway")
    assert result["step_scores"] == []
    assert result["errors"] == ["Malformed pathway: expected dict or list"]


# ---------------------------------------------------------------------------
# empty_pathway_scores
# ---------------------------------------------------------------------------


def test_empty_pathway_scores_schema():
    """The empty payload uses the standard pathway schema with no steps."""
    result = score.empty_pathway_scores()
    assert set(result) == _PATHWAY_KEYS
    assert result["step_scores"] == []
    assert result["errors"] == []


def test_empty_pathway_scores_carries_error():
    """A supplied error message is surfaced in the empty payload."""
    result = score.empty_pathway_scores("boom")
    assert result["errors"] == ["boom"]
