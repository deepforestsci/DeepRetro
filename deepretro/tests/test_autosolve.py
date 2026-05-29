"""Tests for the package AutoSolver orchestration layer."""

from __future__ import annotations

from typing import Any

import pytest

import deepretro.algorithms.autosolve as autosolve
from deepretro.algorithms.autosolve import (
    AutoSolver,
    canonicalize,
    reaction_tree,
    unsolved_leaf,
)

TARGET = "CCCO"
ETHANOL = "CCO"
METHANOL = "CO"


def solved_route(smiles: str) -> dict[str, Any]:
    """Return a minimal solved route node for tests."""
    return {
        "type": "mol",
        "smiles": smiles,
        "is_chemical": True,
        "in_stock": True,
    }


def test_public_lazy_exports_import_autosolver() -> None:
    """Package-level and algorithms-level lazy exports expose AutoSolver."""
    import deepretro
    import deepretro.algorithms

    assert deepretro.AutoSolver is AutoSolver
    assert deepretro.algorithms.AutoSolver is AutoSolver


def test_solve_returns_route_tree_and_solved_flag(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """solve() returns the raw route tree and solved flag on AZ success."""
    first_route = solved_route(TARGET)

    monkeypatch.setattr(
        autosolve,
        "run_az",
        lambda smiles, az_model: (True, [first_route]),
    )

    route, solved = AutoSolver(hallucination_mode="none").solve(TARGET)

    assert solved is True
    assert route == first_route


def test_solve_returns_first_az_route_without_calling_llm(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """AiZynthFinder success short-circuits LLM fallback."""
    llm_calls: list[str] = []
    first_route = solved_route(TARGET)
    second_route = solved_route("CCCC")

    monkeypatch.setattr(
        autosolve,
        "run_az",
        lambda smiles, az_model: (True, [first_route, second_route]),
    )
    monkeypatch.setattr(
        autosolve,
        "llm_pipeline",
        lambda **kwargs: llm_calls.append(kwargs["molecule"]),
    )

    route, solved = AutoSolver(hallucination_mode="none").solve(TARGET)

    assert solved is True
    assert route == first_route
    assert llm_calls == []


def test_solve_calls_llm_with_package_pipeline_arguments(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """LLM fallback uses the target branch ``deepretro.utils.llm`` contract."""
    captured: dict[str, Any] = {}

    def fake_run_az(smiles: str, az_model: str) -> tuple[bool, list[dict[str, Any]]]:
        return (True, [solved_route(smiles)]) if smiles == METHANOL else (False, [])

    def fake_llm_pipeline(
        **kwargs: Any,
    ) -> tuple[list[list[str]], list[str], list[float]]:
        captured.update(kwargs)
        return [[METHANOL]], ["reduce chain"], [0.7]

    monkeypatch.setattr(autosolve, "run_az", fake_run_az)
    monkeypatch.setattr(autosolve, "llm_pipeline", fake_llm_pipeline)

    route, solved = AutoSolver(
        llm="openai/gpt-4o-mini",
        az_model="USPTO",
        stability_check=False,
        hallucination_mode="heuristic",
        enable_thinking=False,
        max_output_tokens=64,
    ).solve(TARGET)

    assert solved is True
    assert route["children"][0]["children"] == [solved_route(METHANOL)]
    assert captured["molecule"] == TARGET
    assert captured["model"] == "openai/gpt-4o-mini"
    assert captured["stability_check"] is False
    assert captured["hallucination_check"] is True
    assert captured["enable_thinking"] is False
    assert captured["max_output_tokens"] == 64


def test_solve_selects_first_fully_solved_candidate(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A failed candidate does not pollute the selected successful route."""

    def fake_run_az(smiles: str, az_model: str) -> tuple[bool, list[dict[str, Any]]]:
        return (True, [solved_route(smiles)]) if smiles == METHANOL else (False, [])

    def fake_llm_pipeline(
        **kwargs: Any,
    ) -> tuple[list[list[str]], list[str], list[float]]:
        if kwargs["molecule"] == TARGET:
            return [[ETHANOL], [METHANOL]], ["bad", "good"], [0.2, 0.9]
        return [], [], []

    monkeypatch.setattr(autosolve, "run_az", fake_run_az)
    monkeypatch.setattr(autosolve, "llm_pipeline", fake_llm_pipeline)

    route, solved = AutoSolver(hallucination_mode="none").solve(TARGET)

    assert solved is True
    assert route["children"][0]["children"] == [solved_route(METHANOL)]
    assert unsolved_leaf(ETHANOL) not in route["children"][0]["children"]


def test_solve_returns_first_partial_candidate_when_none_fully_solve(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """If all candidates fail, the first attempted route remains inspectable."""

    def fake_llm_pipeline(
        **kwargs: Any,
    ) -> tuple[list[list[str]], list[str], list[float]]:
        if kwargs["molecule"] != TARGET:
            return [], [], []
        return [[ETHANOL], [METHANOL]], ["first", "second"], [0.2, 0.9]

    monkeypatch.setattr(autosolve, "run_az", lambda smiles, az_model: (False, []))
    monkeypatch.setattr(autosolve, "llm_pipeline", fake_llm_pipeline)

    route, solved = AutoSolver(hallucination_mode="none").solve(TARGET)

    assert solved is False
    assert route["children"][0]["children"] == [unsolved_leaf(ETHANOL)]


def test_solve_solves_multi_reactant_pathway(monkeypatch: pytest.MonkeyPatch) -> None:
    """Every reactant in a candidate pathway must solve for the route to solve."""

    def fake_run_az(smiles: str, az_model: str) -> tuple[bool, list[dict[str, Any]]]:
        return (
            (True, [solved_route(smiles)])
            if smiles in {ETHANOL, METHANOL}
            else (False, [])
        )

    monkeypatch.setattr(autosolve, "run_az", fake_run_az)
    monkeypatch.setattr(
        autosolve,
        "llm_pipeline",
        lambda **kwargs: ([[ETHANOL, METHANOL]], ["two reactants"], [0.8]),
    )

    route, solved = AutoSolver(hallucination_mode="none").solve(TARGET)

    assert solved is True
    assert route["children"][0]["children"] == [
        solved_route(ETHANOL),
        solved_route(METHANOL),
    ]


def test_solve_returns_unsolved_tree_when_llm_has_no_pathways(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """No LLM candidates yields a stable unsolved leaf."""
    monkeypatch.setattr(autosolve, "run_az", lambda smiles, az_model: (False, []))
    monkeypatch.setattr(autosolve, "llm_pipeline", lambda **kwargs: ([], [], []))

    route, solved = AutoSolver(hallucination_mode="none").solve(TARGET)

    assert solved is False
    assert route == unsolved_leaf(TARGET)


def test_solve_detects_canonical_cycles(monkeypatch: pytest.MonkeyPatch) -> None:
    """Canonicalized revisits stop recursion."""
    monkeypatch.setattr(autosolve, "run_az", lambda smiles, az_model: (False, []))
    monkeypatch.setattr(
        autosolve,
        "llm_pipeline",
        lambda **kwargs: ([["C(O)C"]], ["same molecule"], [0.1]),
    )

    route, solved = AutoSolver(hallucination_mode="none").solve("CCO")

    assert solved is False
    assert route["children"][0]["children"] == [unsolved_leaf("C(O)C")]


def test_solve_stops_at_max_depth(monkeypatch: pytest.MonkeyPatch) -> None:
    """The recursion limit returns an unsolved leaf before external calls."""
    run_az_calls: list[str] = []

    def fake_run_az(smiles: str, az_model: str) -> tuple[bool, list[dict[str, Any]]]:
        run_az_calls.append(smiles)
        return False, []

    monkeypatch.setattr(autosolve, "run_az", fake_run_az)

    route, solved = AutoSolver(hallucination_mode="none", max_depth=0).solve(TARGET)

    assert (route, solved) == (unsolved_leaf(TARGET), False)
    assert run_az_calls == []


def test_constructor_preserves_public_options() -> None:
    """Constructor options are stored without invoking external services."""
    solver = AutoSolver(
        llm="openai/gpt-4o-mini:adv",
        az_model="USPTO",
        stability_check=False,
        hallucination_mode="none",
        max_depth=3,
        enable_thinking=False,
        max_output_tokens=128,
        metadata_model="gpt-4o",
    )

    assert solver.llm == "openai/gpt-4o-mini:adv"
    assert solver.az_model == "USPTO"
    assert solver.stability_check is False
    assert solver.hallucination_mode == "none"
    assert solver.hallucination_checker is None
    assert solver.max_depth == 3
    assert solver.enable_thinking is False
    assert solver.max_output_tokens == 128
    assert solver.metadata_model == "gpt-4o"


def test_constructor_defaults() -> None:
    """Default constructor values match the documented API."""
    solver = AutoSolver(hallucination_mode="none")

    assert solver.llm == "anthropic/claude-opus-4-6:adv"
    assert solver.az_model == "Pistachio_100+"
    assert solver.stability_check is True
    assert solver.hallucination_mode == "none"
    assert solver.max_depth == 50
    assert solver.enable_thinking is True
    assert solver.max_output_tokens is None
    assert solver.metadata_model == "claude-opus-4-20250514"


def test_canonicalize_and_unsolved_leaf_helpers() -> None:
    """Small helpers expose stable deterministic behavior."""
    assert canonicalize("C(O)C") == "CCO"
    assert canonicalize("not_a_smiles") == "not_a_smiles"
    assert unsolved_leaf("CCO") == {
        "type": "mol",
        "smiles": "CCO",
        "is_chemical": True,
        "in_stock": False,
        "children": [],
    }


# ---------------------------------------------------------------------------
# single_step tests
# ---------------------------------------------------------------------------


def test_single_step_returns_az_route_when_solved(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """single_step() returns AZ result without LLM or recursion."""
    az_route = solved_route(TARGET)
    monkeypatch.setattr(
        autosolve,
        "run_az",
        lambda smiles, az_model: (True, [az_route]),
    )

    route, solved = AutoSolver(hallucination_mode="none").single_step(TARGET)

    assert solved is True
    assert route == az_route


def test_single_step_returns_llm_route_without_recursion(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """single_step() builds a one-level tree from LLM output without recursing."""
    monkeypatch.setattr(autosolve, "run_az", lambda smiles, az_model: (False, []))
    monkeypatch.setattr(
        autosolve,
        "llm_pipeline",
        lambda **kwargs: ([[ETHANOL, METHANOL]], ["split"], [0.9]),
    )

    route, solved = AutoSolver(hallucination_mode="none").single_step(TARGET)

    assert solved is False
    children = route["children"][0]["children"]
    assert len(children) == 2
    assert children[0] == unsolved_leaf(ETHANOL)
    assert children[1] == unsolved_leaf(METHANOL)


def test_single_step_returns_unsolved_leaf_when_no_pathways(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """No AZ or LLM result produces an unsolved leaf."""
    monkeypatch.setattr(autosolve, "run_az", lambda smiles, az_model: (False, []))
    monkeypatch.setattr(autosolve, "llm_pipeline", lambda **kwargs: ([], [], []))

    route, solved = AutoSolver(hallucination_mode="none").single_step(TARGET)

    assert solved is False
    assert route == unsolved_leaf(TARGET)


# ---------------------------------------------------------------------------
# parse tests
# ---------------------------------------------------------------------------


def test_parse_formats_route_tree_into_steps(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """parse() delegates to format_output and adds the solved flag."""
    fake_output: dict[str, Any] = {"steps": [{"step": "1"}], "dependencies": {"1": []}}
    monkeypatch.setattr(autosolve, "format_output", lambda tree: dict(fake_output))

    route_tree = reaction_tree(TARGET, [solved_route(ETHANOL)], [0.8])
    result = AutoSolver(hallucination_mode="none").parse(route_tree, solved=True)

    assert result == {**fake_output, "solved": True}


def test_parse_propagates_solved_false(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """parse() correctly propagates unsolved status."""
    monkeypatch.setattr(
        autosolve,
        "format_output",
        lambda tree: {"steps": [], "dependencies": {}},
    )

    result = AutoSolver(hallucination_mode="none").parse(
        unsolved_leaf(TARGET), solved=False
    )

    assert result["solved"] is False


# ---------------------------------------------------------------------------
# add_metadata tests
# ---------------------------------------------------------------------------


def test_add_metadata_enriches_steps_with_recommendation(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """add_metadata() calls recommend_reaction_metadata per step and merges."""
    captured_calls: list[str] = []

    def fake_recommend(
        reaction_smiles: str, **kwargs: Any
    ) -> tuple[int, dict[str, Any]]:
        captured_calls.append(reaction_smiles)
        return 200, {
            "reaction_smiles": reaction_smiles,
            "reactants": [],
            "product": [],
            "reagents": [{"smiles": "O", "reagent_metadata": {"name": "water"}}],
            "conditions": {
                "temperature": "25 C",
                "pressure": "1 atm",
                "solvent": "water",
                "time": "1 h",
            },
            "literature": {"doi": "10.1000/example"},
        }

    monkeypatch.setattr(autosolve, "recommend_reaction_metadata", fake_recommend)

    parsed: dict[str, Any] = {
        "steps": [
            {
                "step": "1",
                "reactants": [{"smiles": "CC", "reactant_metadata": {}}],
                "reagents": [],
                "products": [{"smiles": "CCO", "product_metadata": {}}],
                "conditions": [],
                "reactionmetrics": [{"scalabilityindex": 0, "closestliterature": ""}],
            },
        ],
        "dependencies": {"1": []},
        "solved": True,
    }

    result = AutoSolver(hallucination_mode="none").add_metadata(parsed)

    assert captured_calls == ["CC>>CCO"]
    step = result["steps"][0]
    assert step["reagents"] == [{"smiles": "O", "reagent_metadata": {"name": "water"}}]
    assert step["conditions"] == {
        "temperature": "25 C",
        "pressure": "1 atm",
        "solvent": "water",
        "time": "1 h",
    }
    assert step["reactionmetrics"][0]["closestliterature"] == {"doi": "10.1000/example"}


def test_add_metadata_skips_step_on_recommendation_failure(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Metadata failure for one step does not block the pipeline."""
    monkeypatch.setattr(
        autosolve,
        "recommend_reaction_metadata",
        lambda reaction_smiles, **kwargs: (
            404,
            {"stage": "reagents", "error": "fail"},
        ),
    )

    parsed: dict[str, Any] = {
        "steps": [
            {
                "step": "1",
                "reactants": [{"smiles": "CC", "reactant_metadata": {}}],
                "reagents": [],
                "products": [{"smiles": "CCO", "product_metadata": {}}],
                "conditions": [],
                "reactionmetrics": [{"scalabilityindex": 0, "closestliterature": ""}],
            },
        ],
        "dependencies": {"1": []},
        "solved": True,
    }

    result = AutoSolver(hallucination_mode="none").add_metadata(parsed)

    step = result["steps"][0]
    assert step["reagents"] == []
    assert step["conditions"] == []
    assert step["reactionmetrics"][0]["closestliterature"] == ""


def test_add_metadata_uses_metadata_model_param(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """add_metadata passes the configured metadata_model."""
    captured_kwargs: dict[str, Any] = {}

    def fake_recommend(
        reaction_smiles: str, **kwargs: Any
    ) -> tuple[int, dict[str, Any]]:
        captured_kwargs.update(kwargs)
        return 404, {"stage": "reagents", "error": "fail"}

    monkeypatch.setattr(autosolve, "recommend_reaction_metadata", fake_recommend)

    parsed: dict[str, Any] = {
        "steps": [
            {
                "step": "1",
                "reactants": [{"smiles": "CC", "reactant_metadata": {}}],
                "reagents": [],
                "products": [{"smiles": "CCO", "product_metadata": {}}],
                "conditions": [],
                "reactionmetrics": [{"scalabilityindex": 0, "closestliterature": ""}],
            },
        ],
        "dependencies": {"1": []},
        "solved": True,
    }

    AutoSolver(hallucination_mode="none", metadata_model="gpt-4o").add_metadata(parsed)

    assert captured_kwargs["model"] == "gpt-4o"


def test_add_metadata_handles_multi_reactant_smiles(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Multiple reactants are joined with '.' in the reaction SMILES."""
    captured_calls: list[str] = []

    def fake_recommend(
        reaction_smiles: str, **kwargs: Any
    ) -> tuple[int, dict[str, Any]]:
        captured_calls.append(reaction_smiles)
        return 404, {"stage": "reagents", "error": "fail"}

    monkeypatch.setattr(autosolve, "recommend_reaction_metadata", fake_recommend)

    parsed: dict[str, Any] = {
        "steps": [
            {
                "step": "1",
                "reactants": [
                    {"smiles": "CC", "reactant_metadata": {}},
                    {"smiles": "O", "reactant_metadata": {}},
                ],
                "reagents": [],
                "products": [{"smiles": "CCO", "product_metadata": {}}],
                "conditions": [],
                "reactionmetrics": [{"scalabilityindex": 0, "closestliterature": ""}],
            },
        ],
        "dependencies": {"1": []},
        "solved": True,
    }

    AutoSolver(hallucination_mode="none").add_metadata(parsed)

    assert captured_calls == ["CC.O>>CCO"]


# ---------------------------------------------------------------------------
# autosolve tests
# ---------------------------------------------------------------------------


def test_autosolve_chains_solve_parse_add_metadata(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """autosolve() runs the full pipeline and returns enriched output."""
    call_order: list[str] = []

    def tracking_solve(
        self: AutoSolver, smiles: str, **kw: Any
    ) -> tuple[dict[str, Any], bool]:
        call_order.append("solve")
        return solved_route(smiles), True

    def tracking_parse(
        self: AutoSolver, tree: dict[str, Any], *, solved: bool
    ) -> dict[str, Any]:
        call_order.append("parse")
        return {"steps": [], "dependencies": {}, "solved": solved}

    def tracking_add_metadata(
        self: AutoSolver, parsed: dict[str, Any]
    ) -> dict[str, Any]:
        call_order.append("add_metadata")
        return parsed

    monkeypatch.setattr(AutoSolver, "solve", tracking_solve)
    monkeypatch.setattr(AutoSolver, "parse", tracking_parse)
    monkeypatch.setattr(AutoSolver, "add_metadata", tracking_add_metadata)

    result = AutoSolver(hallucination_mode="none").autosolve(TARGET)

    assert call_order == ["solve", "parse", "add_metadata"]
    assert result["solved"] is True


def test_autosolve_propagates_unsolved_status(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """autosolve() passes solved=False through the pipeline."""
    monkeypatch.setattr(
        AutoSolver,
        "solve",
        lambda self, smiles, **kw: (unsolved_leaf(smiles), False),
    )
    monkeypatch.setattr(
        autosolve,
        "format_output",
        lambda tree: {"steps": [], "dependencies": {}},
    )
    monkeypatch.setattr(
        AutoSolver,
        "add_metadata",
        lambda self, parsed: parsed,
    )

    result = AutoSolver(hallucination_mode="none").autosolve(TARGET)

    assert result["solved"] is False
