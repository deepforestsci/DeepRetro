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


def test_solve_formats_recurse_tree_and_adds_solved_flag(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """solve() delegates to recurse() and formats the returned route tree."""
    fake_tree = reaction_tree(TARGET, [solved_route(ETHANOL)], [0.8])
    fake_output = {"steps": [{"step": "1"}], "dependencies": {"1": []}}
    calls: dict[str, Any] = {}

    def fake_recurse(self: AutoSolver, smiles: str) -> tuple[dict[str, Any], bool]:
        calls["recurse"] = smiles
        return fake_tree, True

    def fake_format(tree: dict[str, Any]) -> dict[str, Any]:
        calls["format"] = tree
        return dict(fake_output)

    monkeypatch.setattr(AutoSolver, "recurse", fake_recurse)
    monkeypatch.setattr(autosolve, "format_output", fake_format)

    result = AutoSolver(hallucination_mode="none").solve(TARGET)

    assert result == {**fake_output, "solved": True}
    assert calls == {"recurse": TARGET, "format": fake_tree}


def test_solve_can_append_rendered_image(monkeypatch: pytest.MonkeyPatch) -> None:
    """return_image=True calls visualization with the formatted output."""
    fake_tree = reaction_tree(TARGET, [solved_route(ETHANOL)], [0.8])
    image = object()
    visualized: dict[str, Any] = {}

    monkeypatch.setattr(
        AutoSolver,
        "recurse",
        lambda self, smiles: (fake_tree, True),
    )
    monkeypatch.setattr(
        autosolve,
        "format_output",
        lambda tree: {"steps": [], "dependencies": {}},
    )

    import deepretro.utils.visualize as visualize

    def fake_visualize(output: dict[str, Any]) -> object:
        visualized["output"] = dict(output)
        return image

    monkeypatch.setattr(visualize, "visualize_pathway", fake_visualize)

    result = AutoSolver(hallucination_mode="none").solve(TARGET, return_image=True)

    assert result["image"] is image
    assert visualized["output"] == {"steps": [], "dependencies": {}, "solved": True}


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
    )

    assert solver.llm == "openai/gpt-4o-mini:adv"
    assert solver.az_model == "USPTO"
    assert solver.stability_check is False
    assert solver.hallucination_mode == "none"
    assert solver.hallucination_checker is None
    assert solver.max_depth == 3
    assert solver.enable_thinking is False
    assert solver.max_output_tokens == 128


def test_recurse_returns_first_az_route_without_calling_llm(
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

    route, solved = AutoSolver(hallucination_mode="none").recurse(TARGET)

    assert solved is True
    assert route == first_route
    assert llm_calls == []


def test_recurse_calls_llm_with_package_pipeline_arguments(
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
    ).recurse(TARGET)

    assert solved is True
    assert route["children"][0]["children"] == [solved_route(METHANOL)]
    assert captured["molecule"] == TARGET
    assert captured["model"] == "openai/gpt-4o-mini"
    assert captured["stability_check"] is False
    assert captured["hallucination_check"] is True
    assert captured["enable_thinking"] is False
    assert captured["max_output_tokens"] == 64


def test_recurse_selects_first_fully_solved_candidate(
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

    route, solved = AutoSolver(hallucination_mode="none").recurse(TARGET)

    assert solved is True
    assert route["children"][0]["children"] == [solved_route(METHANOL)]
    assert unsolved_leaf(ETHANOL) not in route["children"][0]["children"]


def test_recurse_returns_first_partial_candidate_when_none_fully_solve(
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

    route, solved = AutoSolver(hallucination_mode="none").recurse(TARGET)

    assert solved is False
    assert route["children"][0]["children"] == [unsolved_leaf(ETHANOL)]


def test_recurse_solves_multi_reactant_pathway(monkeypatch: pytest.MonkeyPatch) -> None:
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

    route, solved = AutoSolver(hallucination_mode="none").recurse(TARGET)

    assert solved is True
    assert route["children"][0]["children"] == [
        solved_route(ETHANOL),
        solved_route(METHANOL),
    ]


def test_recurse_returns_unsolved_tree_when_llm_has_no_pathways(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """No LLM candidates yields a stable unsolved leaf."""
    monkeypatch.setattr(autosolve, "run_az", lambda smiles, az_model: (False, []))
    monkeypatch.setattr(autosolve, "llm_pipeline", lambda **kwargs: ([], [], []))

    route, solved = AutoSolver(hallucination_mode="none").recurse(TARGET)

    assert solved is False
    assert route == unsolved_leaf(TARGET)


def test_recurse_detects_canonical_cycles(monkeypatch: pytest.MonkeyPatch) -> None:
    """Canonicalized revisits stop recursion."""
    monkeypatch.setattr(autosolve, "run_az", lambda smiles, az_model: (False, []))
    monkeypatch.setattr(
        autosolve,
        "llm_pipeline",
        lambda **kwargs: ([["C(O)C"]], ["same molecule"], [0.1]),
    )

    route, solved = AutoSolver(hallucination_mode="none").recurse("CCO")

    assert solved is False
    assert route["children"][0]["children"] == [unsolved_leaf("C(O)C")]


def test_recurse_stops_at_max_depth(monkeypatch: pytest.MonkeyPatch) -> None:
    """The recursion limit returns an unsolved leaf before external calls."""
    run_az_calls: list[str] = []

    def fake_run_az(smiles: str, az_model: str) -> tuple[bool, list[dict[str, Any]]]:
        run_az_calls.append(smiles)
        return False, []

    monkeypatch.setattr(autosolve, "run_az", fake_run_az)

    route, solved = AutoSolver(hallucination_mode="none", max_depth=0).recurse(TARGET)

    assert (route, solved) == (unsolved_leaf(TARGET), False)
    assert run_az_calls == []


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
