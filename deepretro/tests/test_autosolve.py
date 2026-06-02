"""Tests for the AutoSolver retrosynthesis pipeline."""

from __future__ import annotations

from typing import Any

import pytest

from deepretro.algorithms.autosolve import (
    AutoSolver,
    reaction_tree,
    unsolved_leaf,
)
from deepretro.tests.conftest import (
    ACETIC_ACID,
    ACETIC_ANHYDRIDE,
    ASPIRIN,
    CAFFEINE,
    IBUPROFEN,
    PARA_AMINOPHENOL,
    PARACETAMOL,
    SALICYLIC_ACID,
    aspirin_hydrolysis_tree,
    az_always_fails,
    llm_returns_nothing,
    solved_route,
)
from deepretro.utils.utils_molecule import canonicalize


# ---------------------------------------------------------------------------
# Helper function tests — pure, no DI needed
# ---------------------------------------------------------------------------


class TestCanonicalizeRealMolecules:
    def test_aspirin_non_canonical_form_normalizes(self) -> None:
        non_canonical = "OC(=O)c1ccccc1OC(=O)C"
        result = canonicalize(non_canonical)
        assert result == canonicalize(ASPIRIN)
        assert result != non_canonical

    def test_ibuprofen_idempotent(self) -> None:
        assert canonicalize(IBUPROFEN) == canonicalize(canonicalize(IBUPROFEN))

    def test_caffeine_aromatic_normalization(self) -> None:
        result = canonicalize(CAFFEINE)
        assert result == canonicalize(result)

    def test_invalid_smiles_passthrough(self) -> None:
        assert canonicalize("not_a_smiles") == "not_a_smiles"


class TestUnsolvedLeaf:
    def test_aspirin_leaf_structure(self) -> None:
        node = unsolved_leaf(ASPIRIN)
        assert node["smiles"] == ASPIRIN
        assert node["type"] == "mol"
        assert node["is_chemical"] is True
        assert node["in_stock"] is False
        assert node["children"] == []

    def test_caffeine_leaf_structure(self) -> None:
        node = unsolved_leaf(CAFFEINE)
        assert node["smiles"] == CAFFEINE
        assert node["in_stock"] is False


class TestReactionTree:
    def test_aspirin_hydrolysis_structure(self) -> None:
        children = [unsolved_leaf(SALICYLIC_ACID), unsolved_leaf(ACETIC_ANHYDRIDE)]
        tree = reaction_tree(ASPIRIN, children, [0.85])
        assert tree["smiles"] == ASPIRIN
        assert tree["type"] == "mol"
        rxn = tree["children"][0]
        assert rxn["type"] == "reaction"
        assert rxn["is_reaction"] is True
        assert rxn["metadata"]["policy_probability"] == pytest.approx(0.85)
        assert rxn["children"] == children

    def test_empty_confidence_defaults_to_zero(self) -> None:
        tree = reaction_tree(ASPIRIN, [], [])
        assert tree["children"][0]["metadata"]["policy_probability"] == 0.0


# ---------------------------------------------------------------------------
# Constructor tests
# ---------------------------------------------------------------------------


class TestConstructor:
    def test_preserves_all_options(self) -> None:
        solver = AutoSolver(
            llm="anthropic/claude-sonnet-4-6",
            az_model="USPTO",
            stability_check=False,
            hallucination_mode="none",
            max_depth=3,
            enable_thinking=False,
            max_output_tokens=128,
            metadata_model="anthropic/claude-sonnet-4-6",
        )
        assert solver.llm == "anthropic/claude-sonnet-4-6"
        assert solver.az_model == "USPTO"
        assert solver.stability_check is False
        assert solver.hallucination_mode == "none"
        assert solver.hallucination_checker is None
        assert solver.max_depth == 3
        assert solver.enable_thinking is False
        assert solver.max_output_tokens == 128
        assert solver.metadata_model == "anthropic/claude-sonnet-4-6"

    def test_defaults(self) -> None:
        solver = AutoSolver(hallucination_mode="none")
        assert solver.llm == "anthropic/claude-sonnet-4-6"
        assert solver.az_model == "Pistachio_100+"
        assert solver.stability_check is True
        assert solver.max_depth == 50
        assert solver.enable_thinking is True
        assert solver.max_output_tokens is None
        assert solver.metadata_model == "anthropic/claude-sonnet-4-6"

    def test_normalizes_hallucination_mode_to_lowercase(self) -> None:
        solver = AutoSolver(hallucination_mode="NONE")
        assert solver.hallucination_mode == "none"

    def test_az_runner_injection(self) -> None:
        solver = AutoSolver(az_runner=az_always_fails, hallucination_mode="none")
        assert solver._az_runner is az_always_fails

    def test_lazy_exports(self) -> None:
        import deepretro
        import deepretro.algorithms

        assert deepretro.AutoSolver is AutoSolver
        assert deepretro.algorithms.AutoSolver is AutoSolver


# ---------------------------------------------------------------------------
# solve() tests — DI via constructor
# ---------------------------------------------------------------------------


class TestSolve:
    def test_az_success_returns_route_for_aspirin(self) -> None:
        def az_solves_aspirin(smiles: str, az_model: str) -> tuple[bool, list]:
            return True, [solved_route(smiles)]

        solver = AutoSolver(az_runner=az_solves_aspirin, hallucination_mode="none")
        route, solved = solver.solve(ASPIRIN)

        assert solved is True
        assert route["smiles"] == ASPIRIN
        assert route["in_stock"] is True

    def test_az_success_skips_llm(self) -> None:
        llm_calls: list[str] = []

        def az_solves(smiles: str, az_model: str) -> tuple[bool, list]:
            return True, [solved_route(smiles)]

        def llm_tracks(molecule: str, **kwargs: Any) -> tuple[list, list, list]:
            llm_calls.append(molecule)
            return [], [], []

        solver = AutoSolver(
            az_runner=az_solves,
            llm_runner=llm_tracks,
            hallucination_mode="none",
        )
        solver.solve(ASPIRIN)

        assert llm_calls == []

    def test_llm_fallback_with_solvable_precursors(self) -> None:
        """Aspirin -> salicylic acid + acetic acid, both solved by AZ."""

        def az_solves_precursors(smiles: str, az_model: str) -> tuple[bool, list]:
            canonical = canonicalize(smiles)
            if canonical in {
                canonicalize(SALICYLIC_ACID),
                canonicalize(ACETIC_ACID),
            }:
                return True, [solved_route(smiles)]
            return False, []

        def llm_proposes_hydrolysis(
            molecule: str, **kwargs: Any
        ) -> tuple[list, list, list]:
            return (
                [[SALICYLIC_ACID, ACETIC_ACID]],
                ["aspirin hydrolysis"],
                [0.9],
            )

        solver = AutoSolver(
            az_runner=az_solves_precursors,
            llm_runner=llm_proposes_hydrolysis,
            hallucination_mode="none",
        )
        route, solved = solver.solve(ASPIRIN)

        assert solved is True
        children = route["children"][0]["children"]
        assert len(children) == 2

    def test_cycle_detection_with_non_canonical_aspirin(self) -> None:
        non_canonical = "OC(=O)c1ccccc1OC(=O)C"

        def llm_returns_self(molecule: str, **kwargs: Any) -> tuple[list, list, list]:
            return [[non_canonical]], ["self-reference"], [0.5]

        solver = AutoSolver(
            az_runner=az_always_fails,
            llm_runner=llm_returns_self,
            hallucination_mode="none",
        )
        route, solved = solver.solve(ASPIRIN)

        assert solved is False
        child = route["children"][0]["children"][0]
        assert canonicalize(child["smiles"]) == canonicalize(ASPIRIN)

    def test_max_depth_zero_returns_unsolved_without_az_call(self) -> None:
        az_calls: list[str] = []

        def az_tracks(smiles: str, az_model: str) -> tuple[bool, list]:
            az_calls.append(smiles)
            return False, []

        solver = AutoSolver(az_runner=az_tracks, hallucination_mode="none", max_depth=0)
        route, solved = solver.solve(ASPIRIN)

        assert (route, solved) == (unsolved_leaf(ASPIRIN), False)
        assert az_calls == []

    def test_no_pathways_returns_unsolved_leaf(self) -> None:
        solver = AutoSolver(
            az_runner=az_always_fails,
            llm_runner=llm_returns_nothing,
            hallucination_mode="none",
        )
        route, solved = solver.solve(CAFFEINE)

        assert solved is False
        assert route == unsolved_leaf(CAFFEINE)

    def test_selects_first_fully_solved_pathway(self) -> None:
        """Two pathways: first unsolvable, second solvable."""

        def az_selective(smiles: str, az_model: str) -> tuple[bool, list]:
            if canonicalize(smiles) == canonicalize(ACETIC_ACID):
                return True, [solved_route(smiles)]
            return False, []

        def llm_two_pathways(molecule: str, **kwargs: Any) -> tuple[list, list, list]:
            if canonicalize(molecule) == canonicalize(ASPIRIN):
                return (
                    [[IBUPROFEN], [ACETIC_ACID]],
                    ["wrong", "right"],
                    [0.3, 0.8],
                )
            return [], [], []

        solver = AutoSolver(
            az_runner=az_selective,
            llm_runner=llm_two_pathways,
            hallucination_mode="none",
        )
        route, solved = solver.solve(ASPIRIN)

        assert solved is True
        assert route["children"][0]["children"][0]["smiles"] == ACETIC_ACID

    def test_partial_result_when_no_pathway_fully_solves(self) -> None:
        def llm_one_unsolvable(molecule: str, **kwargs: Any) -> tuple[list, list, list]:
            if canonicalize(molecule) == canonicalize(ASPIRIN):
                return [[IBUPROFEN]], ["try ibuprofen"], [0.4]
            return [], [], []

        solver = AutoSolver(
            az_runner=az_always_fails,
            llm_runner=llm_one_unsolvable,
            hallucination_mode="none",
        )
        route, solved = solver.solve(ASPIRIN)

        assert solved is False
        assert route["children"][0]["children"] == [unsolved_leaf(IBUPROFEN)]

    def test_multi_reactant_pathway_all_must_solve(self) -> None:
        def az_solves_both(smiles: str, az_model: str) -> tuple[bool, list]:
            canonical = canonicalize(smiles)
            if canonical in {
                canonicalize(PARA_AMINOPHENOL),
                canonicalize(ACETIC_ANHYDRIDE),
            }:
                return True, [solved_route(smiles)]
            return False, []

        def llm_proposes(molecule: str, **kwargs: Any) -> tuple[list, list, list]:
            return (
                [[PARA_AMINOPHENOL, ACETIC_ANHYDRIDE]],
                ["acetylation"],
                [0.8],
            )

        solver = AutoSolver(
            az_runner=az_solves_both,
            llm_runner=llm_proposes,
            hallucination_mode="none",
        )
        route, solved = solver.solve(PARACETAMOL)

        assert solved is True
        children = route["children"][0]["children"]
        assert len(children) == 2

    def test_llm_runner_receives_correct_kwargs(self) -> None:
        captured: dict[str, Any] = {}

        def llm_captures(molecule: str, **kwargs: Any) -> tuple[list, list, list]:
            captured["molecule"] = molecule
            captured.update(kwargs)
            return [[SALICYLIC_ACID]], ["hydrolysis"], [0.7]

        def az_solves_precursor(smiles: str, az_model: str) -> tuple[bool, list]:
            if canonicalize(smiles) == canonicalize(SALICYLIC_ACID):
                return True, [solved_route(smiles)]
            return False, []

        solver = AutoSolver(
            llm="anthropic/claude-sonnet-4-6",
            az_runner=az_solves_precursor,
            llm_runner=llm_captures,
            stability_check=False,
            hallucination_mode="heuristic",
            enable_thinking=False,
            max_output_tokens=64,
        )
        solver.solve(ASPIRIN)

        assert captured["molecule"] == ASPIRIN
        assert captured["model"] == "anthropic/claude-sonnet-4-6"
        assert captured["stability_check"] is False
        assert captured["hallucination_check"] is True
        assert captured["enable_thinking"] is False
        assert captured["max_output_tokens"] == 64


# ---------------------------------------------------------------------------
# run_llm() tests
# ---------------------------------------------------------------------------


class TestRunLlm:
    def test_returns_pathways_from_injected_runner(self) -> None:
        def llm_stub(molecule: str, **kwargs: Any) -> tuple[list, list, list]:
            return [["CCO"]], ["reduce"], [0.8]

        solver = AutoSolver(llm_runner=llm_stub, hallucination_mode="none")
        pathways, explanations, confidence = solver.run_llm(ASPIRIN)

        assert pathways == [["CCO"]]
        assert explanations == ["reduce"]
        assert confidence == [0.8]

    def test_returns_empty_when_runner_returns_nothing(self) -> None:
        solver = AutoSolver(llm_runner=llm_returns_nothing, hallucination_mode="none")
        pathways, explanations, confidence = solver.run_llm(ASPIRIN)

        assert pathways == []
        assert explanations == []
        assert confidence == []


# ---------------------------------------------------------------------------
# single_step() tests
# ---------------------------------------------------------------------------


class TestSingleStep:
    def test_az_success_returns_route(self) -> None:
        def az_solves(smiles: str, az_model: str) -> tuple[bool, list]:
            return True, [solved_route(smiles)]

        solver = AutoSolver(az_runner=az_solves, hallucination_mode="none")
        route, solved = solver.single_step(ASPIRIN)

        assert solved is True
        assert route["smiles"] == ASPIRIN

    def test_llm_builds_one_level_tree(self) -> None:
        def llm_proposes(molecule: str, **kwargs: Any) -> tuple[list, list, list]:
            return (
                [[SALICYLIC_ACID, ACETIC_ACID]],
                ["hydrolysis"],
                [0.9],
            )

        solver = AutoSolver(
            az_runner=az_always_fails,
            llm_runner=llm_proposes,
            hallucination_mode="none",
        )
        route, solved = solver.single_step(ASPIRIN)

        assert solved is False
        children = route["children"][0]["children"]
        assert len(children) == 2
        assert children[0] == unsolved_leaf(SALICYLIC_ACID)
        assert children[1] == unsolved_leaf(ACETIC_ACID)

    def test_no_pathways_returns_unsolved(self) -> None:
        solver = AutoSolver(
            az_runner=az_always_fails,
            llm_runner=llm_returns_nothing,
            hallucination_mode="none",
        )
        route, solved = solver.single_step(CAFFEINE)

        assert solved is False
        assert route == unsolved_leaf(CAFFEINE)


# ---------------------------------------------------------------------------
# parse() tests — exercises real RetrosynthesisRouteParser
# ---------------------------------------------------------------------------


class TestParse:
    def test_aspirin_hydrolysis_has_correct_product(self) -> None:
        solver = AutoSolver(hallucination_mode="none")
        result = solver.parse(aspirin_hydrolysis_tree(), solved=True)

        assert result["solved"] is True
        assert len(result["steps"]) >= 1
        top_step = result["steps"][0]
        assert top_step["products"][0]["smiles"] == ASPIRIN
        reactant_smiles = {r["smiles"] for r in top_step["reactants"]}
        assert len(reactant_smiles) == 2

    def test_solved_flag_propagates_false(self) -> None:
        solver = AutoSolver(hallucination_mode="none")
        tree = reaction_tree(ASPIRIN, [unsolved_leaf(SALICYLIC_ACID)], [0.5])
        result = solver.parse(tree, solved=False)

        assert result["solved"] is False

    def test_two_step_route_has_dependencies(self) -> None:
        inner_tree = reaction_tree(
            PARA_AMINOPHENOL,
            [unsolved_leaf("Nc1ccccc1"), unsolved_leaf("O")],
            [0.6],
        )
        outer_tree = reaction_tree(
            PARACETAMOL,
            [inner_tree, unsolved_leaf(ACETIC_ANHYDRIDE)],
            [0.8],
        )
        solver = AutoSolver(hallucination_mode="none")
        result = solver.parse(outer_tree, solved=False)

        assert len(result["steps"]) >= 2
        deps = result["dependencies"]
        assert any(len(v) > 0 for v in deps.values())


# ---------------------------------------------------------------------------
# add_metadata() tests — injectable recommenders, no patching
# ---------------------------------------------------------------------------

try:
    import litellm  # noqa: F401

    _has_litellm = True
except ModuleNotFoundError:
    _has_litellm = False


@pytest.mark.skipif(not _has_litellm, reason="litellm not installed")
class TestAddMetadata:
    def test_enriches_aspirin_step(
        self,
        spy_reagent_recommender,
        spy_conditions_recommender,
        spy_literature_recommender,
        recommender_calls,
    ) -> None:
        solver = AutoSolver(
            hallucination_mode="none", metadata_model="anthropic/claude-sonnet-4-6"
        )
        parsed = solver.parse(aspirin_hydrolysis_tree(), solved=True)
        result = solver.add_metadata(
            parsed,
            reagent_recommender=spy_reagent_recommender,
            conditions_recommender=spy_conditions_recommender,
            literature_recommender=spy_literature_recommender,
        )

        assert any("reagent" in c for c in recommender_calls)
        step = result["steps"][0]
        assert len(step["reagents"]) > 0
        assert step["conditions"]["temperature"] == "25 C"

    def test_skips_step_when_reagent_fails(
        self,
        failing_reagent_recommender,
    ) -> None:
        solver = AutoSolver(hallucination_mode="none")
        parsed = solver.parse(aspirin_hydrolysis_tree(), solved=True)
        original_reagents = list(parsed["steps"][0]["reagents"])

        result = solver.add_metadata(
            parsed,
            reagent_recommender=failing_reagent_recommender,
        )

        assert result["steps"][0]["reagents"] == original_reagents

    def test_passes_metadata_model(
        self,
        spy_reagent_recommender,
        spy_conditions_recommender,
        spy_literature_recommender,
        recommender_calls,
    ) -> None:
        solver = AutoSolver(
            hallucination_mode="none",
            metadata_model="anthropic/claude-sonnet-4-6",
        )
        parsed = solver.parse(aspirin_hydrolysis_tree(), solved=True)
        solver.add_metadata(
            parsed,
            reagent_recommender=spy_reagent_recommender,
            conditions_recommender=spy_conditions_recommender,
            literature_recommender=spy_literature_recommender,
        )

        assert all("anthropic/claude-sonnet-4-6" in c for c in recommender_calls)

    def test_skips_step_with_no_reactants(
        self,
        spy_reagent_recommender,
        recommender_calls,
    ) -> None:
        parsed: dict[str, Any] = {
            "steps": [
                {
                    "step": "1",
                    "reactants": [],
                    "reagents": [],
                    "products": [{"smiles": ASPIRIN, "product_metadata": {}}],
                    "conditions": [],
                    "reactionmetrics": [
                        {"scalabilityindex": 0, "closestliterature": ""}
                    ],
                }
            ],
            "dependencies": {"1": []},
            "solved": True,
        }
        solver = AutoSolver(hallucination_mode="none")
        solver.add_metadata(parsed, reagent_recommender=spy_reagent_recommender)

        assert recommender_calls == []


# ---------------------------------------------------------------------------
# autosolve() tests — subclass override, no patching
# ---------------------------------------------------------------------------


class _TrackingAutoSolver(AutoSolver):
    """Subclass that records call order without patching."""

    def __init__(self, **kwargs: Any) -> None:
        super().__init__(**kwargs)
        self.call_log: list[str] = []

    def solve(
        self,
        smiles: str,
        visited: set[str] | None = None,
        depth: int = 0,
    ) -> tuple[dict[str, Any], bool]:
        self.call_log.append("solve")
        return unsolved_leaf(smiles), False

    def parse(self, route_tree: dict[str, Any], *, solved: bool) -> dict[str, Any]:
        self.call_log.append("parse")
        return {"steps": [], "dependencies": {}, "solved": solved}

    def add_metadata(self, parsed_output: Any, **kwargs: Any) -> Any:
        self.call_log.append("add_metadata")
        return parsed_output


class TestAutosolve:
    def test_chains_solve_parse_add_metadata(self) -> None:
        solver = _TrackingAutoSolver(hallucination_mode="none")
        result = solver.autosolve(ASPIRIN)

        assert solver.call_log == ["solve", "parse", "add_metadata"]
        assert result["solved"] is False

    def test_propagates_solved_true(self) -> None:
        class _SolvedTracker(_TrackingAutoSolver):
            def solve(self, smiles, visited=None, depth=0):
                self.call_log.append("solve")
                return solved_route(smiles), True

        solver = _SolvedTracker(hallucination_mode="none")
        result = solver.autosolve(ASPIRIN)

        assert result["solved"] is True
