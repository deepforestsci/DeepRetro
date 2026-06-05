"""Tests for the AutoSolver retrosynthesis pipeline."""

from __future__ import annotations

import importlib.util
from typing import Any

import pytest

from deepretro.algorithms.autosolve import (
    AutoSolver,
    reaction_tree,
    unsolved_leaf,
)
from deepretro.utils.utils_molecule import canonicalize

ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"
SALICYLIC_ACID = "OC(=O)c1ccccc1O"
ACETIC_ANHYDRIDE = "CC(=O)OC(C)=O"
PARACETAMOL = "CC(=O)Nc1ccc(O)cc1"
PARA_AMINOPHENOL = "Nc1ccc(O)cc1"
IBUPROFEN = "CC(C)Cc1ccc(cc1)C(C)C(=O)O"
CAFFEINE = "Cn1c(=O)c2c(ncn2C)n(C)c1=O"
ACETIC_ACID = "CC(=O)O"


def solved_route(smiles: str) -> dict[str, Any]:
    """Minimal in-stock route node for a solved molecule."""
    return {
        "type": "mol",
        "smiles": smiles,
        "is_chemical": True,
        "in_stock": True,
    }


def az_always_fails(smiles: str, az_model: str) -> tuple[bool, list[Any]]:
    """AZ runner that never solves any molecule."""
    return False, []


def az_solves_all(smiles: str, az_model: str) -> tuple[bool, list[Any]]:
    """AZ runner that solves every molecule."""
    return True, [solved_route(smiles)]


class SelectiveAzRunner:
    """AZ runner that solves only configured canonical molecules."""

    def __init__(self, solvable: set[str]) -> None:
        self.solvable = {canonicalize(smiles) for smiles in solvable}

    def __call__(self, smiles: str, az_model: str) -> tuple[bool, list[Any]]:
        if canonicalize(smiles) in self.solvable:
            return True, [solved_route(smiles)]
        return False, []


class RecordingAzRunner:
    """AZ runner that records calls and never solves."""

    def __init__(self) -> None:
        self.calls: list[str] = []

    def __call__(self, smiles: str, az_model: str) -> tuple[bool, list[Any]]:
        self.calls.append(smiles)
        return False, []


class RecordingLlmRunner:
    """LLM runner that records calls and returns no pathways."""

    def __init__(self) -> None:
        self.calls: list[str] = []

    def __call__(
        self, molecule: str, **kwargs: Any
    ) -> tuple[list[list[str]], list[str], list[float]]:
        self.calls.append(molecule)
        return [], [], []


class CapturingLlmRunner:
    """LLM runner that records call arguments and returns a solvable pathway."""

    def __init__(self) -> None:
        self.captured: dict[str, Any] = {}

    def __call__(
        self, molecule: str, **kwargs: Any
    ) -> tuple[list[list[str]], list[str], list[float]]:
        self.captured["molecule"] = molecule
        self.captured.update(kwargs)
        return [[SALICYLIC_ACID]], ["hydrolysis"], [0.7]


def llm_returns_nothing(
    molecule: str, **kwargs: Any
) -> tuple[list[list[str]], list[str], list[float]]:
    """LLM runner that returns no pathways."""
    return [], [], []


def llm_proposes_hydrolysis(
    molecule: str, **kwargs: Any
) -> tuple[list[list[str]], list[str], list[float]]:
    """LLM runner proposing aspirin hydrolysis."""
    return [[SALICYLIC_ACID, ACETIC_ACID]], ["hydrolysis"], [0.9]


def llm_returns_aspirin_self_reference(
    molecule: str, **kwargs: Any
) -> tuple[list[list[str]], list[str], list[float]]:
    """LLM runner returning aspirin in a non-canonical equivalent form."""
    return [["OC(=O)c1ccccc1OC(=O)C"]], ["self-reference"], [0.5]


def llm_two_pathways(
    molecule: str, **kwargs: Any
) -> tuple[list[list[str]], list[str], list[float]]:
    """LLM runner where only the second pathway can fully solve."""
    if canonicalize(molecule) == canonicalize(ASPIRIN):
        return [[IBUPROFEN], [ACETIC_ACID]], ["wrong", "right"], [0.3, 0.8]
    return [], [], []


def llm_unsolvable_ibuprofen(
    molecule: str, **kwargs: Any
) -> tuple[list[list[str]], list[str], list[float]]:
    """LLM runner proposing an unsolved ibuprofen branch for aspirin."""
    if canonicalize(molecule) == canonicalize(ASPIRIN):
        return [[IBUPROFEN]], ["try ibuprofen"], [0.4]
    return [], [], []


def llm_proposes_paracetamol_acetylation(
    molecule: str, **kwargs: Any
) -> tuple[list[list[str]], list[str], list[float]]:
    """LLM runner proposing paracetamol acetylation."""
    return [[PARA_AMINOPHENOL, ACETIC_ANHYDRIDE]], ["acetylation"], [0.8]


def llm_empty_then_acetic_acid(
    molecule: str, **kwargs: Any
) -> tuple[list[list[str]], list[str], list[float]]:
    """LLM runner returning an empty pathway before a valid pathway."""
    return [[], [ACETIC_ACID]], ["empty", "real"], [0.1, 0.9]


def llm_empty_then_hydrolysis(
    molecule: str, **kwargs: Any
) -> tuple[list[list[str]], list[str], list[float]]:
    """LLM runner returning an empty pathway before aspirin hydrolysis."""
    return [[], [SALICYLIC_ACID, ACETIC_ACID]], ["empty", "hydrolysis"], [0.1, 0.9]


def llm_stub_reduction(
    molecule: str, **kwargs: Any
) -> tuple[list[list[str]], list[str], list[float]]:
    """LLM runner returning one reduction pathway."""
    return [["CCO"]], ["reduce"], [0.8]


def aspirin_hydrolysis_tree() -> dict[str, Any]:
    """One-step route: aspirin -> salicylic acid + acetic anhydride."""
    return reaction_tree(
        ASPIRIN,
        [unsolved_leaf(SALICYLIC_ACID), unsolved_leaf(ACETIC_ANHYDRIDE)],
        [0.9],
    )


class TestUnsolvedLeaf:
    def test_aspirin(self) -> None:
        """Unsolved leaf for aspirin has the expected mol-node fields."""
        node = unsolved_leaf(ASPIRIN)
        assert node["smiles"] == ASPIRIN
        assert node["type"] == "mol"
        assert node["is_chemical"] is True
        assert node["in_stock"] is False
        assert node["children"] == []

    def test_caffeine(self) -> None:
        """Unsolved leaf for caffeine is marked not in stock."""
        node = unsolved_leaf(CAFFEINE)
        assert node["smiles"] == CAFFEINE
        assert node["in_stock"] is False


class TestReactionTree:
    def test_aspirin_hydrolysis(self) -> None:
        """Reaction tree nests children and stores confidence as policy_probability."""
        children = [unsolved_leaf(SALICYLIC_ACID), unsolved_leaf(ACETIC_ANHYDRIDE)]
        tree = reaction_tree(ASPIRIN, children, [0.85])
        assert tree["smiles"] == ASPIRIN
        rxn = tree["children"][0]
        assert rxn["type"] == "reaction"
        assert rxn["metadata"]["policy_probability"] == pytest.approx(0.85)
        assert rxn["children"] == children

    def test_empty_confidence(self) -> None:
        """Empty confidence defaults policy_probability to 0.0."""
        tree = reaction_tree(ASPIRIN, [], [])
        assert tree["children"][0]["metadata"]["policy_probability"] == 0.0


class TestConstructor:
    def test_preserves_options(self) -> None:
        """Constructor stores all explicitly provided options."""
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
        """Constructor applies the documented default options."""
        solver = AutoSolver(hallucination_mode="none")
        assert solver.llm == "anthropic/claude-sonnet-4-6"
        assert solver.az_model == "Pistachio_100+"
        assert solver.stability_check is True
        assert solver.max_depth == 50
        assert solver.enable_thinking is True
        assert solver.max_output_tokens is None
        assert solver.metadata_model == "anthropic/claude-sonnet-4-6"

    def test_normalizes_hallucination_mode(self) -> None:
        """hallucination_mode is lower-cased on construction."""
        solver = AutoSolver(hallucination_mode="NONE")
        assert solver.hallucination_mode == "none"

    def test_rejects_negative_max_depth(self) -> None:
        """A negative max_depth raises ValueError."""
        with pytest.raises(ValueError, match="max_depth"):
            AutoSolver(hallucination_mode="none", max_depth=-1)


class TestSolve:
    def test_az_success_returns_route(self) -> None:
        """AZ success returns the AZ route marked solved."""
        solver = AutoSolver(az_runner=az_solves_all, hallucination_mode="none")
        route, solved = solver.solve(ASPIRIN)

        assert solved is True
        assert route["smiles"] == ASPIRIN
        assert route["in_stock"] is True

    def test_az_success_skips_llm(self) -> None:
        """LLM fallback is skipped when AZ already solves the target."""
        llm_runner = RecordingLlmRunner()
        solver = AutoSolver(
            az_runner=az_solves_all,
            llm_runner=llm_runner,
            hallucination_mode="none",
        )
        solver.solve(ASPIRIN)
        assert llm_runner.calls == []

    def test_llm_fallback_with_solvable_precursors(self) -> None:
        """LLM fallback solves when its precursors are AZ-solvable."""
        solver = AutoSolver(
            az_runner=SelectiveAzRunner({SALICYLIC_ACID, ACETIC_ACID}),
            llm_runner=llm_proposes_hydrolysis,
            hallucination_mode="none",
        )
        route, solved = solver.solve(ASPIRIN)

        assert solved is True
        assert len(route["children"][0]["children"]) == 2

    def test_cycle_detection_with_non_canonical_form(self) -> None:
        """A precursor equal to the target (non-canonical) is detected as a cycle."""
        solver = AutoSolver(
            az_runner=az_always_fails,
            llm_runner=llm_returns_aspirin_self_reference,
            hallucination_mode="none",
        )
        route, solved = solver.solve(ASPIRIN)

        assert solved is False
        child = route["children"][0]["children"][0]
        assert canonicalize(child["smiles"]) == canonicalize(ASPIRIN)

    def test_max_depth_zero_skips_az(self) -> None:
        """max_depth=0 returns an unsolved leaf without calling AZ."""
        az_runner = RecordingAzRunner()
        solver = AutoSolver(az_runner=az_runner, hallucination_mode="none", max_depth=0)
        route, solved = solver.solve(ASPIRIN)

        assert (route, solved) == (unsolved_leaf(ASPIRIN), False)
        assert az_runner.calls == []

    def test_empty_input_skips_az_and_llm(self) -> None:
        """Blank input short-circuits before AZ and LLM are called."""
        az_runner = RecordingAzRunner()
        llm_runner = RecordingLlmRunner()
        solver = AutoSolver(
            az_runner=az_runner,
            llm_runner=llm_runner,
            hallucination_mode="none",
        )
        route, solved = solver.solve("   ")

        assert (route, solved) == (unsolved_leaf(""), False)
        assert az_runner.calls == []
        assert llm_runner.calls == []

    def test_no_pathways_returns_unsolved(self) -> None:
        """No LLM pathways yields an unsolved leaf."""
        solver = AutoSolver(
            az_runner=az_always_fails,
            llm_runner=llm_returns_nothing,
            hallucination_mode="none",
        )
        route, solved = solver.solve(CAFFEINE)

        assert solved is False
        assert route == unsolved_leaf(CAFFEINE)

    def test_selects_first_fully_solved_pathway_with_correct_confidence(self) -> None:
        """The first fully solved pathway is chosen with its own confidence."""
        solver = AutoSolver(
            az_runner=SelectiveAzRunner({ACETIC_ACID}),
            llm_runner=llm_two_pathways,
            hallucination_mode="none",
        )
        route, solved = solver.solve(ASPIRIN)

        assert solved is True
        assert route["children"][0]["children"][0]["smiles"] == ACETIC_ACID
        assert route["children"][0]["metadata"]["policy_probability"] == pytest.approx(
            0.8
        )

    def test_partial_result_when_nothing_fully_solves(self) -> None:
        """When nothing fully solves, the first attempted pathway is returned unsolved."""
        solver = AutoSolver(
            az_runner=az_always_fails,
            llm_runner=llm_unsolvable_ibuprofen,
            hallucination_mode="none",
        )
        route, solved = solver.solve(ASPIRIN)

        assert solved is False
        assert route["children"][0]["children"] == [unsolved_leaf(IBUPROFEN)]

    def test_multi_reactant_pathway(self) -> None:
        """A multi-reactant pathway solves when all reactants are AZ-solvable."""
        solver = AutoSolver(
            az_runner=SelectiveAzRunner({PARA_AMINOPHENOL, ACETIC_ANHYDRIDE}),
            llm_runner=llm_proposes_paracetamol_acetylation,
            hallucination_mode="none",
        )
        route, solved = solver.solve(PARACETAMOL)

        assert solved is True
        assert len(route["children"][0]["children"]) == 2

    def test_empty_pathway_is_not_treated_as_solved(self) -> None:
        """An empty pathway is skipped in favour of the next candidate."""
        solver = AutoSolver(
            az_runner=SelectiveAzRunner({ACETIC_ACID}),
            llm_runner=llm_empty_then_acetic_acid,
            hallucination_mode="none",
        )
        route, solved = solver.solve(ASPIRIN)

        assert solved is True
        assert route["children"][0]["children"][0]["smiles"] == ACETIC_ACID
        assert route["children"][0]["metadata"]["policy_probability"] == pytest.approx(
            0.9
        )

    def test_llm_runner_receives_correct_kwargs(self) -> None:
        """Solver forwards its configured options to the LLM runner."""
        llm_runner = CapturingLlmRunner()
        solver = AutoSolver(
            llm="anthropic/claude-sonnet-4-6",
            az_runner=SelectiveAzRunner({SALICYLIC_ACID}),
            llm_runner=llm_runner,
            stability_check=False,
            hallucination_mode="heuristic",
            enable_thinking=False,
            max_output_tokens=64,
        )
        solver.solve(ASPIRIN)

        captured = llm_runner.captured
        assert captured["molecule"] == ASPIRIN
        assert captured["model"] == "anthropic/claude-sonnet-4-6"
        assert captured["stability_check"] is False
        assert captured["hallucination_check"] is True
        assert captured["enable_thinking"] is False
        assert captured["max_output_tokens"] == 64


class TestRunLlm:
    def test_returns_pathways_from_injected_runner(self) -> None:
        """run_llm returns the injected runner's pathways, explanations, and confidence."""
        solver = AutoSolver(llm_runner=llm_stub_reduction, hallucination_mode="none")
        pathways, explanations, confidence = solver.run_llm(ASPIRIN)

        assert pathways == [["CCO"]]
        assert explanations == ["reduce"]
        assert confidence == [0.8]

    def test_empty_runner(self) -> None:
        """run_llm returns empty results when the runner yields nothing."""
        solver = AutoSolver(llm_runner=llm_returns_nothing, hallucination_mode="none")
        pathways, explanations, confidence = solver.run_llm(ASPIRIN)

        assert pathways == []
        assert explanations == []
        assert confidence == []


class TestSingleStep:
    def test_az_success(self) -> None:
        """single_step returns the AZ route when AZ solves the target."""
        solver = AutoSolver(az_runner=az_solves_all, hallucination_mode="none")
        route, solved = solver.single_step(ASPIRIN)

        assert solved is True
        assert route["smiles"] == ASPIRIN

    def test_llm_builds_one_level_tree(self) -> None:
        """single_step builds a one-level tree of unsolved precursors."""
        solver = AutoSolver(
            az_runner=az_always_fails,
            llm_runner=llm_proposes_hydrolysis,
            hallucination_mode="none",
        )
        route, solved = solver.single_step(ASPIRIN)

        assert solved is False
        children = route["children"][0]["children"]
        assert children[0] == unsolved_leaf(SALICYLIC_ACID)
        assert children[1] == unsolved_leaf(ACETIC_ACID)

    def test_empty_first_pathway_uses_first_non_empty_candidate(self) -> None:
        """single_step skips an empty pathway and uses the next candidate."""
        solver = AutoSolver(
            az_runner=az_always_fails,
            llm_runner=llm_empty_then_hydrolysis,
            hallucination_mode="none",
        )
        route, solved = solver.single_step(ASPIRIN)

        assert solved is False
        assert route["children"][0]["metadata"]["policy_probability"] == pytest.approx(
            0.9
        )
        assert route["children"][0]["children"] == [
            unsolved_leaf(SALICYLIC_ACID),
            unsolved_leaf(ACETIC_ACID),
        ]

    def test_no_pathways(self) -> None:
        """single_step returns an unsolved leaf when no pathways are produced."""
        solver = AutoSolver(
            az_runner=az_always_fails,
            llm_runner=llm_returns_nothing,
            hallucination_mode="none",
        )
        route, solved = solver.single_step(CAFFEINE)

        assert solved is False
        assert route == unsolved_leaf(CAFFEINE)


class TestParse:
    def test_aspirin_hydrolysis_product(self) -> None:
        """parse exposes the product SMILES and the solved flag."""
        solver = AutoSolver(hallucination_mode="none")
        result = solver.parse(aspirin_hydrolysis_tree(), solved=True)

        assert result["solved"] is True
        assert len(result["steps"]) >= 1
        assert result["steps"][0]["products"][0]["smiles"] == ASPIRIN

    def test_solved_false_propagates(self) -> None:
        """parse propagates solved=False into the output."""
        solver = AutoSolver(hallucination_mode="none")
        tree = reaction_tree(ASPIRIN, [unsolved_leaf(SALICYLIC_ACID)], [0.5])
        result = solver.parse(tree, solved=False)

        assert result["solved"] is False

    def test_two_step_route_dependencies(self) -> None:
        """parse emits multiple steps with dependency links for nested routes."""
        inner = reaction_tree(
            PARA_AMINOPHENOL,
            [unsolved_leaf("Nc1ccccc1"), unsolved_leaf("O")],
            [0.6],
        )
        outer = reaction_tree(
            PARACETAMOL,
            [inner, unsolved_leaf(ACETIC_ANHYDRIDE)],
            [0.8],
        )
        solver = AutoSolver(hallucination_mode="none")
        result = solver.parse(outer, solved=False)

        assert len(result["steps"]) >= 2
        assert any(len(v) > 0 for v in result["dependencies"].values())


HAS_LITELLM = importlib.util.find_spec("litellm") is not None


@pytest.mark.skipif(
    not HAS_LITELLM, reason="litellm required for metadata recommendation"
)
class TestAddMetadata:
    def test_enriches_step(
        self,
        spy_reagent_recommender,
        spy_conditions_recommender,
        spy_literature_recommender,
        recommender_calls,
    ) -> None:
        """add_metadata fills reagents and conditions from the recommenders."""
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
        """A failing reagent recommender leaves the step unchanged."""
        solver = AutoSolver(hallucination_mode="none")
        parsed = solver.parse(aspirin_hydrolysis_tree(), solved=True)
        original_reagents = list(parsed["steps"][0]["reagents"])

        result = solver.add_metadata(
            parsed,
            reagent_recommender=failing_reagent_recommender,
        )

        assert result["steps"][0]["reagents"] == original_reagents

    def test_forwards_metadata_model(
        self,
        spy_reagent_recommender,
        spy_conditions_recommender,
        spy_literature_recommender,
        recommender_calls,
    ) -> None:
        """The configured metadata_model is forwarded to the recommenders."""
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
        """Steps without reactants are skipped (no recommender calls)."""
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

    def test_adds_missing_reactionmetrics_for_enriched_step(
        self,
        spy_reagent_recommender,
        spy_conditions_recommender,
        spy_literature_recommender,
    ) -> None:
        """Missing reactionmetrics are created and closestliterature holds the DOI string."""
        parsed = {
            "steps": [
                {
                    "step": "1",
                    "reactants": [{"smiles": SALICYLIC_ACID}],
                    "reagents": [],
                    "products": [{"smiles": ASPIRIN, "product_metadata": {}}],
                    "conditions": {},
                    "reactionmetrics": [],
                }
            ],
            "dependencies": {"1": []},
            "solved": True,
        }
        solver = AutoSolver(hallucination_mode="none")

        result = solver.add_metadata(
            parsed,
            reagent_recommender=spy_reagent_recommender,
            conditions_recommender=spy_conditions_recommender,
            literature_recommender=spy_literature_recommender,
        )

        metrics = result["steps"][0]["reactionmetrics"]
        assert metrics[0]["closestliterature"] == "10.1000/example"


class TestAutosolve:
    def test_returns_parsed_result_for_az_solved_target(self) -> None:
        """autosolve returns a parsed, solved result for an AZ-solved target."""
        solver = AutoSolver(az_runner=az_solves_all, hallucination_mode="none")
        result = solver.autosolve(ASPIRIN)

        assert result["solved"] is True
        assert result["steps"] == []
        assert result["dependencies"] == {}
