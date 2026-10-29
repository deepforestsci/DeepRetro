"""Offline regression tests for heuristic gating and local route repair."""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Any

import pytest

from deepretro.algorithms.guarded import GuardedSolver


def _leaf(smiles: str, in_stock: bool = True) -> dict[str, Any]:
    return {"type": "mol", "smiles": smiles, "in_stock": in_stock, "children": []}


def _route(product: str, *children: dict[str, Any]) -> dict[str, Any]:
    return {
        "type": "mol",
        "smiles": product,
        "children": [{"type": "reaction", "children": list(children)}],
    }


class Checker:
    def __init__(
        self, rejected: set[tuple[str, tuple[str, ...]]] | None = None
    ) -> None:
        self.rejected = rejected or set()
        self.calls: list[tuple[str, tuple[str, ...]]] = []

    def assess(self, product: str, reactants: list[str]) -> dict[str, Any]:
        identity = (product, tuple(reactants))
        self.calls.append(identity)
        return {
            "source": "heuristic",
            "flagged": identity in self.rejected,
            "message": "missing precursor",
        }


class Backend:
    hallucination_mode = "heuristic"
    solve_mode = "pipeline"

    def __init__(
        self,
        routes: dict[str, list[dict[str, Any]]] | None = None,
        proposals: dict[str, list[list[str]]] | None = None,
    ) -> None:
        self.routes = routes or {}
        self.proposals = proposals or {}
        self.hallucination_checker: Any = Checker()
        self.az_calls: list[str] = []
        self.llm_calls: list[tuple[str, str]] = []
        self.error = False

    def _run_az(self, smiles: str) -> tuple[bool, list[dict[str, Any]], bool]:
        self.az_calls.append(smiles)
        return True, self.routes.get(smiles, []), self.error

    def run_llm(
        self, smiles: str, *, depth: int, rejection_feedback: str
    ) -> tuple[list[list[str]], list[str], list[float]]:
        self.llm_calls.append((smiles, rejection_feedback))
        return self.proposals.get(smiles, []), [], []


def test_rejected_az_step_repaired_with_feedback(tmp_path: Path) -> None:
    backend = Backend({"CC=O": [_route("CC=O", _leaf("C"))]}, {"CC=O": [["CCO"]]})
    backend.hallucination_checker.rejected.add(("CC=O", ("C",)))
    solver = GuardedSolver(backend, {"CCO"}.__contains__, tmp_path / "run")
    tree, passed = solver.solve("CC=O")
    assert passed
    assert tree["children"][0]["children"][0]["smiles"] == "CCO"
    assert json.loads(backend.llm_calls[0][1])[0]["reason"] == "missing precursor"
    assert backend.hallucination_checker.calls == [("CC=O", ("C",)), ("CC=O", ("CCO",))]
    assert "C" not in backend.az_calls


def test_stock_flags_do_not_establish_availability(tmp_path: Path) -> None:
    backend = Backend({"C": [_leaf("C")]})
    solver = GuardedSolver(backend, lambda _: False, tmp_path / "run", rounds=0)
    tree, passed = solver.solve("C")
    assert not passed and tree["in_stock"] is False


def test_cycle_rejection_does_not_poison_an_independent_branch(tmp_path: Path) -> None:
    """A previously cyclic disconnection is retried outside its ancestor chain."""

    class ChangingBackend(Backend):
        def run_llm(
            self, smiles: str, *, depth: int, rejection_feedback: str
        ) -> tuple[list[list[str]], list[str], list[float]]:
            count = sum(molecule == smiles for molecule, _ in self.llm_calls)
            self.llm_calls.append((smiles, rejection_feedback))
            if smiles == "CCC":
                return ([["CC"]] if count == 0 else [["C"]]), [], []
            return self.proposals.get(smiles, []), [], []

    backend = ChangingBackend(proposals={"CCCC": [["CCC"], ["CC"]], "CC": [["CCC"]]})
    solver = GuardedSolver(backend, {"C"}.__contains__, tmp_path / "run", rounds=1)
    tree, passed = solver.solve("CCCC")
    assert passed
    assert tree["children"][0]["children"][0]["smiles"] == "CC"
    assert ("CC", ("CCC",)) not in solver.rejected
    assert (
        json.loads([feedback for mol, feedback in backend.llm_calls if mol == "CC"][-1])
        == []
    )


def test_explicit_stock_bypasses_provider_calls(tmp_path: Path) -> None:
    backend = Backend()
    solver = GuardedSolver(backend, {"CCO"}.__contains__, tmp_path / "run", max_depth=0)
    tree, passed = solver.solve("OCC")
    assert passed and tree["stock_source"] == "catalog"
    assert backend.az_calls == backend.llm_calls == []


def test_depth_limit_does_not_expand_reaction(tmp_path: Path) -> None:
    backend = Backend({"CC=O": [_route("CC=O", _leaf("CCO"))]})
    solver = GuardedSolver(backend, {"CCO"}.__contains__, tmp_path / "run", max_depth=0)
    assert not solver.solve("CC=O")[1]
    assert backend.hallucination_checker.calls == []


def test_all_nested_az_edges_are_gated(tmp_path: Path) -> None:
    backend = Backend({"CC(=O)O": [_route("CC(=O)O", _route("CC=O", _leaf("CCO")))]})
    solver = GuardedSolver(backend, {"CCO"}.__contains__, tmp_path / "run")
    assert solver.solve("CC(=O)O")[1]
    assert backend.hallucination_checker.calls == [
        ("CC(=O)O", ("CC=O",)),
        ("CC=O", ("CCO",)),
    ]


def test_repair_keeps_completed_sibling(tmp_path: Path) -> None:
    backend = Backend(
        {"CCC": [_route("CCC", _leaf("C"), _route("CC", _leaf("N")))]}, {"CC": [["CO"]]}
    )
    backend.hallucination_checker.rejected.add(("CC", ("N",)))
    solver = GuardedSolver(backend, {"C", "CO"}.__contains__, tmp_path / "run")
    tree, passed = solver.solve("CCC")
    assert passed
    children = tree["children"][0]["children"]
    assert children[0]["smiles"] == "C" and children[0]["in_stock"]
    assert children[1]["children"][0]["children"][0]["smiles"] == "CO"
    assert [p for p, _ in backend.llm_calls] == ["CC"]
    snapshots = list((tmp_path / "run").glob("checkpoint_*.json"))
    assert any(json.loads(p.read_text()).get("children") for p in snapshots)


def test_cycle_is_rejected_before_checker(tmp_path: Path) -> None:
    backend = Backend({"CCO": [_route("CCO", _leaf("OCC"))]})
    solver = GuardedSolver(backend, lambda _: False, tmp_path / "run", rounds=0)
    assert not solver.solve("CCO")[1]
    assert backend.hallucination_checker.calls == []
    assert "ancestor cycle" in (tmp_path / "run" / "repair_events.jsonl").read_text()


@pytest.mark.parametrize(
    "candidate",
    [
        {"smiles": "CCO", "children": "bad"},
        {"smiles": "CCO", "children": [{"children": [None]}]},
        {
            "smiles": "CCO",
            "children": [{"children": [{"type": "reaction", "smiles": "C"}]}],
        },
        {"smiles": "CC", "in_stock": True},
        {"smiles": "CCO", "children": [{"children": []}]},
    ],
)
def test_malformed_az_routes_fail_closed(
    tmp_path: Path, candidate: dict[str, Any]
) -> None:
    backend = Backend({"CCO": [candidate]})
    solver = GuardedSolver(backend, lambda _: False, tmp_path / "run", rounds=0)
    assert not solver.solve("CCO")[1]


def test_backend_failure_is_not_chemical_rejection(tmp_path: Path) -> None:
    backend = Backend()
    backend.error = True
    solver = GuardedSolver(backend, lambda _: False, tmp_path / "run")
    with pytest.raises(RuntimeError, match="operationally incomplete"):
        solver.solve("CCO")
    assert solver.rejections == {}
    assert list((tmp_path / "run").glob("checkpoint_*.json"))


def test_deadline_stops_before_external_calls(tmp_path: Path) -> None:
    backend = Backend()
    solver = GuardedSolver(
        backend,
        lambda _: pytest.fail("stock call after deadline"),
        tmp_path / "run",
        deadline=time.monotonic() - 1,
    )
    assert not solver.solve("CCO")[1]
    assert backend.az_calls == backend.llm_calls == []


def test_repeated_solve_rechecks_stock_and_does_not_overwrite(tmp_path: Path) -> None:
    backend = Backend()
    stock = {"CCO"}
    solver = GuardedSolver(backend, stock.__contains__, tmp_path / "run", rounds=0)
    assert solver.solve("CCO")[1]
    before = {
        p.name: p.read_bytes() for p in (tmp_path / "run").glob("checkpoint_*.json")
    }
    stock.clear()
    assert not solver.solve("CCO")[1]
    assert all(
        (tmp_path / "run" / name).read_bytes() == data for name, data in before.items()
    )


@pytest.mark.parametrize(
    "option,value", [("max_depth", -1), ("rounds", -1), ("rounds", True)]
)
def test_invalid_budgets(tmp_path: Path, option: str, value: int) -> None:
    with pytest.raises(ValueError):
        GuardedSolver(Backend(), lambda _: False, tmp_path / "run", **{option: value})


def test_requires_heuristic_pipeline(tmp_path: Path) -> None:
    backend = Backend()
    backend.hallucination_mode = "ml"
    with pytest.raises(ValueError, match="heuristic pipeline"):
        GuardedSolver(backend, lambda _: False, tmp_path / "run")


def test_existing_output_directory_is_preserved(tmp_path: Path) -> None:
    marker = tmp_path / "marker"
    marker.write_text("keep")
    with pytest.raises(FileExistsError):
        GuardedSolver(Backend(), lambda _: False, tmp_path)
    assert marker.read_text() == "keep"


def test_real_heuristic_rejects_missing_carbon_skeleton(tmp_path: Path) -> None:
    from deepretro.models.hallucination_checker import HallucinationChecker

    backend = Backend({"c1ccccc1": [_route("c1ccccc1", _leaf("CC"))]})
    backend.hallucination_checker = HallucinationChecker(checker_type="heuristic")
    solver = GuardedSolver(backend, {"CC"}.__contains__, tmp_path / "run", rounds=0)
    assert not solver.solve("c1ccccc1")[1]
    assert solver.rejections["c1ccccc1"]


def test_unassessable_step_is_not_accepted(tmp_path: Path) -> None:
    backend = Backend({"CCO": [_route("CCO", _leaf("CC"))]})

    class Unassessable:
        def assess(self, product: str, reactants: list[str]) -> dict[str, Any]:
            return {"source": "heuristic", "flagged": None, "message": "unsupported"}

    backend.hallucination_checker = Unassessable()
    solver = GuardedSolver(backend, {"CC"}.__contains__, tmp_path / "run", rounds=0)
    assert not solver.solve("CCO")[1]


def test_output_integrates_with_native_parser(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from deepretro.algorithms.autosolve import AutoSolver

    backend = AutoSolver(
        hallucination_mode="heuristic", az_runner=lambda *_: (False, [])
    )
    monkeypatch.setattr(
        backend, "run_llm", lambda *_args, **_kwargs: ([["CCO"]], ["oxidation"], [0.8])
    )
    solver = GuardedSolver(backend, {"CCO"}.__contains__, tmp_path / "run")
    tree, passed = solver.solve("CC=O")
    assert passed
    parsed = backend.parse(tree, solved=passed)
    assert parsed["solved"]
    assert parsed["steps"][0]["products"][0]["smiles"] == "CC=O"
    assert tree["children"][0]["metadata"]["hallucination"]["flagged"] is False
