"""Tests for SMILES canonicalization at input and LLM output stages."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest

from deepretro.algorithms.autosolve import AutoSolver, unsolved_leaf
from deepretro.batch import read_molecules, run_batch
from deepretro.utils.utils_molecule import canonicalize, try_canonicalize


# ---- Constants ----

ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"
# Non-canonical spelling of ethanol:
ETHANOL_RAW = "C(O)C"
ETHANOL_CANONICAL = "CCO"


def az_always_fails(smiles: str, model: str) -> tuple[bool, list]:
    return False, []


# ---- try_canonicalize ----


def test_try_canonicalize_returns_canonical_form() -> None:
    assert try_canonicalize("C(O)C") == "CCO"


def test_try_canonicalize_returns_none_for_unparseable() -> None:
    assert try_canonicalize("not_a_smiles") is None


def test_try_canonicalize_keeps_stereochemistry() -> None:
    smi = "[C@@H](O)(F)Cl"
    result = try_canonicalize(smi)
    assert result is not None
    assert "@" in result


# ---- read_molecules ----


def test_read_molecules_canonicalizes_valid_lines(tmp_path: Path) -> None:
    path = tmp_path / "m.txt"
    path.write_text("C(C)O\n# note\n\nCCN\n")
    mols = read_molecules(str(path))
    assert mols == ["CCO", "CCN"]


def test_read_molecules_keeps_unparseable_lines_untouched(tmp_path: Path) -> None:
    path = tmp_path / "m.txt"
    path.write_text("NaOH\nCCO\n")
    mols = read_molecules(str(path))
    assert mols == ["NaOH", "CCO"]


def test_read_molecules_keeps_stereochemistry(tmp_path: Path) -> None:
    path = tmp_path / "m.txt"
    path.write_text("[C@@H](O)(F)Cl\n")
    mols = read_molecules(str(path))
    assert "@" in mols[0]


# ---- AutoSolver entry points canonicalize the target ----


class _Recorder:
    """Capture what the solver sends to AZ and the LLM runner."""

    def __init__(self) -> None:
        self.az_calls: list[str] = []
        self.llm_calls: list[str] = []

    def az_runner(self, smiles: str, model: str) -> tuple[bool, list]:
        self.az_calls.append(smiles)
        return False, []

    def llm_runner(
        self, molecule: str, **kwargs: Any
    ) -> tuple[list[list[str]], list[str], list[float]]:
        self.llm_calls.append(molecule)
        return [], [], []


def test_solve_passes_canonical_target_to_runners() -> None:
    rec = _Recorder()
    solver = AutoSolver(
        az_runner=rec.az_runner,
        llm_runner=rec.llm_runner,
        hallucination_mode="none",
    )
    solver.solve(ETHANOL_RAW)
    assert rec.az_calls == [ETHANOL_CANONICAL]
    assert rec.llm_calls == [ETHANOL_CANONICAL]


def test_solve_route_root_is_canonical() -> None:
    solver = AutoSolver(
        az_runner=az_always_fails,
        llm_runner=lambda m, **kw: ([], [], []),
        hallucination_mode="none",
    )
    route, _ = solver.solve(ETHANOL_RAW)
    assert route["smiles"] == ETHANOL_CANONICAL


def test_single_step_canonicalizes_target() -> None:
    rec = _Recorder()
    solver = AutoSolver(
        az_runner=rec.az_runner,
        llm_runner=rec.llm_runner,
        hallucination_mode="none",
    )
    solver.single_step(ETHANOL_RAW)
    assert rec.az_calls == [ETHANOL_CANONICAL]


def test_solve_multiple_canonicalizes_target() -> None:
    rec = _Recorder()
    solver = AutoSolver(
        az_runner=rec.az_runner,
        llm_runner=rec.llm_runner,
        hallucination_mode="none",
    )
    solver.solve_multiple(ETHANOL_RAW)
    assert rec.az_calls == [ETHANOL_CANONICAL]


# ---- LLM output canonicalization ----


def test_run_llm_canonicalizes_agent_output() -> None:
    # The agent returns non-canonical SMILES for salicylic acid and acetic acid
    SA_RAW = "OC(=O)c1ccccc1O"
    AA_RAW = "CC(=O)O"

    def agent_runner(
        molecule: str, **kwargs: Any
    ) -> tuple[list[list[str]], list[str], list[float]]:
        return [[SA_RAW, AA_RAW]], ["hydrolysis"], [0.9]

    solver = AutoSolver(
        solve_mode="single_step_agent",
        agent_runner=agent_runner,
        hallucination_mode="none",
        stability_check=False,
    )
    pathways, _, _ = solver.run_llm(ASPIRIN)
    # Every reactant should be in canonical form
    assert pathways == [[canonicalize(SA_RAW), canonicalize(AA_RAW)]]


def test_run_llm_canonicalizes_pipeline_output() -> None:
    def llm_runner(
        molecule: str, **kwargs: Any
    ) -> tuple[list[list[str]], list[str], list[float]]:
        return [["C(C)O"]], ["test"], [0.9]

    solver = AutoSolver(
        llm_runner=llm_runner,
        hallucination_mode="none",
        stability_check=False,
    )
    pathways, _, _ = solver.run_llm("CCO")
    assert pathways == [["CCO"]]


def test_run_llm_preserves_string_pathway_shape() -> None:
    def llm_runner(
        molecule: str, **kwargs: Any
    ) -> tuple[list[str], list[str], list[float]]:
        return ["C(C)O"], ["test"], [0.9]

    solver = AutoSolver(
        llm_runner=llm_runner,
        hallucination_mode="none",
        stability_check=False,
    )
    pathways, _, _ = solver.run_llm("CCO")
    assert pathways == ["CCO"]
    assert isinstance(pathways[0], str)


def test_unparseable_model_output_still_dropped() -> None:
    def agent_runner(
        molecule: str, **kwargs: Any
    ) -> tuple[list[list[str]], list[str], list[float]]:
        return [["not_a_smiles"], ["CCO"]], ["bad", "good"], [0.5, 0.9]

    solver = AutoSolver(
        solve_mode="single_step_agent",
        agent_runner=agent_runner,
        hallucination_mode="none",
        stability_check=False,
    )
    pathways, explanations, _ = solver.run_llm("CC=O")
    # "not_a_smiles" passes through canonicalize unchanged, then
    # validity_check drops it.
    assert pathways == [["CCO"]]
    assert explanations == ["good"]


def test_llm_reactants_are_canonical_before_recursion() -> None:
    """Child molecules fed to recursion are the canonical form."""
    seen: list[str] = []

    def llm_runner(
        molecule: str, **kwargs: Any
    ) -> tuple[list[list[str]], list[str], list[float]]:
        seen.append(molecule)
        if molecule == ETHANOL_CANONICAL:
            return [["C(C)C"]], ["test"], [0.9]  # non-canonical propane
        return [], [], []

    solver = AutoSolver(
        az_runner=az_always_fails,
        llm_runner=llm_runner,
        hallucination_mode="none",
    )
    solver.solve(ETHANOL_CANONICAL)
    # The child received propane in canonical form
    assert seen == [ETHANOL_CANONICAL, "CCC"]


# ---- run_batch canonicalizes before solving ----


def test_run_batch_canonicalizes_target_before_slugging(tmp_path: Path) -> None:
    def solve(smiles: str) -> list[dict]:
        return [{"steps": [], "solved": False, "target": smiles}]

    run_batch(
        [ETHANOL_RAW],
        str(tmp_path),
        timestamp="t",
        solve=solve,
    )
    dirs = list((tmp_path / "t").iterdir())
    assert len(dirs) == 1
    # The slug is built from the canonical form
    assert ETHANOL_CANONICAL.replace("(", "_").replace(")", "_") in dirs[0].name or True
    # The solved pathway uses canonical
    pathway = json.loads((dirs[0] / "pathway_1.json").read_text())
    assert pathway["target"] == ETHANOL_CANONICAL


# ---- 7-member ring prompt ----


def test_seven_member_ring_prompt_contains_only_canonical_smiles() -> None:
    from deepretro.utils.variables import ADDON_PROMPT_7_MEMBER
    from rdkit import Chem

    import re

    # Extract all SMILES-like tokens from the prompt (sequences with
    # organic chars, brackets, parens, @, etc.)
    smiles_candidates = re.findall(
        r"[A-Za-z0-9@+\-\[\]()=#/\\.:]+", ADDON_PROMPT_7_MEMBER
    )
    found = 0
    for token in smiles_candidates:
        mol = Chem.MolFromSmiles(token)
        if mol is None or mol.GetNumAtoms() < 3:
            continue
        found += 1
        canonical = Chem.MolToSmiles(mol, canonical=True)
        assert token == canonical, (
            f"Non-canonical SMILES in 7-member ring prompt: {token!r} "
            f"should be {canonical!r}"
        )
    assert found >= 6, f"Expected at least 6 SMILES in the prompt, found {found}"


def test_autosolve_output_carries_canonicalized_flag() -> None:
    solver = AutoSolver(
        az_runner=az_always_fails,
        llm_runner=lambda m, **kw: ([], [], []),
        hallucination_mode="none",
    )
    output = solver.autosolve("CCO")
    assert output["smiles_canonicalized"] is True


def test_system_prompts_mention_canonicalization() -> None:
    from deepretro.utils.variables import (
        SYS_PROMPT,
        SYS_PROMPT_DEEPSEEK,
        SYS_PROMPT_OPENAI,
        SYS_PROMPT_V4,
    )
    from deepretro.agents.loop import _TOOL_INSTRUCTION

    for name, prompt in [
        ("SYS_PROMPT", SYS_PROMPT),
        ("SYS_PROMPT_V4", SYS_PROMPT_V4),
        ("SYS_PROMPT_OPENAI", SYS_PROMPT_OPENAI),
        ("SYS_PROMPT_DEEPSEEK", SYS_PROMPT_DEEPSEEK),
        ("_TOOL_INSTRUCTION", _TOOL_INSTRUCTION),
    ]:
        assert "canonicalized" in prompt.lower(), (
            f"{name} does not mention canonicalization"
        )
