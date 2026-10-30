"""Guarded retrosynthesis with local, evidence-preserving reaction repair."""

from __future__ import annotations

import copy
import json
import time
from pathlib import Path
from collections.abc import Callable
from typing import Any

from rdkit import Chem


def _canonical(smiles: str) -> str:
    """Return canonical isomeric identity or raise for unspecified atoms."""
    if not isinstance(smiles, str) or not smiles.strip():
        raise ValueError("Empty SMILES")
    mol = Chem.MolFromSmiles(smiles)
    if (
        mol is None
        or not mol.GetNumAtoms()
        or any(a.GetAtomicNum() == 0 for a in mol.GetAtoms())
    ):
        raise ValueError("Invalid SMILES or wildcard atoms")
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    return Chem.MolToSmiles(mol, isomericSmiles=True)


def _leaf(smiles: str) -> dict[str, Any]:
    """Return a fail-closed unsolved molecule node."""
    return dict(
        type="mol", smiles=smiles, is_chemical=True, in_stock=False, children=[]
    )


class GuardedSolver:
    """Check every edge and repair failures at their exact product.

    Passing the heuristic is not evidence of chemical feasibility. Search
    budgets are cooperative and do not interrupt an in-flight provider call.

    Parameters
    ----------
    backend : object
        Pipeline AutoSolver configured with the heuristic checker. Every AZ
        and LLM reaction must pass that checker; flagged fallbacks are rejected.
    stock_checker : callable
        Exact stock membership callback, receiving canonical isomeric SMILES.
        This is required even when an AZ leaf claims to be in stock.
    output_dir : Path
        New event and checkpoint directory. Must not already exist.
    max_depth : int
        Maximum number of reaction edges along a branch.
    rounds : int
        Maximum LLM revision rounds per molecule.
    deadline : float, optional
        Absolute monotonic-clock deadline.

    Examples
    --------
    This example uses a small in-memory stock. It makes no model calls.

    >>> from tempfile import TemporaryDirectory
    >>> from deepretro import AutoSolver
    >>> backend = AutoSolver(
    ...     az_runner=lambda smiles, model: (False, []),
    ...     llm_runner=lambda molecule, **kwargs: ([], [], []),
    ... )
    >>> with TemporaryDirectory() as directory:
    ...     solver = GuardedSolver(
    ...         backend, {"CCO"}.__contains__, Path(directory) / "run"
    ...     )
    ...     tree, passed = solver.solve("CCO")
    >>> passed, tree["in_stock"]
    (True, True)
    """

    def __init__(
        self,
        backend: Any,
        stock_checker: Callable[[str], bool],
        output_dir: Path,
        max_depth: int = 50,
        rounds: int = 3,
        deadline: float | None = None,
    ) -> None:
        """Validate the backend and create a new run directory."""
        if (
            backend.hallucination_mode != "heuristic"
            or backend.solve_mode != "pipeline"
        ):
            raise ValueError("GuardedSolver requires a heuristic pipeline backend")
        if backend.hallucination_checker is None:
            raise ValueError("A heuristic checker is required")
        if (
            isinstance(max_depth, bool)
            or not isinstance(max_depth, int)
            or max_depth < 0
        ):
            raise ValueError("max_depth must be a nonnegative integer")
        if isinstance(rounds, bool) or not isinstance(rounds, int) or rounds < 0:
            raise ValueError("rounds must be a nonnegative integer")
        if not callable(stock_checker):
            raise TypeError("stock_checker must be callable")
        self.backend = backend
        self.stock_checker = stock_checker
        self.output_dir = Path(output_dir)
        # Each solver owns a fresh directory, so events/checkpoints cannot be
        # confused with a previous run or another writer.
        self.output_dir.mkdir(parents=True, exist_ok=False)
        self.max_depth = max_depth
        self.rounds = rounds
        self.deadline = deadline
        self.az_cache: dict[str, list[dict[str, Any]]] = {}
        self.solved_memo: dict[str, dict[str, Any]] = {}
        self.rejections: dict[str, list[dict[str, Any]]] = {}
        self.rejected: set[tuple[str, tuple[str, ...]]] = set()
        self.sequence = 0
        self.root: dict[str, Any] | None = None

    def _event(self, kind: str, **details: Any) -> None:
        """Append an event without overwriting prior evidence."""
        with (self.output_dir / "repair_events.jsonl").open("a") as handle:
            handle.write(
                json.dumps(dict(event=kind, time=time.time(), **details), default=str)
                + "\n"
            )

    def _checkpoint(self) -> None:
        """Preserve the currently connected root including completed siblings."""
        self.sequence += 1
        path = self.output_dir / f"checkpoint_{self.sequence:06d}.json"
        with path.open("x") as handle:
            json.dump(self.root, handle, indent=2)

    def _stopped(self, depth: int) -> bool:
        """Check depth and wall-clock budgets before external calls."""
        return depth > self.max_depth or (
            self.deadline is not None and time.monotonic() >= self.deadline
        )

    def _reject(
        self,
        product: str,
        reactants: list[str],
        reason: str,
        source: str,
        *,
        branch_cycle: bool = False,
    ) -> None:
        """Record exact repair location and feed the reason to subsequent proposals."""
        if not branch_cycle:
            self.rejected.add((product, tuple(sorted(reactants))))
        item: dict[str, Any] = dict(
            product=product, reactants=reactants, reason=reason, source=source
        )
        if branch_cycle:
            item["branch_cycle"] = True
        self.rejections.setdefault(product, []).append(item)
        self._event("reaction_rejected", **item)
        self._checkpoint()

    def _memo_safe(self, node: dict[str, Any], ancestors: set[str], depth: int) -> bool:
        """Ensure reusing a solved subtree cannot cross an ancestor or depth bound."""
        if self._stopped(depth) or _canonical(node["smiles"]) in ancestors:
            return False
        lineage = ancestors | {_canonical(node["smiles"])}
        return all(
            self._memo_safe(child, lineage, depth + 1)
            for reaction in node.get("children", [])
            for child in reaction.get("children", [])
        )

    def solve(
        self, smiles: str, visited: set[str] | None = None, depth: int = 0
    ) -> tuple[dict[str, Any], bool]:
        """Return a fully guarded route or a checkpointed partial route.

        Parameters
        ----------
        smiles : str
            Target molecular structure.
        visited : set of str, optional
            Ancestors on the current branch, never a global explored set.
        depth : int
            Current reaction depth.

        Returns
        -------
        tuple
            Raw molecule tree and heuristic-passing, stocked-route status.
            Each call starts a fresh search; checkpoint numbering continues.

        Raises
        ------
        ValueError
            If the molecular identity or depth is invalid.
        RuntimeError
            If the AZ backend reports an operational error.

        See Also
        --------
        GuardedSolver : Construction and a complete offline example.
        """
        if isinstance(depth, bool) or not isinstance(depth, int) or depth < 0:
            raise ValueError("depth must be a nonnegative integer")
        target = _canonical(smiles)
        self.az_cache.clear()
        self.solved_memo.clear()
        self.rejections.clear()
        self.rejected.clear()
        node = _leaf(target)
        self.root = node
        self._checkpoint()
        solved = self._solve_into(node, {_canonical(s) for s in (visited or ())}, depth)
        self._checkpoint()
        return node, solved

    def _solve_into(
        self,
        node: dict[str, Any],
        ancestors: set[str],
        depth: int,
        az_hint: dict[str, Any] | None = None,
    ) -> bool:
        """Solve a mutable attached node so local repair retains its ancestors."""
        product = _canonical(node["smiles"])
        if self._stopped(depth) or product in ancestors:
            self._event(
                "branch_blocked", product=product, depth=depth, reason="cycle_or_budget"
            )
            return False
        if self.stock_checker(product) is True:
            node.update(_leaf(product), in_stock=True, stock_source="catalog")
            return True
        memo = self.solved_memo.get(product)
        if memo is not None and self._memo_safe(memo, ancestors, depth):
            node.clear()
            node.update(copy.deepcopy(memo))
            return True
        candidates: list[dict[str, Any]] = []
        if az_hint is not None:
            candidates.append(az_hint)
        if product not in self.az_cache:
            if az_hint is None:
                az_solved, routes, az_errored = self.backend._run_az(product)
                if az_errored:
                    self._event("backend_error", product=product, source="az")
                    raise RuntimeError(
                        "AZ backend failed; search is operationally incomplete"
                    )
                self.az_cache[product] = [copy.deepcopy(r) for r in routes[:3]]
                self._event(
                    "az_candidates",
                    product=product,
                    reported_solved=az_solved,
                    errored=az_errored,
                    routes=self.az_cache[product],
                )
            else:
                # The inherited route is already the first AZ search result for this node.
                self.az_cache[product] = [copy.deepcopy(az_hint)]
        if az_hint is None:
            candidates.extend(self.az_cache[product])
        first_partial: dict[str, Any] | None = None
        for candidate in candidates:
            if self._stopped(depth):
                break
            if self._try_az(node, candidate, ancestors, depth):
                self.solved_memo[product] = copy.deepcopy(node)
                self._checkpoint()
                return True
            if node.get("children") and first_partial is None:
                first_partial = copy.deepcopy(node)
        for revision in range(self.rounds):
            if self._stopped(depth) or depth >= self.max_depth:
                break
            # A reaction that cycles on one branch can be valid on another.
            feedback = [
                item
                for item in self.rejections.get(product, [])
                if not item.get("branch_cycle")
                or set(item["reactants"]) & (ancestors | {product})
            ]
            self._event(
                "generation",
                product=product,
                revision=revision,
                rejection_feedback=feedback,
            )
            paths, _explanations, _confidence = self.backend.run_llm(
                product,
                depth=depth,
                rejection_feedback=json.dumps(feedback),
            )
            for path in paths:
                if self._stopped(depth):
                    break
                raw = path.split(".") if isinstance(path, str) else path
                try:
                    if not isinstance(raw, (list, tuple)) or not raw:
                        raise ValueError("Empty or malformed reactant list")
                    reactants = [_canonical(s) for s in raw]
                except (ValueError, TypeError) as exc:
                    self._reject(product, [], str(exc), "llm")
                    continue
                if self._try_reaction(node, reactants, "llm", ancestors, depth):
                    self.solved_memo[product] = copy.deepcopy(node)
                    self._checkpoint()
                    return True
                if node.get("children") and first_partial is None:
                    first_partial = copy.deepcopy(node)
        if first_partial is not None:
            node.clear()
            node.update(first_partial)
        self._checkpoint()
        return False

    def _try_az(
        self,
        node: dict[str, Any],
        candidate: dict[str, Any],
        ancestors: set[str],
        depth: int,
    ) -> bool:
        """Review an inherited AZ route without trusting its solved flag."""
        product = node["smiles"]
        try:
            if not isinstance(candidate, dict) or candidate.get("type", "mol") != "mol":
                raise ValueError("AZ candidate must be a molecule node")
            if _canonical(candidate["smiles"]) != product:
                raise ValueError("AZ route root does not match requested product")
            reactions = candidate.get("children", [])
            if not reactions:
                # Exact membership was checked before reading the AZ hint.
                return False
            if not isinstance(reactions, list) or len(reactions) != 1:
                raise ValueError("Molecule must have exactly one reaction child")
            if (
                not isinstance(reactions[0], dict)
                or reactions[0].get("type", "reaction") != "reaction"
            ):
                raise ValueError("AZ molecule child must be a reaction node")
            hints = reactions[0].get("children", [])
            if not isinstance(hints, list) or not all(
                isinstance(child, dict) and child.get("type", "mol") == "mol"
                for child in hints
            ):
                raise ValueError("Reaction children must be molecule nodes")
            reactants = [_canonical(child["smiles"]) for child in hints]
        except (ValueError, KeyError, TypeError) as exc:
            self._reject(product, [], str(exc), "az")
            return False
        return self._try_reaction(node, reactants, "az", ancestors, depth, hints=hints)

    def _try_reaction(
        self,
        node: dict[str, Any],
        reactants: list[str],
        source: str,
        ancestors: set[str],
        depth: int,
        hints: list[dict[str, Any]] | None = None,
    ) -> bool:
        """Gate an edge before recursively expanding any of its precursors."""
        product = node["smiles"]
        signature = (product, tuple(sorted(reactants)))
        if signature in self.rejected:
            return False
        if not reactants or any(s in ancestors | {product} for s in reactants):
            self._reject(
                product,
                reactants,
                "Empty reaction or ancestor cycle",
                source,
                branch_cycle=bool(reactants) and product not in reactants,
            )
            return False
        if self._stopped(depth) or depth >= self.max_depth:
            return False
        review = self.backend.hallucination_checker.assess(product, reactants)
        if not isinstance(review, dict) or review.get("source") != "heuristic":
            raise RuntimeError("Heuristic checker returned a malformed assessment")
        if review.get("flagged") is not False:
            self._reject(
                product,
                reactants,
                str(
                    review.get("message")
                    or review.get("explanation")
                    or "Heuristic check did not pass"
                ),
                source,
            )
            return False
        if self._stopped(depth):
            return False
        children = [_leaf(s) for s in reactants]
        node.update(
            in_stock=False,
            children=[
                dict(
                    type="reaction",
                    is_reaction=True,
                    solved_by=source,
                    metadata=dict(hallucination=review, guarded=True),
                    children=children,
                )
            ],
        )
        self._checkpoint()
        solved = True
        for index, child in enumerate(children):
            child_solved = self._solve_into(
                child,
                ancestors | {product},
                depth + 1,
                az_hint=hints[index] if hints is not None else None,
            )
            solved = child_solved and solved
            self._checkpoint()
        return solved
