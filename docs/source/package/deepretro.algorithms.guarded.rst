Guarded heuristic search
========================

``GuardedSolver`` wraps a heuristic pipeline ``AutoSolver``. It checks every
reaction from both AiZynthFinder and the proposal pipeline before expanding
precursors. Unlike ordinary search, it does not retain heuristic-flagged
fallbacks. Rejections are fed back at the failed molecule; accepted ancestors
and completed sibling branches remain in the current tree.

Usage
-----

This example runs offline with injected proposals and a small example stock:

.. code-block:: python

   from pathlib import Path
   from tempfile import TemporaryDirectory
   from deepretro import AutoSolver
   from deepretro.algorithms.guarded import GuardedSolver

   backend = AutoSolver(
       hallucination_mode="heuristic",
       az_runner=lambda smiles, model: (False, []),
       llm_runner=lambda molecule, **kwargs: ([["CCO"]], ["oxidation"], [0.9]),
   )
   with TemporaryDirectory() as directory:
       solver = GuardedSolver(
           backend, {"CCO"}.__contains__, Path(directory) / "run", rounds=3,
       )
       tree, passed = solver.solve("CC=O")
       parsed = backend.parse(tree, solved=passed)

The required ``stock_checker`` receives canonical isomeric SMILES. It must
return a Python boolean based on actual stock membership. Template ``in_stock``
flags and the ordinary small-molecule shortcut do not satisfy this check.
Catalog hits are labelled ``stock_source="catalog"``; they are not attributed
to template search. Consequently the parser's AZ-only completion flags may be
false for a complete guarded route; use the returned ``passed`` value for
guarded completion.

The output directory must be new. It receives ``repair_events.jsonl`` and
numbered root-tree checkpoints, including incomplete searches. Reusing a
solver for another target resets search caches and continues checkpoint
numbering. Checkpoints are evidence snapshots, not restart files.

``rounds`` bounds proposal calls per molecule; each built-in feedback proposal
uses one transport attempt with nested retries disabled. Custom runners must
honor the forwarded ``messages`` and ``max_attempts`` arguments. ``max_depth``
bounds reaction edges. An optional absolute ``time.monotonic()`` deadline is
cooperative: it does not interrupt an in-flight call. Provider request timeouts
still apply. AZ operational errors raise rather than becoming chemical
rejections. An empty proposal response can also reflect provider or parsing
failure; inspect the normal LLM logs before interpreting an incomplete route.

Passing the heuristic is not a chemical-feasibility judgment. These checks
do not replace expert review or experimental validation. Existing ordinary
``AutoSolver`` callers keep their current behavior.

Disk-backed stock
-----------------

``SortedInchiKeyQuery`` performs exact binary-search membership against a
read-only memory-mapped ``S27`` NumPy array. The stock producer must verify
sorting, key validity and unique cardinality. Opening the query validates
array metadata but does not scan a large stock or repair its ordering.

.. code-block:: python

   from types import SimpleNamespace
   from rdkit import Chem
   from deepretro.utils.disk_stock import SortedInchiKeyQuery

   with SortedInchiKeyQuery("stock.npy", length=stock_cardinality) as stock:
       def in_stock(smiles: str) -> bool:
           molecule = Chem.MolFromSmiles(smiles)
           return SimpleNamespace(inchi_key=Chem.MolToInchiKey(molecule)) in stock

       # Keep the mapping open for the entire solve.
       solver = GuardedSolver(backend, in_stock, Path("new-run"))
       tree, passed = solver.solve(target_smiles)

AiZynthFinder installations can also reference the query in their stock
configuration as ``deepretro.utils.disk_stock.SortedInchiKeyQuery``.

API
---

.. automodule:: deepretro.algorithms.guarded
   :members:

.. automodule:: deepretro.utils.disk_stock
   :members:
