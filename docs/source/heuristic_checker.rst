Heuristic reaction checker
==========================

The heuristic checker compares proposed reactants with a target product and
assigns a structural consistency score from 0 to 100. Higher scores mean fewer
penalties. The score is not a probability of chemical correctness: reaction
conditions, mechanisms, stereoselectivity and experimental feasibility are
outside its scope.

Scoring a reaction
------------------

The API takes reactants first and product second::

   from deepretro.algorithms import calculate_hallucination_score

   report = calculate_hallucination_score("CCO", "CC=O")
   print(report["score"], report["severity"])
   print(report["penalties"])

The result contains ``score``, ``severity``, ``penalties``, ``penalty_total`` and
``message``. ``penalty_total`` records the sum before the score is clamped at zero.
Invalid or empty structures, unchanged reactions and carbon-radical inputs return
``unassessable=True`` with score zero. Radical chemistry is an unsupported case;
its exclusion does not establish that the reaction is chemically impossible.
Atom-map labels are ignored; isotope and stereo differences remain meaningful
when detecting unchanged reactions. Isotopes are normalized for skeleton scoring.

Candidate ranking and warnings
------------------------------

The pipeline ranks assessable candidates by descending score. Scores below
``reject_below`` remain fallback candidates, so a successful solver route can
contain flagged steps. Retention is not a clean verdict::

   from deepretro.algorithms.pipeline_checks import hallucination_checker

   status, candidates = hallucination_checker("CC=O", [["CCO"], ["CC=O"]])
   assert status == 200
   assert candidates == [["CCO"]]  # the unchanged reaction is excluded

Pass ``rank_only=False`` for hard filtering. The wrapper's
``check_single_pathway(product, reactants)`` returns 1 for a flag or an unsupported
input, and 0 otherwise; it uses the threshold rather than candidate retention.
The LLM tool uses the same individual verdict.

Enable heuristic checking in the solver with
``AutoSolver(hallucination_mode="heuristic")``. Use
``hallucination_mode="none"`` to disable it. Exported routes include
``hallucination_summary`` with the number of annotated steps, flagged steps and
minimum score. This summary covers recorded solver verdicts; it is not a complete
chemical audit of all AiZynthFinder or ML-generated steps.

Configuration
-------------

Weights are immutable and validated. Each mismatch family has its own weight and
cap. Severity display cutoffs (``cut_high``, ``cut_medium``, ``cut_low``) are
independent of the ``reject_below`` threshold::

   from deepretro.algorithms import DEFAULT_WEIGHTS, HallucinationWeights

   weights = DEFAULT_WEIGHTS.replace(reject_below=50, w_ring=20)
   weights.to_json("weights.json")
   restored = HallucinationWeights.from_json("weights.json")
   assert restored == weights

Pass ``weights`` to the scoring function or wrapper, or
``hallucination_weights=weights`` to ``AutoSolver``. The batch CLI accepts
``--hallucination-weights weights.json``. Heuristic and disabled modes do not
require a training sheet::

   python -m deepretro.batch --molecules molecules.txt \
       --hallucination-mode heuristic --hallucination-weights weights.json

JSON files must contain the weight fields; unknown and missing fields are
rejected. Legacy files without ``reject_below`` inherit their own ``cut_medium``.
Defaults are heuristic policy choices, not fitted chemical probabilities.

Behavior and limitations
------------------------

V2 sums per-element atom penalties, counts ring-size multiplicity differences,
compares aromatic substitution patterns and ranks rather than discards flagged
fallbacks. Comparison and scoring are separate functions, allowing a cached
comparison to be rescored without repeating RDKit operations.

Spectator detection uses bounded maximum-common-substructure searches. A timed
out comparison retains the fragment conservatively. It removes candidate
spectators only when the remaining heavy-atom count can cover the product; this
is not atom mapping or proof that every required atom has a valid source.
Ring matching and relative substituent positions are heuristics, especially for
multiple similar rings and complex ring systems. Legitimate ring formation,
protection, redox chemistry or omitted reagents can receive penalties. Review
flagged steps in their reaction context, and do not interpret an unflagged step
as validated chemistry.
