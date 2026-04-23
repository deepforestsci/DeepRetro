deepretro.utils.parse_metrics
=============================

Reaction metric helpers for route formatting.

This module is used by ``deepretro.utils.parse`` when a precursor is attached
to a parsed reaction step. It is intentionally small: it predicts an optional
reaction class and maps that class to a scalability label. When the classifier
is not configured or cannot make a prediction, route formatting continues and
the scalability value is ``"N/A"``.

What It Produces
----------------

``ReactionMetricCalculator`` exposes two operations:

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Method
     - Behavior
   * - ``reaction_type(mol1, mol2)``
     - Returns ``(reaction_name, reaction_index)`` for a reactant/product pair.
       If the model file is missing, fingerprints cannot be computed, or
       prediction fails, it returns ``("Unknown Reaction", -1)``.
   * - ``scalability_index(mol1, mol2)``
     - Returns the configured scalability label for the predicted reaction
       index. If there is no configured model path or prediction fails, it
       returns ``"N/A"``.

Confidence estimates are not part of this module and are not emitted by the
package route parser.

Model Configuration
-------------------

By default, ``ReactionMetricCalculator()`` resolves the classifier path from
``RXN_CLASSIFICATION_MODEL_PATH``.

- If the environment variable is unset or empty, scalability is disabled and
  ``scalability_index`` returns ``"N/A"``.
- If the value is absolute, it is used directly.
- If the value is relative, it is resolved from the repository root. The root
  is identified by ``config/langfuse_config.json``.

Use an explicit path when the caller owns model configuration:

.. code-block:: python

   from deepretro.utils.parse_metrics import ReactionMetricCalculator

   calculator = ReactionMetricCalculator(
       model_path="model_out/model.joblib",
   )
   scalability = calculator.scalability_index("CC", "CCO")

Use an empty path when tests or lightweight callers should avoid loading a
classifier:

.. code-block:: python

   calculator = ReactionMetricCalculator(model_path="")
   assert calculator.scalability_index("CC", "CCO") == "N/A"

Testing and Dependency Injection
--------------------------------

The calculator accepts injected dependencies so tests do not need to read a
joblib model or compute RDKit fingerprints:

.. code-block:: python

   class FakeClassifier:
       def predict(self, fingerprints):
           return [0]

   calculator = ReactionMetricCalculator(
       model_path="unused-in-test.joblib",
       model_loader=lambda path: FakeClassifier(),
       fingerprint_calculator=lambda smiles: [1, 0],
       reaction_encoding_names={0: "Example reaction"},
       scalability_encoding={0: "medium"},
   )

   assert calculator.reaction_type("CC", "CCO") == ("Example reaction", 0)
   assert calculator.scalability_index("CC", "CCO") == "medium"

Logging and Failure Behavior
----------------------------

Recoverable model and fingerprint failures do not raise out of the public
methods. They are logged with the module-level ``structlog`` logger using
``logger.error`` and structured fields:

.. code-block:: text

   event="Error in metric calculation"
   function="get_reaction_type"
   error="<exception message>"

Pass a custom logger to ``ReactionMetricCalculator(logger=...)`` when tests or
applications need to capture those recoverable errors.

Compatibility Functions
-----------------------

The module-level functions remain for existing callers:

- ``get_reaction_type(mol1, mol2, model_path)`` delegates to
  ``ReactionMetricCalculator(model_path=model_path).reaction_type(...)``.
- ``calc_scalability_index(mol1, mol2)`` delegates to
  ``ReactionMetricCalculator().scalability_index(...)`` and therefore uses
  ``RXN_CLASSIFICATION_MODEL_PATH``.

API Reference
-------------

.. automodule:: deepretro.utils.parse_metrics
   :members:
   :undoc-members:
