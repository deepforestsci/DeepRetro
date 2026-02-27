DeepRetro
=========

``deepretro`` is a chemistry ML utility package for retrosynthesis workflows.
It focuses on robust reaction-step featurization and practical integrations
that can be dropped into DeepChem training pipelines or custom research code.

Overview
--------

The package currently provides:

- Reaction-step vectorization using product/reactant fingerprints plus
  handcrafted chemistry descriptors.
- Domain-feature extraction helpers for product/reactant SMILES pairs.
- AiZynthFinder wrappers for template-based route search.

Primary Use Cases
-----------------

- Train scoring/ranking models on reaction-step data.
- Build hybrid systems that mix template search with ML/LLM decision layers.
- Generate reproducible numerical features from retrosynthesis candidates.

Input and Output Conventions
----------------------------

Reaction steps are represented as:

.. code-block:: python

   (product_smiles, reactants_smiles)

where ``reactants_smiles`` may contain multiple molecules separated by ``.``.

For featurization, output dimensionality is:

.. math::

   \text{feature\_dim} = 2 * \text{size} + (15 \text{ if use\_domain\_features else } 0)

Quickstart
----------

.. code-block:: python

   from deepretro import ReactionStepFeaturizer

   featurizer = ReactionStepFeaturizer(radius=2, size=2048, use_domain_features=True)
   X = featurizer.featurize([
       ("CCO", "CC.O"),
       ("c1ccccc1", "c1ccccc1.Cl"),
   ])
   print(X.shape)  # (2, 4111)

Reliability Contract
--------------------

- Deterministic vectors for identical inputs.
- Invalid SMILES produce all-``NaN`` rows rather than silent coercion.
- Domain feature order is fixed and documented in :doc:`deepretro.utils_pkg`.

Top-Level API
-------------

.. automodule:: deepretro
   :members:
   :undoc-members:

Subpackages
-----------

.. toctree::
   :maxdepth: 2

   deepretro.featurizers
   deepretro.utils_pkg
