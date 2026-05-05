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
- Heuristic hallucination detection and scoring for retrosynthetic steps.
- ML-based hallucination classification (XGBoost via DeepChem ``GBDTModel``).
- Dataset loading with DeepChem ``DiskDataset`` sharding and stratified splitting.


Input and Output Conventions
----------------------------

Reaction steps are represented as:

.. code-block:: python

   (product_smiles, reactants_smiles)

where ``reactants_smiles`` may contain multiple molecules separated by ``.``.


Quickstart
----------

.. code-block:: python

   from deepretro import ReactionStepFeaturizer

   featurizer = ReactionStepFeaturizer(radius=3, size=2048, use_domain_features=True)
   X = featurizer.featurize([
       ("CCO", "CC.O"),
       ("c1ccccc1", "c1ccccc1.Cl"),
   ])
   print(X.shape)  # (2, 4135)


Top-Level API
-------------

.. automodule:: deepretro
   :members:
   :undoc-members:

Subpackages
-----------

.. toctree::
   :maxdepth: 2

   deepretro.algorithms.hallucination_checker
   deepretro.data.loader
   deepretro.models.hallucination_classifier
   deepretro.algorithms.stability_checker
   deepretro.featurizers
   deepretro.utils_pkg
