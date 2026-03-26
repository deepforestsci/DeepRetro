deepretro.models.hallucination\_classifier
==========================================

XGBoost-based binary classifier for detecting hallucinated retrosynthesis
reactions.  Built on DeepChem's :class:`~deepchem.models.GBDTModel`, which
wraps an ``XGBClassifier`` and adds automatic early-stopping via an
internal 80/20 train/validation split.

Training a new model
--------------------

Prepare a CSV with columns ``product``, ``reactants``, and ``label``
(1 = hallucinated, 0 = valid).  Then:

.. code-block:: python

   from deepretro.data import ReactionDataLoader, stratified_split
   from deepretro.models import HallucinationClassifier

   # Load and featurize
   loader = ReactionDataLoader()
   dataset = loader.create_dataset("data/hallucination_dataset.csv")
   train, valid, test = stratified_split(dataset)

   # Train
   clf = HallucinationClassifier(model_dir="my_models/")
   clf.fit(train)

   # Evaluate (also sets the optimal probability threshold)
   scores = clf.evaluate(test)
   print(scores)

Saving and loading
------------------

The model is auto-saved to ``model_dir`` after training.  To reload:

.. code-block:: python

   clf = HallucinationClassifier(model_dir="my_models/")
   clf.load("my_models/")

The saved artifacts include the XGBoost model weights and the optimal
classification threshold.

Configuration
-------------

No environment variables are required.  All paths are passed as
arguments to the constructor and ``load()`` / ``save()`` methods.

API
---

.. automodule:: deepretro.models.hallucination_classifier
   :members:
   :undoc-members:
