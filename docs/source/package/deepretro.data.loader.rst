deepretro.data.loader
=====================

Dataset loading pipeline for reaction-step data, built on DeepChem's
:class:`~deepchem.data.DataLoader` base class.

``ReactionDataLoader`` reads a CSV with product SMILES, reactant SMILES,
and a binary label column, featurizes each row with a reaction featurizer
(by default :class:`~deepretro.featurizers.ReactionStepFeaturizer`), and
writes the result to a :class:`~deepchem.data.DiskDataset` with automatic
sharding for memory-efficient handling of large files.

A convenience function ``stratified_split`` wraps DeepChem's
``SingletaskStratifiedSplitter`` to split any ``Dataset`` into
train / valid / test sets while preserving class balance.

Usage
-----

.. code-block:: python

   from deepretro.data import ReactionDataLoader, stratified_split

   # Load and featurize a reaction CSV into a DiskDataset
   loader = ReactionDataLoader(
       product_col="product",
       reactants_col="reactants",
       label_col="label",
   )
   dataset = loader.create_dataset("data/hallucination_dataset.csv", shard_size=1000)

   # Stratified train/valid/test split (70/15/15)
   train, valid, test = stratified_split(dataset)

API
---

.. automodule:: deepretro.data.loader
   :members:
   :undoc-members:

