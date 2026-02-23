deepretro Package
=================

``deepretro`` is a standalone Python package that provides algorithms, models,featurizers and utilities for retrosynthesis step modelling.

Installation
------------

.. code-block:: bash

   pip install -e "deepretro/[dev]"

Featurizers
-----------

.. autoclass:: deepretro.featurizers.reactionstep.ReactionStepFeaturizer
   :members:
   :undoc-members:
   :show-inheritance:

Utilities
---------

.. automodule:: deepretro.utils.domain_features
   :members:
   :undoc-members:

Running Tests
-------------

.. code-block:: bash

   pytest deepretro/tests/ -v

