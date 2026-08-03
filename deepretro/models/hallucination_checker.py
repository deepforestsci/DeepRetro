"""Inference-time hallucination checking for retrosynthetic steps.

This module exposes :class:`HallucinationChecker`, the *inference* entry point
that decides whether a proposed ``product <- reactants`` step is chemically
plausible or a model hallucination. It supports two interchangeable backends:

* ``"heuristic"`` -- deterministic rule-based checks from
  :mod:`deepretro.algorithms.pipeline_checks`; needs no trained model.
* ``"ml"`` -- a trained DeepChem/XGBoost classifier loaded from a directory
  produced by :class:`~deepretro.models.hallucination_trainer.HallucinationTrainer`.

For the offline *training* side of this pipeline see
:mod:`deepretro.models.hallucination_trainer`.
"""

import os
import json
import numpy as np
import deepchem as dc
from typing import Tuple, List
from deepretro.models.hallucination_utils import create_model_instance
from deepretro.featurizers.reactionstep import FEATURIZER_MAP
from deepretro.algorithms.pipeline_checks import (
    hallucination_checker as heuristic_checker,
)


class HallucinationChecker:
    """
    Inference infrastructure for detecting pathway hallucinations using either rule-based heuristics or unified
    machine learning models.

    Examples
    --------
    The heuristic backend needs no trained model and classifies a single step,
    returning ``0`` for a plausible step and ``1`` for a hallucination:

    >>> from deepretro.models.hallucination_checker import HallucinationChecker
    >>> checker = HallucinationChecker(checker_type="heuristic")
    >>> checker.check_single_pathway(
    ...     "Cn1nccc1[C@]1(O)CCCC[C@H]1O", "Cn1nccc1Br.O=C1CCCC[C@H]1O"
    ... )
    0

    The ML backend loads a directory previously written by
    :class:`~deepretro.models.hallucination_trainer.HallucinationTrainer`::

        checker = HallucinationChecker(checker_type="ml", model_path="my_model")
        status, valid_pathways = checker(product, candidate_pathways)
    """

    def __init__(
        self,
        checker_type: str = "ml",
        model_path: str | None = None,
    ) -> None:
        """
        Initializes a hallucination checker.

        Parameters
        ----------
        checker_type : str
            Backend used for hallucination detection. Supported values are ``"heuristic"`` and ``"ml"``, default="ml"
        model_path : str, optional
            Path to a saved model directory. Required when ``checker_type="ml"``.
        """
        self.checker_type: str = checker_type.lower()
        self.model_path: str | None = model_path
        self.model: dc.models.Model | None = None
        self.featurizer: dc.feat.Featurizer | None = None
        self.threshold: float = 0.5
        self.model_type: str | None = None

        if self.checker_type not in ["heuristic", "ml"]:
            raise ValueError(
                f"Unknown checker_type: {self.checker_type}. Allowed: 'heuristic', 'ml'"
            )

        if self.checker_type == "ml":
            if self.model_path is None:
                raise ValueError(
                    "A valid 'model_path' must be provided when checker_type is 'ml'."
                )
            self.load_model(self.model_path)

    def load_model(self, model_path: str) -> None:
        """
        Loads the saved trainer configuration, instantiates the required underlying featurizers, and restores
        model weights from disk.

        Parameters
        ----------
        model_path : str
            Path to the saved model directory.
        """
        config_path = os.path.join(model_path, "config.json")
        if not os.path.exists(config_path):
            raise FileNotFoundError(
                f"No config.json file found in model path: {model_path}"
            )

        with open(config_path, "r") as f:
            config = json.load(f)

        # Extract metadata directly from the saved JSON configuration
        self.model_type = config.get("model_type")
        model_params = config.get("model_params", {})
        feat_name = config.get("feat_name")
        feat_params = config.get("feat_params", {})
        self.threshold = config.get("threshold", 0.5)
        n_tasks = config.get("n_tasks", 1)

        print(f"Restoring saved {self.model_type.upper()} pipeline from disk...")
        print(f"Loaded Decision Threshold: {self.threshold:.4f}")

        # 1. Dynamically initialize the custom operational featurizer wrapper
        if feat_name.lower() not in FEATURIZER_MAP:
            raise ValueError(
                f"Saved featurizer '{feat_name}' is not registered in FEATURIZER_MAP."
            )

        feat_config = FEATURIZER_MAP[feat_name.lower()]
        merged_feat_params = feat_config["default_params"].copy()
        merged_feat_params.update(feat_params)
        self.featurizer = feat_config["class"](**merged_feat_params)

        # 2. Package parameters and construct the model instance using our standalone factory
        params = model_params.copy()
        params["model_dir"] = model_path
        params["n_tasks"] = n_tasks

        self.model = create_model_instance(self.model_type, **params)
        try:
            self.model.restore()
        except Exception:
            self.model.reload()

    def check_single_pathway(self, target: str, reactants: str) -> int:
        """
        Public inference API route. Accepts a single pathway definition step
        and routes it to the designated heuristic or ML backend processor.

        Parameters
        ----------
        target : str
            SMILES representation of the target molecule (product).
        reactants : str
            SMILES representation of the reactant molecules (dot-separated).

        Returns
        -------
        int
            Binary prediction: 1 if it's a hallucination, 0 if it's a valid step.
        """
        if self.checker_type == "heuristic":
            _, valid_pathways = self._check_pathway_heur(target, [reactants])
            return 0 if len(valid_pathways) > 0 else 1
        else:
            return self._check_pathway_ml(target, reactants)

    def __call__(self, target: str, pathways: List[str] | List[List[str]]):
        """Filter candidate pathways -- a convenience alias for :meth:`check_pathways`.

        Parameters
        ----------
        target : str
            Product molecule SMILES.
        pathways : List[str] or List[List[str]]
            Candidate pathways as reactant SMILES strings or lists thereof.

        Returns
        -------
        Tuple[int, List]
            Status code and the subset of pathways classified as valid.
        """
        return self.check_pathways(target, pathways)

    def check_pathways(
        self,
        target: str,
        pathways: List[str] | List[List[str]],
    ) -> Tuple[int, List]:
        """
        Filters a collection of candidate pathways.

        Parameters
        ----------
        target: str
            Product molecule SMILES.
        pathways: List[str]
            Candidate pathways represented as reactant SMILES strings or lists of reactant SMILES.

        Returns
        -------
        Tuple[int, List]
            Status code and the subset of pathways classified as valid.
        """
        if self.checker_type == "heuristic":
            return self._check_pathway_heur(target, pathways)
        else:
            from deepretro.utils.utils_molecule import is_valid_smiles

            valid_pathways = []
            for pathway in pathways:
                if isinstance(pathway, list):
                    reactants_smi = ".".join(pathway)
                else:
                    reactants_smi = pathway

                if not is_valid_smiles(reactants_smi):
                    continue

                pred = self._check_pathway_ml(target, reactants_smi)
                if pred == 0:
                    valid_pathways.append(pathway)

            return 200, valid_pathways

    def _check_pathway_heur(
        self,
        target: str,
        reactants: List[str] | str,
    ) -> Tuple[int, List]:
        """Evaluates pathways using deterministic chemical heuristics.

        Parameters
        ----------
        target : str
            Product molecule SMILES.
        reactants : str or List[str]
            Reactant representation passed to the heuristic checker.

        Returns
        -------
        Tuple[int, List]
            Status code and valid pathways.
        """
        status_code, valid_pathways = heuristic_checker(target, reactants)
        return status_code, valid_pathways

    def _check_pathway_ml(self, target: str, reactants: str) -> int:
        """
        Evaluates a reaction step using the trained machine learning model.

        Parameters
        ----------
        target : str
            Product molecule SMILES.
        reactants : str
            Reactant molecule SMILES.

        Returns
        -------
        int
            Binary prediction where ``1`` indicates a hallucination and
            ``0`` indicates a valid reaction step.
        """
        datapoint = (target, reactants)
        features = self.featurizer.featurize([datapoint])

        # Check if featurization completely failed natively inside DeepChem
        if features.size == 0 or (
            isinstance(features[0], np.ndarray) and features[0].size == 0
        ):
            print(
                "Warning: Featurization failed on input strings. Marking as hallucination by default."
            )
            return 1

        dataset = dc.data.NumpyDataset(X=features)
        y_pred = self.model.predict(dataset)

        # Isolate raw positive-class scalar probabilities across Sklearn vs. GNN shape spaces
        if len(y_pred.shape) == 3 and y_pred.shape[2] == 2:
            prob = y_pred[0, 0, 1]
        elif len(y_pred.shape) == 2 and y_pred.shape[1] == 2:
            prob = y_pred[0, 1]
        else:
            prob = y_pred[0].flatten()[0]

        return 1 if prob >= self.threshold else 0
