"""DeepChem-compatible featurizer for reaction-step (product + reactants) pairs."""

from __future__ import annotations

import logging
from collections.abc import Iterable

import numpy as np
import deepchem as dc
import deepchem.feat  # noqa: F401  # ensure ``dc.feat`` submodule is registered
from deepchem.feat.graph_data import GraphData
from deepretro.utils import NUM_DOMAIN_FEATURES, extract_domain_features_single
from typing import Any, Callable, Dict, Tuple, TypedDict

logger = logging.getLogger(__name__)

# ``class`` is a Python keyword, so this TypedDict uses the functional syntax.
FeaturizerConfig = TypedDict(
    "FeaturizerConfig",
    {"class": Callable[..., Any], "type": str, "default_params": Dict[str, Any]},
)


class ReactionStepFeaturizer(dc.feat.Featurizer):
    """Represent a reaction as fingerprints or a disconnected molecular graph.

    Input pairs contain ``(product_smiles, reactants_smiles)``. Join multiple
    reactants with a dot. Fingerprints concatenate product features first,
    then reactant features. Graphs combine both sides into one disconnected
    graph; they do not encode an atom mapping or reaction direction.

    Parameters
    ----------
    base_featurizer : str, optional
        ``"circularfingerprint"`` (default) or ``"molgraphconv"``.
    use_domain_features : bool, optional
        Append 15 reaction descriptors to circular fingerprints. Default False.
        This option has no effect on graph features.
    **kwargs : Any
        Options for the DeepChem featurizer, such as ``radius`` and ``size``.
        Sparse fingerprints are not supported.

    Examples
    --------
    >>> feat = ReactionStepFeaturizer(size=32)
    >>> feat.featurize([("CCO", "CC.O")]).shape
    (1, 64)
    >>> ReactionStepFeaturizer(size=32, use_domain_features=True).feature_dim
    79
    """

    def __init__(
        self,
        base_featurizer: str = "circularfingerprint",
        use_domain_features: bool = False,
        **kwargs: Any,
    ) -> None:
        """Resolve the registered base featurizer and its configuration."""
        self.use_domain_features = use_domain_features
        self._vector_feature_dim: int | None = None

        feat_lower = base_featurizer.lower()
        if feat_lower not in FEATURIZER_MAP:
            raise ValueError(
                f"Inner base featurizer '{base_featurizer}' not found in FEATURIZER_MAP."
            )

        config = FEATURIZER_MAP[feat_lower]
        self.base_feat_type = config["type"]
        if self.base_feat_type not in {"1d_vector", "graph"}:
            raise ValueError(
                f"'{base_featurizer}' cannot be used as a base featurizer."
            )
        if kwargs.get("sparse", False):
            raise ValueError("ReactionStepFeaturizer requires dense fingerprints.")

        # Combine default config parameters with remaining explicit user overrides
        merged_params = config["default_params"].copy()
        merged_params.update(kwargs)
        self.base_featurizer = config["class"](**merged_params)

    @property
    def feature_dim(self) -> int:
        """Return the number of columns in a vector feature row.

        Returns
        -------
        int
            Twice the fingerprint size, plus 15 when domain features are enabled.
            For custom vector featurizers, the width of the first successful
            reaction pair.

        Raises
        ------
        ValueError
            If this instance produces graphs, or a custom vector featurizer
            has not yet processed a valid reaction pair.
        """
        if self.base_feat_type != "1d_vector":
            raise ValueError("Graph features do not have a fixed vector dimension.")
        if not isinstance(self.base_featurizer, dc.feat.CircularFingerprint):
            if self._vector_feature_dim is None:
                raise ValueError(
                    "Custom vector dimension requires a valid reaction pair."
                )
            return self._vector_feature_dim
        return 2 * self.base_featurizer.size + (
            NUM_DOMAIN_FEATURES if self.use_domain_features else 0
        )

    def featurize(
        self,
        datapoints: Iterable[Tuple[str, str]],
        log_every_n: int = 1000,
        **kwargs: Any,
    ) -> np.ndarray:
        """Featurize a batch while preserving invalid input positions.

        Parameters
        ----------
        datapoints : iterable of tuple of str and str
            Product and dot-separated reactant SMILES pairs.
        log_every_n : int, optional
            Log progress every this many rows. Must be positive.
        **kwargs : Any
            Options forwarded to the base featurizer.

        Returns
        -------
        numpy.ndarray
            A numeric matrix for fingerprints, including all-NaN failure rows.
            For graphs, a one-dimensional object array with empty arrays at
            failed positions. Empty fingerprint input has shape
            ``(0, feature_dim)``. Custom vector widths are learned from valid
            rows. Before any valid row is seen, failed custom batches have
            shape ``(n_rows, 0)``.

        Raises
        ------
        ValueError
            If ``log_every_n`` is not positive.
        """
        if log_every_n <= 0:
            raise ValueError("log_every_n must be positive.")
        features: list[np.ndarray | GraphData] = []
        for index, datapoint in enumerate(datapoints):
            if index % log_every_n == 0:
                logger.info("Featurizing reaction %d", index)
            try:
                features.append(self._featurize(datapoint, **kwargs))
            except Exception:
                logger.warning("Failed to featurize reaction %d", index, exc_info=True)
                features.append(self._failure_features())
        if self.base_feat_type == "graph":
            result = np.empty(len(features), dtype=object)
            result[:] = features
            return result
        try:
            width = self.feature_dim
        except ValueError:
            return np.empty((len(features), 0))
        vectors = [
            np.full(width, np.nan) if np.asarray(feature).size == 0 else feature
            for feature in features
        ]
        return np.asarray(vectors).reshape(len(vectors), width)

    def _failure_features(self) -> np.ndarray:
        """Return a NaN row when its width is known, otherwise an empty sentinel."""
        try:
            return np.full(self.feature_dim, np.nan)
        except ValueError:
            return np.array([])

    def _featurize(
        self, datapoint: Tuple[str, str], **kwargs: Any
    ) -> np.ndarray | GraphData:
        """
        Featurize a single reaction step.

        Parameters
        ----------
        datapoint : Tuple[str, str]
            Tuple containing:
                - ``product_smiles``: SMILES string of the product molecule.
                - ``reactants_smiles``: SMILES string of the reactant set.

        **kwargs
            Additional keyword arguments passed to the underlying DeepChem
            featurizer.

        Returns
        -------
        numpy.ndarray or deepchem.feat.graph_data.GraphData
            Feature representation of the reaction step.

            For vector-based featurizers, returns a concatenated feature vector consisting of:
            ``[product_features, reactant_features, domain_features]``, where domain features are included only
            when enabled and supported.

            For graph-based featurizers, returns a graph representation generated from the combined product and
            reactant SMILES.

        Notes
        -----
        If either side fails vector featurization, return an all-NaN row.
        Graph failures use DeepChem's empty-array sentinel. The batch method
        preserves these failed positions so callers can keep labels aligned.
        """
        product_smiles, reactants_smiles = datapoint

        if self.base_feat_type == "1d_vector":
            prod_feats = self.base_featurizer.featurize(product_smiles, **kwargs)[0]

            if prod_feats.size == 0:
                return self._failure_features()

            reac_feats = self.base_featurizer.featurize(reactants_smiles, **kwargs)[0]
            if reac_feats.size == 0:
                return self._failure_features()

            parts = [prod_feats, reac_feats]

            # Append domain chemistry metrics only for ECFP calculations
            if self.use_domain_features and isinstance(
                self.base_featurizer, dc.feat.CircularFingerprint
            ):
                parts.append(
                    extract_domain_features_single(product_smiles, reactants_smiles)
                )
            result = np.concatenate(parts)
            self._vector_feature_dim = result.size
            return result

        elif self.base_feat_type == "graph":
            combined_smiles = f"{product_smiles}.{reactants_smiles}"
            return self.base_featurizer.featurize(combined_smiles, **kwargs)[0]

        raise ValueError(f"Unsupported base featurizer type: '{self.base_feat_type}'.")


FEATURIZER_MAP: Dict[str, FeaturizerConfig] = {
    "molgraphconv": {
        "class": dc.feat.MolGraphConvFeaturizer,
        "type": "graph",
        "default_params": {},
    },
    "circularfingerprint": {
        "class": dc.feat.CircularFingerprint,
        "type": "1d_vector",
        "default_params": {"radius": 2, "size": 2048},
    },
    "reactionfeaturizer": {
        "class": ReactionStepFeaturizer,
        "type": "composite",
        "default_params": {},
    },
}
