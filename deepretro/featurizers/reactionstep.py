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

    Raises
    ------
    ValueError
        If ``base_featurizer`` is not a key of ``FEATURIZER_MAP``.
    ValueError
        If the base featurizer does not make vectors or graphs.
    ValueError
        If ``sparse=True`` is passed.
    ValueError
        If a vector featurizer fails on the probe molecule ``"C"``.

    Notes
    -----
    For vector featurizers, the constructor featurizes the probe molecule
    ``"C"`` once. The size of that result sets the row width.

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
        """Build the base featurizer and find the vector row width."""
        self.use_domain_features = use_domain_features

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

        self._feature_dim: int | None = None
        if self.base_feat_type == "1d_vector":
            probe_size = self.base_featurizer.featurize(["C"])[0].size
            if probe_size == 0:
                raise ValueError(
                    "Cannot find the vector width: the probe molecule failed."
                )
            self._feature_dim = 2 * probe_size + (
                NUM_DOMAIN_FEATURES if self._adds_domain_features else 0
            )

    @property
    def _adds_domain_features(self) -> bool:
        """Return True when rows end with the domain descriptors."""
        return self.use_domain_features and isinstance(
            self.base_featurizer, dc.feat.CircularFingerprint
        )

    @property
    def feature_dim(self) -> int:
        """Return the number of columns in a vector feature row.

        Returns
        -------
        int
            Twice the base vector size.
            Add 15 when domain features are on for circular fingerprints.

        Raises
        ------
        ValueError
            If this instance makes graphs.

        Examples
        --------
        >>> ReactionStepFeaturizer(size=32).feature_dim
        64
        """
        if self._feature_dim is None:
            raise ValueError("Graph features do not have a fixed vector dimension.")
        return self._feature_dim

    def featurize(
        self,
        datapoints: Iterable[Tuple[str, str]],
        log_every_n: int = 1000,
        **kwargs: Any,
    ) -> np.ndarray:
        """Featurize a batch and keep one output row per input pair.

        Parameters
        ----------
        datapoints : iterable of tuple of str and str
            Product and dot-separated reactant SMILES pairs.
        log_every_n : int, optional
            Log progress every this many rows. Default 1000.
        **kwargs : Any
            Options forwarded to the base featurizer.

        Returns
        -------
        numpy.ndarray
            For vectors: a float matrix of shape ``(n_rows, feature_dim)``.
            A failed row is all NaN.

            For graphs: an object array of shape ``(n_rows,)``.
            A failed row is an empty array.

        Notes
        -----
        This method does not raise for a bad row.
        It logs a warning with the traceback and keeps a failed row.

        Examples
        --------
        >>> import numpy as np
        >>> X = ReactionStepFeaturizer(size=8).featurize(
        ...     [("CCO", "CC.O"), ("bad", "CC.O")]
        ... )
        >>> X.shape, bool(np.isnan(X[1]).all())
        ((2, 16), True)
        """
        features: list[np.ndarray | GraphData] = []
        for index, datapoint in enumerate(datapoints):
            if index % log_every_n == 0:
                logger.info("Featurizing reaction %d", index)
            try:
                features.append(self._featurize(datapoint, **kwargs))
            except Exception:
                logger.warning("Failed to featurize reaction %d", index, exc_info=True)
                features.append(
                    np.array([])
                    if self.base_feat_type == "graph"
                    else np.full(self.feature_dim, np.nan)
                )
        if self.base_feat_type == "graph":
            result = np.empty(len(features), dtype=object)
            result[:] = features
            return result
        return np.asarray(features).reshape(len(features), self.feature_dim)

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
        A failed graph is an empty array, as in DeepChem. The batch method
        keeps these failed rows so callers can keep labels aligned.
        """
        product_smiles, reactants_smiles = datapoint

        if self.base_feat_type == "graph":
            combined_smiles = f"{product_smiles}.{reactants_smiles}"
            return self.base_featurizer.featurize(combined_smiles, **kwargs)[0]

        prod_feats = self.base_featurizer.featurize(product_smiles, **kwargs)[0]
        if prod_feats.size == 0:
            return np.full(self.feature_dim, np.nan)

        reac_feats = self.base_featurizer.featurize(reactants_smiles, **kwargs)[0]
        if reac_feats.size == 0:
            return np.full(self.feature_dim, np.nan)

        parts = [prod_feats, reac_feats]
        # Append domain chemistry metrics only for ECFP calculations
        if self._adds_domain_features:
            parts.append(
                extract_domain_features_single(product_smiles, reactants_smiles)
            )
        return np.concatenate(parts)


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
