"""DeepChem-compatible featurizer for reaction-step (product + reactants) pairs."""

import numpy as np
from deepchem.feat import Featurizer, CircularFingerprint

from deepretro.utils import extract_domain_features_single, NUM_DOMAIN_FEATURES


class ReactionStepFeaturizer(Featurizer):
    """
    Featurize a reaction step (product + reactants) into a flat numeric vector.

    Concatenates three parts:

    1. CircularFingerprint (Morgan/ECFP) for the product     — ``size`` bits
    2. CircularFingerprint (Morgan/ECFP) for the reactants   — ``size`` bits
    3. 39 hand-crafted domain features (optional)

    Parameters
    ----------
    radius : int, optional (default 3)
        Morgan fingerprint radius.  ``radius=3`` corresponds to ECFP6,
        which gave the best held-out ROC-AUC during the autoresearch
        sweep (see ``autoresearch/results.tsv``).  Pass ``radius=2``
        explicitly to recover ECFP4 behaviour.
    size : int, optional (default 2048)
        Fingerprint bit length for each molecule.
    use_domain_features : bool, optional (default True)
        If True, appends 39 domain features (atom/bond/ring/MW deltas,
        Tanimoto similarity, macrocycle flags, etc.).

    Notes
    -----
    This class requires RDKit to be installed.

    Examples
    --------
    >>> from deepretro.featurizers.reactionstep import ReactionStepFeaturizer
    >>> featurizer = ReactionStepFeaturizer(radius=3, size=2048)
    >>> reactions = [("CCO", "CC.O"), ("c1ccccc1", "c1ccccc1.Cl")]
    >>> X = featurizer.featurize(reactions)
    >>> X.shape
    (2, 4135)
    """

    def __init__(
        self, radius: int = 3, size: int = 2048, use_domain_features: bool = True
    ) -> None:
        self.radius = radius
        self.size = size
        self.use_domain_features = use_domain_features
        self._fp = CircularFingerprint(radius=radius, size=size)

    @property
    def feature_dim(self) -> int:
        """
        Total length of one feature vector.

        Returns
        -------
        dim : int
            ``2 * size + NUM_DOMAIN_FEATURES`` when
            ``use_domain_features=True``, ``2 * size`` otherwise.
        """
        return 2 * self.size + (NUM_DOMAIN_FEATURES if self.use_domain_features else 0)

    def _featurize(self, datapoint: tuple) -> np.ndarray:
        """
        Featurize a single reaction step.

        Parameters
        ----------
        datapoint : tuple of (str, str)
            ``(product_smiles, reactants_smiles)`` where reactants may be
            dot-separated when there are multiple reactants.

        Returns
        -------
        features : np.ndarray, shape (feature_dim,)
            Flat feature vector.  Returns a NaN vector if the *product*
            SMILES cannot be parsed.  When the joined dot-separated
            reactant string fails to parse but at least one fragment is
            valid, falls back to fingerprinting the first valid fragment
            (instead of dropping the whole row).
        """
        try:
            product_smiles, reactants_smiles = datapoint

            prod_fp = self._fp.featurize([product_smiles])[0]
            # CircularFingerprint returns an empty array (shape (0,)) rather
            # than raising on invalid SMILES.  Treat unparseable product as a
            # full-row NaN — without a valid product we have no anchor.
            if prod_fp.shape != (self.size,):
                return np.full(self.feature_dim, np.nan)

            reac_fp = self._fp.featurize([reactants_smiles])[0]
            if reac_fp.shape != (self.size,):
                # Graceful dot-SMILES fallback: try each fragment individually
                # and use the first valid one.  This recovers rows where a
                # single bad fragment otherwise loses the whole reaction
                # (autoresearch wins ~+0.005 ROC-AUC from this).
                fallback_fp = None
                for frag in reactants_smiles.split("."):
                    if not frag:
                        continue
                    fp = self._fp.featurize([frag])[0]
                    if fp.shape == (self.size,):
                        fallback_fp = fp
                        break
                if fallback_fp is None:
                    return np.full(self.feature_dim, np.nan)
                reac_fp = fallback_fp

            parts = [prod_fp, reac_fp]
            if self.use_domain_features:
                parts.append(
                    extract_domain_features_single(product_smiles, reactants_smiles)
                )
            return np.concatenate(parts)
        except Exception:
            return np.full(self.feature_dim, np.nan)
