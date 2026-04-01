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
    3. 27 domain features (optional)

    Parameters
    ----------
    radius : int, optional (default 2)
        Morgan fingerprint radius.  radius=2 corresponds to ECFP4.
    size : int, optional (default 2048)
        Fingerprint bit length for each molecule.
    use_domain_features : bool, optional (default True)
        If True, appends 27 domain features.

    Notes
    -----
    This class requires RDKit to be installed.

    Examples
    --------
    >>> from deepretro.featurizers.reactionstep import ReactionStepFeaturizer
    >>> featurizer = ReactionStepFeaturizer(radius=2, size=2048)
    >>> reactions = [("CCO", "CC.O"), ("c1ccccc1", "c1ccccc1.Cl")]
    >>> X = featurizer.featurize(reactions)
    >>> X.shape
    (2, 4123)
    """

    def __init__(
        self, radius: int = 2, size: int = 2048, use_domain_features: bool = True
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
            ``2 * size + 27`` when ``use_domain_features=True``,
            ``2 * size`` otherwise.
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
            Flat feature vector.  Returns a NaN vector if either SMILES
            cannot be parsed, so invalid rows are distinguishable from
            real data downstream.
        """
        try:
            product_smiles, reactants_smiles = datapoint
            prod_fp = self._fp.featurize([product_smiles])[0]
            reac_fp = self._fp.featurize([reactants_smiles])[0]
            # CircularFingerprint returns an empty array (shape (0,)) rather than raising an error when SMILES is invalid. We detect that here so we always return a well-shaped NaN vector on bad input.
            if prod_fp.shape != (self.size,) or reac_fp.shape != (self.size,):
                return np.full(self.feature_dim, np.nan)
            parts = [prod_fp, reac_fp]
            if self.use_domain_features:
                parts.append(
                    extract_domain_features_single(product_smiles, reactants_smiles)
                )
            return np.concatenate(parts)
        except Exception:
            return np.full(self.feature_dim, np.nan)
