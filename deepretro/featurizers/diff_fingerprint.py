"""DeepChem-compatible featurizer using dual fingerprints and a difference vector."""

import numpy as np
from deepchem.feat import Featurizer, CircularFingerprint

from deepretro.featurizers.rdkit_fingerprint import RDKTopologicalFingerprint


class DiffFingerprint(Featurizer):
    """
    Featurize a reaction step into a flat vector with a difference fingerprint.

    For each molecule (product and reactants) two fingerprints are computed
    and concatenated: RDKit topological (path-based) + Morgan (circular/ECFP).
    A pointwise difference ``FP(product) − FP(reactant)`` is appended to
    capture the transformation pattern independent of the specific molecules.

    Final layout: ``[FP(reactant), FP(product), FP(diff)]``
    where each FP block is ``RDKit(size) + Morgan(size)`` and
    ``FP(diff) = FP(product) - FP(reactant)``.

    Parameters
    ----------
    radius : int, optional (default 2)
        Morgan fingerprint radius.  radius=2 corresponds to ECFP4.
    size : int, optional (default 2048)
        Bit-vector length for each fingerprint type per molecule.
    min_path : int, optional (default 1)
        Minimum bond-path length for the RDKit topological fingerprint.
    max_path : int, optional (default 7)
        Maximum bond-path length for the RDKit topological fingerprint.

    Notes
    -----
    This class requires RDKit to be installed.

    Examples
    --------
    >>> from deepretro.featurizers.diff_fingerprint import DiffFingerprint
    >>> featurizer = DiffFingerprint(radius=2, size=2048)
    >>> reactions = [("CCO", "CC.O"), ("c1ccccc1", "c1ccccc1.Cl")]
    >>> X = featurizer.featurize(reactions)
    >>> X.shape
    (2, 12288)
    """

    def __init__(
        self,
        radius: int = 2,
        size: int = 2048,
        min_path: int = 1,
        max_path: int = 7,
    ) -> None:
        self.radius = radius
        self.size = size
        self.min_path = min_path
        self.max_path = max_path
        self._morgan = CircularFingerprint(radius=radius, size=size)
        self._rdk = RDKTopologicalFingerprint(
            fpSize=size, minPath=min_path, maxPath=max_path
        )

    @property
    def feature_dim(self) -> int:
        """
        Total length of one feature vector.

        Returns
        -------
        dim : int
            ``6 * size``.
        """
        return 6 * self.size

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

            prod_rdk = self._rdk.featurize([product_smiles])[0]
            prod_morgan = self._morgan.featurize([product_smiles])[0]
            reac_rdk = self._rdk.featurize([reactants_smiles])[0]
            reac_morgan = self._morgan.featurize([reactants_smiles])[0]

            if (
                prod_rdk.shape != (self.size,)
                or prod_morgan.shape != (self.size,)
                or reac_rdk.shape != (self.size,)
                or reac_morgan.shape != (self.size,)
            ):
                return np.full(self.feature_dim, np.nan)

            prod_fp = np.concatenate([prod_rdk, prod_morgan])
            reac_fp = np.concatenate([reac_rdk, reac_morgan])
            diff_fp = prod_fp - reac_fp

            return np.concatenate([reac_fp, prod_fp, diff_fp])
        except Exception:
            return np.full(self.feature_dim, np.nan)
