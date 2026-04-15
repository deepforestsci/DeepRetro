"""DeepChem-compatible featurizer for RDKit topological path fingerprints."""

import numpy as np
from rdkit.Chem import RDKFingerprint

from deepchem.feat.base_classes import MolecularFeaturizer
from deepchem.utils.typing import RDKitMol


class RDKTopologicalFingerprint(MolecularFeaturizer):
    """RDKit topological (path-based) fingerprint.

    Encodes all linear paths of length ``minPath`` to ``maxPath`` in the
    molecular graph and hashes them into a bit vector.  This captures
    overall molecular shape and substructure connectivity — complementary
    to circular (Morgan/ECFP) fingerprints which focus on local atom
    neighbourhoods.

    Parameters
    ----------
    fpSize : int, optional (default 2048)
        Length of the generated bit vector.
    minPath : int, optional (default 1)
        Minimum path length (in bonds).
    maxPath : int, optional (default 7)
        Maximum path length (in bonds).

    Notes
    -----
    This class requires RDKit to be installed.

    Examples
    --------
    >>> from deepretro.featurizers.rdkit_fingerprint import RDKTopologicalFingerprint
    >>> fp = RDKTopologicalFingerprint(fpSize=2048)
    >>> X = fp.featurize(["CCO", "c1ccccc1"])
    >>> X.shape
    (2, 2048)
    """

    def __init__(
        self,
        fpSize: int = 2048,
        minPath: int = 1,
        maxPath: int = 7,
    ) -> None:
        self.fpSize = fpSize
        self.minPath = minPath
        self.maxPath = maxPath

    def _featurize(self, mol: RDKitMol) -> np.ndarray:
        """Compute the RDKit topological fingerprint for one molecule.

        Parameters
        ----------
        mol : rdkit.Chem.rdchem.Mol
            RDKit Mol object.

        Returns
        -------
        np.ndarray, shape (fpSize,)
            Bit vector as a float array.
        """
        fp = RDKFingerprint(
            mol, fpSize=self.fpSize, minPath=self.minPath, maxPath=self.maxPath
        )
        return np.asarray(fp, dtype=float)
