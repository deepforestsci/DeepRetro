"""DeepChem-compatible featurizer for reaction-step (product + reactants) pairs."""

import numpy as np
import deepchem as dc
from deepchem.feat.graph_data import GraphData
from deepretro.utils import extract_domain_features_single
from typing import Any, Tuple


class ReactionStepFeaturizer(dc.feat.Featurizer):
    """
    Featurizer for a single retrosynthesis reaction step.

    Each datapoint is represented as a tuple of ``(product_smiles, reactants_smiles)``.

    The featurizer wraps a DeepChem molecular featurizer and applies it to both the product and reactant molecules. 
    Depending on the selected base featurizer, the output is either:

    * **Vector-based features** (e.g. CircularFingerprint): product and reactant feature vectors are concatenated. 
    Optional reaction-level domain features can also be appended.

    * **Graph-based features** (e.g. MolGraphConvFeaturizer): product and reactant SMILES are combined into a single
    molecular graph representation.

    Parameters
    ----------
    base_featurizer : str
        Name of the base featurizer to use. Must correspond to a key in ``FEATURIZER_MAP``. default="circularfingerprint"

    use_domain_features : bool
        Whether to append handcrafted reaction-level domain features when using circular fingerprints. default=True

    **kwargs
        Additional keyword arguments passed to the underlying DeepChem featurizer. These override any defaults defined 
        in ``FEATURIZER_MAP``.

    Raises
    ------
    ValueError
        If ``base_featurizer`` is not present in ``FEATURIZER_MAP``.
    """

    def __init__(self, 
                base_featurizer: str ="circularfingerprint", 
                use_domain_features: bool = False, 
                **kwargs: Any) -> None:
        """
        Initialize the reaction-step featurizer.

        Parameters
        ----------
        base_featurizer: str
            Base DeepChem featurizer to wrap. default="circularfingerprint"
        use_domain_features: bool
            controls whether reaction-level domain features are appended for fingerprint-based features. default = False
        **kwargs
            Parameters passed to the selected base featurizer. The special argument ``use_domain_features`` 
        """
        # Extract use_domain_features safely from kwargs with a standard True fallback
        self.use_domain_features = use_domain_features
        
        feat_lower = base_featurizer.lower()
        if feat_lower not in FEATURIZER_MAP:
            raise ValueError(f"Inner base featurizer '{base_featurizer}' not found in FEATURIZER_MAP.")
            
        config = FEATURIZER_MAP[feat_lower]
        self.base_feat_type = config["type"]
        
        # Combine default config parameters with remaining explicit user overrides
        merged_params = config["default_params"].copy()
        merged_params.update(kwargs)
        self.base_featurizer = config["class"](**merged_params)

    def _featurize(self, datapoint: Tuple[str, str], **kwargs: Any) -> np.ndarray | GraphData:
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
        Invalid molecules are handled by DeepChem's parent ``Featurizer.featurize()`` implementation. If either 
        the product or reactant featurization fails for vector-based featurizers, an empty NumPy array is returned.
        """
        product_smiles, reactants_smiles = datapoint
        
        if self.base_feat_type == "1d_vector":
            # --- 1D VECTOR ROUTE (Fingerprints, Descriptors) ---
            prod_feats = self.base_featurizer.featurize(product_smiles, **kwargs)[0]
            
            if prod_feats.size == 0:
                return np.array([])
                
            reac_feats = self.base_featurizer.featurize(reactants_smiles, **kwargs)[0]
            if reac_feats.size == 0:
                return np.array([])
                
            parts = [prod_feats, reac_feats]
            
            # Append domain chemistry metrics only for ECFP calculations
            if self.use_domain_features and isinstance(self.base_featurizer, dc.feat.CircularFingerprint):
                parts.append(
                    extract_domain_features_single(product_smiles, reactants_smiles)
                )
            return np.concatenate(parts)
            
        elif self.base_feat_type == "graph":
            # --- GRAPH ROUTE (MolGraphConv, DMPNN) ---
            combined_smiles = f"{product_smiles}.{reactants_smiles}"
            return self.base_featurizer.featurize(combined_smiles, **kwargs)[0]
