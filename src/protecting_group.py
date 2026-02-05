"""
Protecting Group Module for Retrosynthesis

This module provides comprehensive functionality for:
1. Identifying functional groups that may need protection
2. Suggesting appropriate protecting groups for each site
3. Supporting auto and human-in-the-loop (HITL) modes
4. Tracking active protecting groups across synthesis steps

References:
    Greene, T.W. and Wuts, P.G.M. (2006) Protective Groups in Organic Synthesis.
    4th Edition, John Wiley & Sons, Inc., Hoboken.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
import re
import hashlib
from dataclasses import dataclass, field
from typing import Literal, Optional

# =============================================================================
# FUNCTIONAL GROUPS DATABASE (55+ SMARTS patterns in 10 categories)
# =============================================================================

FUNCTIONAL_GROUPS = {
    # -------------------------------------------------------------------------
    # Category 1: Alcohols (12 patterns)
    # -------------------------------------------------------------------------
    "primary_alcohol": {
        "smarts":
        "[CX4H2][OX2H]",
        "name":
        "Primary Alcohol",
        "category":
        "alcohol",
        "reactivity":
        "high",
        "compatible_pgs":
        ["TBS", "TBDPS", "Bn", "Ac", "PMB", "MOM", "THP", "TES", "TIPS"],
    },
    "secondary_alcohol": {
        "smarts": "[CX4H1]([#6])([#6])[OX2H]",
        "name": "Secondary Alcohol",
        "category": "alcohol",
        "reactivity": "medium",
        "compatible_pgs": ["TBS", "TBDPS", "Bn", "Ac", "PMB", "THP", "TES"],
    },
    "tertiary_alcohol": {
        "smarts": "[CX4H0]([#6])([#6])([#6])[OX2H]",
        "name": "Tertiary Alcohol",
        "category": "alcohol",
        "reactivity": "low",
        "compatible_pgs": ["TBS", "Ac", "TMS"],
    },
    "allylic_alcohol": {
        "smarts": "[CX3]=[CX3][CX4H][OX2H]",
        "name": "Allylic Alcohol",
        "category": "alcohol",
        "reactivity": "high",
        "compatible_pgs": ["TBS", "Bn", "PMB", "Ac", "THP"],
    },
    "propargylic_alcohol": {
        "smarts": "[CX2]#[CX2][CX4H][OX2H]",
        "name": "Propargylic Alcohol",
        "category": "alcohol",
        "reactivity": "high",
        "compatible_pgs": ["TBS", "Bn", "THP", "TBDPS"],
    },
    "benzylic_alcohol": {
        "smarts": "[cX3][CX4H][OX2H]",
        "name": "Benzylic Alcohol",
        "category": "alcohol",
        "reactivity": "high",
        "compatible_pgs": ["TBS", "Ac", "Piv", "THP"],
    },
    "homoallylic_alcohol": {
        "smarts": "[CX3]=[CX3][CX4][CX4H][OX2H]",
        "name": "Homoallylic Alcohol",
        "category": "alcohol",
        "reactivity": "medium",
        "compatible_pgs": ["TBS", "Bn", "PMB", "THP"],
    },
    "alpha_hydroxy_ketone": {
        "smarts": "[CX4H]([OX2H])[CX3](=[OX1])",
        "name": "α-Hydroxy Ketone",
        "category": "alcohol",
        "reactivity": "high",
        "compatible_pgs": ["TBS", "Bn", "acetonide"],
    },
    "beta_hydroxy_ketone": {
        "smarts": "[CX4H]([OX2H])[CX4][CX3](=[OX1])",
        "name": "β-Hydroxy Ketone",
        "category": "alcohol",
        "reactivity": "medium",
        "compatible_pgs": ["TBS", "Bn", "Ac"],
    },
    "alpha_hydroxy_ester": {
        "smarts": "[CX4H]([OX2H])[CX3](=[OX1])[OX2]",
        "name": "α-Hydroxy Ester",
        "category": "alcohol",
        "reactivity": "high",
        "compatible_pgs": ["TBS", "Bn", "acetonide", "MOM"],
    },
    "hemiacetal": {
        "smarts": "[CX4H]([OX2H])([OX2][#6])",
        "name": "Hemiacetal",
        "category": "alcohol",
        "reactivity": "high",
        "compatible_pgs": ["Ac", "Bn", "TBS"],
    },
    "enol": {
        "smarts": "[CX3]=[CX3][OX2H]",
        "name": "Enol",
        "category": "alcohol",
        "reactivity": "very_high",
        "compatible_pgs": ["TBS", "TMS", "Ac", "TES"],
    },

    # -------------------------------------------------------------------------
    # Category 2: Phenols (6 patterns)
    # -------------------------------------------------------------------------
    "phenol": {
        "smarts": "[cX3][OX2H]",
        "name": "Phenol",
        "category": "phenol",
        "reactivity": "medium",
        "compatible_pgs": ["Bn", "Ac", "MOM", "MEM", "TBS", "PMB"],
    },
    "catechol": {
        "smarts": "[cX3]([OX2H])[cX3][OX2H]",
        "name": "Catechol (1,2-dihydroxybenzene)",
        "category": "phenol",
        "reactivity": "medium",
        "compatible_pgs": ["methylenedioxy", "acetonide", "Bn"],
    },
    "hydroquinone": {
        "smarts": "[cX3]1([OX2H])[cX3][cX3][cX3]([OX2H])[cX3][cX3]1",
        "name": "Hydroquinone (1,4-dihydroxybenzene)",
        "category": "phenol",
        "reactivity": "medium",
        "compatible_pgs": ["Bn", "Ac", "TBS"],
    },
    "resorcinol": {
        "smarts": "[cX3]1([OX2H])[cX3][cX3]([OX2H])[cX3][cX3][cX3]1",
        "name": "Resorcinol (1,3-dihydroxybenzene)",
        "category": "phenol",
        "reactivity": "medium",
        "compatible_pgs": ["Bn", "MOM", "Ac"],
    },
    "naphthol": {
        "smarts": "[cX3]1[cX3][cX3][cX3]2[cX3]([OX2H])[cX3][cX3][cX3][cX3]12",
        "name": "Naphthol",
        "category": "phenol",
        "reactivity": "medium",
        "compatible_pgs": ["Bn", "MOM", "TBS", "Ac"],
    },
    "hydroxypyridine": {
        "smarts": "[nX2][cX3][cX3][OX2H]",
        "name": "Hydroxypyridine",
        "category": "phenol",
        "reactivity": "medium",
        "compatible_pgs": ["Bn", "MOM", "Ac"],
    },

    # -------------------------------------------------------------------------
    # Category 3: Diols and Polyols (8 patterns)
    # -------------------------------------------------------------------------
    "vicinal_diol": {
        "smarts": "[CX4H]([OX2H])[CX4H]([OX2H])",
        "name": "1,2-Diol (vicinal)",
        "category": "diol",
        "reactivity": "high",
        "compatible_pgs": ["acetonide", "benzylidene", "cyclohexylidene"],
    },
    "1_3_diol": {
        "smarts": "[CX4H]([OX2H])[CX4H2][CX4H]([OX2H])",
        "name": "1,3-Diol",
        "category": "diol",
        "reactivity": "high",
        "compatible_pgs": ["benzylidene", "acetonide", "cyclohexylidene"],
    },
    "1_2_amino_alcohol": {
        "smarts": "[CX4H]([OX2H])[CX4H]([NX3H2])",
        "name": "1,2-Amino Alcohol",
        "category": "diol",
        "reactivity": "high",
        "compatible_pgs": ["oxazolidinone", "Boc", "acetonide"],
    },
    "triol": {
        "smarts": "[CX4H]([OX2H])[CX4H]([OX2H])[CX4H]([OX2H])",
        "name": "1,2,3-Triol",
        "category": "diol",
        "reactivity": "high",
        "compatible_pgs": ["acetonide", "benzylidene"],
    },
    "cis_diol": {
        "smarts": "[CX4H]1([OX2H])[CX4H]([OX2H])[#6][#6]1",
        "name": "cis-1,2-Diol (cyclic)",
        "category": "diol",
        "reactivity": "high",
        "compatible_pgs": ["acetonide", "cyclohexylidene"],
    },
    "glycol": {
        "smarts": "[OX2H][CX4H2][CX4H2][OX2H]",
        "name": "Ethylene Glycol derivative",
        "category": "diol",
        "reactivity": "high",
        "compatible_pgs": ["acetonide", "benzylidene"],
    },
    "pinacol": {
        "smarts": "[CX4]([OX2H])([#6])[CX4]([OX2H])([#6])",
        "name": "Pinacol (tertiary 1,2-diol)",
        "category": "diol",
        "reactivity": "low",
        "compatible_pgs": ["acetonide"],
    },
    "sugar_diol": {
        "smarts": "[OX2][CX4H]([OX2H])[CX4H]([OX2H])[OX2]",
        "name": "Sugar vicinal diol",
        "category": "diol",
        "reactivity": "high",
        "compatible_pgs": ["acetonide", "benzylidene"],
    },

    # -------------------------------------------------------------------------
    # Category 4: Amines (10 patterns)
    # -------------------------------------------------------------------------
    "primary_amine": {
        "smarts": "[CX4][NX3H2]",
        "name": "Primary Amine",
        "category": "amine",
        "reactivity": "high",
        "compatible_pgs": ["Boc", "Cbz", "Fmoc", "Ac", "Ts"],
    },
    "secondary_amine": {
        "smarts": "[CX4][NX3H1][CX4]",
        "name": "Secondary Amine",
        "category": "amine",
        "reactivity": "medium",
        "compatible_pgs": ["Boc", "Cbz", "Ts"],
    },
    "aniline": {
        "smarts": "[cX3][NX3H2]",
        "name": "Aniline",
        "category": "amine",
        "reactivity": "medium",
        "compatible_pgs": ["Ac", "Boc", "Cbz", "Ts"],
    },
    "n_methyl_aniline": {
        "smarts": "[cX3][NX3H1][CH3]",
        "name": "N-Methyl Aniline",
        "category": "amine",
        "reactivity": "low",
        "compatible_pgs": ["Boc", "Ts"],
    },
    "alpha_amino_acid": {
        "smarts": "[NX3H2][CX4H]([#6])[CX3](=[OX1])[OX2H]",
        "name": "α-Amino Acid",
        "category": "amine",
        "reactivity": "high",
        "compatible_pgs": ["Fmoc", "Boc", "Cbz"],
    },
    "beta_amino_alcohol": {
        "smarts": "[NX3H2][CX4H][CX4H]([OX2H])",
        "name": "β-Amino Alcohol",
        "category": "amine",
        "reactivity": "high",
        "compatible_pgs": ["Boc", "Cbz", "oxazolidinone"],
    },
    "amino_acid_ester": {
        "smarts": "[NX3H2][CX4H]([#6])[CX3](=[OX1])[OX2][#6]",
        "name": "Amino Acid Ester",
        "category": "amine",
        "reactivity": "high",
        "compatible_pgs": ["Fmoc", "Boc", "Cbz"],
    },
    "diamino": {
        "smarts": "[NX3H2][CX4H][CX4H][NX3H2]",
        "name": "1,2-Diamine",
        "category": "amine",
        "reactivity": "high",
        "compatible_pgs": ["Boc", "Cbz"],
    },
    "hydrazine": {
        "smarts": "[NX3H2][NX3H]",
        "name": "Hydrazine",
        "category": "amine",
        "reactivity": "very_high",
        "compatible_pgs": ["Boc", "Cbz"],
    },
    "hydroxylamine": {
        "smarts": "[NX3H2][OX2H]",
        "name": "Hydroxylamine",
        "category": "amine",
        "reactivity": "very_high",
        "compatible_pgs": ["Boc", "Bn"],
    },

    # -------------------------------------------------------------------------
    # Category 5: Heterocyclic N-H (8 patterns)
    # -------------------------------------------------------------------------
    "indole_nh": {
        "smarts": "[nX3H1]1[cX3][cX3]2[cX3][cX3][cX3][cX3][cX3]12",
        "name": "Indole N-H",
        "category": "heterocyclic_nh",
        "reactivity": "medium",
        "compatible_pgs": ["Boc", "SEM", "Ts"],
    },
    "pyrrole_nh": {
        "smarts": "[nX3H1]1[cX3][cX3][cX3][cX3]1",
        "name": "Pyrrole N-H",
        "category": "heterocyclic_nh",
        "reactivity": "medium",
        "compatible_pgs": ["Boc", "SEM", "Ts", "TIPS"],
    },
    "imidazole_nh": {
        "smarts": "[nX3H1]1[cX3][nX2][cX3][cX3]1",
        "name": "Imidazole N-H",
        "category": "heterocyclic_nh",
        "reactivity": "medium",
        "compatible_pgs": ["Boc", "SEM", "Ts", "Bn", "Trt"],
    },
    "pyrazole_nh": {
        "smarts": "[nX3H1]1[nX2][cX3][cX3][cX3]1",
        "name": "Pyrazole N-H",
        "category": "heterocyclic_nh",
        "reactivity": "medium",
        "compatible_pgs": ["Boc", "SEM", "Ts"],
    },
    "benzimidazole_nh": {
        "smarts": "[nX3H1]1[cX3][nX2][cX3]2[cX3][cX3][cX3][cX3][cX3]12",
        "name": "Benzimidazole N-H",
        "category": "heterocyclic_nh",
        "reactivity": "medium",
        "compatible_pgs": ["Boc", "SEM"],
    },
    "carbazole_nh": {
        "smarts":
        "[nX3H1]1[cX3]2[cX3][cX3][cX3][cX3][cX3]2[cX3]2[cX3]1[cX3][cX3][cX3][cX3]2",
        "name": "Carbazole N-H",
        "category": "heterocyclic_nh",
        "reactivity": "low",
        "compatible_pgs": ["Boc", "Bn"],
    },
    "triazole_nh": {
        "smarts": "[nX3H1]1[nX2][nX2][cX3][cX3]1",
        "name": "Triazole N-H",
        "category": "heterocyclic_nh",
        "reactivity": "medium",
        "compatible_pgs": ["Boc", "SEM"],
    },
    "tetrazole_nh": {
        "smarts": "[nX3H1]1[nX2][nX2][nX2][cX3]1",
        "name": "Tetrazole N-H",
        "category": "heterocyclic_nh",
        "reactivity": "medium",
        "compatible_pgs": ["Bn", "Trt", "PMB"],
    },

    # -------------------------------------------------------------------------
    # Category 6: Carbonyls (6 patterns)
    # -------------------------------------------------------------------------
    "aldehyde": {
        "smarts": "[CX3H1](=[OX1])[#6]",
        "name": "Aldehyde",
        "category": "carbonyl",
        "reactivity": "very_high",
        "compatible_pgs": ["acetal", "dithiane", "dithiolane"],
    },
    "ketone": {
        "smarts": "[CX3](=[OX1])([#6])[#6]",
        "name": "Ketone",
        "category": "carbonyl",
        "reactivity": "high",
        "compatible_pgs": ["ketal", "dithiane", "dithiolane"],
    },
    "alpha_keto_ester": {
        "smarts": "[CX3](=[OX1])[CX3](=[OX1])[OX2]",
        "name": "α-Keto Ester",
        "category": "carbonyl",
        "reactivity": "very_high",
        "compatible_pgs": ["ketal"],
    },
    "beta_keto_ester": {
        "smarts": "[CX3](=[OX1])[CX4H2][CX3](=[OX1])[OX2]",
        "name": "β-Keto Ester",
        "category": "carbonyl",
        "reactivity": "high",
        "compatible_pgs": ["enol_silyl_ether", "ketal"],
    },
    "1_3_dicarbonyl": {
        "smarts": "[CX3](=[OX1])[CX4H2][CX3](=[OX1])",
        "name": "1,3-Dicarbonyl",
        "category": "carbonyl",
        "reactivity": "high",
        "compatible_pgs": ["enol_silyl_ether"],
    },
    "formyl": {
        "smarts": "[cX3][CX3H1](=[OX1])",
        "name": "Formyl (on aromatic)",
        "category": "carbonyl",
        "reactivity": "very_high",
        "compatible_pgs": ["acetal", "dithiane"],
    },

    # -------------------------------------------------------------------------
    # Category 7: Carboxylic Acids and Derivatives (5 patterns)
    # -------------------------------------------------------------------------
    "carboxylic_acid": {
        "smarts":
        "[CX3](=[OX1])[OX2H]",
        "name":
        "Carboxylic Acid",
        "category":
        "carboxylic_acid",
        "reactivity":
        "high",
        "compatible_pgs": [
            "Me_ester", "Et_ester", "tBu_ester", "Bn_ester", "allyl_ester",
            "TMS_ester"
        ],
    },
    "amide_nh": {
        "smarts": "[CX3](=[OX1])[NX3H2]",
        "name": "Primary Amide",
        "category": "carboxylic_acid",
        "reactivity": "low",
        "compatible_pgs": ["Boc", "Bn"],
    },
    "carbamate_nh": {
        "smarts": "[OX2][CX3](=[OX1])[NX3H]",
        "name": "Carbamate N-H",
        "category": "carboxylic_acid",
        "reactivity": "low",
        "compatible_pgs": ["Bn"],
    },
    "urea_nh": {
        "smarts": "[NX3H][CX3](=[OX1])[NX3H]",
        "name": "Urea N-H",
        "category": "carboxylic_acid",
        "reactivity": "low",
        "compatible_pgs": ["Boc"],
    },
    "sulfonamide_nh": {
        "smarts": "[SX4](=[OX1])(=[OX1])[NX3H]",
        "name": "Sulfonamide N-H",
        "category": "carboxylic_acid",
        "reactivity": "low",
        "compatible_pgs": ["Boc", "Bn"],
    },

    # -------------------------------------------------------------------------
    # Category 8: Thiols (4 patterns)
    # -------------------------------------------------------------------------
    "thiol": {
        "smarts": "[CX4][SX2H]",
        "name": "Thiol",
        "category": "thiol",
        "reactivity": "very_high",
        "compatible_pgs": ["Trt", "Acm", "STbu", "Bn", "PMB"],
    },
    "thiophenol": {
        "smarts": "[cX3][SX2H]",
        "name": "Thiophenol",
        "category": "thiol",
        "reactivity": "very_high",
        "compatible_pgs": ["Ac", "Bn"],
    },
    "cysteine": {
        "smarts": "[NX3][CX4H]([CX4H2][SX2H])[CX3](=[OX1])",
        "name": "Cysteine (thiol)",
        "category": "thiol",
        "reactivity": "very_high",
        "compatible_pgs": ["Trt", "Acm", "STbu", "Mmt"],
    },
    "dithiol": {
        "smarts": "[CX4H]([SX2H])[CX4H]([SX2H])",
        "name": "1,2-Dithiol",
        "category": "thiol",
        "reactivity": "very_high",
        "compatible_pgs": ["Acm"],
    },

    # -------------------------------------------------------------------------
    # Category 9: Terminal Alkynes (2 patterns)
    # -------------------------------------------------------------------------
    "terminal_alkyne": {
        "smarts": "[CX2]#[CX2H]",
        "name": "Terminal Alkyne",
        "category": "alkyne",
        "reactivity": "medium",
        "compatible_pgs": ["TMS", "TIPS", "TES"],
    },
    "propargyl": {
        "smarts": "[CX4H2][CX2]#[CX2H]",
        "name": "Propargyl (with alkyne)",
        "category": "alkyne",
        "reactivity": "medium",
        "compatible_pgs": ["TMS", "TIPS"],
    },

    # -------------------------------------------------------------------------
    # Category 10: Phosphorus (4 patterns)
    # -------------------------------------------------------------------------
    "phosphate": {
        "smarts": "[PX4](=[OX1])([OX2H])([OX2])[OX2]",
        "name": "Phosphate",
        "category": "phosphorus",
        "reactivity": "medium",
        "compatible_pgs": ["Bn", "cyanoethyl", "methyl_phosphate"],
    },
    "phosphonate": {
        "smarts": "[PX4](=[OX1])([OX2H])([OX2])[#6]",
        "name": "Phosphonate",
        "category": "phosphorus",
        "reactivity": "medium",
        "compatible_pgs": ["Bn", "methyl_phosphate", "cyanoethyl"],
    },
    "phosphoric_acid": {
        "smarts": "[PX4](=[OX1])([OX2H])([OX2H])[OX2]",
        "name": "Phosphoric Acid",
        "category": "phosphorus",
        "reactivity": "high",
        "compatible_pgs": ["Bn", "cyanoethyl"],
    },
    "phosphoramidite": {
        "smarts": "[PX3]([OX2])([OX2])[NX3]",
        "name": "Phosphoramidite",
        "category": "phosphorus",
        "reactivity": "medium",
        "compatible_pgs": ["cyanoethyl"],
    },

    # -------------------------------------------------------------------------
    # Category 11: Esters and Lactones (3 patterns)
    # -------------------------------------------------------------------------
    "ester": {
        "smarts": "[CX3](=[OX1])[OX2][#6]",
        "name": "Ester",
        "category": "ester",
        "reactivity": "low",
        "compatible_pgs": [],
    },
    "lactone": {
        "smarts": "[CX3R](=[OX1])[OX2R]",
        "name": "Lactone (cyclic ester)",
        "category": "ester",
        "reactivity": "low",
        "compatible_pgs": [],
    },
    "thioester": {
        "smarts": "[CX3](=[OX1])[SX2][#6]",
        "name": "Thioester",
        "category": "ester",
        "reactivity": "medium",
        "compatible_pgs": [],
    },

    # -------------------------------------------------------------------------
    # Category 12: Tertiary Amines (2 patterns)
    # -------------------------------------------------------------------------
    "tertiary_amine": {
        "smarts": "[NX3;H0]([#6])([#6])[#6]",
        "name": "Tertiary Amine",
        "category": "amine",
        "reactivity": "low",
        "compatible_pgs": [],
    },
    "tertiary_amine_aromatic": {
        "smarts": "[nX3;H0]([#6])[#6]",
        "name": "Aromatic Tertiary Amine",
        "category": "amine",
        "reactivity": "low",
        "compatible_pgs": [],
    },

    # -------------------------------------------------------------------------
    # Category 13: Amides (3 patterns)
    # -------------------------------------------------------------------------
    "secondary_amide": {
        "smarts": "[CX3](=[OX1])[NX3H1][#6]",
        "name": "Secondary Amide",
        "category": "amide",
        "reactivity": "low",
        "compatible_pgs": ["Boc", "Bn"],
    },
    "tertiary_amide": {
        "smarts": "[CX3](=[OX1])[NX3;H0]([#6])[#6]",
        "name": "Tertiary Amide",
        "category": "amide",
        "reactivity": "low",
        "compatible_pgs": [],
    },
    "lactam": {
        "smarts": "[CX3R](=[OX1])[NX3R]",
        "name": "Lactam (cyclic amide)",
        "category": "amide",
        "reactivity": "low",
        "compatible_pgs": ["Boc"],
    },

    # -------------------------------------------------------------------------
    # Category 14: Epoxides (1 pattern)
    # -------------------------------------------------------------------------
    "epoxide": {
        "smarts": "[OX2r3]1[CX4r3][CX4r3]1",
        "name": "Epoxide",
        "category": "epoxide",
        "reactivity": "very_high",
        "compatible_pgs": [],
    },
}

# =============================================================================
# PROTECTING GROUPS DATABASE (40+ groups with full metadata)
# =============================================================================

PROTECTING_GROUPS = {
    # -------------------------------------------------------------------------
    # Silyl Ethers (for alcohols/phenols)
    # -------------------------------------------------------------------------
    "TMS": {
        "name": "Trimethylsilyl",
        "abbreviation": "TMS",
        "protects": ["alcohol", "phenol", "alkyne", "carboxylic_acid"],
        "protection_smarts": "[Si](C)(C)C",
        "protection_reagent": "TMSCl, Et3N or TMSOTf, 2,6-lutidine",
        "deprotection": ["TBAF", "K2CO3/MeOH", "dilute HCl", "aqueous workup"],
        "stability": {
            "acid": "very_low",
            "base": "medium",
            "hydrogenation": "high",
            "oxidation": "high"
        },
        "orthogonality_score": 0.5,
        "cost_score": 0.95,
        "selectivity": "low",
    },
    "TES": {
        "name": "Triethylsilyl",
        "abbreviation": "TES",
        "protects": ["alcohol", "phenol", "alkyne"],
        "protection_smarts": "[Si](CC)(CC)CC",
        "protection_reagent": "TESCl, imidazole, DMF",
        "deprotection": ["TBAF", "K2CO3/MeOH", "AcOH"],
        "stability": {
            "acid": "low",
            "base": "high",
            "hydrogenation": "high",
            "oxidation": "high"
        },
        "orthogonality_score": 0.6,
        "cost_score": 0.85,
        "selectivity": "medium",
    },
    "TBS": {
        "name": "tert-Butyldimethylsilyl",
        "abbreviation": "TBS",
        "protects": ["alcohol", "phenol"],
        "protection_smarts": "[Si](C)(C)C(C)(C)C",
        "protection_reagent": "TBSCl, imidazole, DMF or TBSOTf, 2,6-lutidine",
        "deprotection": ["TBAF", "HF-pyridine", "AcOH", "PPTS/MeOH"],
        "stability": {
            "acid": "low",
            "base": "high",
            "hydrogenation": "high",
            "oxidation": "high"
        },
        "orthogonality_score": 0.85,
        "cost_score": 0.7,
        "selectivity": "medium",
    },
    "TIPS": {
        "name": "Triisopropylsilyl",
        "abbreviation": "TIPS",
        "protects": ["alcohol", "phenol", "alkyne", "heterocyclic_nh"],
        "protection_smarts": "[Si](C(C)C)(C(C)C)C(C)C",
        "protection_reagent":
        "TIPSCl, imidazole, DMF or TIPSOTf, 2,6-lutidine",
        "deprotection": ["TBAF", "HF-pyridine"],
        "stability": {
            "acid": "medium",
            "base": "high",
            "hydrogenation": "high",
            "oxidation": "high"
        },
        "orthogonality_score": 0.9,
        "cost_score": 0.5,
        "selectivity": "high",
    },
    "TBDPS": {
        "name": "tert-Butyldiphenylsilyl",
        "abbreviation": "TBDPS",
        "protects": ["alcohol", "phenol"],
        "protection_smarts": "[Si](c1ccccc1)(c2ccccc2)C(C)(C)C",
        "protection_reagent": "TBDPSCl, imidazole, DMF",
        "deprotection": ["TBAF", "HF-pyridine"],
        "stability": {
            "acid": "medium",
            "base": "high",
            "hydrogenation": "high",
            "oxidation": "high"
        },
        "orthogonality_score": 0.95,
        "cost_score": 0.4,
        "selectivity": "very_high_primary",
    },
    "TPS": {
        "name": "Triphenylsilyl",
        "abbreviation": "TPS",
        "protects": ["alcohol"],
        "protection_smarts": "[Si](c1ccccc1)(c2ccccc2)c3ccccc3",
        "protection_reagent": "TPSCl, imidazole, DMF",
        "deprotection": ["TBAF", "HF-pyridine"],
        "stability": {
            "acid": "medium",
            "base": "high",
            "hydrogenation": "high",
            "oxidation": "high"
        },
        "orthogonality_score": 0.9,
        "cost_score": 0.3,
        "selectivity": "very_high_primary",
    },

    # -------------------------------------------------------------------------
    # Benzyl Ethers (for alcohols/phenols)
    # -------------------------------------------------------------------------
    "Bn": {
        "name":
        "Benzyl",
        "abbreviation":
        "Bn",
        "protects": [
            "alcohol", "phenol", "amine", "thiol", "carboxylic_acid",
            "heterocyclic_nh"
        ],
        "protection_smarts":
        "Cc1ccccc1",
        "protection_reagent":
        "BnBr, NaH, DMF or BnBr, Ag2O",
        "deprotection": ["H2/Pd-C", "Na/NH3 (Birch)", "DDQ"],
        "stability": {
            "acid": "high",
            "base": "high",
            "hydrogenation": "low",
            "oxidation": "medium"
        },
        "orthogonality_score":
        0.8,
        "cost_score":
        0.8,
        "selectivity":
        "medium",
    },
    "PMB": {
        "name": "p-Methoxybenzyl",
        "abbreviation": "PMB",
        "protects": ["alcohol", "phenol", "amine", "thiol"],
        "protection_smarts": "Cc1ccc(OC)cc1",
        "protection_reagent": "PMBCl, NaH, DMF or PMBOH, CSA",
        "deprotection": ["DDQ", "CAN", "TFA"],
        "stability": {
            "acid": "low",
            "base": "high",
            "hydrogenation": "low",
            "oxidation": "low"
        },
        "orthogonality_score": 0.75,
        "cost_score": 0.65,
        "selectivity": "medium",
    },
    "DMB": {
        "name": "3,4-Dimethoxybenzyl",
        "abbreviation": "DMB",
        "protects": ["alcohol", "amine"],
        "protection_smarts": "Cc1ccc(OC)c(OC)c1",
        "protection_reagent": "DMBCl, NaH, DMF",
        "deprotection": ["DDQ", "CAN"],
        "stability": {
            "acid": "low",
            "base": "high",
            "hydrogenation": "low",
            "oxidation": "low"
        },
        "orthogonality_score": 0.7,
        "cost_score": 0.55,
        "selectivity": "medium",
    },
    "Trt": {
        "name": "Trityl (Triphenylmethyl)",
        "abbreviation": "Trt",
        "protects": ["alcohol", "amine", "thiol", "heterocyclic_nh"],
        "protection_smarts": "C(c1ccccc1)(c2ccccc2)c3ccccc3",
        "protection_reagent": "TrtCl, Et3N, DMAP",
        "deprotection": ["TFA", "dilute AcOH", "BF3·Et2O"],
        "stability": {
            "acid": "very_low",
            "base": "high",
            "hydrogenation": "medium",
            "oxidation": "high"
        },
        "orthogonality_score": 0.85,
        "cost_score": 0.5,
        "selectivity": "very_high_primary",
    },

    # -------------------------------------------------------------------------
    # Acetals (for alcohols)
    # -------------------------------------------------------------------------
    "MOM": {
        "name": "Methoxymethyl",
        "abbreviation": "MOM",
        "protects": ["alcohol", "phenol"],
        "protection_smarts": "COC",
        "protection_reagent": "MOMCl, i-Pr2NEt, CH2Cl2",
        "deprotection": ["dilute HCl", "TFA", "PPTS/MeOH", "ZnBr2"],
        "stability": {
            "acid": "low",
            "base": "high",
            "hydrogenation": "high",
            "oxidation": "high"
        },
        "orthogonality_score": 0.7,
        "cost_score": 0.85,
        "selectivity": "medium",
    },
    "MEM": {
        "name": "Methoxyethoxymethyl",
        "abbreviation": "MEM",
        "protects": ["alcohol", "phenol"],
        "protection_smarts": "COCCOC",
        "protection_reagent": "MEMCl, i-Pr2NEt, CH2Cl2",
        "deprotection": ["ZnBr2", "TiCl4", "dilute HCl"],
        "stability": {
            "acid": "low",
            "base": "high",
            "hydrogenation": "high",
            "oxidation": "high"
        },
        "orthogonality_score": 0.7,
        "cost_score": 0.75,
        "selectivity": "medium",
    },
    "BOM": {
        "name": "Benzyloxymethyl",
        "abbreviation": "BOM",
        "protects": ["alcohol"],
        "protection_smarts": "COCc1ccccc1",
        "protection_reagent": "BOMCl, i-Pr2NEt",
        "deprotection": ["H2/Pd-C", "Na/NH3"],
        "stability": {
            "acid": "medium",
            "base": "high",
            "hydrogenation": "low",
            "oxidation": "medium"
        },
        "orthogonality_score": 0.75,
        "cost_score": 0.7,
        "selectivity": "medium",
    },
    "THP": {
        "name": "Tetrahydropyranyl",
        "abbreviation": "THP",
        "protects": ["alcohol", "phenol"],
        "protection_smarts": "C1CCCCO1",
        "protection_reagent": "DHP, PPTS, CH2Cl2",
        "deprotection": ["PPTS/MeOH", "dilute HCl", "AcOH"],
        "stability": {
            "acid": "low",
            "base": "high",
            "hydrogenation": "high",
            "oxidation": "high"
        },
        "orthogonality_score": 0.65,
        "cost_score": 0.9,
        "selectivity": "low",
    },
    "EE": {
        "name": "Ethoxyethyl",
        "abbreviation": "EE",
        "protects": ["alcohol"],
        "protection_smarts": "CCOC(C)",
        "protection_reagent": "Ethyl vinyl ether, PPTS",
        "deprotection": ["PPTS/MeOH", "dilute HCl"],
        "stability": {
            "acid": "very_low",
            "base": "high",
            "hydrogenation": "high",
            "oxidation": "high"
        },
        "orthogonality_score": 0.55,
        "cost_score": 0.9,
        "selectivity": "low",
    },
    "SEM": {
        "name": "2-(Trimethylsilyl)ethoxymethyl",
        "abbreviation": "SEM",
        "protects": ["alcohol", "heterocyclic_nh"],
        "protection_smarts": "COC[Si](C)(C)C",
        "protection_reagent": "SEMCl, i-Pr2NEt, CH2Cl2",
        "deprotection": ["TBAF", "MgBr2"],
        "stability": {
            "acid": "medium",
            "base": "high",
            "hydrogenation": "high",
            "oxidation": "high"
        },
        "orthogonality_score": 0.8,
        "cost_score": 0.6,
        "selectivity": "medium",
    },

    # -------------------------------------------------------------------------
    # Esters (for alcohols)
    # -------------------------------------------------------------------------
    "Ac": {
        "name": "Acetyl",
        "abbreviation": "Ac",
        "protects": ["alcohol", "phenol", "amine", "thiol"],
        "protection_smarts": "CC(=O)",
        "protection_reagent": "Ac2O, pyridine or AcCl, Et3N",
        "deprotection": ["K2CO3/MeOH", "NH3/MeOH", "LiOH"],
        "stability": {
            "acid": "high",
            "base": "low",
            "hydrogenation": "high",
            "oxidation": "high"
        },
        "orthogonality_score": 0.6,
        "cost_score": 0.95,
        "selectivity": "low",
    },
    "Piv": {
        "name": "Pivaloyl",
        "abbreviation": "Piv",
        "protects": ["alcohol", "amine"],
        "protection_smarts": "CC(C)(C)C(=O)",
        "protection_reagent": "PivCl, pyridine",
        "deprotection": ["LiAlH4", "DIBAL", "LiOH (slow)"],
        "stability": {
            "acid": "high",
            "base": "medium",
            "hydrogenation": "high",
            "oxidation": "high"
        },
        "orthogonality_score": 0.75,
        "cost_score": 0.8,
        "selectivity": "high_primary",
    },
    "Bz": {
        "name": "Benzoyl",
        "abbreviation": "Bz",
        "protects": ["alcohol", "amine"],
        "protection_smarts": "c1ccccc1C(=O)",
        "protection_reagent": "BzCl, pyridine",
        "deprotection": ["K2CO3/MeOH", "NH3/MeOH", "LiOH"],
        "stability": {
            "acid": "high",
            "base": "low",
            "hydrogenation": "high",
            "oxidation": "high"
        },
        "orthogonality_score": 0.65,
        "cost_score": 0.85,
        "selectivity": "medium",
    },
    "Lev": {
        "name": "Levulinoyl",
        "abbreviation": "Lev",
        "protects": ["alcohol"],
        "protection_smarts": "CC(=O)CCC(=O)",
        "protection_reagent": "Levulinic acid, DCC, DMAP",
        "deprotection": ["H2NNH2·AcOH"],
        "stability": {
            "acid": "high",
            "base": "low",
            "hydrogenation": "high",
            "oxidation": "high"
        },
        "orthogonality_score": 0.8,
        "cost_score": 0.7,
        "selectivity": "medium",
    },

    # -------------------------------------------------------------------------
    # Carbamates (for amines)
    # -------------------------------------------------------------------------
    "Boc": {
        "name": "tert-Butoxycarbonyl",
        "abbreviation": "Boc",
        "protects": ["amine", "heterocyclic_nh"],
        "protection_smarts": "CC(C)(C)OC(=O)",
        "protection_reagent": "Boc2O, Et3N or Boc2O, DMAP",
        "deprotection": ["TFA", "HCl/dioxane", "TMSOTf"],
        "stability": {
            "acid": "low",
            "base": "high",
            "hydrogenation": "high",
            "oxidation": "high"
        },
        "orthogonality_score": 0.9,
        "cost_score": 0.75,
        "selectivity": "high",
    },
    "Cbz": {
        "name": "Benzyloxycarbonyl (Carboxybenzyl)",
        "abbreviation": "Cbz",
        "protects": ["amine"],
        "protection_smarts": "c1ccccc1COC(=O)",
        "protection_reagent": "CbzCl, NaHCO3",
        "deprotection": ["H2/Pd-C", "HBr/AcOH", "Na/NH3"],
        "stability": {
            "acid": "high",
            "base": "high",
            "hydrogenation": "low",
            "oxidation": "medium"
        },
        "orthogonality_score": 0.85,
        "cost_score": 0.8,
        "selectivity": "high",
    },
    "Fmoc": {
        "name": "9-Fluorenylmethoxycarbonyl",
        "abbreviation": "Fmoc",
        "protects": ["amine"],
        "protection_smarts": "C1c2ccccc2c3ccccc13COC(=O)",
        "protection_reagent": "Fmoc-Cl, NaHCO3 or Fmoc-OSu",
        "deprotection": ["piperidine/DMF", "DBU", "morpholine"],
        "stability": {
            "acid": "high",
            "base": "low",
            "hydrogenation": "high",
            "oxidation": "high"
        },
        "orthogonality_score": 0.95,
        "cost_score": 0.5,
        "selectivity": "high",
    },
    "Alloc": {
        "name": "Allyloxycarbonyl",
        "abbreviation": "Alloc",
        "protects": ["amine", "alcohol"],
        "protection_smarts": "C=CCOC(=O)",
        "protection_reagent": "Alloc-Cl, NaHCO3",
        "deprotection": ["Pd(PPh3)4, morpholine", "Pd(0), Bu3SnH"],
        "stability": {
            "acid": "high",
            "base": "high",
            "hydrogenation": "medium",
            "oxidation": "medium"
        },
        "orthogonality_score": 0.85,
        "cost_score": 0.7,
        "selectivity": "high",
    },
    "Troc": {
        "name": "2,2,2-Trichloroethoxycarbonyl",
        "abbreviation": "Troc",
        "protects": ["amine", "alcohol"],
        "protection_smarts": "ClC(Cl)(Cl)COC(=O)",
        "protection_reagent": "Troc-Cl, pyridine",
        "deprotection": ["Zn/AcOH", "Cd/AcOH"],
        "stability": {
            "acid": "high",
            "base": "medium",
            "hydrogenation": "medium",
            "oxidation": "high"
        },
        "orthogonality_score": 0.8,
        "cost_score": 0.7,
        "selectivity": "high",
    },

    # -------------------------------------------------------------------------
    # Sulfonamides (for amines)
    # -------------------------------------------------------------------------
    "Ts": {
        "name": "Tosyl (p-Toluenesulfonyl)",
        "abbreviation": "Ts",
        "protects": ["amine", "heterocyclic_nh"],
        "protection_smarts": "Cc1ccc(S(=O)(=O))cc1",
        "protection_reagent": "TsCl, pyridine",
        "deprotection": ["Na/NH3", "Mg/MeOH", "SmI2"],
        "stability": {
            "acid": "high",
            "base": "high",
            "hydrogenation": "medium",
            "oxidation": "high"
        },
        "orthogonality_score": 0.7,
        "cost_score": 0.85,
        "selectivity": "medium",
    },
    "Ns": {
        "name": "Nosyl (4-Nitrobenzenesulfonyl)",
        "abbreviation": "Ns",
        "protects": ["amine"],
        "protection_smarts": "[O-][N+](=O)c1ccc(S(=O)(=O))cc1",
        "protection_reagent": "NsCl, Et3N",
        "deprotection": ["PhSH, K2CO3", "HSCH2CO2H, base"],
        "stability": {
            "acid": "high",
            "base": "medium",
            "hydrogenation": "low",
            "oxidation": "high"
        },
        "orthogonality_score": 0.8,
        "cost_score": 0.75,
        "selectivity": "high",
    },
    "SES": {
        "name": "2-(Trimethylsilyl)ethanesulfonyl",
        "abbreviation": "SES",
        "protects": ["amine"],
        "protection_smarts": "[Si](C)(C)(C)CCS(=O)(=O)",
        "protection_reagent": "SESCl, Et3N",
        "deprotection": ["TBAF", "CsF"],
        "stability": {
            "acid": "high",
            "base": "high",
            "hydrogenation": "high",
            "oxidation": "high"
        },
        "orthogonality_score": 0.85,
        "cost_score": 0.6,
        "selectivity": "high",
    },

    # -------------------------------------------------------------------------
    # Thiol Protecting Groups
    # -------------------------------------------------------------------------
    "Acm": {
        "name": "Acetamidomethyl",
        "abbreviation": "Acm",
        "protects": ["thiol"],
        "protection_smarts": "CC(=O)NC",
        "protection_reagent": "AcmCl, base",
        "deprotection": ["I2", "Hg(OAc)2", "Ag(I)", "Tl(TFA)3"],
        "stability": {
            "acid": "high",
            "base": "high",
            "hydrogenation": "high",
            "oxidation": "medium"
        },
        "orthogonality_score": 0.85,
        "cost_score": 0.7,
        "selectivity": "high",
    },
    "STbu": {
        "name": "S-tert-Butyl",
        "abbreviation": "STbu",
        "protects": ["thiol"],
        "protection_smarts": "CC(C)(C)S",
        "protection_reagent": "tBuSH, I2",
        "deprotection": ["Bu3P", "DTT", "TCEP"],
        "stability": {
            "acid": "high",
            "base": "high",
            "hydrogenation": "high",
            "oxidation": "low"
        },
        "orthogonality_score": 0.8,
        "cost_score": 0.75,
        "selectivity": "high",
    },
    "Mmt": {
        "name": "4-Methoxytrityl",
        "abbreviation": "Mmt",
        "protects": ["thiol"],
        "protection_smarts": "COc1ccc(C(c2ccccc2)c3ccccc3)cc1",
        "protection_reagent": "MmtCl, Et3N",
        "deprotection": ["1% TFA", "AgNO3"],
        "stability": {
            "acid": "very_low",
            "base": "high",
            "hydrogenation": "medium",
            "oxidation": "high"
        },
        "orthogonality_score": 0.9,
        "cost_score": 0.5,
        "selectivity": "very_high",
    },
    "Npys": {
        "name": "3-Nitro-2-pyridinesulfenyl",
        "abbreviation": "Npys",
        "protects": ["thiol"],
        "protection_smarts": "[O-][N+](=O)c1cccnc1S",
        "protection_reagent": "NpysCl, base",
        "deprotection": ["PPh3", "Bu3P", "thiols"],
        "stability": {
            "acid": "high",
            "base": "medium",
            "hydrogenation": "medium",
            "oxidation": "medium"
        },
        "orthogonality_score": 0.85,
        "cost_score": 0.55,
        "selectivity": "high",
    },

    # -------------------------------------------------------------------------
    # Carbonyl Protecting Groups
    # -------------------------------------------------------------------------
    "acetal": {
        "name": "Dimethyl Acetal",
        "abbreviation": "acetal",
        "protects": ["carbonyl"],
        "protection_smarts": "C(OC)(OC)",
        "protection_reagent": "HC(OMe)3, PPTS or MeOH, PPTS",
        "deprotection": ["dilute HCl", "PPTS/acetone-H2O", "Amberlyst-15"],
        "stability": {
            "acid": "low",
            "base": "high",
            "hydrogenation": "high",
            "oxidation": "high"
        },
        "orthogonality_score": 0.7,
        "cost_score": 0.9,
        "selectivity": "medium",
    },
    "ketal": {
        "name": "Dimethyl Ketal",
        "abbreviation": "ketal",
        "protects": ["carbonyl"],
        "protection_smarts": "C(OC)(OC)([#6])[#6]",
        "protection_reagent": "HC(OMe)3, PPTS or MeOH, PPTS",
        "deprotection": ["dilute HCl", "PPTS/acetone-H2O"],
        "stability": {
            "acid": "low",
            "base": "high",
            "hydrogenation": "high",
            "oxidation": "high"
        },
        "orthogonality_score": 0.7,
        "cost_score": 0.9,
        "selectivity": "medium",
    },
    "dioxolane": {
        "name": "1,3-Dioxolane (Ethylene Acetal)",
        "abbreviation": "dioxolane",
        "protects": ["carbonyl"],
        "protection_smarts": "C1OCCO1",
        "protection_reagent": "Ethylene glycol, PPTS, Dean-Stark",
        "deprotection": ["dilute HCl", "PPTS/acetone-H2O", "I2/acetone"],
        "stability": {
            "acid": "low",
            "base": "high",
            "hydrogenation": "high",
            "oxidation": "high"
        },
        "orthogonality_score": 0.75,
        "cost_score": 0.85,
        "selectivity": "medium",
    },
    "dioxane": {
        "name": "1,3-Dioxane",
        "abbreviation": "dioxane",
        "protects": ["carbonyl"],
        "protection_smarts": "C1OCCCO1",
        "protection_reagent": "1,3-Propanediol, PPTS, Dean-Stark",
        "deprotection": ["dilute HCl", "PPTS/acetone-H2O"],
        "stability": {
            "acid": "low",
            "base": "high",
            "hydrogenation": "high",
            "oxidation": "high"
        },
        "orthogonality_score": 0.75,
        "cost_score": 0.85,
        "selectivity": "medium",
    },
    "dithiane": {
        "name": "1,3-Dithiane",
        "abbreviation": "dithiane",
        "protects": ["carbonyl"],
        "protection_smarts": "C1SCCCS1",
        "protection_reagent": "1,3-Propanedithiol, BF3·Et2O",
        "deprotection": ["Hg(ClO4)2", "NBS", "I2", "MeI then hydrolysis"],
        "stability": {
            "acid": "high",
            "base": "high",
            "hydrogenation": "high",
            "oxidation": "low"
        },
        "orthogonality_score": 0.85,
        "cost_score": 0.7,
        "selectivity": "high",
    },
    "dithiolane": {
        "name": "1,3-Dithiolane",
        "abbreviation": "dithiolane",
        "protects": ["carbonyl"],
        "protection_smarts": "C1SCCS1",
        "protection_reagent": "1,2-Ethanedithiol, BF3·Et2O",
        "deprotection": ["Hg(ClO4)2", "NBS", "I2"],
        "stability": {
            "acid": "high",
            "base": "high",
            "hydrogenation": "high",
            "oxidation": "low"
        },
        "orthogonality_score": 0.85,
        "cost_score": 0.7,
        "selectivity": "high",
    },
    "enol_silyl_ether": {
        "name": "Enol Silyl Ether",
        "abbreviation": "enol_silyl_ether",
        "protects": ["carbonyl"],
        "protection_smarts": "[Si]OC=C",
        "protection_reagent": "TMSCl, Et3N or LDA, TMSCl",
        "deprotection": ["TBAF", "dilute HCl", "aqueous workup"],
        "stability": {
            "acid": "low",
            "base": "medium",
            "hydrogenation": "high",
            "oxidation": "medium"
        },
        "orthogonality_score": 0.6,
        "cost_score": 0.8,
        "selectivity": "medium",
    },

    # -------------------------------------------------------------------------
    # Diol Protecting Groups
    # -------------------------------------------------------------------------
    "acetonide": {
        "name": "Acetonide (Isopropylidene)",
        "abbreviation": "acetonide",
        "protects": ["diol"],
        "protection_smarts": "CC1(C)OC([#6])C([#6])O1",
        "protection_reagent": "2,2-Dimethoxypropane, PPTS or acetone, PPTS",
        "deprotection": ["80% AcOH", "PPTS/MeOH", "Amberlyst-15"],
        "stability": {
            "acid": "low",
            "base": "high",
            "hydrogenation": "high",
            "oxidation": "high"
        },
        "orthogonality_score": 0.8,
        "cost_score": 0.85,
        "selectivity": "high_cis_1_2",
    },
    "benzylidene": {
        "name": "Benzylidene Acetal",
        "abbreviation": "benzylidene",
        "protects": ["diol"],
        "protection_smarts": "c1ccccc1C2OC([#6])C([#6])O2",
        "protection_reagent": "PhCH(OMe)2, CSA or PhCHO, PPTS",
        "deprotection": ["H2/Pd-C", "DIBAL", "80% AcOH"],
        "stability": {
            "acid": "low",
            "base": "high",
            "hydrogenation": "low",
            "oxidation": "medium"
        },
        "orthogonality_score": 0.8,
        "cost_score": 0.75,
        "selectivity": "high_1_3",
    },
    "cyclohexylidene": {
        "name": "Cyclohexylidene Acetal",
        "abbreviation": "cyclohexylidene",
        "protects": ["diol"],
        "protection_smarts": "C1(CCCCC1)OC([#6])C([#6])O1",
        "protection_reagent": "Cyclohexanone, PPTS, Dean-Stark",
        "deprotection": ["80% AcOH (slow)", "TFA"],
        "stability": {
            "acid": "medium",
            "base": "high",
            "hydrogenation": "high",
            "oxidation": "high"
        },
        "orthogonality_score": 0.85,
        "cost_score": 0.8,
        "selectivity": "high_cis_1_2",
    },
    "methylenedioxy": {
        "name": "Methylenedioxy",
        "abbreviation": "methylenedioxy",
        "protects": ["diol", "phenol"],
        "protection_smarts": "OCO",
        "protection_reagent": "CH2Br2, NaH or CH2I2, Ag2O",
        "deprotection": ["BBr3", "strong acid"],
        "stability": {
            "acid": "high",
            "base": "high",
            "hydrogenation": "high",
            "oxidation": "high"
        },
        "orthogonality_score": 0.9,
        "cost_score": 0.75,
        "selectivity": "high_catechol",
    },
    "oxazolidinone": {
        "name": "Oxazolidinone",
        "abbreviation": "oxazolidinone",
        "protects": ["diol", "amine"],
        "protection_smarts": "O=C1OC([#6])C([#6])N1",
        "protection_reagent": "CDI, then base",
        "deprotection": ["LiOH/H2O2", "Ba(OH)2"],
        "stability": {
            "acid": "medium",
            "base": "low",
            "hydrogenation": "high",
            "oxidation": "high"
        },
        "orthogonality_score": 0.85,
        "cost_score": 0.7,
        "selectivity": "high_amino_alcohol",
    },

    # -------------------------------------------------------------------------
    # Carboxylic Acid Protecting Groups
    # -------------------------------------------------------------------------
    "Me_ester": {
        "name": "Methyl Ester",
        "abbreviation": "Me_ester",
        "protects": ["carboxylic_acid"],
        "protection_smarts": "COC(=O)",
        "protection_reagent": "MeOH, H2SO4 or CH2N2",
        "deprotection": ["LiOH", "NaOH", "LiI"],
        "stability": {
            "acid": "high",
            "base": "low",
            "hydrogenation": "high",
            "oxidation": "high"
        },
        "orthogonality_score": 0.5,
        "cost_score": 0.95,
        "selectivity": "low",
    },
    "Et_ester": {
        "name": "Ethyl Ester",
        "abbreviation": "Et_ester",
        "protects": ["carboxylic_acid"],
        "protection_smarts": "CCOC(=O)",
        "protection_reagent": "EtOH, H2SO4",
        "deprotection": ["LiOH", "NaOH"],
        "stability": {
            "acid": "high",
            "base": "low",
            "hydrogenation": "high",
            "oxidation": "high"
        },
        "orthogonality_score": 0.5,
        "cost_score": 0.95,
        "selectivity": "low",
    },
    "tBu_ester": {
        "name": "tert-Butyl Ester",
        "abbreviation": "tBu_ester",
        "protects": ["carboxylic_acid"],
        "protection_smarts": "CC(C)(C)OC(=O)",
        "protection_reagent": "Isobutene, H2SO4 or tBuOH, DCC",
        "deprotection": ["TFA", "HCl/dioxane", "ZnBr2"],
        "stability": {
            "acid": "low",
            "base": "high",
            "hydrogenation": "high",
            "oxidation": "high"
        },
        "orthogonality_score": 0.85,
        "cost_score": 0.75,
        "selectivity": "high",
    },
    "Bn_ester": {
        "name": "Benzyl Ester",
        "abbreviation": "Bn_ester",
        "protects": ["carboxylic_acid"],
        "protection_smarts": "c1ccccc1COC(=O)",
        "protection_reagent": "BnBr, K2CO3 or BnOH, DCC",
        "deprotection": ["H2/Pd-C", "Na/NH3"],
        "stability": {
            "acid": "high",
            "base": "high",
            "hydrogenation": "low",
            "oxidation": "medium"
        },
        "orthogonality_score": 0.8,
        "cost_score": 0.8,
        "selectivity": "medium",
    },
    "allyl_ester": {
        "name": "Allyl Ester",
        "abbreviation": "allyl_ester",
        "protects": ["carboxylic_acid"],
        "protection_smarts": "C=CCOC(=O)",
        "protection_reagent": "Allyl alcohol, DCC or allyl bromide, K2CO3",
        "deprotection": ["Pd(PPh3)4, morpholine", "Pd(0), Bu3SnH"],
        "stability": {
            "acid": "high",
            "base": "high",
            "hydrogenation": "low",
            "oxidation": "medium"
        },
        "orthogonality_score": 0.85,
        "cost_score": 0.75,
        "selectivity": "high",
    },
    "TMS_ester": {
        "name": "Trimethylsilyl Ester",
        "abbreviation": "TMS_ester",
        "protects": ["carboxylic_acid"],
        "protection_smarts": "[Si](C)(C)(C)OC(=O)",
        "protection_reagent": "TMSCl, Et3N or BSA",
        "deprotection": ["aqueous workup", "MeOH"],
        "stability": {
            "acid": "very_low",
            "base": "medium",
            "hydrogenation": "high",
            "oxidation": "high"
        },
        "orthogonality_score": 0.4,
        "cost_score": 0.9,
        "selectivity": "low",
    },

    # -------------------------------------------------------------------------
    # Phosphorus Protecting Groups
    # -------------------------------------------------------------------------
    "cyanoethyl": {
        "name": "2-Cyanoethyl",
        "abbreviation": "cyanoethyl",
        "protects": ["phosphorus"],
        "protection_smarts": "N#CCCO",
        "protection_reagent": "2-Cyanoethanol, coupling reagent",
        "deprotection": ["DBU", "Et3N", "NH3"],
        "stability": {
            "acid": "high",
            "base": "low",
            "hydrogenation": "high",
            "oxidation": "high"
        },
        "orthogonality_score": 0.8,
        "cost_score": 0.8,
        "selectivity": "high",
    },
    "methyl_phosphate": {
        "name": "Methyl Phosphate",
        "abbreviation": "methyl_phosphate",
        "protects": ["phosphorus"],
        "protection_smarts": "COP(=O)",
        "protection_reagent": "MeI, base",
        "deprotection": ["LiI", "TMSBr"],
        "stability": {
            "acid": "medium",
            "base": "medium",
            "hydrogenation": "high",
            "oxidation": "high"
        },
        "orthogonality_score": 0.6,
        "cost_score": 0.9,
        "selectivity": "low",
    },
}

# =============================================================================
# ATOM MAPPING
# =============================================================================


def _normalize_implicit_hs(mol: object) -> object:
    """
    Fix implicit hydrogen counts on bracket stereo atoms.

    SMILES bracket notation (e.g. ``[C@]``) sets ``NoImplicit=True`` on the
    atom, causing RDKit to assign 0 implicit Hs even when the atom is
    undervalent. This is chemically incorrect for common organic atoms like
    ``[C@]`` with only 3 bonds (should have 1 H).

    This function detects such atoms, adds the missing explicit Hs, and
    returns a clean molecule with correct H counts.

    Parameters
    ----------
    mol : Mol
        RDKit Mol object that may have incorrect implicit H counts.

    Returns
    -------
    Mol
        A new Mol with correct hydrogen counts, or the original if
        no fix was needed.
    """
    # Default valences for common organic atoms
    default_valence = {6: 4, 7: 3, 8: 2, 16: 2}

    rw = Chem.RWMol(mol)
    modified = False
    for atom in rw.GetAtoms():
        if not atom.GetNoImplicit():
            continue
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in default_valence:
            continue
        expected = default_valence[atomic_num]
        bond_order_sum = sum(b.GetBondTypeAsDouble() for b in atom.GetBonds())
        current_total = bond_order_sum + atom.GetNumExplicitHs()
        if current_total < expected:
            missing = int(round(expected - current_total))
            atom.SetNoImplicit(False)
            atom.SetNumExplicitHs(missing)
            modified = True

    if not modified:
        return mol

    fixed = rw.GetMol()
    # Use partial sanitization to avoid strict valence check on transition
    Chem.SanitizeMol(
        fixed,
        Chem.SanitizeFlags.SANITIZE_ALL
        ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES)

    # Round-trip through SMILES to get a fully clean molecule
    clean_smiles = Chem.MolToSmiles(fixed, canonical=True)
    clean_mol = Chem.MolFromSmiles(clean_smiles)
    return clean_mol if clean_mol is not None else mol


def canonicalize_and_map_atoms(smiles: str) -> tuple[str, dict, object]:
    """
    Canonicalize a SMILES and assign stable atom map numbers.

    Parses the SMILES with RDKit, canonicalizes it, normalizes implicit
    hydrogen counts (fixing bracket stereo atoms like ``[C@]``), and assigns
    atom map numbers based on the canonical atom ordering. The resulting
    mapped SMILES uses notation like ``[CH3:1][C:2](=[O:3])[OH:4]`` so that
    atom map numbers survive round-trips and can be used as stable site
    identifiers.

    Parameters
    ----------
    smiles : str
        SMILES string of the molecule.

    Returns
    -------
    tuple[str, dict, Mol]
        - mapped_smiles: SMILES with atom map numbers embedded.
        - atom_map_dict: ``{map_num: {"symbol": str, "idx": int}}`` mapping
          each atom map number to its element symbol and canonical index.
        - mol: the RDKit Mol object (canonical, with map numbers set).
        Returns ``("", {}, None)`` for invalid SMILES.

    Examples
    --------
    >>> mapped, amap, mol = canonicalize_and_map_atoms("CCO")
    >>> len(amap) > 0
    True
    >>> all(isinstance(k, int) for k in amap)
    True
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return "", {}, None

    # Normalize implicit Hs on bracket stereo atoms (e.g. [C@] -> [C@@H])
    mol = _normalize_implicit_hs(mol)

    # Canonicalize so atom ordering is deterministic
    canonical = Chem.MolToSmiles(mol, canonical=True)
    mol = Chem.MolFromSmiles(canonical)
    if mol is None:
        return "", {}, None

    atom_map_dict: dict[int, dict] = {}
    for atom in mol.GetAtoms():
        map_num = atom.GetIdx() + 1  # 1-based
        atom.SetAtomMapNum(map_num)
        atom_map_dict[map_num] = {
            "symbol": atom.GetSymbol(),
            "idx": atom.GetIdx(),
        }

    mapped_smiles = Chem.MolToSmiles(mol, canonical=True)
    return mapped_smiles, atom_map_dict, mol


def _generate_site_id(functional_group_id: str,
                      atom_map_numbers: tuple) -> str:
    """Generate a deterministic site ID from functional group and atom maps."""
    key = f"{functional_group_id}:{','.join(str(n) for n in sorted(atom_map_numbers))}"
    return hashlib.sha256(key.encode()).hexdigest()[:12]


# =============================================================================
# DATA CLASSES
# =============================================================================


@dataclass
class ProtectionSite:
    """Represents a functional group site that can be protected."""

    site_id: str
    """Deterministic hash key derived from functional_group_id + atom_map_numbers."""

    atom_map_numbers: tuple
    """Canonical atom map numbers involved in this functional group match."""

    functional_group_id: str
    """ID key from FUNCTIONAL_GROUPS dictionary."""

    functional_group_name: str
    """Human-readable name of the functional group."""

    category: str
    """Category of the functional group (e.g., 'alcohol', 'amine')."""

    reactivity: str
    """Reactivity level: 'low', 'medium', 'high', 'very_high'."""

    compatible_pgs: list = field(default_factory=list)
    """List of compatible protecting group IDs."""


@dataclass
class PGSuggestion:
    """Represents a protecting group suggestion for a site."""

    protecting_group_id: str
    """ID key from PROTECTING_GROUPS dictionary."""

    name: str
    """Full name of the protecting group."""

    abbreviation: str
    """Common abbreviation."""

    score: float
    """Ranking score (higher is better)."""

    protection_reagent: str
    """Reagents needed for protection."""

    deprotection_conditions: list = field(default_factory=list)
    """List of deprotection methods."""

    stability: dict = field(default_factory=dict)
    """Stability profile under various conditions."""

    orthogonality_score: float = 0.0
    """Score indicating orthogonality to other PGs."""


@dataclass
class ProtectionRecommendation:
    """Complete recommendation for a protection site."""

    site: ProtectionSite
    """The protection site."""

    suggestions: list = field(default_factory=list)
    """List of PGSuggestion objects, ranked by score."""


@dataclass
class ActiveProtection:
    """Represents a single active protecting group in the synthesis tree."""

    site_id: str
    """Links to ProtectionSite.site_id."""

    protecting_group_id: str
    """PG ID from PROTECTING_GROUPS dictionary."""

    atom_map_numbers: tuple
    """Which atoms are protected."""

    step_added: str
    """Step number where PG was added."""

    step_removed: Optional[str] = None
    """Step number where PG was removed (None = still active)."""


@dataclass
class ProtectionState:
    """Tracks all active protecting groups across the synthesis tree."""

    active_protections: list = field(default_factory=list)
    """List of ActiveProtection objects."""

    def get_active_at_step(self, step: str) -> list:
        """Return protections that are active at a given step."""
        return [
            p for p in self.active_protections
            if p.step_removed is None or (step < p.step_removed)
        ]

    def get_currently_active(self) -> list:
        """Return protections that have not been removed."""
        return [p for p in self.active_protections if p.step_removed is None]

    def add_protection(self, site_id: str, pg_id: str, atom_map_numbers: tuple,
                       step: str) -> None:
        """Record that a protecting group was added at a step."""
        self.active_protections.append(
            ActiveProtection(
                site_id=site_id,
                protecting_group_id=pg_id,
                atom_map_numbers=atom_map_numbers,
                step_added=step,
            ))

    def remove_protection(self, site_id: str, step: str) -> bool:
        """Mark a protecting group as removed at a step. Returns True if found."""
        for p in self.active_protections:
            if p.site_id == site_id and p.step_removed is None:
                p.step_removed = step
                return True
        return False

    def get_active_pg_ids(self) -> list[str]:
        """Get list of currently active protecting group IDs (for orthogonality scoring)."""
        return [
            p.protecting_group_id for p in self.active_protections
            if p.step_removed is None
        ]

    def to_dict(self) -> dict:
        """Serialize to a JSON-compatible dictionary."""
        return {
            "active_protections": [{
                "site_id": p.site_id,
                "protecting_group_id": p.protecting_group_id,
                "atom_map_numbers": list(p.atom_map_numbers),
                "step_added": p.step_added,
                "step_removed": p.step_removed,
            } for p in self.active_protections]
        }

    @classmethod
    def from_dict(cls, data: dict) -> 'ProtectionState':
        """Deserialize from a dictionary."""
        if not data:
            return cls()
        protections = []
        for p in data.get("active_protections", []):
            protections.append(
                ActiveProtection(
                    site_id=p["site_id"],
                    protecting_group_id=p["protecting_group_id"],
                    atom_map_numbers=tuple(p["atom_map_numbers"]),
                    step_added=p["step_added"],
                    step_removed=p.get("step_removed"),
                ))
        return cls(active_protections=protections)


# =============================================================================
# CORE FUNCTIONS
# =============================================================================


def identify_protection_sites(
        smiles: str,
        mapped_mol: object = None,
        atom_map: dict = None,
        allow_overlap: bool = False) -> list[ProtectionSite]:
    """
    Identify all functional groups in a molecule that may need protection.

    Uses RDKit substructure matching with SMARTS patterns to identify
    functional groups that commonly require protection in organic synthesis.
    Uses canonical atom map numbers for stable site identification.

    Parameters
    ----------
    smiles : str
        SMILES string of the molecule to analyze.
    mapped_mol : Mol, optional
        Pre-mapped RDKit Mol object from ``canonicalize_and_map_atoms``.
        If not provided, atom mapping is performed internally.
    atom_map : dict, optional
        Atom map dictionary from ``canonicalize_and_map_atoms``.
        Required when ``mapped_mol`` is provided.
    allow_overlap : bool, optional
        If True, allow overlapping matches even within the same category.
        This shows all possible interpretations of each site (e.g. a site
        matching both ``secondary_alcohol`` and ``beta_hydroxy_ketone``).
        Useful for HITL mode where the chemist picks what to protect.
        If False (default), within-category overlaps are deduplicated,
        keeping the more specific (longer SMARTS) match first.

    Returns
    -------
    list[ProtectionSite]
        List of ProtectionSite objects, each representing a site that
        could potentially be protected. Returns empty list for invalid SMILES.

    Examples
    --------
    >>> sites = identify_protection_sites("CC(O)CC(=O)C")
    >>> len(sites) > 0
    True
    >>> any(s.category == 'alcohol' for s in sites)
    True
    """
    # Perform atom mapping if not provided
    if mapped_mol is None or atom_map is None:
        _, atom_map, mapped_mol = canonicalize_and_map_atoms(smiles)
        if mapped_mol is None:
            return []

    mol = mapped_mol

    # Build index-to-map-number lookup from the mapped molecule
    idx_to_map = {}
    for atom in mol.GetAtoms():
        map_num = atom.GetAtomMapNum()
        if map_num > 0:
            idx_to_map[atom.GetIdx()] = map_num

    sites = []
    # Track exact (fg_id, atoms) pairs to avoid true duplicates
    seen_exact: set = set()
    # Category-aware deduplication (only used when allow_overlap=False):
    # block overlapping matches within the same category, allowing
    # cross-category overlap.
    seen_atoms_by_category: dict[str, set] = {}

    # Sort functional groups by specificity (more specific patterns first)
    # This is approximated by pattern length
    sorted_fg = sorted(FUNCTIONAL_GROUPS.items(),
                       key=lambda x: len(x[1]["smarts"]),
                       reverse=True)

    for fg_id, fg_data in sorted_fg:
        smarts = fg_data["smarts"]
        pattern = Chem.MolFromSmarts(smarts)

        if pattern is None:
            continue

        matches = mol.GetSubstructMatches(pattern)
        category = fg_data["category"]

        for match in matches:
            match_set = frozenset(match)

            # Always skip true duplicates (exact same pattern + atoms)
            exact_key = (fg_id, match_set)
            if exact_key in seen_exact:
                continue

            if not allow_overlap:
                # Default mode: block overlapping within same category
                category_seen = seen_atoms_by_category.get(category, set())
                already_matched = any(match_set & seen
                                      for seen in category_seen)
                if already_matched:
                    continue

            # Convert RDKit indices to canonical atom map numbers
            atom_map_numbers = tuple(
                idx_to_map.get(idx, idx + 1) for idx in match)

            site_id = _generate_site_id(fg_id, atom_map_numbers)

            site = ProtectionSite(
                site_id=site_id,
                atom_map_numbers=atom_map_numbers,
                functional_group_id=fg_id,
                functional_group_name=fg_data["name"],
                category=fg_data["category"],
                reactivity=fg_data["reactivity"],
                compatible_pgs=fg_data["compatible_pgs"].copy())
            sites.append(site)
            seen_exact.add(exact_key)
            seen_atoms_by_category.setdefault(category, set()).add(match_set)

    return sites


def _calculate_pg_score(pg_id: str,
                        site: ProtectionSite,
                        existing_pgs: list[str] = None) -> float:
    """
    Calculate a score for a protecting group based on compatibility and context.
    
    Parameters
    ----------
    pg_id : str
        Protecting group ID.
    site : ProtectionSite
        The site being protected.
    existing_pgs : list[str], optional
        List of protecting groups already present/planned.
        
    Returns
    -------
    float
        Score between 0 and 1, higher is better.
    """
    if pg_id not in PROTECTING_GROUPS:
        return 0.0

    pg_data = PROTECTING_GROUPS[pg_id]

    # Base score from orthogonality
    score = pg_data.get("orthogonality_score", 0.5)

    # Boost for cost efficiency
    cost_score = pg_data.get("cost_score", 0.5)
    score = score * 0.7 + cost_score * 0.3

    # Reactivity matching - high reactivity sites need more robust PGs
    reactivity_map = {"low": 0.2, "medium": 0.5, "high": 0.8, "very_high": 1.0}
    site_reactivity = reactivity_map.get(site.reactivity, 0.5)

    # Stability consideration - more reactive sites need more stable PGs
    stability = pg_data.get("stability", {})
    avg_stability = sum(1.0 if v == "high" else 0.5 if v == "medium" else 0.2
                        for v in stability.values()) / max(len(stability), 1)

    # Adjust score based on reactivity/stability match
    if site_reactivity > 0.7 and avg_stability < 0.5:
        score *= 0.7  # Penalize unstable PGs for reactive sites

    # Check orthogonality with existing PGs
    if existing_pgs:
        for existing_pg in existing_pgs:
            if existing_pg in PROTECTING_GROUPS:
                existing_data = PROTECTING_GROUPS[existing_pg]
                # Check for overlapping deprotection conditions
                existing_deprot = set(existing_data.get("deprotection", []))
                new_deprot = set(pg_data.get("deprotection", []))
                overlap = existing_deprot & new_deprot
                if overlap:
                    score *= 0.8  # Slight penalty for non-orthogonal PGs

    return min(1.0, max(0.0, score))


def suggest_protecting_groups(site: ProtectionSite,
                              existing_pgs: list[str] = None,
                              max_suggestions: int = 5) -> list[PGSuggestion]:
    """
    Suggest protecting groups for a given site, ranked by compatibility.
    
    Parameters
    ----------
    site : ProtectionSite
        The protection site to find PGs for.
    existing_pgs : list[str], optional
        List of protecting groups already present/planned.
    max_suggestions : int
        Maximum number of suggestions to return.
        
    Returns
    -------
    list[PGSuggestion]
        List of PGSuggestion objects, ranked by score (highest first).
    """
    if existing_pgs is None:
        existing_pgs = []

    suggestions = []

    for pg_id in site.compatible_pgs:
        if pg_id not in PROTECTING_GROUPS:
            continue

        pg_data = PROTECTING_GROUPS[pg_id]
        score = _calculate_pg_score(pg_id, site, existing_pgs)

        suggestion = PGSuggestion(
            protecting_group_id=pg_id,
            name=pg_data["name"],
            abbreviation=pg_data["abbreviation"],
            score=score,
            protection_reagent=pg_data["protection_reagent"],
            deprotection_conditions=pg_data.get("deprotection", []).copy(),
            stability=pg_data.get("stability", {}).copy(),
            orthogonality_score=pg_data.get("orthogonality_score", 0.5))
        suggestions.append(suggestion)

    # Sort by score (highest first)
    suggestions.sort(key=lambda x: x.score, reverse=True)

    return suggestions[:max_suggestions]


def get_protection_recommendations(
        smiles: str,
        mode: Literal["auto", "hitl"] = "auto",
        existing_pgs: list[str] = None,
        max_suggestions_per_site: int = 5,
        mapped_mol: object = None,
        atom_map: dict = None) -> list[ProtectionRecommendation]:
    """
    Get protection recommendations for a molecule.

    This is the main entry point for the protecting group recommendation system.

    Parameters
    ----------
    smiles : str
        SMILES string of the molecule.
    mode : Literal["auto", "hitl"]
        - "auto": Returns one representative site per functional group type
          with only the top PG suggestion. When multiple instances of the same
          functional group exist, the highest-reactivity instance is kept.
        - "hitl": Returns ALL identified sites (including every instance of
          the same functional group) with a ranked list of compatible PGs.
    existing_pgs : list[str], optional
        List of protecting groups already present/planned (for orthogonality scoring).
    max_suggestions_per_site : int
        Maximum suggestions per site (only used in "hitl" mode).
    mapped_mol : Mol, optional
        Pre-mapped RDKit Mol from ``canonicalize_and_map_atoms``.
    atom_map : dict, optional
        Atom map dictionary from ``canonicalize_and_map_atoms``.

    Returns
    -------
    list[ProtectionRecommendation]
        List of recommendations. In "hitl" mode, one per identified site
        (including duplicate functional group types). In "auto" mode, one
        per functional group type.

    Examples
    --------
    >>> recs = get_protection_recommendations("CC(O)CC(=O)C", mode="auto")
    >>> for rec in recs:
    ...     print(f"{rec.site.functional_group_name}: {rec.suggestions[0].abbreviation}")

    >>> recs = get_protection_recommendations("CC(O)CC(=O)C", mode="hitl")
    >>> for rec in recs:
    ...     print(f"{rec.site.functional_group_name}:")
    ...     for sug in rec.suggestions:
    ...         print(f"  - {sug.abbreviation} (score: {sug.score:.2f})")
    """
    if existing_pgs is None:
        existing_pgs = []

    # In HITL mode, allow overlapping matches so the chemist sees every
    # possible interpretation (e.g. a site matching both secondary_alcohol
    # and beta_hydroxy_ketone). In auto mode, deduplicate within categories.
    sites = identify_protection_sites(smiles,
                                      mapped_mol=mapped_mol,
                                      atom_map=atom_map,
                                      allow_overlap=(mode == "hitl"))
    recommendations = []

    for site in sites:
        suggestions = suggest_protecting_groups(
            site,
            existing_pgs=existing_pgs,
            max_suggestions=max_suggestions_per_site if mode == "hitl" else 1)

        rec = ProtectionRecommendation(site=site, suggestions=suggestions)
        recommendations.append(rec)

    # In auto mode, collapse to one representative per functional group type.
    # Keep the instance with the highest reactivity (most important to protect).
    if mode == "auto":
        reactivity_rank = {"very_high": 4, "high": 3, "medium": 2, "low": 1}
        best_per_fg: dict[str, ProtectionRecommendation] = {}
        for rec in recommendations:
            fg_id = rec.site.functional_group_id
            if fg_id not in best_per_fg:
                best_per_fg[fg_id] = rec
            else:
                current_rank = reactivity_rank.get(rec.site.reactivity, 0)
                best_rank = reactivity_rank.get(
                    best_per_fg[fg_id].site.reactivity, 0)
                if current_rank > best_rank:
                    best_per_fg[fg_id] = rec
        recommendations = list(best_per_fg.values())

    return recommendations


def format_recommendations_for_prompt(
        recommendations: list[ProtectionRecommendation],
        mode: Literal["auto", "hitl"] = "auto") -> str:
    """
    Format protection recommendations as a string for LLM prompt context.

    Includes canonical atom map numbers so the LLM can reference specific atoms.

    Parameters
    ----------
    recommendations : list[ProtectionRecommendation]
        List of recommendations from get_protection_recommendations.
    mode : Literal["auto", "hitl"]
        Display mode for formatting.

    Returns
    -------
    str
        Formatted string suitable for inclusion in LLM prompts.
    """
    if not recommendations:
        return "No protection sites identified in this molecule."

    lines = ["IDENTIFIED PROTECTION SITES AND RECOMMENDATIONS:", ""]

    for i, rec in enumerate(recommendations, 1):
        site = rec.site
        lines.append(
            f"Site {i}: {site.functional_group_name} (id: {site.site_id})")
        lines.append(f"  Category: {site.category}")
        lines.append(f"  Reactivity: {site.reactivity}")
        lines.append(f"  Atom map numbers: {site.atom_map_numbers}")

        if mode == "auto" and rec.suggestions:
            sug = rec.suggestions[0]
            lines.append(f"  Recommended PG: {sug.name} ({sug.abbreviation})")
            lines.append(f"    Protection: {sug.protection_reagent}")
            lines.append(
                f"    Deprotection: {', '.join(sug.deprotection_conditions[:3])}"
            )
        elif mode == "hitl":
            lines.append("  Available protecting groups (ranked):")
            for j, sug in enumerate(rec.suggestions, 1):
                lines.append(
                    f"    {j}. {sug.name} ({sug.abbreviation}) - Score: {sug.score:.2f}"
                )
                lines.append(f"       Protection: {sug.protection_reagent}")
                lines.append(
                    f"       Deprotection: {', '.join(sug.deprotection_conditions[:2])}"
                )

        lines.append("")

    return "\n".join(lines)


# =============================================================================
# HITL (Human-in-the-Loop) API FUNCTIONS
# =============================================================================


@dataclass
class ProtectionSelection:
    """Represents a user's selection of a site and protecting group."""

    site_index: int
    """Index of the selected site (0-based, from identify_protection_sites)."""

    site: ProtectionSite
    """The selected protection site."""

    selected_pg: PGSuggestion
    """The selected protecting group."""


@dataclass
class ProtectionPlan:
    """Complete protection plan with all user selections."""

    smiles: str
    """Original molecule SMILES."""

    selections: list  # list[ProtectionSelection]
    """List of site/PG selections made by the user."""

    def get_all_pgs(self) -> list[str]:
        """Get list of all selected protecting group IDs."""
        return [sel.selected_pg.protecting_group_id for sel in self.selections]

    def format_summary(self) -> str:
        """Format the protection plan as a readable summary."""
        if not self.selections:
            return "No protection selections made."

        lines = [f"Protection Plan for: {self.smiles}", "=" * 50, ""]

        for i, sel in enumerate(self.selections, 1):
            lines.append(f"{i}. {sel.site.functional_group_name}")
            lines.append(
                f"   Protecting Group: {sel.selected_pg.name} ({sel.selected_pg.abbreviation})"
            )
            lines.append(
                f"   Protection: {sel.selected_pg.protection_reagent}")
            lines.append(
                f"   Deprotection: {', '.join(sel.selected_pg.deprotection_conditions[:2])}"
            )
            lines.append("")

        return "\n".join(lines)


def get_sites_for_selection(smiles: str) -> list[dict]:
    """
    Get all protection sites in a format suitable for user selection.

    This is the first step in the HITL workflow. Returns a list of sites
    with their indices that the user can select from.

    Parameters
    ----------
    smiles : str
        SMILES string of the molecule.

    Returns
    -------
    list[dict]
        List of dictionaries, each containing:
        - 'index': int, the site index for selection
        - 'site': ProtectionSite object
        - 'site_id': str, deterministic hash key for stable identification
        - 'functional_group': str, human-readable name
        - 'category': str, category of functional group
        - 'reactivity': str, reactivity level
        - 'atom_map_numbers': tuple, canonical atom map numbers

    Examples
    --------
    >>> sites = get_sites_for_selection("CCO")
    >>> for s in sites:
    ...     print(f"{s['index']}: {s['functional_group']} ({s['category']})")
    0: Primary Alcohol (alcohol)
    """
    sites = identify_protection_sites(smiles)

    return [{
        'index': i,
        'site': site,
        'site_id': site.site_id,
        'functional_group': site.functional_group_name,
        'category': site.category,
        'reactivity': site.reactivity,
        'atom_map_numbers': site.atom_map_numbers,
    } for i, site in enumerate(sites)]


def get_pg_options_for_site(site: ProtectionSite,
                            existing_pgs: list[str] = None,
                            max_options: int = 10) -> list[dict]:
    """
    Get protecting group options for a selected site.
    
    This is the second step in the HITL workflow. After the user selects
    a site, this returns the PG options they can choose from.
    
    Parameters
    ----------
    site : ProtectionSite
        The selected protection site.
    existing_pgs : list[str], optional
        List of protecting groups already selected/present (for orthogonality).
    max_options : int
        Maximum number of options to return.
        
    Returns
    -------
    list[dict]
        List of dictionaries, each containing:
        - 'index': int, the PG index for selection
        - 'suggestion': PGSuggestion object
        - 'name': str, full name
        - 'abbreviation': str, common abbreviation
        - 'score': float, ranking score
        - 'protection_reagent': str
        - 'deprotection': list[str]
        
    Examples
    --------
    >>> sites = get_sites_for_selection("CCO")
    >>> pg_options = get_pg_options_for_site(sites[0]['site'])
    >>> for opt in pg_options:
    ...     print(f"{opt['index']}: {opt['abbreviation']} (score: {opt['score']:.2f})")
    """
    suggestions = suggest_protecting_groups(site,
                                            existing_pgs=existing_pgs,
                                            max_suggestions=max_options)

    return [{
        'index': i,
        'suggestion': sug,
        'name': sug.name,
        'abbreviation': sug.abbreviation,
        'score': sug.score,
        'protection_reagent': sug.protection_reagent,
        'deprotection': sug.deprotection_conditions,
    } for i, sug in enumerate(suggestions)]


def create_protection_selection(site_info: dict,
                                pg_info: dict) -> ProtectionSelection:
    """
    Create a ProtectionSelection from user choices.
    
    Parameters
    ----------
    site_info : dict
        Site dictionary from get_sites_for_selection.
    pg_info : dict
        PG dictionary from get_pg_options_for_site.
        
    Returns
    -------
    ProtectionSelection
        The complete selection object.
    """
    return ProtectionSelection(site_index=site_info['index'],
                               site=site_info['site'],
                               selected_pg=pg_info['suggestion'])


def create_protection_plan(
        smiles: str, selections: list[ProtectionSelection]) -> ProtectionPlan:
    """
    Create a complete protection plan from user selections.
    
    Parameters
    ----------
    smiles : str
        Original molecule SMILES.
    selections : list[ProtectionSelection]
        List of user selections.
        
    Returns
    -------
    ProtectionPlan
        The complete protection plan.
    """
    return ProtectionPlan(smiles=smiles, selections=selections)


def hitl_workflow_example(smiles: str) -> str:
    """
    Example of the complete HITL workflow.

    This demonstrates the API usage pattern for human-in-the-loop
    protecting group selection. Shows ALL identified sites (including
    multiple instances of the same functional group and overlapping
    interpretations) with ranked PG options for each.

    Parameters
    ----------
    smiles : str
        SMILES string of the molecule.

    Returns
    -------
    str
        Description of the workflow with example data.

    Examples
    --------
    >>> print(hitl_workflow_example("CC(O)C=O"))
    """
    lines = ["HITL Workflow Example", "=" * 50, ""]
    lines.append(f"Molecule: {smiles}")
    lines.append("")

    # Step 1: Get all recommendations in HITL mode (shows every site)
    recs = get_protection_recommendations(smiles, mode="hitl")

    lines.append(
        f"STEP 1: All Identified Protection Sites ({len(recs)} found)")
    lines.append("-" * 50)

    if not recs:
        lines.append("No protection sites found.")
        return "\n".join(lines)

    for i, rec in enumerate(recs):
        site = rec.site
        has_pgs = "yes" if rec.suggestions else "no"
        lines.append(
            f"  [{i}] {site.functional_group_name} ({site.category})"
            f" - {site.reactivity} reactivity"
            f" - atoms {site.atom_map_numbers}"
            f" - PGs available: {has_pgs}")

    lines.append("")
    lines.append("→ Chemist selects which site(s) to protect")
    lines.append("")

    # Step 2: Show PG options for EACH site that has suggestions
    lines.append("STEP 2: Protecting Group Options Per Site")
    lines.append("-" * 50)

    selections = []
    for i, rec in enumerate(recs):
        site = rec.site
        if not rec.suggestions:
            lines.append(
                f"  Site [{i}] {site.functional_group_name}: No compatible PGs (info-only site)"
            )
            lines.append("")
            continue

        lines.append(
            f"  Site [{i}] {site.functional_group_name} (atoms {site.atom_map_numbers}):"
        )
        for j, sug in enumerate(rec.suggestions[:5]):
            lines.append(
                f"    [{j}] {sug.abbreviation} - Score: {sug.score:.2f}"
                f"  |  {sug.protection_reagent}")
        lines.append("")

        # For the example, auto-select top PG for the first protectable site
        if not selections:
            site_info = {
                'index': i,
                'site': site,
                'site_id': site.site_id,
                'functional_group': site.functional_group_name,
                'category': site.category,
                'reactivity': site.reactivity,
                'atom_map_numbers': site.atom_map_numbers,
            }
            pg_info = {
                'index': 0,
                'suggestion': rec.suggestions[0],
                'name': rec.suggestions[0].name,
                'abbreviation': rec.suggestions[0].abbreviation,
                'score': rec.suggestions[0].score,
                'protection_reagent': rec.suggestions[0].protection_reagent,
                'deprotection': rec.suggestions[0].deprotection_conditions,
            }
            selections.append(
                create_protection_selection(site_info, pg_info))

    lines.append("→ Chemist selects protecting group for each chosen site")
    lines.append("")

    # Step 3: Create plan from selections
    if selections:
        plan = create_protection_plan(smiles, selections)

        lines.append("STEP 3: Example Protection Plan (auto-selected top PG for first site)")
        lines.append("-" * 50)
        lines.append(plan.format_summary())

    return "\n".join(lines)


def format_sites_for_display(sites: list[dict]) -> str:
    """
    Format sites list as a readable string for display.
    
    Parameters
    ----------
    sites : list[dict]
        Sites from get_sites_for_selection.
        
    Returns
    -------
    str
        Formatted string.
    """
    if not sites:
        return "No protection sites identified."

    lines = ["Available Protection Sites:", ""]
    for s in sites:
        lines.append(f"[{s['index']}] {s['functional_group']}")
        lines.append(f"    Site ID: {s['site_id']}")
        lines.append(f"    Category: {s['category']}")
        lines.append(f"    Reactivity: {s['reactivity']}")
        lines.append(f"    Atom map numbers: {s['atom_map_numbers']}")
        lines.append("")

    return "\n".join(lines)


def format_pg_options_for_display(options: list[dict]) -> str:
    """
    Format PG options list as a readable string for display.
    
    Parameters
    ----------
    options : list[dict]
        Options from get_pg_options_for_site.
        
    Returns
    -------
    str
        Formatted string.
    """
    if not options:
        return "No protecting group options available."

    lines = ["Available Protecting Groups (ranked by score):", ""]
    for opt in options:
        lines.append(f"[{opt['index']}] {opt['name']} ({opt['abbreviation']})")
        lines.append(f"    Score: {opt['score']:.2f}")
        lines.append(f"    Protection: {opt['protection_reagent']}")
        lines.append(f"    Deprotection: {', '.join(opt['deprotection'][:2])}")
        lines.append("")

    return "\n".join(lines)
