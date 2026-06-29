from collections import Counter

from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold


def count_functional_groups(mol):
    """Count common functional groups using SMARTS patterns."""

    smarts = {
        "Alcohols": "[OX2H][CX4]",
        "Phenols": "[OX2H][c]",
        "Amines": "[NX3;H2,H1,H0;!$(NC=O)]",
        "Amides": "[NX3][CX3](=[OX1])",
        "Esters": "[CX3](=[OX1])[OX2][#6]",
        "CarboxylicAcids": "[CX3](=O)[OX2H1]",
        "Ketones": "[#6][CX3](=O)[#6]",
        "Aldehydes": "[CX3H1](=O)[#6]",
        "Ethers": "[OD2]([#6])[#6]",
        "Sulfonamides": "[SX4](=[OX1])(=[OX1])[NX3]",
        "Nitriles": "[CX2]#N",
        "Halides": "[F,Cl,Br,I]",
        "BoronicAcids": "[BX3]([OX2H])[OX2H]",
        "Alkenes": "[CX3]=[CX3]",
        "Alkynes": "[CX2]#[CX2]",
    }

    counts = {}
    for name, pattern in smarts.items():
        patt = Chem.MolFromSmarts(pattern)
        counts[name] = len(mol.GetSubstructMatches(patt))

    return counts


def count_atom_types(mol):
    """Count heteroatoms and halogens."""

    atoms = Counter(atom.GetSymbol() for atom in mol.GetAtoms())

    return {
        "NitrogenCount": atoms["N"],
        "OxygenCount": atoms["O"],
        "SulfurCount": atoms["S"],
        "PhosphorusCount": atoms["P"],
        "FluorineCount": atoms["F"],
        "ChlorineCount": atoms["Cl"],
        "BromineCount": atoms["Br"],
        "IodineCount": atoms["I"],
    }


def count_bond_types(mol):
    """Count unsaturation-related bond types."""

    double_bonds = 0
    triple_bonds = 0
    aromatic_bonds = 0

    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            double_bonds += 1
        elif bond.GetBondType() == Chem.rdchem.BondType.TRIPLE:
            triple_bonds += 1

        if bond.GetIsAromatic():
            aromatic_bonds += 1

    return {
        "DoubleBondCount": double_bonds,
        "TripleBondCount": triple_bonds,
        "AromaticBondCount": aromatic_bonds,
    }


def count_stereochemistry(mol):
    """Count stereocenters."""

    chiral_centers = Chem.FindMolChiralCenters(
        mol,
        includeUnassigned=True,
    )

    defined = sum(1 for _, c in chiral_centers if c != "?")
    undefined = sum(1 for _, c in chiral_centers if c == "?")

    ez_bonds = sum(
        1
        for bond in mol.GetBonds()
        if bond.GetStereo()
        != Chem.rdchem.BondStereo.STEREONONE
    )

    return {
        "ChiralCenters": len(chiral_centers),
        "DefinedChiralCenters": defined,
        "UndefinedChiralCenters": undefined,
        "EZDoubleBonds": ez_bonds,
    }


def murcko_stats(mol):
    """Murcko scaffold statistics."""

    scaffold = MurckoScaffold.GetScaffoldForMol(mol)

    return {
        "MurckoScaffoldAtoms": scaffold.GetNumAtoms(),
        "MurckoScaffoldBonds": scaffold.GetNumBonds(),
    }


def get_retrosynthesis_descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")

    descriptors = {
        # Size
        "MolecularWeight": round(Descriptors.MolWt(mol), 3),
        "ExactMolWt": round(Descriptors.ExactMolWt(mol), 3),
        "HeavyAtomCount": Lipinski.HeavyAtomCount(mol),

        # Physchem
        "LogP": round(Descriptors.MolLogP(mol), 3),
        "TPSA": round(rdMolDescriptors.CalcTPSA(mol), 3),

        # Complexity
        "BertzCT": round(Descriptors.BertzCT(mol), 3),

        # Rings
        "RingCount": Lipinski.RingCount(mol),
        "NumAromaticRings": Lipinski.NumAromaticRings(mol),
        "NumSaturatedRings": Lipinski.NumSaturatedRings(mol),
        "NumAliphaticRings": Lipinski.NumAliphaticRings(mol),
        "NumHeterocycles": rdMolDescriptors.CalcNumHeterocycles(mol),
        "NumAromaticHeterocycles": Lipinski.NumAromaticHeterocycles(mol),
        "NumBridgeheadAtoms": rdMolDescriptors.CalcNumBridgeheadAtoms(mol),
        "NumSpiroAtoms": rdMolDescriptors.CalcNumSpiroAtoms(mol),

        # Flexibility
        "NumRotatableBonds": Lipinski.NumRotatableBonds(mol),
        "FractionCSP3": round(rdMolDescriptors.CalcFractionCSP3(mol), 3),

        # H-bonding
        "NumHDonors": Lipinski.NumHDonors(mol),
        "NumHAcceptors": Lipinski.NumHAcceptors(mol),

        # Charges / radicals
        "FormalCharge": Chem.GetFormalCharge(mol),
        "NumRadicalElectrons": Descriptors.NumRadicalElectrons(mol),
    }

    descriptors.update(count_atom_types(mol))
    descriptors.update(count_bond_types(mol))
    descriptors.update(count_stereochemistry(mol))
    descriptors.update(murcko_stats(mol))
    descriptors.update(count_functional_groups(mol))

    return descriptors


def build_rdkit_user_prompt(smiles):
    descriptors = get_retrosynthesis_descriptors(smiles)

    prompt = f"""Target Molecule

SMILES:
{smiles}

Retrosynthesis-Relevant Properties

"""

    for name, value in descriptors.items():
        prompt += f"- {name}: {value}\n"

    prompt += """

Use these descriptors together with the molecular structure to:
1. Identify plausible retrosynthetic disconnections.
2. Suggest key intermediates and synthons.
3. Identify likely reaction classes.
4. Assess synthetic complexity and major synthetic challenges.
5. Highlight stereochemical, ring-system, and functional-group considerations.
"""

    return prompt
