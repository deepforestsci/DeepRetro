"""LLM-assisted reaction metadata helpers."""

from __future__ import annotations

import ast
import json
from typing import TypeGuard, cast

from deepretro.utils.llm_helpers import ChatMessage
from deepretro.utils.utils_molecule import (
    calc_chemical_formula,
    calc_mol_wt,
    is_valid_smiles,
)
from deepretro.utils.variables import (
    CONDITIONS_SYS_PROMPT,
    CONDITIONS_USER_PROMPT,
    LITERATURE_SYS_PROMPT,
    LITERATURE_USER_PROMPT,
    REAGENT_SYS_PROMPT,
    REAGENT_USER_PROMPT,
)
from deepretro.metadata_types import (
    ConditionsPayload,
    ConditionsRecommender,
    ConditionsStatusPayload,
    ConditionValue,
    JSONValue,
    LiteratureRecommender,
    LiteratureStatusPayload,
    MetadataError,
    MetadataRecommendation,
    MetadataResponse,
    MetadataResponseStatusPayload,
    MoleculeRecord,
    ReactionParticipants,
    RecommendationStatusPayload,
    ReagentRecommender,
    ReagentStatusPayload,
    StatusPayload,
)

DEFAULT_METADATA_MODEL = "claude-opus-4-20250514"
MAX_COMPLETION_TOKENS = 4096
TOP_P = 0.9

__all__ = [
    "ConditionsRecommender",
    "ConditionsPayload",
    "ConditionsStatusPayload",
    "ConditionValue",
    "JSONValue",
    "LiteratureRecommender",
    "LiteratureStatusPayload",
    "MetadataError",
    "MetadataRecommendation",
    "MetadataResponse",
    "MetadataResponseStatusPayload",
    "MoleculeRecord",
    "RecommendationStatusPayload",
    "ReactionParticipants",
    "ReagentRecommender",
    "ReagentStatusPayload",
    "StatusPayload",
    "build_reagent_records",
    "build_conditions_messages",
    "build_literature_messages",
    "build_reagent_messages",
    "extract_smiles",
    "parse_metadata_response",
    "parse_literature_reaction",
    "parse_reaction_smiles",
    "valid_conditions_payload",
]


def extract_smiles(records: list[MoleculeRecord]) -> list[str]:
    """Extract SMILES strings from molecule records.

    Parameters
    ----------
    records : list[MoleculeRecord]
        Molecule records containing a ``"smiles"`` key.

    Returns
    -------
    list[str]
        SMILES strings in input order.

    Examples
    --------
    >>> extract_smiles([{"smiles": "CCO"}])
    ['CCO']
    """
    return [str(record["smiles"]) for record in records]


def molecule_records_from_smiles(section: str) -> list[MoleculeRecord]:
    """Convert a dot-separated SMILES section into molecule records.

    Parameters
    ----------
    section : str
        Dot-separated SMILES string section.

    Returns
    -------
    list[MoleculeRecord]
        Molecule records with ``"smiles"`` keys.

    Examples
    --------
    >>> molecule_records_from_smiles("CCO.O")
    [{'smiles': 'CCO'}, {'smiles': 'O'}]
    """
    return [
        {"smiles": smiles}
        for smiles in (part.strip() for part in section.split("."))
        if smiles
    ]


def parse_reaction_smiles(reaction_smiles: str) -> ReactionParticipants:
    """Parse a reaction SMILES string into reactant, reagent, and product records.

    Parameters
    ----------
    reaction_smiles : str
        Reaction SMILES in ``reactants>agents>products`` form. The common
        ``reactants>>products`` form is accepted and yields no input reagents.

    Returns
    -------
    ReactionParticipants
        Parsed molecule records.

    Raises
    ------
    ValueError
        If the string does not contain exactly three reaction sections or lacks
        reactants/products.

    Examples
    --------
    >>> parse_reaction_smiles("CCO>>CC=O").reactants
    [{'smiles': 'CCO'}]
    """
    sections = reaction_smiles.split(">")
    if len(sections) != 3:
        raise ValueError(
            "reaction SMILES must use 'reactants>agents>products' or "
            "'reactants>>products' form"
        )

    reactants = molecule_records_from_smiles(sections[0])
    reagents = molecule_records_from_smiles(sections[1])
    product = molecule_records_from_smiles(sections[2])
    if not reactants or not product:
        raise ValueError("reaction SMILES must include reactants and products")
    return ReactionParticipants(reactants=reactants, reagents=reagents, product=product)


def build_reagent_messages(reactants: list[str], product: str) -> list[ChatMessage]:
    """Build the chat messages for reagent prediction.

    Parameters
    ----------
    reactants : list[str]
        Reactant SMILES strings.
    product : str
        Product SMILES string.

    Returns
    -------
    list[ChatMessage]
        System and user messages for the metadata LLM.

    Examples
    --------
    >>> messages = build_reagent_messages(["CCO"], "CC=O")
    >>> [message["role"] for message in messages]
    ['system', 'user']
    """
    user_prompt = REAGENT_USER_PROMPT.replace("{product}", product)
    user_prompt = user_prompt.replace("{reactants}", ", ".join(reactants))
    return [
        {"role": "system", "content": REAGENT_SYS_PROMPT},
        {"role": "user", "content": user_prompt},
    ]


def build_conditions_messages(
    reactants: list[str],
    product: str,
    reagents: list[str],
) -> list[ChatMessage]:
    """Build the chat messages for reaction-condition prediction.

    Parameters
    ----------
    reactants : list[str]
        Reactant SMILES strings.
    product : str
        Product SMILES string.
    reagents : list[str]
        Reagent SMILES strings.

    Returns
    -------
    list[ChatMessage]
        System and user messages for the metadata LLM.

    Examples
    --------
    >>> messages = build_conditions_messages(["CCO"], "CC=O", ["O"])
    >>> "CCO" in messages[1]["content"]
    True
    """
    user_prompt = CONDITIONS_USER_PROMPT.replace("{product}", product)
    user_prompt = user_prompt.replace("{reactants}", ", ".join(reactants))
    user_prompt = user_prompt.replace("{reagents}", ", ".join(reagents))
    return [
        {"role": "system", "content": CONDITIONS_SYS_PROMPT},
        {"role": "user", "content": user_prompt},
    ]


def build_literature_messages(
    reactants: list[str],
    product: str,
    reagents: list[str],
    conditions: ConditionsPayload,
) -> list[ChatMessage]:
    """Build the chat messages for nearest-literature prediction.

    Parameters
    ----------
    reactants : list[str]
        Reactant SMILES strings.
    product : str
        Product SMILES string.
    reagents : list[str]
        Reagent SMILES strings.
    conditions : ConditionsPayload
        Reaction conditions payload.

    Returns
    -------
    list[ChatMessage]
        System and user messages for the metadata LLM.

    Examples
    --------
    >>> messages = build_literature_messages(["CCO"], "CC=O", ["O"], {"temperature": "25 C"})
    >>> "25 C" in messages[1]["content"]
    True
    """
    user_prompt = LITERATURE_USER_PROMPT.replace("{product}", product)
    user_prompt = user_prompt.replace("{reactants}", ", ".join(reactants))
    user_prompt = user_prompt.replace("{reagents}", ", ".join(reagents))
    user_prompt = user_prompt.replace("{conditions}", str(conditions))
    return [
        {"role": "system", "content": LITERATURE_SYS_PROMPT},
        {"role": "user", "content": user_prompt},
    ]


def parse_metadata_response(response_text: str) -> MetadataResponse:
    """Parse an LLM response containing a JSON or Python-literal dictionary.

    Parameters
    ----------
    response_text : str
        Raw LLM response text.

    Returns
    -------
    MetadataResponse
        Parsed dictionary payload.

    Raises
    ------
    ValueError
        If the parsed payload is not a dictionary.

    Examples
    --------
    >>> parse_metadata_response("{'data': ['O']}")
    {'data': ['O']}
    >>> parse_metadata_response("```json\\n{\\"data\\": [\\"O\\"]}\\n```")
    {'data': ['O']}
    """
    text = response_text.strip()
    candidates = [text]
    if text.startswith("```"):
        lines = text.splitlines()
        if lines and lines[0].startswith("```"):
            lines = lines[1:]
        if lines and lines[-1].strip() == "```":
            lines = lines[:-1]
        fenced_text = "\n".join(lines).strip()
        if fenced_text:
            candidates.append(fenced_text)

    start = text.find("{")
    end = text.rfind("}")
    if start != -1 and end > start:
        candidates.append(text[start : end + 1])

    errors: list[str] = []
    for candidate in dict.fromkeys(candidates):
        for parser in (json.loads, ast.literal_eval):
            try:
                parsed: object = parser(candidate)
            except Exception as exc:
                errors.append(str(exc))
                continue
            if isinstance(parsed, dict):
                return cast(MetadataResponse, parsed)
            errors.append("metadata response must be a dictionary")

    raise ValueError("; ".join(errors) or "metadata response must be a dictionary")


def build_reagent_records(reagent_smiles: list[str]) -> list[MoleculeRecord]:
    """Return reagent records for valid SMILES with formula and mass metadata.

    Parameters
    ----------
    reagent_smiles : list[str]
        Candidate reagent SMILES strings from metadata prediction.

    Returns
    -------
    list[MoleculeRecord]
        Reagent records for valid SMILES only.

    Examples
    --------
    >>> records = build_reagent_records(["O", "not-smiles"])
    >>> records[0]["reagent_metadata"]["chemical_formula"]
    'H2O'
    """
    return [
        {
            "smiles": smiles,
            "reagent_metadata": {
                "name": "",
                "chemical_formula": calc_chemical_formula(smiles),
                "mass": calc_mol_wt(smiles),
            },
        }
        for smiles in reagent_smiles
        if is_valid_smiles(smiles)
    ]


def valid_conditions_payload(payload: object) -> TypeGuard[ConditionsPayload]:
    """Return whether a payload contains the required condition fields.

    Parameters
    ----------
    payload : object
        Candidate conditions payload parsed from LLM output.

    Returns
    -------
    bool
        ``True`` when all required condition fields are present.

    Examples
    --------
    >>> valid_conditions_payload({"temperature": "25 C", "pressure": "1 atm", "solvent": "water", "time": "1 h"})
    True
    >>> valid_conditions_payload({"temperature": "25 C"})
    False
    >>> valid_conditions_payload({"temperature": "", "pressure": "1 atm", "solvent": "water", "time": "1 h"})
    False
    """
    required_fields = {"temperature", "pressure", "solvent", "time"}
    if not isinstance(payload, dict) or not required_fields.issubset(payload):
        return False

    condition_fields = cast(dict[str, object], payload)
    for field in required_fields:
        value = condition_fields[field]
        if not isinstance(value, str | int | float):
            return False
        if isinstance(value, str) and not value.strip():
            return False
    return True


def parse_literature_reaction(response_text: str) -> JSONValue:
    """Parse a literature reaction payload from LLM response text.

    Parameters
    ----------
    response_text : str
        Response text containing a Python-literal dictionary with a
        ``"literature_reaction"`` key.

    Returns
    -------
    JSONValue
        Parsed literature reaction payload.

    Examples
    --------
    >>> parse_literature_reaction("{'literature_reaction': 'example'}")
    'example'
    """
    return parse_metadata_response(response_text)["literature_reaction"]
