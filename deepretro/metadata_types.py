"""Type definitions for reaction metadata helpers."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from typing import Literal, TypeAlias, TypedDict

from deepretro.types import ConditionValue, JSONValue


class ReagentMetadata(TypedDict, total=False):
    """Calculated metadata attached to a recommended reagent."""

    name: str
    chemical_formula: str
    mass: float


class MoleculeIdentity(TypedDict):
    """Required molecule identity fields used by metadata helpers."""

    smiles: str


class MoleculeRecord(MoleculeIdentity, total=False):
    """Molecule record used by metadata helpers."""

    name: str
    reagent_metadata: ReagentMetadata


class ConditionsPayload(TypedDict):
    """Reaction condition fields required by the metadata workflow."""

    temperature: ConditionValue
    pressure: ConditionValue
    solvent: ConditionValue
    time: ConditionValue


class MetadataError(TypedDict):
    """Structured error payload returned by metadata orchestration."""

    stage: Literal["parse", "reagents", "conditions", "literature"]
    error: str


class MetadataRecommendation(TypedDict):
    """Successful three-stage metadata recommendation payload."""

    reaction_smiles: str
    reactants: list[MoleculeRecord]
    product: list[MoleculeRecord]
    reagents: list[MoleculeRecord]
    conditions: ConditionsPayload
    literature: JSONValue


MetadataResponse: TypeAlias = dict[str, JSONValue]
StatusPayload: TypeAlias = tuple[int, object]
MetadataResponseStatusPayload: TypeAlias = tuple[int, MetadataResponse | str]
ReagentStatusPayload: TypeAlias = tuple[int, list[MoleculeRecord] | str]
ConditionsStatusPayload: TypeAlias = tuple[int, ConditionsPayload | str]
LiteratureStatusPayload: TypeAlias = tuple[int, JSONValue | str]
RecommendationStatusPayload: TypeAlias = tuple[
    int, MetadataRecommendation | MetadataError
]

ReagentRecommender = Callable[
    [list[MoleculeRecord], list[MoleculeRecord], str, float],
    ReagentStatusPayload,
]
ConditionsRecommender = Callable[
    [list[MoleculeRecord], list[MoleculeRecord], list[MoleculeRecord], str, float],
    ConditionsStatusPayload,
]
LiteratureRecommender = Callable[
    [
        list[MoleculeRecord],
        list[MoleculeRecord],
        list[MoleculeRecord],
        ConditionsPayload,
        str,
        float,
    ],
    LiteratureStatusPayload,
]


@dataclass(frozen=True)
class ReactionParticipants:
    """Structured participants parsed from a reaction SMILES string."""

    reactants: list[MoleculeRecord]
    reagents: list[MoleculeRecord]
    product: list[MoleculeRecord]
