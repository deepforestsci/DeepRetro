"""LLM-assisted reaction metadata helpers."""

from __future__ import annotations

import ast
import json
from typing import TypeGuard, cast

import litellm
import structlog
from litellm import completion

from deepretro.utils.cache import CacheManager, make_cache_key
from deepretro.utils.langfuse_config import get_langfuse_metadata
from deepretro.utils.llm_helpers import ChatMessage, build_completion_params
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

logger = structlog.get_logger(__name__)

_litellm_initialized = False


def _ensure_litellm_configured() -> None:
    # Deferred so importing this module doesn't trigger side-effects at load time.
    global _litellm_initialized
    if _litellm_initialized:
        return
    from dotenv import load_dotenv

    load_dotenv()
    litellm.success_callback = ["langfuse"]
    litellm.drop_params = True
    _litellm_initialized = True


DEFAULT_METADATA_MODEL = "claude-opus-4-20250514"
MAX_COMPLETION_TOKENS = 4096
TOP_P = 0.9
metadata_cache = CacheManager()

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
    "call_metadata_llm",
    "conditions_agent",
    "conditions_llm_call",
    "extract_smiles",
    "literature_agent",
    "metadata_cache",
    "parse_metadata_response",
    "parse_literature_reaction",
    "parse_reaction_smiles",
    "reagent_agent",
    "reagent_llm_call",
    "recommend_reaction_metadata",
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


def call_metadata_llm(
    messages: list[ChatMessage],
    model: str,
    temperature: float,
) -> tuple[int, str]:
    """Call LiteLLM with Langfuse metadata, then retry without it on failure.

    Parameters
    ----------
    messages : list[ChatMessage]
        Chat messages sent to LiteLLM.
    model : str
        LiteLLM model identifier.
    temperature : float
        Sampling temperature.

    Returns
    -------
    tuple[int, str]
        ``(200, response_text)`` on success or ``(404, "")`` after both
        attempts fail.

    Examples
    --------
    >>> call_metadata_llm([{"role": "user", "content": "Return {}"}], "claude-opus-4-20250514", 0.0)  # doctest: +SKIP
    (200, '{}')
    """
    _ensure_litellm_configured()

    metadata: dict[str, JSONValue] | None
    try:
        metadata = get_langfuse_metadata("metadata")
    except Exception as exc:
        logger.info("metadata.langfuse_metadata_unavailable", error=str(exc))
        metadata = None

    params = build_completion_params(
        model=model,
        messages=messages,
        max_completion_tokens=MAX_COMPLETION_TOKENS,
        temperature=temperature,
        enable_thinking=False,
        metadata=metadata,
    )
    params["top_p"] = TOP_P

    try:
        response = completion(**params)
        return 200, str(response.choices[0].message.content)
    except (litellm.AuthenticationError, litellm.PermissionDeniedError) as exc:
        logger.warning("metadata.llm_auth_failed", model=model, error=str(exc))
        return 404, ""
    except litellm.APIError as exc:
        logger.info("metadata.llm_call_failed", model=model, error=str(exc))

    retry_params = dict(params)
    retry_params.pop("metadata", None)
    try:
        response = completion(**retry_params)
        return 200, str(response.choices[0].message.content)
    except litellm.APIError as exc:
        logger.info("metadata.llm_retry_failed", model=model, error=str(exc))
        return 404, ""


def reagent_agent(
    reactants: list[MoleculeRecord],
    product: list[MoleculeRecord],
    model: str = DEFAULT_METADATA_MODEL,
    temperature: float = 0.0,
) -> ReagentStatusPayload:
    """Predict reaction reagents and attach lightweight molecule metadata.

    Parameters
    ----------
    reactants : list[MoleculeRecord]
        Reactant molecule records containing a ``"smiles"`` key.
    product : list[MoleculeRecord]
        Product molecule records containing a ``"smiles"`` key.
    model : str, optional
        LiteLLM model identifier.
    temperature : float, optional
        Sampling temperature.

    Returns
    -------
    ReagentStatusPayload
        ``(200, reagents)`` on success, otherwise ``(404, "")`` or the LLM
        status code and an empty payload.
    """
    status, result = reagent_llm_call(
        extract_smiles(reactants), extract_smiles(product)[0], model, temperature
    )
    if status != 200:
        return status, ""
    if not isinstance(result, dict):
        return 404, ""

    try:
        data = result["data"]
        if not isinstance(data, list):
            raise ValueError("reagent response data must be a list")
        reagent_smiles = [
            smiles
            for smiles in data
            if isinstance(smiles, str) and is_valid_smiles(smiles)
        ]
    except Exception as exc:
        logger.info("metadata.reagent_parse_failed", error=str(exc))
        return 404, ""

    if not reagent_smiles:
        return 404, ""

    try:
        reagents = build_reagent_records(reagent_smiles)
    except Exception as exc:
        logger.info("metadata.reagent_enrichment_failed", error=str(exc))
        return 404, ""

    return 200, reagents


def reagent_llm_call(
    reactants: list[str],
    product: str,
    model: str = DEFAULT_METADATA_MODEL,
    temperature: float = 0.0,
) -> MetadataResponseStatusPayload:
    """Call an LLM to predict reaction reagent SMILES.

    Parameters
    ----------
    reactants : list[str]
        Reactant SMILES strings.
    product : str
        Product SMILES string.
    model : str, optional
        LiteLLM model identifier.
    temperature : float, optional
        Sampling temperature.

    Returns
    -------
    MetadataResponseStatusPayload
        ``(200, parsed_response)`` on success or ``(404, "")`` on failure.

    Examples
    --------
    >>> reagent_llm_call(["CCO"], "CC=O")  # doctest: +SKIP
    (200, {'data': [...]})
    """
    messages = build_reagent_messages(reactants, product)

    status, response_text = call_metadata_llm(messages, model, temperature)
    if status != 200:
        return status, ""

    try:
        return 200, parse_metadata_response(response_text)
    except Exception as exc:
        logger.info("metadata.reagent_literal_parse_failed", error=str(exc))
        return 404, ""


def conditions_agent(
    reactants: list[MoleculeRecord],
    product: list[MoleculeRecord],
    reagents: list[MoleculeRecord],
    model: str = DEFAULT_METADATA_MODEL,
    temperature: float = 0.0,
) -> ConditionsStatusPayload:
    """Predict reaction conditions for a retrosynthetic step.

    Parameters
    ----------
    reactants : list[MoleculeRecord]
        Reactant molecule records containing a ``"smiles"`` key.
    product : list[MoleculeRecord]
        Product molecule records containing a ``"smiles"`` key.
    reagents : list[MoleculeRecord]
        Reagent molecule records containing a ``"smiles"`` key.
    model : str, optional
        LiteLLM model identifier.
    temperature : float, optional
        Sampling temperature.

    Returns
    -------
    ConditionsStatusPayload
        ``(200, conditions)`` on success or ``(404, "")`` on failure.

    Examples
    --------
    >>> conditions_agent([{"smiles": "CCO"}], [{"smiles": "CC=O"}], [{"smiles": "O"}])  # doctest: +SKIP
    (200, {'temperature': ...})
    """
    status, result = conditions_llm_call(
        extract_smiles(reactants),
        extract_smiles(product)[0],
        extract_smiles(reagents),
        model,
        temperature,
    )
    if status != 200:
        return status, ""

    if not valid_conditions_payload(result):
        fields = sorted(result) if isinstance(result, dict) else []
        logger.info("metadata.conditions_incomplete", fields=fields)
        return 404, ""
    return 200, result


def conditions_llm_call(
    reactants: list[str],
    product: str,
    reagents: list[str],
    model: str = DEFAULT_METADATA_MODEL,
    temperature: float = 0.0,
) -> MetadataResponseStatusPayload:
    """Call an LLM to predict reaction conditions.

    Parameters
    ----------
    reactants : list[str]
        Reactant SMILES strings.
    product : str
        Product SMILES string.
    reagents : list[str]
        Reagent SMILES strings.
    model : str, optional
        LiteLLM model identifier.
    temperature : float, optional
        Sampling temperature.

    Returns
    -------
    MetadataResponseStatusPayload
        ``(200, parsed_response)`` on success or ``(404, "")`` on failure.

    Examples
    --------
    >>> conditions_llm_call(["CCO"], "CC=O", ["O"])  # doctest: +SKIP
    (200, {'temperature': ...})
    """
    messages = build_conditions_messages(reactants, product, reagents)

    status, response_text = call_metadata_llm(messages, model, temperature)
    if status != 200:
        return status, ""

    try:
        return 200, parse_metadata_response(response_text)
    except Exception as exc:
        logger.info("metadata.conditions_literal_parse_failed", error=str(exc))
        return 404, ""


def literature_agent(
    reactants: list[MoleculeRecord],
    product: list[MoleculeRecord],
    reagents: list[MoleculeRecord],
    conditions: ConditionsPayload,
    model: str = DEFAULT_METADATA_MODEL,
    temperature: float = 0.0,
) -> LiteratureStatusPayload:
    """Find the closest literature reaction for a retrosynthetic step.

    Parameters
    ----------
    reactants : list[MoleculeRecord]
        Reactant molecule records containing a ``"smiles"`` key.
    product : list[MoleculeRecord]
        Product molecule records containing a ``"smiles"`` key.
    reagents : list[MoleculeRecord]
        Reagent molecule records containing a ``"smiles"`` key.
    conditions : ConditionsPayload
        Reaction conditions payload.
    model : str, optional
        LiteLLM model identifier.
    temperature : float, optional
        Sampling temperature.

    Returns
    -------
    LiteratureStatusPayload
        ``(200, literature_reaction)`` on success or ``(404, "")`` on failure.

    Examples
    --------
    >>> literature_agent([{"smiles": "CCO"}], [{"smiles": "CC=O"}], [{"smiles": "O"}], {"temperature": "25 C"})  # doctest: +SKIP
    (200, ...)
    """
    messages = build_literature_messages(
        extract_smiles(reactants),
        extract_smiles(product)[0],
        extract_smiles(reagents),
        conditions,
    )

    status, response_text = call_metadata_llm(messages, model, temperature)
    if status != 200:
        return status, ""

    try:
        return 200, parse_literature_reaction(response_text)
    except Exception as exc:
        logger.info("metadata.literature_parse_failed", error=str(exc))
        return 404, ""


def recommend_reaction_metadata(
    reaction_smiles: str,
    model: str = DEFAULT_METADATA_MODEL,
    temperature: float = 0.0,
    reagent_recommender: ReagentRecommender | None = None,
    conditions_recommender: ConditionsRecommender | None = None,
    literature_recommender: LiteratureRecommender | None = None,
    cache: CacheManager | None = metadata_cache,
) -> RecommendationStatusPayload:
    """Recommend reagents, conditions, and literature for one reaction string.

    Parameters
    ----------
    reaction_smiles : str
        Reaction SMILES in ``reactants>agents>products`` or
        ``reactants>>products`` form.
    model : str, optional
        LiteLLM model identifier used by the default recommenders.
    temperature : float, optional
        Sampling temperature.
    reagent_recommender : callable, optional
        Alternate reagent recommendation function. Intended for tests and
        controlled integrations.
    conditions_recommender : callable, optional
        Alternate conditions recommendation function.
    literature_recommender : callable, optional
        Alternate literature recommendation function.
    cache : CacheManager | None, optional
        Public cache instance used for the complete recommendation. Pass
        ``None`` to disable this top-level cache.

    Returns
    -------
    RecommendationStatusPayload
        ``(200, recommendation)`` on success. The recommendation contains the
        parsed reaction participants plus ``"reagents"``, ``"conditions"``, and
        ``"literature"``. Failure returns the failing status and a small error
        payload identifying the stage.

    Examples
    --------
    >>> status, result = recommend_reaction_metadata("CCO>>CC=O")  # doctest: +SKIP
    >>> status  # doctest: +SKIP
    200
    """
    active_cache = cache
    use_cache = (
        active_cache is not None
        and reagent_recommender is None
        and conditions_recommender is None
        and literature_recommender is None
    )
    if use_cache:
        assert active_cache is not None
        cache_key = make_cache_key(
            "recommend_reaction_metadata",
            reaction_smiles,
            model=model,
            temperature=temperature,
        )
        miss = object()
        cached = active_cache.get(cache_key, default=miss)
        if cached is not miss:
            return cast(RecommendationStatusPayload, cached)

    try:
        participants = parse_reaction_smiles(reaction_smiles)
    except ValueError as exc:
        return 400, {"stage": "parse", "error": str(exc)}

    reagent_recommender = (
        reagent_agent if reagent_recommender is None else reagent_recommender
    )
    conditions_recommender = (
        conditions_agent if conditions_recommender is None else conditions_recommender
    )
    literature_recommender = (
        literature_agent if literature_recommender is None else literature_recommender
    )

    reagent_status, reagents = reagent_recommender(
        participants.reactants,
        participants.product,
        model,
        temperature,
    )
    if reagent_status != 200:
        return reagent_status, {
            "stage": "reagents",
            "error": "metadata recommendation failed",
        }
    if not isinstance(reagents, list):
        return 404, {
            "stage": "reagents",
            "error": "metadata recommendation failed",
        }

    conditions_status, conditions = conditions_recommender(
        participants.reactants,
        participants.product,
        reagents,
        model,
        temperature,
    )
    if conditions_status != 200:
        return conditions_status, {
            "stage": "conditions",
            "error": "metadata recommendation failed",
        }
    if not valid_conditions_payload(conditions):
        return 404, {
            "stage": "conditions",
            "error": "metadata recommendation failed",
        }

    literature_status, literature = literature_recommender(
        participants.reactants,
        participants.product,
        reagents,
        conditions,
        model,
        temperature,
    )
    if literature_status != 200:
        return literature_status, {
            "stage": "literature",
            "error": "metadata recommendation failed",
        }

    recommendation: MetadataRecommendation = {
        "reaction_smiles": reaction_smiles,
        "reactants": participants.reactants,
        "product": participants.product,
        "reagents": reagents,
        "conditions": conditions,
        "literature": literature,
    }
    result: RecommendationStatusPayload = (200, recommendation)
    if use_cache:
        assert active_cache is not None
        active_cache.set(cache_key, result, tag="metadata")
    return result
