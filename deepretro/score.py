from __future__ import annotations

import contextlib
import importlib.util
import io
import os
import pickle
import sys
import warnings
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any

from deepretro._pathway import (
    _extract_step_smiles,
    _get_pathway_steps,
    _split_reactants,
)
from deepretro.ursa import UrsaAdapter, UrsaConfig


@dataclass(frozen=True)
class SCScoreConfig:
    module_path: str | None = None
    weight_path: str | None = None
    fp_len: int = 1024
    fp_rad: int = 2


def canonicalize_smiles(smiles: str) -> str | None:
    """Return canonical isomeric SMILES, or None for invalid input/unavailable RDKit."""
    text = str(smiles).strip() if smiles is not None else ""
    if not text:
        return None

    chem = _get_chem()
    if chem is None:
        return None

    try:
        with _quiet_rdkit_logs():
            mol = chem.MolFromSmiles(text)
    except Exception:
        return None
    if mol is None:
        return None
    return str(chem.MolToSmiles(mol, isomericSmiles=True))


def sa_score(smiles: str) -> float | None:
    """Return RDKit Contrib SA_Score for a SMILES string, or None if unavailable/invalid."""
    value, _error = _sa_score_with_error(smiles)
    return value


def sc_score(smiles: str, config: SCScoreConfig | None = None) -> float | None:
    """Return optional SCScore for a SMILES string, or None if unavailable/invalid."""
    value, _error = _sc_score_with_error(smiles, config)
    return value


def score_molecule(
    smiles: str, sc_config: SCScoreConfig | None = None
) -> dict[str, Any]:
    """Score one molecule; unavailable optional metrics are None with warnings."""
    warning_messages: list[str] = []
    error_messages: list[str] = []
    canonical = canonicalize_smiles(smiles)

    if canonical is None:
        if _get_chem() is None:
            warning_messages.extend(
                ["RDKit unavailable", "SA_Score unavailable", "SCScore unavailable"]
            )
            error_messages.append("RDKit unavailable")
        else:
            error_messages.append("Invalid SMILES")
        return {
            "input_smiles": smiles,
            "canonical_smiles": None,
            "valid": False,
            "sa_score": None,
            "sc_score": None,
            "warnings": _dedupe(warning_messages),
            "errors": _dedupe(error_messages),
        }

    sa_value, sa_error = _sa_score_with_error(canonical)
    if sa_value is None:
        warning_messages.append(sa_error or "SA_Score unavailable")

    sc_value, sc_error = _sc_score_with_error(canonical, sc_config)
    if sc_value is None:
        warning_messages.append(sc_error or "SCScore unavailable")

    return {
        "input_smiles": smiles,
        "canonical_smiles": canonical,
        "valid": True,
        "sa_score": sa_value,
        "sc_score": sc_value,
        "warnings": _dedupe(warning_messages),
        "errors": error_messages,
    }


def score_step(
    product_smiles: str,
    reactants_smiles: str | list[str],
    sc_config: SCScoreConfig | None = None,
) -> dict[str, Any]:
    """Score a retrosynthesis step with dot-separated or list reactant SMILES."""
    reactant_inputs = _split_reactants(reactants_smiles)
    if not reactant_inputs:
        return _malformed_step_score("Malformed step: missing reactant SMILES")

    product = score_molecule(product_smiles, sc_config)
    reactants = [score_molecule(smiles, sc_config) for smiles in reactant_inputs]

    reactant_sa = [item["sa_score"] for item in reactants]
    reactant_sc = [item["sc_score"] for item in reactants]
    mean_reactant_sa = _mean(reactant_sa)
    max_reactant_sa = _max_or_none(reactant_sa)
    min_reactant_sa = _min_or_none(reactant_sa)
    mean_reactant_sc = _mean(reactant_sc)
    max_reactant_sc = _max_or_none(reactant_sc)
    min_reactant_sc = _min_or_none(reactant_sc)

    delta_sa_mean = _delta(product["sa_score"], mean_reactant_sa)
    delta_sa_max = _delta(product["sa_score"], max_reactant_sa)
    delta_sc_mean = _delta(product["sc_score"], mean_reactant_sc)
    delta_sc_max = _delta(product["sc_score"], max_reactant_sc)

    warning_messages = list(product["warnings"])
    error_messages = list(product["errors"])
    for reactant in reactants:
        warning_messages.extend(reactant["warnings"])
        error_messages.extend(reactant["errors"])

    return {
        "product": product,
        "reactants": reactants,
        "mean_reactant_sa": mean_reactant_sa,
        "max_reactant_sa": max_reactant_sa,
        "min_reactant_sa": min_reactant_sa,
        "mean_reactant_sc": mean_reactant_sc,
        "max_reactant_sc": max_reactant_sc,
        "min_reactant_sc": min_reactant_sc,
        "delta_sa_mean": delta_sa_mean,
        "delta_sa_max": delta_sa_max,
        "delta_sc_mean": delta_sc_mean,
        "delta_sc_max": delta_sc_max,
        "sa_simplifies": delta_sa_mean > 0 if delta_sa_mean is not None else False,
        "sc_simplifies": delta_sc_mean > 0 if delta_sc_mean is not None else False,
        "warnings": _dedupe(warning_messages),
        "errors": _dedupe(error_messages),
    }


def score_pathway(
    pathway: dict[str, Any] | list[dict[str, Any]],
    sc_config: SCScoreConfig | None = None,
    *,
    ursa_config: UrsaConfig | None = None,
) -> dict[str, Any]:
    """Score viewer-style {'steps': [...]} pathways or a raw list of step dicts."""
    steps, pathway_error = _get_pathway_steps(pathway)
    if pathway_error is not None:
        return _pathway_result([], [pathway_error])

    step_scores: list[dict[str, Any]] = []
    warning_messages: list[str] = []
    error_messages: list[str] = []

    for index, step in enumerate(steps):
        try:
            product_smiles, reactant_smiles = _extract_step_smiles(step)
            step_score = score_step(product_smiles, reactant_smiles, sc_config)
        except ValueError as exc:
            step_score = _malformed_step_score(f"Malformed step {index}: {exc}")
        step_scores.append(step_score)
        warning_messages.extend(step_score["warnings"])
        error_messages.extend(step_score["errors"])

    result = {
        "step_scores": step_scores,
        "route_summary": _summarize_steps(step_scores),
        "warnings": _dedupe(warning_messages),
        "errors": _dedupe(error_messages),
    }
    if (
        ursa_config is not None
        and ursa_config.enabled
        and steps
        and not result["errors"]
    ):
        result["ursa"] = UrsaAdapter(ursa_config).annotate(pathway)
    return result


def empty_pathway_scores(error: str | None = None) -> dict[str, Any]:
    """Return an empty pathway score payload using the standard schema."""
    errors = [error] if error else []
    return _pathway_result([], errors)


@lru_cache(maxsize=1)
def _get_chem() -> Any | None:
    try:
        with _quiet_output():
            from rdkit import Chem
    except Exception:
        return None
    return Chem


@lru_cache(maxsize=1)
def _get_sascorer() -> Any | None:
    try:
        from rdkit.Contrib.SA_Score import sascorer

        return sascorer
    except Exception:
        pass

    try:
        contrib_path = os.path.join(
            os.environ["CONDA_PREFIX"], "share", "RDKit", "Contrib"
        )
        if contrib_path not in sys.path:
            sys.path.append(contrib_path)
        from SA_Score import sascorer

        return sascorer
    except Exception:
        return None


def _sa_warning() -> str:
    if _get_chem() is None:
        return "RDKit unavailable"
    if _get_sascorer() is None:
        return "SA_Score unavailable"
    return "SA_Score failed"


def _sa_score_with_error(smiles: str) -> tuple[float | None, str | None]:
    canonical = canonicalize_smiles(smiles)
    if canonical is None:
        return (
            None,
            "Invalid SMILES" if _get_chem() is not None else "RDKit unavailable",
        )

    chem = _get_chem()
    if chem is None:
        return None, "RDKit unavailable"

    scorer = _get_sascorer()
    if scorer is None:
        return None, "SA_Score unavailable"

    try:
        mol = chem.MolFromSmiles(canonical)
        if mol is None:
            return None, "Invalid SMILES"
        return float(scorer.calculateScore(mol)), None
    except Exception as exc:
        return None, _format_error("SA_Score failed", exc)


def _sc_score_with_error(
    smiles: str, config: SCScoreConfig | None
) -> tuple[float | None, str | None]:
    canonical = canonicalize_smiles(smiles)
    if canonical is None:
        return (
            None,
            "Invalid SMILES" if _get_chem() is not None else "RDKit unavailable",
        )

    cfg = config or SCScoreConfig()
    module_path, weight_path = _sc_paths(cfg)
    scorer, load_error = _get_scscorer(module_path, weight_path, cfg.fp_len, cfg.fp_rad)
    if scorer is None:
        return None, load_error or "SCScore unavailable"

    try:
        value = scorer.get_score_from_smi(canonical)
    except Exception as exc:
        return None, _format_error("SCScore scoring failed", exc)
    if isinstance(value, (list, tuple)):
        value = value[1] if len(value) > 1 else value[0]
    try:
        return float(value), None
    except (TypeError, ValueError) as exc:
        return None, _format_error("SCScore scoring failed", exc)


@lru_cache(maxsize=8)
def _get_scscorer(
    module_path: str | None,
    weight_path: str | None,
    fp_len: int,
    fp_rad: int,
) -> tuple[Any | None, str | None]:
    _install_scscore_shims()
    try:
        scorer_class = _load_scscorer_class(module_path)
        scorer = scorer_class()
        with _quiet_output():
            if weight_path is None:
                scorer.restore(FP_rad=fp_rad, FP_len=fp_len)
            else:
                scorer.restore(weight_path, FP_rad=fp_rad, FP_len=fp_len)
        return scorer, None
    except Exception as exc:
        return None, _format_error("SCScore unavailable", exc)


def _install_scscore_shims() -> None:
    sys.modules.setdefault("cPickle", pickle)
    try:
        import numpy
    except Exception:
        return
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FutureWarning)
        has_numpy_bool = hasattr(numpy, "bool")
    if not has_numpy_bool:
        setattr(numpy, "bool", numpy.bool_)


def _load_scscorer_class(module_path: str | None) -> Any:
    if module_path is None:
        from scscore.standalone_model_numpy import SCScorer

        return SCScorer

    spec = importlib.util.spec_from_file_location(
        "_deepretro_scscore_standalone", module_path
    )
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not load SCScore module from {module_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.SCScorer


def _sc_paths(config: SCScoreConfig) -> tuple[str | None, str | None]:
    # Fall back to env vars so the SCScore model can be pointed at without
    # changing call sites: SCSCORE_MODULE_PATH (standalone_model_numpy.py) and
    # SCSCORE_WEIGHT_PATH (model.ckpt-*.as_numpy.json.gz).
    raw_module = config.module_path or os.getenv("SCSCORE_MODULE_PATH")
    raw_weight = config.weight_path or os.getenv("SCSCORE_WEIGHT_PATH")
    module_path = str(Path(raw_module).expanduser().resolve()) if raw_module else None
    weight_path = str(Path(raw_weight).expanduser().resolve()) if raw_weight else None
    return module_path, weight_path


@contextlib.contextmanager
def _quiet_output() -> Any:
    with (
        contextlib.redirect_stdout(io.StringIO()),
        contextlib.redirect_stderr(io.StringIO()),
    ):
        yield


@contextlib.contextmanager
def _quiet_rdkit_logs() -> Any:
    try:
        from rdkit import RDLogger

        RDLogger.DisableLog("rdApp.*")
    except Exception:
        RDLogger = None
    try:
        yield
    finally:
        if RDLogger is not None:
            RDLogger.EnableLog("rdApp.*")


def _pathway_result(
    step_scores: list[dict[str, Any]], error_messages: list[str]
) -> dict[str, Any]:
    return {
        "step_scores": step_scores,
        "route_summary": _summarize_steps(step_scores),
        "warnings": [],
        "errors": _dedupe(error_messages),
    }


def _malformed_step_score(message: str) -> dict[str, Any]:
    return {
        "product": None,
        "reactants": [],
        "mean_reactant_sa": None,
        "max_reactant_sa": None,
        "min_reactant_sa": None,
        "mean_reactant_sc": None,
        "max_reactant_sc": None,
        "min_reactant_sc": None,
        "delta_sa_mean": None,
        "delta_sa_max": None,
        "delta_sc_mean": None,
        "delta_sc_max": None,
        "sa_simplifies": False,
        "sc_simplifies": False,
        "warnings": [],
        "errors": [message],
    }


def _summarize_steps(step_scores: list[dict[str, Any]]) -> dict[str, Any]:
    delta_sa = [step.get("delta_sa_mean") for step in step_scores]
    delta_sc = [step.get("delta_sc_mean") for step in step_scores]
    complexity_increasing_steps = [
        index
        for index, step in enumerate(step_scores)
        if _is_negative(step.get("delta_sa_mean"))
        or _is_negative(step.get("delta_sc_mean"))
    ]

    return {
        "n_steps": len(step_scores),
        "mean_product_sa": _mean(_product_scores(step_scores, "sa_score")),
        "mean_product_sc": _mean(_product_scores(step_scores, "sc_score")),
        "mean_delta_sa": _mean(delta_sa),
        "mean_delta_sc": _mean(delta_sc),
        "fraction_steps_sa_simplifying": _fraction_simplifying(
            step_scores, "delta_sa_mean", "sa_simplifies"
        ),
        "fraction_steps_sc_simplifying": _fraction_simplifying(
            step_scores, "delta_sc_mean", "sc_simplifies"
        ),
        "complexity_increasing_steps": complexity_increasing_steps,
    }


def _product_scores(step_scores: list[dict[str, Any]], key: str) -> list[float | None]:
    return [
        step["product"][key]
        for step in step_scores
        if isinstance(step.get("product"), dict)
    ]


def _mean(values: list[float | None]) -> float | None:
    real_values = [float(value) for value in values if value is not None]
    return sum(real_values) / len(real_values) if real_values else None


def _max_or_none(values: list[float | None]) -> float | None:
    real_values = [float(value) for value in values if value is not None]
    return max(real_values) if real_values else None


def _min_or_none(values: list[float | None]) -> float | None:
    real_values = [float(value) for value in values if value is not None]
    return min(real_values) if real_values else None


def _delta(left: float | None, right: float | None) -> float | None:
    return None if left is None or right is None else float(left) - float(right)


def _fraction_simplifying(
    step_scores: list[dict[str, Any]], delta_key: str, simplify_key: str
) -> float | None:
    scored_steps = [step for step in step_scores if step.get(delta_key) is not None]
    if not scored_steps:
        return None
    return sum(1 for step in scored_steps if step.get(simplify_key) is True) / len(
        scored_steps
    )


def _is_negative(value: Any) -> bool:
    return value is not None and float(value) < 0


def _format_error(prefix: str, exc: Exception) -> str:
    message = str(exc).strip() or exc.__class__.__name__
    return f"{prefix}: {message}"


def _dedupe(messages: list[str]) -> list[str]:
    seen: set[str] = set()
    deduped: list[str] = []
    for message in messages:
        if message not in seen:
            deduped.append(message)
            seen.add(message)
    return deduped
