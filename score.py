from __future__ import annotations

import importlib.util
import contextlib
import io
import os
import pickle
import sys
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any


_CHEM: Any | None = None
_CHEM_ATTEMPTED = False
_SASCORER: Any | None = None
_SASCORE_ATTEMPTED = False
_SCORER_CACHE: dict[tuple[str | None, str | None, int, int], Any | None] = {}
_SCORER_ERRORS: dict[tuple[str | None, str | None, int, int], str] = {}


@dataclass(frozen=True)
class SCScoreConfig:
    module_path: str | None = None
    weight_path: str | None = None
    fp_len: int = 1024
    fp_rad: int = 2


def canonicalize_smiles(smiles: str) -> str | None:
    text = str(smiles).strip() if smiles is not None else ""
    if not text:
        return None

    chem = _get_chem()
    if chem is None:
        return None

    try:
        with _quiet_optional_output():
            mol = chem.MolFromSmiles(text)
            if mol is None:
                return None
            return str(chem.MolToSmiles(mol, isomericSmiles=True))
    except Exception:
        return None


def sa_score(smiles: str) -> float | None:
    text = canonicalize_smiles(smiles)
    if text is None:
        return None

    chem = _get_chem()
    scorer = _get_sascorer()
    if chem is None or scorer is None:
        return None

    try:
        with _quiet_optional_output():
            mol = chem.MolFromSmiles(text)
            if mol is None:
                return None
            return float(scorer.calculateScore(mol))
    except Exception:
        return None


def sc_score(smiles: str, config: SCScoreConfig | None = None) -> float | None:
    text = canonicalize_smiles(smiles)
    if text is None:
        return None

    cfg = config or SCScoreConfig()
    key = _sc_cache_key(cfg)
    scorer = _get_scscorer(cfg)
    if scorer is None:
        return None

    try:
        with _quiet_optional_output():
            score = scorer.get_score_from_smi(text)
        if isinstance(score, (list, tuple)):
            score = score[1] if len(score) > 1 else score[0]
        _SCORER_ERRORS.pop(key, None)
        return float(score)
    except Exception as exc:
        _SCORER_ERRORS[key] = _short_error("SCScore scoring failed", exc)
        return None


def score_molecule(smiles: str, sc_config: SCScoreConfig | None = None) -> dict[str, Any]:
    canonical = canonicalize_smiles(smiles)
    warnings: list[str] = []
    errors: list[str] = []

    if canonical is None:
        if _get_chem() is None:
            warnings.extend(["RDKit unavailable", "SAscore unavailable", "SCScore unavailable"])
            errors.append("RDKit unavailable")
        else:
            errors.append("Invalid SMILES")
        return {
            "input_smiles": smiles,
            "canonical_smiles": None,
            "valid": False,
            "sa_score": None,
            "sc_score": None,
            "warnings": _dedupe(warnings),
            "errors": _dedupe(errors),
        }

    sa_value = sa_score(canonical)
    if sa_value is None:
        warnings.append("SAscore unavailable")

    sc_value = sc_score(canonical, sc_config)
    if sc_value is None:
        warnings.append(_scscore_warning(sc_config))

    return {
        "input_smiles": smiles,
        "canonical_smiles": canonical,
        "valid": True,
        "sa_score": sa_value,
        "sc_score": sc_value,
        "warnings": _dedupe(warnings),
        "errors": errors,
    }


def score_step(
    product_smiles: str,
    reactants_smiles: str | list[str],
    sc_config: SCScoreConfig | None = None,
) -> dict[str, Any]:
    product = score_molecule(product_smiles, sc_config)
    reactants = [score_molecule(smiles, sc_config) for smiles in _split_reactants(reactants_smiles)]

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

    warnings = list(product["warnings"])
    errors = list(product["errors"])
    for reactant in reactants:
        warnings.extend(reactant["warnings"])
        errors.extend(reactant["errors"])

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
        "warnings": _dedupe(warnings),
        "errors": _dedupe(errors),
    }


def score_pathway(pathway: dict[str, Any] | list[dict[str, Any]], sc_config: SCScoreConfig | None = None) -> dict[str, Any]:
    steps = pathway.get("steps", []) if isinstance(pathway, dict) else pathway
    if not isinstance(steps, list):
        steps = []

    step_scores: list[dict[str, Any]] = []
    warnings: list[str] = []
    errors: list[str] = []

    for index, step in enumerate(steps):
        try:
            product_smiles, reactant_smiles = _extract_step_smiles(step)
            step_score = score_step(product_smiles, reactant_smiles, sc_config)
        except Exception as exc:
            step_score = _malformed_step_score(f"Malformed step {index}: {exc}")
        step_scores.append(step_score)
        warnings.extend(step_score["warnings"])
        errors.extend(step_score["errors"])

    route_summary = _summarize_steps(step_scores)
    return {
        "step_scores": step_scores,
        "route_summary": route_summary,
        "warnings": _dedupe(warnings),
        "errors": _dedupe(errors),
    }


def _get_chem() -> Any | None:
    global _CHEM, _CHEM_ATTEMPTED
    if _CHEM_ATTEMPTED:
        return _CHEM

    _CHEM_ATTEMPTED = True
    try:
        with _quiet_optional_output():
            from rdkit import Chem

        _CHEM = Chem
        return _CHEM
    except Exception:
        _CHEM = None
        return None


def _get_sascorer() -> Any | None:
    global _SASCORER, _SASCORE_ATTEMPTED
    if _SASCORE_ATTEMPTED:
        return _SASCORER

    _SASCORE_ATTEMPTED = True
    try:
        with _quiet_optional_output():
            from rdkit.Contrib.SA_Score import sascorer

        _SASCORER = sascorer
        return _SASCORER
    except Exception:
        pass

    try:
        contrib_path = os.path.join(os.environ["CONDA_PREFIX"], "share", "RDKit", "Contrib")
        if contrib_path not in sys.path:
            sys.path.append(contrib_path)
        with _quiet_optional_output():
            from SA_Score import sascorer

        _SASCORER = sascorer
        return _SASCORER
    except Exception:
        _SASCORER = None
        return None


def _get_scscorer(config: SCScoreConfig) -> Any | None:
    key = _sc_cache_key(config)
    if key in _SCORER_CACHE:
        return _SCORER_CACHE[key]

    try:
        _install_scscore_compatibility_shims()
        with _quiet_optional_output():
            scorer_class = _load_scscorer_class(config)
            scorer = scorer_class()
            if key[1] is not None:
                scorer.restore(key[1], FP_rad=config.fp_rad, FP_len=config.fp_len)
            else:
                scorer.restore(FP_rad=config.fp_rad, FP_len=config.fp_len)
        _SCORER_CACHE[key] = scorer
        _SCORER_ERRORS.pop(key, None)
        return scorer
    except Exception as exc:
        _SCORER_CACHE[key] = None
        _SCORER_ERRORS[key] = _short_error("SCScore unavailable", exc)
        return None


def _install_scscore_compatibility_shims() -> None:
    try:
        import numpy

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", FutureWarning)
            has_numpy_bool = hasattr(numpy, "bool")
        if not has_numpy_bool:
            numpy.bool = numpy.bool_  # type: ignore[attr-defined]
    except Exception:
        pass
    sys.modules.setdefault("cPickle", pickle)


def _load_scscorer_class(config: SCScoreConfig) -> Any:
    if config.module_path:
        module_path = Path(config.module_path).expanduser().resolve()
        spec = importlib.util.spec_from_file_location("_deepretro_scscore_standalone", module_path)
        if spec is None or spec.loader is None:
            raise ImportError(f"Could not load SCScore module from {module_path}")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module.SCScorer

    from scscore.standalone_model_numpy import SCScorer

    return SCScorer


def _scscore_warning(config: SCScoreConfig | None) -> str:
    key = _sc_cache_key(config or SCScoreConfig())
    return _SCORER_ERRORS.get(key, "SCScore unavailable")


def _short_error(prefix: str, exc: Exception) -> str:
    message = str(exc).strip()
    if not message:
        message = exc.__class__.__name__
    return f"{prefix}: {message}"


@contextlib.contextmanager
def _quiet_optional_output() -> Any:
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
        yield


def _sc_cache_key(config: SCScoreConfig) -> tuple[str | None, str | None, int, int]:
    module_path = str(Path(config.module_path).expanduser().resolve()) if config.module_path else None
    weight_path = str(Path(config.weight_path).expanduser().resolve()) if config.weight_path else None
    return module_path, weight_path, int(config.fp_len), int(config.fp_rad)


def _split_reactants(reactants_smiles: str | list[str]) -> list[str]:
    if isinstance(reactants_smiles, str):
        parts = reactants_smiles.split(".")
    else:
        parts = [str(item) for item in reactants_smiles]
    return [part.strip() for part in parts if part.strip()]


def _extract_step_smiles(step: dict[str, Any]) -> tuple[str, list[str]]:
    if not isinstance(step, dict):
        raise ValueError("step is not a dictionary")

    if "products" in step or "reactants" in step:
        products = step.get("products")
        reactants = step.get("reactants")
        if not isinstance(products, list) or not products:
            raise ValueError("missing product SMILES")
        if not isinstance(reactants, list):
            raise ValueError("missing reactant SMILES")
        product_smiles = _smiles_from_item(products[0])
        reactant_smiles = [_smiles_from_item(item) for item in reactants]
        return product_smiles, reactant_smiles

    product_smiles = step.get("product_smiles") or step.get("product")
    reactant_smiles = step.get("reactants_smiles") or step.get("reactants")
    if product_smiles is None or reactant_smiles is None:
        raise ValueError("missing product or reactant SMILES")
    return str(product_smiles), _split_reactants(reactant_smiles)


def _smiles_from_item(item: Any) -> str:
    if isinstance(item, dict):
        value = item.get("smiles")
    else:
        value = item
    if value is None:
        raise ValueError("missing SMILES")
    return str(value)


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
    n_steps = len(step_scores)
    product_sa = [
        step["product"]["sa_score"]
        for step in step_scores
        if isinstance(step.get("product"), dict)
    ]
    product_sc = [
        step["product"]["sc_score"]
        for step in step_scores
        if isinstance(step.get("product"), dict)
    ]
    delta_sa = [step.get("delta_sa_mean") for step in step_scores]
    delta_sc = [step.get("delta_sc_mean") for step in step_scores]
    complexity_increasing_steps = [
        index
        for index, step in enumerate(step_scores)
        if (
            step.get("delta_sa_mean") is not None
            and step.get("sa_simplifies") is False
        )
        or (
            step.get("delta_sc_mean") is not None
            and step.get("sc_simplifies") is False
        )
    ]

    return {
        "n_steps": n_steps,
        "mean_product_sa": _mean(product_sa),
        "mean_product_sc": _mean(product_sc),
        "mean_delta_sa": _mean(delta_sa),
        "mean_delta_sc": _mean(delta_sc),
        "fraction_steps_sa_simplifying": _fraction_simplifying(step_scores, "sa_simplifies") if n_steps else None,
        "fraction_steps_sc_simplifying": _fraction_simplifying(step_scores, "sc_simplifies") if n_steps else None,
        "complexity_increasing_steps": complexity_increasing_steps,
    }


def _mean(values: list[float | None]) -> float | None:
    real_values = [float(value) for value in values if value is not None]
    if not real_values:
        return None
    return sum(real_values) / len(real_values)


def _max_or_none(values: list[float | None]) -> float | None:
    real_values = [float(value) for value in values if value is not None]
    return max(real_values) if real_values else None


def _min_or_none(values: list[float | None]) -> float | None:
    real_values = [float(value) for value in values if value is not None]
    return min(real_values) if real_values else None


def _delta(left: float | None, right: float | None) -> float | None:
    if left is None or right is None:
        return None
    return float(left) - float(right)


def _fraction_simplifying(step_scores: list[dict[str, Any]], key: str) -> float:
    return sum(1 for step in step_scores if step.get(key) is True) / len(step_scores)


def _dedupe(messages: list[str]) -> list[str]:
    seen: set[str] = set()
    result: list[str] = []
    for message in messages:
        if message not in seen:
            result.append(message)
            seen.add(message)
    return result
