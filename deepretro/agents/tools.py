"""Tool registry, executors, and dispatcher for the retrosynthesis agent.

Schemas use the OpenAI function-calling format
(``{"type": "function", "function": {...}}``) because LiteLLM normalizes every
provider (including Anthropic) to that shape. :func:`build_tool_registry`
assembles the tools appropriate for a ``tool_backend`` and binds
``check_hallucination`` to a resolved checker so the ML classifier can be
swapped in without touching tool code.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from typing import Any

import structlog

from deepretro.algorithms.stability_checker import check_molecule_stability
from deepretro.utils.typing import HallucinationChecker
from deepretro.utils.utils_molecule import canonicalize, is_valid_smiles, validity_check

logger = structlog.get_logger(__name__)

ToolExecutor = Callable[..., dict[str, Any]]


@dataclass(frozen=True)
class ToolRegistry:
    """Bundle of LiteLLM tool schemas and their executors.

    Attributes
    ----------
    schemas : list[dict]
        OpenAI function-format tool definitions to pass to ``completion``.
    executors : dict[str, ToolExecutor]
        Mapping from tool name to its executor callable.

    Examples
    --------
    >>> registry = build_tool_registry()
    >>> registry.execute("validate_smiles", {"smiles": "C(O)C"})["canonical_smiles"]
    'CCO'
    """

    schemas: list[dict[str, Any]]
    executors: dict[str, ToolExecutor]

    @property
    def names(self) -> list[str]:
        """Return the registered tool names."""
        return list(self.executors)

    def execute(self, name: str, arguments: dict[str, Any]) -> dict[str, Any]:
        """Dispatch a tool call, capturing unknown-tool and executor errors.

        Parameters
        ----------
        name : str
            Tool name from the model's tool call.
        arguments : dict
            Parsed tool-call arguments.

        Returns
        -------
        dict
            The tool result, or ``{"error": ...}`` on an unknown tool or an
            executor exception.

        Examples
        --------
        >>> build_tool_registry().execute("nope", {})["error"]
        'unknown tool: nope'
        """
        executor = self.executors.get(name)
        if executor is None:
            return {"error": f"unknown tool: {name}"}
        try:
            return executor(**arguments)
        except Exception as exc:  # ponytail: model-supplied args; never crash the loop
            logger.warning("Tool execution failed", tool=name, error=str(exc))
            return {"error": f"tool {name} failed: {exc}"}


def _validate_smiles(smiles: str, context: str | None = None) -> dict[str, Any]:
    """Validate a SMILES string and return its canonical form."""
    if not is_valid_smiles(smiles):
        return {
            "valid": False,
            "canonical_smiles": None,
            "error": "invalid SMILES string — cannot be parsed by RDKit",
            "context": context or "",
        }
    return {
        "valid": True,
        "canonical_smiles": canonicalize(smiles),
        "error": None,
        "context": context or "",
    }


def _check_stability(smiles: str) -> dict[str, Any]:
    """Assess heuristic molecular stability for a SMILES string."""
    report = check_molecule_stability(smiles)
    return {
        "assessment": report.get("assessment"),
        "stability_score": report.get("stability_score"),
        "issues": report.get("issues", []),
    }


def _make_check_hallucination(checker: HallucinationChecker) -> ToolExecutor:
    """Bind a resolved hallucination checker into a tool executor.

    TEMPLATE seam: the executor delegates entirely to *checker*. In heuristic
    mode this is the working package heuristic; the ML wiring is completed by
    a downstream contributor by passing a trained-classifier checker.
    """

    def _check(product: str, reactants: str | list[str]) -> dict[str, Any]:
        raw_pathway = reactants if isinstance(reactants, list) else [reactants]
        pathways, _, _ = validity_check(
            product, [raw_pathway], ["agent tool check"], [1.0]
        )
        if not pathways:
            return {
                "is_hallucination": True,
                "note": "pathway failed validity check",
            }
        status, kept = checker(product, pathways)  # use normalized pathway
        is_hallucination = status != 200 or not kept
        return {
            "is_hallucination": is_hallucination,
            "note": (
                "kept by checker" if not is_hallucination else "flagged by checker"
            ),
        }

    return _check


def _make_run_python(sandbox: Any) -> ToolExecutor:
    """Bind a sandbox into a ``run_python`` tool executor."""

    def _run(code: str) -> dict[str, Any]:
        result = sandbox.run(code)
        return {
            "ok": result.ok,
            "stdout": result.stdout,
            "stderr": result.stderr,
            "error": result.error,
        }

    return _run


_VALIDATE_SMILES_SCHEMA = {
    "type": "function",
    "function": {
        "name": "validate_smiles",
        "description": (
            "Validate a SMILES string and return its canonical form. Call this "
            "before proposing any precursor to ensure it is chemically valid."
        ),
        "parameters": {
            "type": "object",
            "properties": {
                "smiles": {"type": "string", "description": "SMILES to validate"},
                "context": {
                    "type": "string",
                    "description": "Optional note, e.g. 'reactant 1'",
                },
            },
            "required": ["smiles"],
        },
    },
}

_CHECK_STABILITY_SCHEMA = {
    "type": "function",
    "function": {
        "name": "check_stability",
        "description": (
            "Assess whether a molecule is chemically stable using heuristic "
            "structural checks. Returns an assessment and a 0-100 score."
        ),
        "parameters": {
            "type": "object",
            "properties": {
                "smiles": {"type": "string", "description": "SMILES to assess"},
            },
            "required": ["smiles"],
        },
    },
}

_CHECK_HALLUCINATION_SCHEMA = {
    "type": "function",
    "function": {
        "name": "check_hallucination",
        "description": (
            "Check whether a proposed retrosynthesis step is a likely "
            "hallucination given the product and its reactants."
        ),
        "parameters": {
            "type": "object",
            "properties": {
                "product": {"type": "string", "description": "Product SMILES"},
                "reactants": {
                    "type": "array",
                    "items": {"type": "string"},
                    "description": "Reactant SMILES for the proposed step",
                },
            },
            "required": ["product", "reactants"],
        },
    },
}

_RUN_PYTHON_SCHEMA = {
    "type": "function",
    "function": {
        "name": "run_python",
        "description": (
            "Execute Python in a sandbox to compute chemistry (RDKit and the "
            "deepretro package are importable). Print results to stdout. No "
            "network access and a short time/memory budget."
        ),
        "parameters": {
            "type": "object",
            "properties": {
                "code": {"type": "string", "description": "Python source to run"},
            },
            "required": ["code"],
        },
    },
}


def build_tool_registry(
    hallucination_checker: HallucinationChecker | None = None,
    sandbox: Any | None = None,
    tool_backend: str = "structured",
) -> ToolRegistry:
    """Assemble the tool registry for a given backend.

    Parameters
    ----------
    hallucination_checker : callable or None, optional
        Resolved checker ``(product, pathways) -> (status, kept)``. When
        ``None`` the ``check_hallucination`` tool is omitted.
    sandbox : Sandbox or None, optional
        Sandbox for ``run_python``. A default :class:`SubprocessSandbox` is
        created when needed and none is supplied.
    tool_backend : {"structured", "sandbox"}, optional
        ``structured`` exposes the read-only check tools. ``sandbox`` adds the
        ``run_python`` code-execution tool.

    Returns
    -------
    ToolRegistry
        Schemas and executors for the selected backend.

    Examples
    --------
    >>> "run_python" in build_tool_registry(tool_backend="sandbox").names
    True
    >>> "run_python" in build_tool_registry(tool_backend="structured").names
    False
    """
    schemas: list[dict[str, Any]] = [_VALIDATE_SMILES_SCHEMA, _CHECK_STABILITY_SCHEMA]
    executors: dict[str, ToolExecutor] = {
        "validate_smiles": _validate_smiles,
        "check_stability": _check_stability,
    }

    if hallucination_checker is not None:
        schemas.append(_CHECK_HALLUCINATION_SCHEMA)
        executors["check_hallucination"] = _make_check_hallucination(
            hallucination_checker
        )

    if tool_backend == "sandbox":
        if sandbox is None:
            from deepretro.agents.sandbox import SubprocessSandbox

            sandbox = SubprocessSandbox()
        schemas.append(_RUN_PYTHON_SCHEMA)
        executors["run_python"] = _make_run_python(sandbox)

    return ToolRegistry(schemas=schemas, executors=executors)
