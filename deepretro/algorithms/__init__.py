from .hallucination_checker import (
    calculate_hallucination_score,
    hallucination_compare_molecules,
    interpret_score,
)

from .stability_checker import (
    check_molecule_stability,
    is_valid_smiles,
    check_carbocations,
    check_carbenes,
    check_fused_cyclopentane,
)

__all__ = [
    "calculate_hallucination_score",
    "hallucination_compare_molecules",
    "interpret_score",
    "check_molecule_stability",
    "is_valid_smiles",
    "check_carbocations",
    "check_carbenes",
    "check_fused_cyclopentane",
    "autosolve",
    "autosolve_async",
]


def __getattr__(name: str):
    if name in ("autosolve", "autosolve_async"):
        from .autosolve import autosolve, autosolve_async

        globals()["autosolve"] = autosolve
        globals()["autosolve_async"] = autosolve_async
        return globals()[name]
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
