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
    "AutoSolver",
]


def __getattr__(name: str):
    if name == "AutoSolver":
        from .autosolve import AutoSolver

        globals()["AutoSolver"] = AutoSolver
        return AutoSolver
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
