"""Unit tests for :class:`deepretro.models.hallucination_checker.HallucinationChecker`.

This is the ML/heuristic *inference* checker in ``deepretro/models/`` -- not
to be confused with the pure-heuristic checker tested in
``test_hallucination_checker.py`` (which targets
``deepretro/algorithms/hallucination_checker.py``).

The constructor-validation and heuristic-routing tests are fast. The ML
round-trip test trains and saves a small model with the trainer, then reloads
it through the checker, so it is marked ``slow``.
"""

import csv

import pytest
from deepretro.models.hallucination_checker import HallucinationChecker

# A real, atom-balanced retrosynthetic disconnection the heuristic accepts.
VALID_PRODUCT = "Cn1nccc1[C@]1(O)CCCC[C@H]1O"
VALID_REACTANTS = "Cn1nccc1Br.O=C1CCCC[C@H]1O"
INVALID_SMILES = "not_a_smiles!!!"


# ---------------------------------------------------------------------------
# Constructor validation
# ---------------------------------------------------------------------------


def test_unknown_checker_type_raises():
    """An unsupported backend name is rejected up front."""
    with pytest.raises(ValueError, match="Unknown checker_type"):
        HallucinationChecker(checker_type="magic")


def test_ml_checker_without_model_path_raises():
    """The ML backend requires a model path."""
    with pytest.raises(ValueError, match="model_path"):
        HallucinationChecker(checker_type="ml", model_path=None)


def test_checker_type_is_lowercased():
    """The backend selector is normalised to lower case."""
    checker = HallucinationChecker(checker_type="Heuristic")
    assert checker.checker_type == "heuristic"


def test_heuristic_checker_loads_no_model():
    """The heuristic backend never touches the ML model machinery."""
    checker = HallucinationChecker(checker_type="heuristic")
    assert checker.model is None
    assert checker.featurizer is None


# ---------------------------------------------------------------------------
# Heuristic routing
# ---------------------------------------------------------------------------


@pytest.fixture
def heuristic_checker() -> HallucinationChecker:
    """A checker wired to the deterministic heuristic backend."""
    return HallucinationChecker(checker_type="heuristic")


def test_check_single_pathway_accepts_valid_step(heuristic_checker):
    """A balanced disconnection is classified as a valid step (0)."""
    assert heuristic_checker.check_single_pathway(VALID_PRODUCT, VALID_REACTANTS) == 0


def test_check_single_pathway_flags_invalid_smiles(heuristic_checker):
    """Unparseable reactants leave no valid pathway, so it is a hallucination (1)."""
    assert heuristic_checker.check_single_pathway(VALID_PRODUCT, INVALID_SMILES) == 1


def test_check_pathways_returns_valid_subset(heuristic_checker):
    """``check_pathways`` returns a 200 status and keeps the valid pathway."""
    status, valid = heuristic_checker.check_pathways(
        VALID_PRODUCT, [["Cn1nccc1Br", "O=C1CCCC[C@H]1O"]]
    )
    assert status == 200
    assert valid == [["Cn1nccc1Br", "O=C1CCCC[C@H]1O"]]


def test_check_pathways_rejects_all_invalid(heuristic_checker):
    """When nothing is valid the status is 400 and the subset is empty."""
    status, valid = heuristic_checker.check_pathways(VALID_PRODUCT, [INVALID_SMILES])
    assert status == 400
    assert valid == []


def test_call_delegates_to_check_pathways(heuristic_checker):
    """``__call__`` is a thin alias for ``check_pathways``."""
    assert heuristic_checker(VALID_PRODUCT, [INVALID_SMILES]) == (400, [])


# ---------------------------------------------------------------------------
# load_model error handling
# ---------------------------------------------------------------------------


def test_load_model_missing_config_raises(tmp_path):
    """Loading a directory without ``config.json`` fails loudly."""
    checker = HallucinationChecker(checker_type="heuristic")
    with pytest.raises(FileNotFoundError, match="config.json"):
        checker.load_model(str(tmp_path))


# ---------------------------------------------------------------------------
# ML round-trip (train -> save -> load -> predict)
# ---------------------------------------------------------------------------


def _write_reaction_csv(path, rows):
    """Write a product/reactants/label CSV the trainer can load."""
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["product", "reactants", "label"])
        writer.writerows(rows)


@pytest.fixture(scope="module")
def saved_ml_model_dir(tmp_path_factory):
    """Train and persist a tiny xgboost checker model, returning its directory."""
    from deepretro.models.hallucination_trainer import HallucinationTrainer

    rows = [
        (VALID_PRODUCT, VALID_REACTANTS, 0),
        ("CCO", "CC=O", 0),
        ("CC(=O)O", "CC(=O)Cl.O", 0),
        ("c1ccccc1O", "c1ccccc1Br.O", 0),
        ("CCN", "CC#N", 1),
        ("CCCBr", "CCCO", 1),
        ("c1ccccc1N", "c1ccccc1", 1),
        ("CC(C)O", "CC(C)Cl", 1),
    ]
    workdir = tmp_path_factory.mktemp("ml_checker")
    train_csv = workdir / "train.csv"
    test_csv = workdir / "test.csv"
    _write_reaction_csv(train_csv, rows)
    _write_reaction_csv(test_csv, rows)

    trainer = HallucinationTrainer(trainer_dir=str(workdir), model_type="xgboost")
    train_ds, test_ds = trainer.load_dataset(str(train_csv), str(test_csv))
    trainer.train_model(train_ds, test_ds)
    return trainer.model_dir


@pytest.mark.slow
def test_ml_checker_loads_saved_model(saved_ml_model_dir):
    """A saved model directory can be restored by the ML checker."""
    checker = HallucinationChecker(checker_type="ml", model_path=saved_ml_model_dir)
    assert checker.model is not None
    assert checker.featurizer is not None
    assert checker.model_type == "xgboost"


@pytest.mark.slow
def test_ml_checker_predicts_binary_label(saved_ml_model_dir):
    """The ML backend returns a binary 0/1 hallucination verdict."""
    checker = HallucinationChecker(checker_type="ml", model_path=saved_ml_model_dir)
    prediction = checker.check_single_pathway(VALID_PRODUCT, VALID_REACTANTS)
    assert prediction in (0, 1)


@pytest.mark.slow
def test_ml_check_pathways_skips_invalid_smiles(saved_ml_model_dir):
    """Invalid SMILES are dropped before ML inference and never appear as valid."""
    checker = HallucinationChecker(checker_type="ml", model_path=saved_ml_model_dir)
    status, valid = checker.check_pathways(VALID_PRODUCT, [INVALID_SMILES])
    assert status == 200
    assert valid == []
