"""Tests for the inference-time hallucination checker."""

import csv

import numpy as np
import pytest

from deepretro.models.hallucination_checker import HallucinationChecker


VALID_PRODUCT = "Cn1nccc1[C@]1(O)CCCC[C@H]1O"
VALID_REACTANTS = "Cn1nccc1Br.O=C1CCCC[C@H]1O"
INVALID_SMILES = "not_a_smiles!!!"


class _FixedFeaturizer:
    def featurize(self, datapoints):
        return np.ones((len(datapoints), 2))


class _FixedPredictionModel:
    def __init__(self, prediction):
        self.prediction = np.asarray(prediction)

    def predict(self, dataset):
        return self.prediction


def _ml_checker(prediction, threshold=0.5):
    checker = object.__new__(HallucinationChecker)
    checker.checker_type = "ml"
    checker.featurizer = _FixedFeaturizer()
    checker.model = _FixedPredictionModel(prediction)
    checker.threshold = threshold
    return checker


def _write_reaction_csv(path, rows):
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["product", "reactants", "label"])
        writer.writerows(rows)


def test_unknown_checker_type_raises():
    with pytest.raises(ValueError, match="Unknown checker_type"):
        HallucinationChecker(checker_type="magic")


def test_ml_checker_without_model_path_raises():
    with pytest.raises(ValueError, match="model_path"):
        HallucinationChecker(checker_type="ml", model_path=None)


def test_checker_type_is_lowercased():
    checker = HallucinationChecker(checker_type="Heuristic")
    assert checker.checker_type == "heuristic"


def test_heuristic_checker_loads_no_model():
    checker = HallucinationChecker(checker_type="heuristic")
    assert checker.model is None
    assert checker.featurizer is None


@pytest.fixture
def heuristic_checker():
    return HallucinationChecker(checker_type="heuristic")


def test_check_single_pathway_accepts_valid_step(heuristic_checker):
    assert heuristic_checker.check_single_pathway(VALID_PRODUCT, VALID_REACTANTS) == 0


def test_check_single_pathway_flags_invalid_smiles(heuristic_checker):
    assert heuristic_checker.check_single_pathway(VALID_PRODUCT, INVALID_SMILES) == 1


def test_check_pathways_returns_valid_subset(heuristic_checker):
    status, valid = heuristic_checker.check_pathways(
        VALID_PRODUCT, [["Cn1nccc1Br", "O=C1CCCC[C@H]1O"]]
    )
    assert status == 200
    assert valid == [["Cn1nccc1Br", "O=C1CCCC[C@H]1O"]]


def test_check_pathways_rejects_all_invalid(heuristic_checker):
    status, valid = heuristic_checker.check_pathways(VALID_PRODUCT, [INVALID_SMILES])
    assert status == 400
    assert valid == []


def test_call_delegates_to_check_pathways(heuristic_checker):
    assert heuristic_checker(VALID_PRODUCT, [INVALID_SMILES]) == (400, [])


def test_load_model_missing_config_raises(tmp_path):
    checker = HallucinationChecker(checker_type="heuristic")
    with pytest.raises(FileNotFoundError, match="config.json"):
        checker.load_model(str(tmp_path))


@pytest.mark.parametrize(
    ("prediction", "expected"),
    [
        ([[0.8, 0.2]], 0),
        ([[0.2, 0.8]], 1),
    ],
)
def test_ml_prediction_uses_positive_column_for_n_by_2(prediction, expected):
    checker = _ml_checker(prediction)
    assert checker.check_single_pathway("product", "reactants") == expected


def test_ml_prediction_uses_positive_column_for_n_by_1_by_2():
    checker = _ml_checker([[[0.1, 0.9]]])
    assert checker.check_single_pathway("product", "reactants") == 1


def test_ml_prediction_includes_threshold_boundary():
    checker = _ml_checker([[0.5, 0.5]], threshold=0.5)
    assert checker.check_single_pathway("product", "reactants") == 1


def test_ml_check_pathways_skips_invalid_smiles(monkeypatch):
    checker = _ml_checker([[0.9, 0.1]])
    calls = []
    monkeypatch.setattr(
        checker,
        "_check_pathway_ml",
        lambda target, reactants: calls.append((target, reactants)) or 0,
    )

    status, valid = checker.check_pathways(VALID_PRODUCT, [INVALID_SMILES])

    assert status == 200
    assert valid == []
    assert calls == []


@pytest.mark.slow
def test_ml_train_save_reload_predict(tmp_path):
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
    train_csv = tmp_path / "train.csv"
    test_csv = tmp_path / "test.csv"
    _write_reaction_csv(train_csv, rows)
    _write_reaction_csv(test_csv, rows)

    trainer = HallucinationTrainer(trainer_dir=str(tmp_path), model_type="xgboost")
    train_ds, test_ds = trainer.load_dataset(str(train_csv), str(test_csv))
    trainer.train_model(train_ds, test_ds)

    checker = HallucinationChecker(checker_type="ml", model_path=trainer.model_dir)
    prediction = checker.check_single_pathway(VALID_PRODUCT, VALID_REACTANTS)

    assert checker.model_type == "xgboost"
    assert checker.model is not None
    assert checker.featurizer is not None
    assert prediction in (0, 1)
