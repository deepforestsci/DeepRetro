"""
Train XGBoost hallucination classifier using DeepChem featurization.

Pipeline: Load CSV -> Featurize SMILES -> Stratified split -> Train XGBoost -> Evaluate -> Save
"""

import sys
import os
import json
import pickle

import numpy as np
import pandas as pd
from pathlib import Path
from deepchem.feat import CircularFingerprint
from deepchem.data import NumpyDataset
from deepchem.splits import SingletaskStratifiedSplitter
from sklearn.metrics import (roc_auc_score, accuracy_score, classification_report,
                             confusion_matrix, f1_score, precision_recall_curve)
from xgboost import XGBClassifier

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from src.utils.ml_utils import extract_domain_features_batch, NUM_DOMAIN_FEATURES


def load_and_prepare_data(csv_path):
    """
    Load reaction dataset from CSV.

    Parameters
    ----------
    csv_path : str
        Path to CSV with columns: product, reactants, label.

    Returns
    -------
    products : list of str
        Product SMILES strings.
    reactants : list of str
        Reactant SMILES strings (dot-separated when multiple).
    labels : list of int
        Binary labels (1 = hallucination, 0 = real).
    """
    print("Step 1: Loading data")
    df = pd.read_csv(csv_path)
    products = df['product'].tolist()
    reactants = df['reactants'].tolist()
    labels = df['label'].tolist()
    print(f"  {len(products)} reactions | "
          f"{sum(labels)} hallucination / {len(labels)-sum(labels)} not")
    return products, reactants, labels


def featurize_smiles(products, reactants, use_domain_features=True):
    """
    Build concatenated circular fingerprints for product+reactant SMILES.

    Parameters
    ----------
    products : list of str
        Product SMILES strings.
    reactants : list of str
        Reactant SMILES strings.
    use_domain_features : bool, optional
        Append 15 hand-crafted domain features. Default True.

    Returns
    -------
    features : np.ndarray, shape (n_samples, n_features)
        Combined feature matrix.
    featurizer : CircularFingerprint
        Fitted featurizer (saved for inference).
    """
    print("\nStep 2: Featurizing SMILES...")
    featurizer = CircularFingerprint(radius=2, size=2048)
    product_features = featurizer.featurize(products)
    reactant_features = featurizer.featurize(reactants)
    combined = np.concatenate([product_features, reactant_features], axis=1)
    if use_domain_features:
        domain = extract_domain_features_batch(products, reactants)
        combined = np.concatenate([combined, domain], axis=1)
    print(f"  Product: {product_features.shape} | Reactant: {reactant_features.shape} | "
          f"Combined: {combined.shape}")
    return combined, featurizer


def create_datasets(features, labels):
    """
    Stratified split into train (70%), valid (15%), test (15%).

    Uses DeepChem's SingletaskStratifiedSplitter to preserve label balance.

    Parameters
    ----------
    features : np.ndarray, shape (n_samples, n_features)
        Feature matrix.
    labels : list of int
        Binary labels.

    Returns
    -------
    train_ds : NumpyDataset
    valid_ds : NumpyDataset
    test_ds : NumpyDataset
    """
    print("\nStep 3: Stratified splitting")
    dataset = NumpyDataset(X=features, y=np.array(labels).reshape(-1, 1))
    splitter = SingletaskStratifiedSplitter()
    train_inds, valid_inds, test_inds = splitter.split(
        dataset, frac_train=0.7, frac_valid=0.15, frac_test=0.15, seed=42
    )
    train_ds = dataset.select(train_inds.astype(int))
    valid_ds = dataset.select(valid_inds.astype(int))
    test_ds = dataset.select(test_inds.astype(int))

    for name, ds in [('train', train_ds), ('valid', valid_ds), ('test', test_ds)]:
        lbl = ds.y.flatten()
        print(f"  {name:>5}: {len(ds)} samples (0={np.sum(lbl == 0)}, 1={np.sum(lbl == 1)})")
    return train_ds, valid_ds, test_ds


def train_model(train_dataset, valid_dataset):
    """
    Train XGBoost with early stopping and class-imbalance weighting.

    Parameters
    ----------
    train_dataset : NumpyDataset
        Training data.
    valid_dataset : NumpyDataset
        Validation data for early stopping.

    Returns
    -------
    model : XGBClassifier
        Trained model.
    """
    print("\nStep 4: Training XGBoost...")
    X_train, y_train = train_dataset.X, train_dataset.y.flatten()
    X_valid, y_valid = valid_dataset.X, valid_dataset.y.flatten()

    neg, pos = np.sum(y_train == 0), np.sum(y_train == 1)
    scale_pos_weight = neg / pos if pos > 0 else 1.0
    print(f"  {len(X_train)} train / {len(X_valid)} valid | "
          f"{X_train.shape[1]} features | scale_pos_weight={scale_pos_weight:.2f}")

    model = XGBClassifier(
        max_depth=6, learning_rate=0.05, n_estimators=300,
        subsample=0.8, colsample_bytree=0.8, min_child_weight=3,
        gamma=0.1, scale_pos_weight=scale_pos_weight, random_state=42,
        eval_metric='logloss', early_stopping_rounds=20,
    )
    model.fit(X_train, y_train, eval_set=[(X_valid, y_valid)], verbose=True)
    print(f"\n  Best iteration: {model.best_iteration}/{model.n_estimators} | "
          f"Best score: {model.best_score:.6f}")
    return model


def evaluate_model(model, test_dataset):
    """
    Evaluate on held-out test set and find optimal classification threshold.

    Threshold is chosen to maximise F1 on the precision-recall curve.

    Parameters
    ----------
    model : XGBClassifier
        Trained model.
    test_dataset : NumpyDataset
        Held-out test data.

    Returns
    -------
    scores : dict
        Keys: test_auc, test_acc, test_f1, optimal_threshold,
        optimal_f1, optimal_acc.
    """
    print("\nStep 5: Evaluating on test set...")
    X, y = test_dataset.X, test_dataset.y.flatten()
    pred = model.predict(X)
    proba = model.predict_proba(X)

    auc = roc_auc_score(y, proba[:, 1])
    acc = accuracy_score(y, pred)
    f1 = f1_score(y, pred)
    print(f"  AUC={auc:.4f} | Acc={acc:.4f} | F1={f1:.4f}")

    print(f"\n  Classification Report:")
    print(classification_report(y, pred, target_names=['Not Hallucination', 'Hallucination']))
    cm = confusion_matrix(y, pred)
    print(f"  Confusion: TN={cm[0,0]} FP={cm[0,1]} FN={cm[1,0]} TP={cm[1,1]}")

    # Optimal threshold via precision-recall curve
    precision, recall, thresholds = precision_recall_curve(y, proba[:, 1])
    f1_scores = 2 * (precision * recall) / (precision + recall + 1e-10)
    opt_idx = np.argmax(f1_scores)
    opt_threshold = thresholds[opt_idx]
    opt_pred = (proba[:, 1] >= opt_threshold).astype(int)
    opt_f1 = f1_score(y, opt_pred)
    opt_acc = accuracy_score(y, opt_pred)
    print(f"\n  Optimal threshold: {opt_threshold:.3f} -> F1={opt_f1:.4f} Acc={opt_acc:.4f}")

    return {'test_auc': auc, 'test_acc': acc, 'test_f1': f1,
            'optimal_threshold': opt_threshold, 'optimal_f1': opt_f1, 'optimal_acc': opt_acc}


def save_model_artifacts(model, featurizer, optimal_threshold, use_domain_features=True):
    """
    Persist model, featurizer, and metadata to ``models/`` directory.

    Parameters
    ----------
    model : XGBClassifier
        Trained model.
    featurizer : CircularFingerprint
        Featurizer instance.
    optimal_threshold : float
        Classification threshold from evaluation.
    use_domain_features : bool, optional
        Whether domain features were used during training. Default True.
    """
    print("\nStep 6: Saving artifacts...")
    models_dir = Path("models")
    models_dir.mkdir(exist_ok=True)

    for name, obj in [("xgboost_hallucination_model.pkl", model),
                      ("xgboost_featurizer.pkl", featurizer)]:
        with open(models_dir / name, 'wb') as f:
            pickle.dump(obj, f)
        print(f"  ✓ {name}")

    feat_dim = 4096 + (NUM_DOMAIN_FEATURES if use_domain_features else 0)
    metadata = {
        'optimal_threshold': float(optimal_threshold),
        'model_type': 'XGBoost', 'featurizer_type': 'CircularFingerprint',
        'featurizer_radius': 2, 'featurizer_size': 2048,
        'feature_dimension': feat_dim, 'domain_features': use_domain_features,
    }
    with open(models_dir / "xgboost_metadata.json", 'w') as f:
        json.dump(metadata, f, indent=2)
    print(f"  ✓ xgboost_metadata.json (threshold={optimal_threshold:.4f})")


def main():
    """Orchestrate the full training pipeline."""
    csv_path = "data/xgboost_dataset_merged.csv"
    use_domain_features = True

    print("=" * 50)
    print("XGBoost Hallucination Classifier Training")
    print("=" * 50)

    products, reactants, labels = load_and_prepare_data(csv_path)
    features, featurizer = featurize_smiles(products, reactants, use_domain_features)
    train_ds, valid_ds, test_ds = create_datasets(features, labels)
    model = train_model(train_ds, valid_ds)
    scores = evaluate_model(model, test_ds)
    save_model_artifacts(model, featurizer, scores['optimal_threshold'], use_domain_features)

    print("\n" + "=" * 50)
    print("Training Complete!")
    print("=" * 50)
    return model


if __name__ == "__main__":
    main()
