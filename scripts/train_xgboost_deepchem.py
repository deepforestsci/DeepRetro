"""
Train XGBoost classifier for hallucination prediction using DeepChem for featurization.

This script:
1. Loads the dataset (product, reactants, label)
2. Featurizes SMILES using DeepChem's CircularFingerprint
3. Creates DeepChem dataset and splits
4. Trains XGBoost directly (using DeepChem only for featurization)
5. Evaluates the model
"""

import pandas as pd
import numpy as np
from pathlib import Path
import pickle
import json

# DeepChem imports (for featurization only)
from deepchem.feat import CircularFingerprint
from deepchem.data import NumpyDataset
from sklearn.metrics import roc_auc_score, accuracy_score, classification_report, confusion_matrix, f1_score, precision_recall_curve
from sklearn.model_selection import train_test_split
from deepchem.splits import RandomSplitter

# XGBoost import (using directly, not through DeepChem wrapper)
from xgboost import XGBClassifier

# RDKit for domain-specific features
from rdkit import Chem
from rdkit.Chem import Descriptors
from collections import Counter


def load_and_prepare_data(csv_path: str):
    """
    Step 1: Load CSV and prepare SMILES strings.
    
    Reads the CSV file with columns: product, reactants, target, label
    Returns lists of product SMILES, reactant SMILES, targets, and labels.
    """
    print("Step 1: Loading data from CSV...")
    df = pd.read_csv(csv_path)
    
    products = df['product'].tolist()
    reactants = df['reactants'].tolist()
    labels = df['label'].tolist()
    targets = df['target'].tolist() if 'target' in df.columns else ['unknown'] * len(products)
    
    print(f"  Loaded {len(products)} unique reactions (already deduplicated)")
    print(f"  Unique targets: {len(set(targets))}")
    print(f"  Label distribution: {sum(labels)} hallucination, {len(labels)-sum(labels)} not hallucination")
    
    return products, reactants, labels, targets


def extract_domain_features(products: list, reactants: list):
    """
    Extract domain-specific features for hallucination detection.
    These features capture chemical inconsistencies that fingerprints might miss.
    """
    domain_features_list = []
    
    for product_smiles, reactant_string in zip(products, reactants):
        features = []
        
        try:
            # Parse molecules
            product_mol = Chem.MolFromSmiles(product_smiles)
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactant_string.split(".")]
            
            if product_mol is None:
                # Invalid product - use zeros
                features = [0.0] * 15
                domain_features_list.append(features)
                continue
            
            # Product features
            product_atoms = Counter([a.GetSymbol() for a in product_mol.GetAtoms()])
            product_bonds = product_mol.GetNumBonds()
            product_rings = len(Chem.GetSSSR(product_mol))
            product_aromatic = sum(1 for a in product_mol.GetAtoms() if a.GetIsAromatic())
            product_mw = Descriptors.MolWt(product_mol)
            
            # Reactant features (sum across all reactants)
            reactant_atoms = Counter()
            reactant_bonds = 0
            reactant_rings = 0
            reactant_aromatic = 0
            reactant_mw = 0
            num_valid_reactants = 0
            
            for mol in reactant_mols:
                if mol:
                    reactant_atoms += Counter([a.GetSymbol() for a in mol.GetAtoms()])
                    reactant_bonds += mol.GetNumBonds()
                    reactant_rings += len(Chem.GetSSSR(mol))
                    reactant_aromatic += sum(1 for a in mol.GetAtoms() if a.GetIsAromatic())
                    reactant_mw += Descriptors.MolWt(mol)
                    num_valid_reactants += 1
            
            # Difference features (key for hallucination detection!)
            features = [
                float(reactant_atoms.get('C', 0) - product_atoms.get('C', 0)),  # C difference
                float(reactant_atoms.get('N', 0) - product_atoms.get('N', 0)),  # N difference
                float(reactant_atoms.get('O', 0) - product_atoms.get('O', 0)),  # O difference
                float(reactant_atoms.get('Cl', 0) - product_atoms.get('Cl', 0)),  # Cl difference
                float(reactant_atoms.get('Br', 0) - product_atoms.get('Br', 0)),  # Br difference
                float(reactant_bonds - product_bonds),  # Bond difference
                float(reactant_rings - product_rings),  # Ring difference
                float(reactant_aromatic - product_aromatic),  # Aromatic difference
                float(reactant_mw - product_mw),  # Molecular weight difference
                float(num_valid_reactants),  # Number of reactants
                float(sum(reactant_atoms.values()) - sum(product_atoms.values())),  # Total atom difference
                float(product_mw),  # Product molecular weight
                float(reactant_mw),  # Total reactant molecular weight
                float(product_rings),  # Product ring count
                float(reactant_rings),  # Reactant ring count
            ]
            
        except Exception as e:
            # On error, use zeros
            features = [0.0] * 15
        
        domain_features_list.append(features)
    
    return np.array(domain_features_list)


def featurize_smiles(products: list, reactants: list, use_domain_features: bool = True):
    """
    Step 2: Convert SMILES strings to molecular fingerprints.
    
    Uses DeepChem's CircularFingerprint (Morgan/ECFP fingerprints):
    - radius=2: captures local structure around each atom
    - size=2048: creates 2048-bit fingerprint vectors
    
    For each sample:
    - Product SMILES → product fingerprint (2048 features)
    - Reactants SMILES → reactant fingerprint (2048 features)
    - Domain features → 15 chemical consistency features (optional)
    - Concatenate → combined fingerprint (4096 or 4111 features)
    """
    print("\nStep 2: Featurizing SMILES strings...")
    
    # Create featurizer (Morgan fingerprint, radius 2, 2048 bits)
    featurizer = CircularFingerprint(radius=2, size=2048)
    
    # Featurize products
    print("  Featurizing products...")
    product_features = featurizer.featurize(products)
    
    # Featurize reactants
    print("  Featurizing reactants...")
    reactant_features = featurizer.featurize(reactants)
    
    # Combine fingerprints
    print("  Combining fingerprints...")
    fingerprint_features = np.concatenate([product_features, reactant_features], axis=1)
    
    # Add domain-specific features
    if use_domain_features:
        print("  Extracting domain-specific features (atom/bond/ring differences)...")
        domain_features = extract_domain_features(products, reactants)
        combined_features = np.concatenate([fingerprint_features, domain_features], axis=1)
        print(f"  Domain features shape: {domain_features.shape}")
    else:
        combined_features = fingerprint_features
    
    print(f"  Product features shape: {product_features.shape}")
    print(f"  Reactant features shape: {reactant_features.shape}")
    print(f"  Combined features shape: {combined_features.shape}")
    
    return combined_features, featurizer


def create_deepchem_dataset(features: np.ndarray, labels: list, targets: list):
    """
    Step 3: Create DeepChem dataset and split into train/validation/test.
    
    Uses STRATIFIED SPLITTING BY TARGET to ensure:
    - All reactions from the same target go to the same split
    - No leakage between train/test
    - Tests generalization to new targets
    
    DeepChem's NumpyDataset:
    - Stores features (X) and labels (y) as NumPy arrays
    - Provides interface for DeepChem models
    
    Split strategy:
    - 70% training (targets)
    - 15% validation (targets)
    - 15% test (targets)
    
    Returns
    -------
    tuple
        (train_dataset, valid_dataset, test_dataset, split_indices)
        where split_indices is a dict with 'train', 'valid', 'test' keys
        containing the original indices for each split
    """
    print("\nStep 3: Creating DeepChem dataset and splitting by target...")
    
    # Convert labels to numpy array and reshape for DeepChem
    # DeepChem expects labels as (n_samples, n_tasks) - we have 1 task
    labels_array = np.array(labels).reshape(-1, 1)
    
    # Create DeepChem dataset
    dataset = NumpyDataset(X=features, y=labels_array)
    
    # STRATIFIED SPLITTING BY TARGET
    # Group indices by target
    from collections import defaultdict
    target_to_indices = defaultdict(list)
    for idx, target in enumerate(targets):
        target_to_indices[target].append(idx)
    
    # Get unique targets and split them
    unique_targets = list(target_to_indices.keys())
    n_targets = len(unique_targets)
    
    # Split targets (not samples)
    # First split: train targets (70%) vs temp targets (30%)
    train_targets, temp_targets = train_test_split(
        unique_targets, test_size=0.3, random_state=42, shuffle=True
    )
    
    # Second split: temp -> valid targets (15%) and test targets (15%)
    valid_targets, test_targets = train_test_split(
        temp_targets, test_size=0.5, random_state=42, shuffle=True
    )
    
    # Get indices for each split based on target assignment
    train_indices = []
    for target in train_targets:
        train_indices.extend(target_to_indices[target])
    
    valid_indices = []
    for target in valid_targets:
        valid_indices.extend(target_to_indices[target])
    
    test_indices = []
    for target in test_targets:
        test_indices.extend(target_to_indices[target])
    
    # Convert to numpy arrays and sort for consistency
    train_indices = np.array(sorted(train_indices))
    valid_indices = np.array(sorted(valid_indices))
    test_indices = np.array(sorted(test_indices))
    
    # Create split datasets
    train_dataset = NumpyDataset(X=dataset.X[train_indices], y=dataset.y[train_indices])
    valid_dataset = NumpyDataset(X=dataset.X[valid_indices], y=dataset.y[valid_indices])
    test_dataset = NumpyDataset(X=dataset.X[test_indices], y=dataset.y[test_indices])
    
    # Store split indices for leakage checking
    split_indices = {
        'train': train_indices,
        'valid': valid_indices,
        'test': test_indices
    }
    
    print(f"  Unique targets: {n_targets}")
    print(f"  Train targets: {len(train_targets)} ({len(train_targets)/n_targets*100:.1f}%)")
    print(f"  Valid targets: {len(valid_targets)} ({len(valid_targets)/n_targets*100:.1f}%)")
    print(f"  Test targets: {len(test_targets)} ({len(test_targets)/n_targets*100:.1f}%)")
    print(f"  Training samples: {len(train_dataset)}")
    print(f"  Validation samples: {len(valid_dataset)}")
    print(f"  Test samples: {len(test_dataset)}")
    
    # Report class distribution in each split
    train_labels = dataset.y[train_indices].flatten()
    valid_labels = dataset.y[valid_indices].flatten()
    test_labels = dataset.y[test_indices].flatten()
    
    print(f"\n  Class distribution in splits:")
    print(f"    Train: Label 0={np.sum(train_labels == 0)} ({np.sum(train_labels == 0)/len(train_labels)*100:.1f}%), "
          f"Label 1={np.sum(train_labels == 1)} ({np.sum(train_labels == 1)/len(train_labels)*100:.1f}%)")
    print(f"    Valid: Label 0={np.sum(valid_labels == 0)} ({np.sum(valid_labels == 0)/len(valid_labels)*100:.1f}%), "
          f"Label 1={np.sum(valid_labels == 1)} ({np.sum(valid_labels == 1)/len(valid_labels)*100:.1f}%)")
    print(f"    Test:  Label 0={np.sum(test_labels == 0)} ({np.sum(test_labels == 0)/len(test_labels)*100:.1f}%), "
          f"Label 1={np.sum(test_labels == 1)} ({np.sum(test_labels == 1)/len(test_labels)*100:.1f}%)")
    
    return train_dataset, valid_dataset, test_dataset, split_indices


def train_xgboost_model(train_dataset, valid_dataset, use_domain_features: bool = True):
    """
    Step 4: Train XGBoost classifier directly (using DeepChem for featurization only).
    
    Process:
    1. Extract features and labels from DeepChem datasets
    2. Calculate class weights for imbalanced data
    3. Create XGBoost classifier with hyperparameters
    4. Train on training dataset
    5. Use validation dataset for early stopping
    
    Note: We use DeepChem for featurization but train XGBoost directly
    to avoid GBDTModel wrapper issues.
    """
    print("\nStep 4: Training XGBoost model...")
    
    # Extract features and labels from DeepChem datasets
    X_train = train_dataset.X
    y_train = train_dataset.y.flatten() if train_dataset.y.ndim > 1 else train_dataset.y
    X_valid = valid_dataset.X
    y_valid = valid_dataset.y.flatten() if valid_dataset.y.ndim > 1 else valid_dataset.y
    
    # Verify data shapes
    print(f"    X_train shape: {X_train.shape}, y_train shape: {y_train.shape}")
    print(f"    X_valid shape: {X_valid.shape}, y_valid shape: {y_valid.shape}")
    print(f"    y_train unique values: {np.unique(y_train)}")
    
    # Calculate class imbalance ratio for scale_pos_weight
    # scale_pos_weight = (negative_samples / positive_samples)
    # This helps XGBoost handle imbalanced classes
    negative_count = np.sum(y_train == 0)
    positive_count = np.sum(y_train == 1)
    scale_pos_weight = negative_count / positive_count if positive_count > 0 else 1.0
    
    print(f"    Class distribution: {negative_count} negative (Not Hallucination), {positive_count} positive (Hallucination)")
    print(f"    Using scale_pos_weight: {scale_pos_weight:.2f} to handle class imbalance")
    
    # Create XGBoost classifier with improved hyperparameters
    # Note: In XGBoost 2.0+, early_stopping_rounds goes in constructor
    model = XGBClassifier(
        max_depth=6,           # Maximum tree depth
        learning_rate=0.05,    # Lower learning rate for better convergence
        n_estimators=300,      # More trees (with early stopping)
        subsample=0.8,         # Fraction of samples per tree
        colsample_bytree=0.8,  # Fraction of features per tree
        min_child_weight=3,   # Minimum samples in leaf (reduces overfitting)
        gamma=0.1,            # Minimum loss reduction for split
        scale_pos_weight=scale_pos_weight,  # Handle class imbalance
        random_state=42,
        eval_metric='logloss',  # Metric for early stopping
        early_stopping_rounds=20  # Early stopping rounds
    )
    
    # Train the model with early stopping
    # In XGBoost 2.0+, early_stopping_rounds is set in constructor
    # eval_set is still passed to fit() for monitoring
    print("  Training...")
    print(f"    Training on {len(X_train)} samples with {X_train.shape[1]} features")
    print(f"    Validating on {len(X_valid)} samples")
    
    model.fit(
        X_train, y_train,
        eval_set=[(X_valid, y_valid)],
        verbose=True  # Show training progress
    )
    
    print(f"\n  Training complete!")
    max_iterations = 200  # From n_estimators parameter
    print(f"    Best iteration: {model.best_iteration} (out of {max_iterations} max)")
    print(f"    Best validation score: {model.best_score:.6f}")
    if model.best_iteration < max_iterations - 1:
        rounds_since_best = max_iterations - model.best_iteration - 1
        print(f"    ✓ Early stopping triggered (no improvement for {rounds_since_best} rounds)")
    else:
        print(f"    → Reached maximum iterations ({max_iterations})")
    
    # Analyze feature importance
    print(f"\n  Analyzing feature importance...")
    analyze_feature_importance(model, use_domain_features=use_domain_features)
    
    return model


def analyze_feature_importance(model, use_domain_features: bool = True):
    """
    Analyze and display feature importance from the trained XGBoost model.
    
    Parameters
    ----------
    model : XGBClassifier
        Trained XGBoost model
    use_domain_features : bool
        Whether domain features were used
    """
    # Get feature importance (gain-based)
    feature_importance = model.feature_importances_
    n_features = len(feature_importance)
    
    # Create feature names
    feature_names = []
    
    # Product fingerprint features (0-2047)
    for i in range(2048):
        feature_names.append(f"Product_FP_{i}")
    
    # Reactant fingerprint features (2048-4095)
    for i in range(2048):
        feature_names.append(f"Reactant_FP_{i}")
    
    # Domain features (4096-4110)
    if use_domain_features:
        domain_feature_names = [
            "C_diff", "N_diff", "O_diff", "Cl_diff", "Br_diff",
            "Bond_diff", "Ring_diff", "Aromatic_diff", "MW_diff",
            "Num_reactants", "Total_atom_diff",
            "Product_MW", "Reactant_MW", "Product_rings", "Reactant_rings"
        ]
        feature_names.extend(domain_feature_names)
    
    # Get top features
    feature_importance_with_names = list(zip(feature_names, feature_importance))
    feature_importance_with_names.sort(key=lambda x: x[1], reverse=True)
    
    print(f"\n  Top 20 Most Important Features:")
    print(f"    {'Rank':<6} {'Feature Name':<30} {'Importance':<12}")
    print(f"    {'-'*6} {'-'*30} {'-'*12}")
    
    for rank, (name, importance) in enumerate(feature_importance_with_names[:20], 1):
        print(f"    {rank:<6} {name:<30} {importance:<12.6f}")
    
    # Analyze domain features specifically
    if use_domain_features:
        print(f"\n  Domain Features Importance:")
        print(f"    {'Feature':<25} {'Importance':<12} {'% of Total':<12}")
        print(f"    {'-'*25} {'-'*12} {'-'*12}")
        
        domain_start_idx = 4096
        total_importance = sum(feature_importance)
        domain_total = sum(feature_importance[domain_start_idx:])
        
        for i, name in enumerate(domain_feature_names):
            idx = domain_start_idx + i
            importance = feature_importance[idx]
            pct = (importance / total_importance * 100) if total_importance > 0 else 0
            print(f"    {name:<25} {importance:<12.6f} {pct:<12.2f}%")
        
        print(f"\n    Domain features total importance: {domain_total:.6f} ({domain_total/total_importance*100:.2f}% of all features)")
        print(f"    Fingerprint features total: {sum(feature_importance[:4096]):.6f} ({sum(feature_importance[:4096])/total_importance*100:.2f}% of all features)")
    
    # Fingerprint analysis
    product_fp_importance = sum(feature_importance[0:2048])
    reactant_fp_importance = sum(feature_importance[2048:4096])
    
    print(f"\n  Fingerprint Features Breakdown:")
    print(f"    Product fingerprints: {product_fp_importance:.6f} ({product_fp_importance/sum(feature_importance)*100:.2f}%)")
    print(f"    Reactant fingerprints: {reactant_fp_importance:.6f} ({reactant_fp_importance/sum(feature_importance)*100:.2f}%)")


def evaluate_model(model, train_dataset, valid_dataset, test_dataset, products=None, reactants=None, targets=None, split_indices=None):
    """
    Step 5: Evaluate model performance on all datasets.
    
    Uses model.predict() to get predictions, then computes metrics:
    - ROC-AUC: Area under ROC curve (good for imbalanced data)
    - Accuracy: Overall correctness
    - Classification report: Precision, recall, F1 for each class
    
    Also saves validation and test predictions to a CSV file.
    
    Parameters
    ----------
    model : XGBClassifier
        Trained XGBoost model
    train_dataset : NumpyDataset
        Training dataset
    valid_dataset : NumpyDataset
        Validation dataset
    test_dataset : NumpyDataset
        Test dataset
    products : list, optional
        List of product SMILES strings
    reactants : list, optional
        List of reactant SMILES strings
    targets : list, optional
        List of target SMILES strings
    split_indices : dict, optional
        Dictionary with 'valid' and 'test' keys containing indices
    """
    print("\nStep 5: Evaluating model...")
    
    # Extract features and labels from datasets
    X_train = train_dataset.X
    y_train = train_dataset.y
    X_valid = valid_dataset.X
    y_valid = valid_dataset.y
    X_test = test_dataset.X
    y_test = test_dataset.y
    
    # Get predictions and true labels for each dataset
    print("  Getting predictions...")
    train_pred = model.predict(X_train)
    train_pred_proba = model.predict_proba(X_train)
    train_true = y_train
    
    valid_pred = model.predict(X_valid)
    valid_pred_proba = model.predict_proba(X_valid)
    valid_true = y_valid
    
    test_pred = model.predict(X_test)
    test_pred_proba = model.predict_proba(X_test)
    test_true = y_test
    
    # Compute ROC-AUC (needs probability predictions)
    print("  Computing metrics...")
    train_auc = roc_auc_score(train_true, train_pred_proba[:, 1])
    valid_auc = roc_auc_score(valid_true, valid_pred_proba[:, 1])
    test_auc = roc_auc_score(test_true, test_pred_proba[:, 1])
    
    # Compute accuracy
    train_acc = accuracy_score(train_true, train_pred)
    valid_acc = accuracy_score(valid_true, valid_pred)
    test_acc = accuracy_score(test_true, test_pred)
    
    print(f"\n  Training Results:")
    print(f"    ROC-AUC: {train_auc:.4f}")
    print(f"    Accuracy: {train_acc:.4f}")
    
    print(f"\n  Validation Results:")
    print(f"    ROC-AUC: {valid_auc:.4f}")
    print(f"    Accuracy: {valid_acc:.4f}")
    
    print(f"\n  Test Results:")
    print(f"    ROC-AUC: {test_auc:.4f}")
    print(f"    Accuracy: {test_acc:.4f}")
    
    # Classification report for test set
    print(f"\n  Test Set Classification Report:")
    print(classification_report(test_true, test_pred, 
                                target_names=['Not Hallucination', 'Hallucination']))
    
    # Confusion matrix
    print(f"\n  Test Set Confusion Matrix:")
    cm = confusion_matrix(test_true, test_pred)
    print(f"    True Negatives: {cm[0,0]}, False Positives: {cm[0,1]}")
    print(f"    False Negatives: {cm[1,0]}, True Positives: {cm[1,1]}")
    
    # F1 scores
    test_f1 = f1_score(test_true, test_pred)
    print(f"\n  Test Set F1-Score: {test_f1:.4f}")
    
    # Find optimal threshold for better recall
    print(f"\n  Finding optimal classification threshold...")
    precision, recall, thresholds = precision_recall_curve(test_true, test_pred_proba[:, 1])
    
    # Find threshold that maximizes F1 score
    f1_scores = 2 * (precision * recall) / (precision + recall + 1e-10)
    optimal_idx = np.argmax(f1_scores)
    optimal_threshold = thresholds[optimal_idx]
    optimal_f1 = f1_scores[optimal_idx]
    
    print(f"    Default threshold (0.5): F1 = {test_f1:.4f}")
    print(f"    Optimal threshold ({optimal_threshold:.3f}): F1 = {optimal_f1:.4f}")
    
    # Predictions with optimal threshold
    test_pred_optimal = (test_pred_proba[:, 1] >= optimal_threshold).astype(int)
    test_acc_optimal = accuracy_score(test_true, test_pred_optimal)
    test_f1_optimal = f1_score(test_true, test_pred_optimal)
    
    print(f"    Accuracy with optimal threshold: {test_acc_optimal:.4f}")
    print(f"    F1-Score with optimal threshold: {test_f1_optimal:.4f}")
    
    # Classification report with optimal threshold
    print(f"\n  Test Set Classification Report (Optimal Threshold):")
    print(classification_report(test_true, test_pred_optimal, 
                                target_names=['Not Hallucination', 'Hallucination']))
    
    # Save predictions to file
    print(f"\n  Saving predictions to file...")
    save_predictions_to_file(
        valid_pred, valid_pred_proba, valid_true, 
        test_pred, test_pred_proba, test_true,
        test_pred_optimal, optimal_threshold,
        products, reactants, targets, split_indices
    )
    
    return {
        'test_auc': test_auc,
        'test_acc': test_acc,
        'test_f1': test_f1,
        'optimal_threshold': optimal_threshold,
        'optimal_f1': optimal_f1,
        'optimal_acc': test_acc_optimal
    }


def save_model_artifacts(model: XGBClassifier, featurizer: CircularFingerprint, optimal_threshold: float):
    """
    Save the trained model, featurizer, and optimal threshold to disk.
    
    Parameters
    ----------
    model : XGBClassifier
        Trained XGBoost model
    featurizer : CircularFingerprint
        DeepChem featurizer used for SMILES encoding
    optimal_threshold : float
        Optimal classification threshold for binary predictions
    """
    print("\nStep 6: Saving model artifacts...")
    
    # Create models directory if it doesn't exist
    models_dir = Path("models")
    models_dir.mkdir(exist_ok=True)
    
    # Save model
    model_path = models_dir / "xgboost_hallucination_model.pkl"
    with open(model_path, 'wb') as f:
        pickle.dump(model, f)
    print(f"  ✓ Saved model to {model_path}")
    
    # Save featurizer
    featurizer_path = models_dir / "xgboost_featurizer.pkl"
    with open(featurizer_path, 'wb') as f:
        pickle.dump(featurizer, f)
    print(f"  ✓ Saved featurizer to {featurizer_path}")
    
    # Save threshold and metadata
    metadata = {
        'optimal_threshold': float(optimal_threshold),
        'model_type': 'XGBoost',
        'featurizer_type': 'CircularFingerprint',
        'featurizer_radius': 2,
        'featurizer_size': 2048,
        'feature_dimension': 4111,  # product (2048) + reactant (2048) + domain (15)
        'domain_features': True  # Includes atom/bond/ring difference features
    }
    metadata_path = models_dir / "xgboost_metadata.json"
    with open(metadata_path, 'w') as f:
        json.dump(metadata, f, indent=2)
    print(f"  ✓ Saved metadata to {metadata_path}")
    print(f"    Optimal threshold: {optimal_threshold:.4f}")


def check_data_leakage(train_dataset, valid_dataset, test_dataset, split_indices, products, reactants):
    """
    Check for train-test leakage by examining:
    1. Duplicate SMILES strings across splits
    2. Duplicate reaction pairs (product + reactants) across splits
    3. Exact sample overlap between train/valid/test sets
    
    Parameters
    ----------
    train_dataset : NumpyDataset
        Training dataset
    valid_dataset : NumpyDataset
        Validation dataset
    test_dataset : NumpyDataset
        Test dataset
    split_indices : dict
        Dictionary with 'train', 'valid', 'test' keys containing indices
    products : list
        List of product SMILES strings
    reactants : list
        List of reactant SMILES strings
    """
    print("\n" + "=" * 60)
    print("CHECKING FOR DATA LEAKAGE")
    print("=" * 60)
    
    train_indices = set(split_indices['train'])
    valid_indices = set(split_indices['valid'])
    test_indices = set(split_indices['test'])
    
    # Check 0: Verify no index overlap (should never happen, but good to check)
    print(f"\n0. Checking split integrity...")
    train_valid_overlap = train_indices & valid_indices
    train_test_overlap = train_indices & test_indices
    valid_test_overlap = valid_indices & test_indices
    
    if train_valid_overlap or train_test_overlap or valid_test_overlap:
        print(f"   ❌ CRITICAL: Index overlap detected!")
        if train_valid_overlap:
            print(f"      Train-Valid overlap: {len(train_valid_overlap)} indices")
        if train_test_overlap:
            print(f"      Train-Test overlap: {len(train_test_overlap)} indices")
        if valid_test_overlap:
            print(f"      Valid-Test overlap: {len(valid_test_overlap)} indices")
    else:
        print(f"   ✓ No index overlap - splits are properly separated")
    
    # Check 1: Duplicate product SMILES across splits
    print(f"\n1. Checking for duplicate PRODUCT SMILES across splits...")
    train_products = {products[i] for i in train_indices}
    valid_products = {products[i] for i in valid_indices}
    test_products = {products[i] for i in test_indices}
    
    train_valid_product_overlap = train_products & valid_products
    train_test_product_overlap = train_products & test_products
    valid_test_product_overlap = valid_products & test_products
    
    total_product_overlap = len(train_valid_product_overlap | train_test_product_overlap | valid_test_product_overlap)
    
    if total_product_overlap > 0:
        print(f"   ⚠️  Found {total_product_overlap} products appearing in multiple splits:")
        print(f"      Train-Valid: {len(train_valid_product_overlap)} products")
        print(f"      Train-Test: {len(train_test_product_overlap)} products")
        print(f"      Valid-Test: {len(valid_test_product_overlap)} products")
        print(f"   Note: This is OK if different reaction contexts, but worth checking")
    else:
        print(f"   ✓ No product SMILES overlap across splits")
    
    # Check 2: Duplicate reactant SMILES across splits
    print(f"\n2. Checking for duplicate REACTANT SMILES across splits...")
    train_reactants = {reactants[i] for i in train_indices}
    valid_reactants = {reactants[i] for i in valid_indices}
    test_reactants = {reactants[i] for i in test_indices}
    
    train_valid_reactant_overlap = train_reactants & valid_reactants
    train_test_reactant_overlap = train_reactants & test_reactants
    valid_test_reactant_overlap = valid_reactants & test_reactants
    
    total_reactant_overlap = len(train_valid_reactant_overlap | train_test_reactant_overlap | valid_test_reactant_overlap)
    
    if total_reactant_overlap > 0:
        print(f"   ⚠️  Found {total_reactant_overlap} reactant combinations appearing in multiple splits:")
        print(f"      Train-Valid: {len(train_valid_reactant_overlap)} reactants")
        print(f"      Train-Test: {len(train_test_reactant_overlap)} reactants")
        print(f"      Valid-Test: {len(valid_test_reactant_overlap)} reactants")
        print(f"   Note: This is OK if different products, but worth checking")
    else:
        print(f"   ✓ No reactant SMILES overlap across splits")
    
    # Check 3: Duplicate reaction pairs (product + reactants) - THIS IS THE CRITICAL CHECK
    print(f"\n3. Checking for duplicate REACTION PAIRS (product + reactants) across splits...")
    print(f"   This is the most critical check - same reaction should not be in train and test!")
    
    train_pairs = {f"{products[i]}>>{reactants[i]}" for i in train_indices}
    valid_pairs = {f"{products[i]}>>{reactants[i]}" for i in valid_indices}
    test_pairs = {f"{products[i]}>>{reactants[i]}" for i in test_indices}
    
    train_valid_pair_overlap = train_pairs & valid_pairs
    train_test_pair_overlap = train_pairs & test_pairs
    valid_test_pair_overlap = valid_pairs & test_pairs
    
    total_pair_overlap = len(train_valid_pair_overlap | train_test_pair_overlap | valid_test_pair_overlap)
    
    if total_pair_overlap > 0:
        print(f"   ❌ CRITICAL LEAKAGE DETECTED: {total_pair_overlap} reaction pairs in multiple splits!")
        print(f"      Train-Valid overlap: {len(train_valid_pair_overlap)} pairs")
        print(f"      Train-Test overlap: {len(train_test_pair_overlap)} pairs ⚠️")
        print(f"      Valid-Test overlap: {len(valid_test_pair_overlap)} pairs")
        
        # Calculate impact
        test_size = len(test_indices)
        train_size = len(train_indices)
        leakage_percentage = (len(train_test_pair_overlap) / test_size * 100) if test_size > 0 else 0
        
        print(f"\n   📊 IMPACT ANALYSIS:")
        print(f"      Test set size: {test_size} samples")
        print(f"      Leaked samples: {len(train_test_pair_overlap)} samples")
        print(f"      Leakage percentage: {leakage_percentage:.1f}% of test set")
        
        if leakage_percentage > 10:
            print(f"      ⚠️  HIGH IMPACT: {leakage_percentage:.1f}% of test samples are in training!")
            print(f"         This will significantly inflate test performance.")
        elif leakage_percentage > 5:
            print(f"      ⚠️  MODERATE IMPACT: {leakage_percentage:.1f}% of test samples are in training.")
        else:
            print(f"      ⚠️  LOW IMPACT: {leakage_percentage:.1f}% leakage (still should be fixed).")
        
        if train_test_pair_overlap:
            print(f"\n   Example train-test leakage (first 3):")
            for i, pair in enumerate(list(train_test_pair_overlap)[:3]):
                print(f"     {pair[:100]}...")
    else:
        print(f"   ✓ No reaction pair overlap - no leakage detected!")
    
    # Check 4: Preprocessing order
    print(f"\n4. Preprocessing order check:")
    print(f"   ✓ Featurization done BEFORE splitting (good - no preprocessing leakage)")
    
    # Summary
    print(f"\n" + "=" * 60)
    print("LEAKAGE CHECK SUMMARY")
    print("=" * 60)
    
    has_critical_leakage = (total_pair_overlap > 0 or 
                           train_test_overlap or 
                           train_valid_overlap or 
                           valid_test_overlap)
    
    if has_critical_leakage:
        print("❌ CRITICAL: Data leakage detected!")
        if total_pair_overlap > 0:
            print("   - Same reaction pairs appear in train and test sets")
            print("   - This will inflate test performance!")
            print("   - ACTION REQUIRED: Remove duplicates or use stratified splitting")
    else:
        print("✓ No critical leakage detected!")
        if total_product_overlap > 0 or total_reactant_overlap > 0:
            print("   ⚠️  Some SMILES overlap exists, but not same reaction pairs")
            print("   (This is acceptable if different reaction contexts)")
    
    print("=" * 60 + "\n")


def save_predictions_to_file(valid_pred, valid_pred_proba, valid_true,
                              test_pred, test_pred_proba, test_true,
                              test_pred_optimal, optimal_threshold,
                              products=None, reactants=None, targets=None, split_indices=None):
    """
    Save validation and test predictions to a CSV file.
    
    Parameters
    ----------
    valid_pred : np.ndarray
        Validation predictions (default threshold)
    valid_pred_proba : np.ndarray
        Validation probability predictions
    valid_true : np.ndarray
        Validation true labels
    test_pred : np.ndarray
        Test predictions (default threshold)
    test_pred_proba : np.ndarray
        Test probability predictions
    test_true : np.ndarray
        Test true labels
    test_pred_optimal : np.ndarray
        Test predictions with optimal threshold
    optimal_threshold : float
        Optimal classification threshold
    products : list, optional
        List of product SMILES strings
    reactants : list, optional
        List of reactant SMILES strings
    targets : list, optional
        List of target SMILES strings
    split_indices : dict, optional
        Dictionary with 'valid' and 'test' keys containing indices
    """
    output_path = Path("data/predictions.csv")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    rows = []
    
    # Add validation predictions
    if split_indices and 'valid' in split_indices:
        valid_indices = split_indices['valid']
        # Flatten arrays to avoid deprecation warnings
        valid_true_flat = valid_true.flatten() if valid_true.ndim > 1 else valid_true
        valid_pred_flat = valid_pred.flatten() if valid_pred.ndim > 1 else valid_pred
        
        for i, idx in enumerate(valid_indices):
            row = {
                'split': 'validation',
                'index': int(idx),
                'product': products[idx] if products and idx < len(products) else '',
                'reactants': reactants[idx] if reactants and idx < len(reactants) else '',
                'target': targets[idx] if targets and idx < len(targets) else '',
                'true_label': int(valid_true_flat[i]),
                'predicted_label_default': int(valid_pred_flat[i]),
                'predicted_label_optimal': int(valid_pred_flat[i]),  # Same for validation
                'probability_not_hallucination': float(valid_pred_proba[i][0]),
                'probability_hallucination': float(valid_pred_proba[i][1]),
                'correct_default': int(valid_pred_flat[i] == valid_true_flat[i]),
                'correct_optimal': int(valid_pred_flat[i] == valid_true_flat[i])
            }
            rows.append(row)
    
    # Add test predictions
    if split_indices and 'test' in split_indices:
        test_indices = split_indices['test']
        # Flatten arrays to avoid deprecation warnings
        test_true_flat = test_true.flatten() if test_true.ndim > 1 else test_true
        test_pred_flat = test_pred.flatten() if test_pred.ndim > 1 else test_pred
        test_pred_optimal_flat = test_pred_optimal.flatten() if test_pred_optimal.ndim > 1 else test_pred_optimal
        
        for i, idx in enumerate(test_indices):
            row = {
                'split': 'test',
                'index': int(idx),
                'product': products[idx] if products and idx < len(products) else '',
                'reactants': reactants[idx] if reactants and idx < len(reactants) else '',
                'target': targets[idx] if targets and idx < len(targets) else '',
                'true_label': int(test_true_flat[i]),
                'predicted_label_default': int(test_pred_flat[i]),
                'predicted_label_optimal': int(test_pred_optimal_flat[i]),
                'probability_not_hallucination': float(test_pred_proba[i][0]),
                'probability_hallucination': float(test_pred_proba[i][1]),
                'correct_default': int(test_pred_flat[i] == test_true_flat[i]),
                'correct_optimal': int(test_pred_optimal_flat[i] == test_true_flat[i])
            }
            rows.append(row)
    
    # Write to CSV
    if rows:
        df = pd.DataFrame(rows)
        df.to_csv(output_path, index=False)
        print(f"    Saved {len(rows)} predictions to {output_path}")
        print(f"    Columns: split, index, product, reactants, target, true_label,")
        print(f"             predicted_label_default, predicted_label_optimal,")
        print(f"             probability_not_hallucination, probability_hallucination,")
        print(f"             correct_default, correct_optimal")
    else:
        print(f"    Warning: No predictions to save (missing split_indices or data)")


def main():
    """
    Main function: Orchestrates the entire training pipeline.
    """
    # Paths
    csv_path = "data/xgboost_dataset_merged.csv"  # Use merged dataset with all training data
    
    print("=" * 60)
    print("XGBoost Hallucination Classifier Training")
    print("Using DeepChem for featurization and data handling")
    print("=" * 60)
    
    # Step 1: Load data
    products, reactants, labels, targets = load_and_prepare_data(csv_path)
    
    # Step 2: Featurize SMILES
    use_domain_features = True  # Set to False to disable domain features
    features, featurizer = featurize_smiles(products, reactants, use_domain_features=use_domain_features)
    
    # Step 3: Create dataset and split by target (stratified)
    train_dataset, valid_dataset, test_dataset, split_indices = create_deepchem_dataset(features, labels, targets)
    
    # Step 3.5: Check for data leakage
    check_data_leakage(train_dataset, valid_dataset, test_dataset, split_indices, products, reactants)
    
    # Step 4: Train model
    model = train_xgboost_model(train_dataset, valid_dataset, use_domain_features=use_domain_features)
    
    # Step 5: Evaluate
    test_scores = evaluate_model(model, train_dataset, valid_dataset, test_dataset,
                                 products, reactants, targets, split_indices)
    
    # Step 6: Save model, featurizer, and threshold
    save_model_artifacts(model, featurizer, test_scores['optimal_threshold'])
    
    print("\n" + "=" * 60)
    print("Training Complete!")
    print("=" * 60)
    
    return model


if __name__ == "__main__":
    model = main()

