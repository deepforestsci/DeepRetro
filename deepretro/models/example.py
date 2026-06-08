import os
import json
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from deepretro.models.hallucination_trainer import HallucinationTrainer
from deepretro.models.hallucination_checker import HallucinationChecker

# Load data
df = pd.read_csv("/home/riya/Documents/Chiron/DeepRetro/data/hallucination_dataset.csv")

# Stratified split
train_df, test_df = train_test_split(
    df,
    test_size=0.2,
    random_state=42,
    stratify=df["label"]
)

# Save to CSV
train_df.to_csv("train.csv", index=False)
test_df.to_csv("test.csv", index=False)

print("Saved train.csv and test.csv")

WORKSPACE_DIR = "./hallucination_workspace"

# -----------------------------------------------------------------
# STAGE 1: THE TRAINING & OPTIMIZATION PHASE (Trainer)
# -----------------------------------------------------------------
print("====================================================")
print("STAGE 1: INITIALIZING TRAINER & TUNING HYPERPARAMETERS")
print("====================================================")

# Initialize the trainer infrastructure for an XGBoost classification task
# Note: Target columns default to 1 task ('is_hallucination')
trainer = HallucinationTrainer(
    trainer_dir=WORKSPACE_DIR,
    model_type="xgboost",
    n_tasks=1
)

# Load and featurize the explicit training and testing partitions.
# CSVLoader parses the compound tuples through our internal ReactionStepFeaturizer
train_dataset, test_dataset = trainer.load_dataset(
    train_csv="train.csv",
    test_csv="test.csv",
    product_col="product",
    reactants_col="reactants",
    label_col="label"
)

# Execute training. We flag 'tune_params=True' and specify 'k_folds=2' to leverage 
# our custom KFoldRandomHyperparamOpt optimizer child class.
print("\nStarting automated optimization loop...")
final_model, test_performance = trainer.train_model(
    train_dataset=train_dataset,
    test_dataset=test_dataset,
    tune_params=True,
    n_trials=3,   # Kept low for execution speed in this example
    k_folds=2     # Runs cross-validation to isolate parameters and evaluate average optimal threshold
)

print("\n--- Training Stage Complete ---")
print(f"Final Test Partition Scores: {test_performance}")
print(f"Production Artifact Directory: {trainer.model_dir}\n")


# -----------------------------------------------------------------
# STAGE 2: THE INFRASTRUCTURE INFERENCE PHASE (Checker)
# -----------------------------------------------------------------
print("====================================================")
print("STAGE 2: DEPLOYING LIVE INFERENCE CHECKER Pipeline")
print("====================================================")

# Point the deployment Checker strictly to the artifacts directory.
# Setting checker_type="ml" tells it to parse and reconstruct the exact model
# (XGBoost vs GCN) and custom featurizer settings from the generated config.json.
production_model_path = os.path.join(WORKSPACE_DIR, "xgboost_model")

ml_checker = HallucinationChecker(
    checker_type="ml",
    model_path=production_model_path
)

# Validate a completely new out-of-distribution pathway step
product = "CCO"
proposed_reactants = "CC.N"

print(f"\nEvaluating real-time pathway proposal:")
print(f" -> Product:   {product}")
print(f" -> Reactants: {proposed_reactants}")

# Run the live assessment loop
prediction = ml_checker.check_single_pathway(product, proposed_reactants)

print("\n--- Inference Result ---")
if prediction == 1:
    print("RESULT: [!] CHOSEN STEP DETECTED AS AN UNSTABLE/HALLUCINATED PATHWAY.")
else:
    print("RESULT: [✓] CHOSEN STEP CONFIRMED AS A VALID PATHWAY.")
    
# Demonstration of the clean decoupling: Heuristic checking fallback option
print("\nBonus: Quick swap checking via fallback heuristic processor:")
heur_checker = HallucinationChecker(checker_type="heur")
heur_pred = heur_checker.check_single_pathway("C1=CC=CC=C1", "NNN.P") # Will match dummy heuristic criteria
print(f"Heuristic code prediction: {heur_pred}")
