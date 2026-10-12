from pathlib import Path

import pandas as pd
from sklearn.metrics import (
    accuracy_score,
    precision_score,
    recall_score,
    f1_score,
)

ROOT = Path(__file__).resolve().parents[1]
rows = []

for model in ["logistic_regression", "random_forest", "xgboost"]:
    data = pd.read_csv(
        ROOT / f"models/metrics/{model}_oof_predictions.csv"
    )

    for (antibiotic, fold), group in data.groupby(["antibiotic", "fold"]):
        for rule in ["fixed_0.5", "selected_threshold"]:
            predictions = (
                (group.y_prob >= 0.5).astype(int)
                if rule == "fixed_0.5"
                else group.y_pred.astype(int)
            )

            rows.append({
                "antibiotic": antibiotic,
                "model": model,
                "fold": fold,
                "rule": rule,
                "threshold": (
                    0.5 if rule == "fixed_0.5"
                    else float(group.threshold.iloc[0])
                ),
                "accuracy": accuracy_score(group.y_true, predictions),
                "precision": precision_score(
                    group.y_true, predictions, zero_division=0
                ),
                "recall": recall_score(
                    group.y_true, predictions, zero_division=0
                ),
                "f1": f1_score(
                    group.y_true, predictions, zero_division=0
                ),
            })

fold_results = pd.DataFrame(rows)
metrics = ["accuracy", "precision", "recall", "f1"]

means = fold_results.groupby(
    ["antibiotic", "model", "rule"]
)[metrics].mean()

fixed = means.xs("fixed_0.5", level="rule")
selected = means.xs("selected_threshold", level="rule")

comparison = (
    fixed.add_suffix("_fixed")
    .join(selected.add_suffix("_selected"))
)

for metric in metrics:
    comparison[f"{metric}_change"] = (
        comparison[f"{metric}_selected"]
        - comparison[f"{metric}_fixed"]
    )

destination = ROOT / "outputs/tables"
destination.mkdir(parents=True, exist_ok=True)

fold_results.to_csv(
    destination / "threshold_comparison_folds.csv", index=False
)
comparison.reset_index().to_csv(
    destination / "threshold_comparison_by_antibiotic.csv", index=False
)

counts = comparison.groupby(level="model")["f1_change"].agg(
    improved=lambda values: int((values > 1e-12).sum()),
    unchanged=lambda values: int((values.abs() <= 1e-12).sum()),
    declined=lambda values: int((values < -1e-12).sum()),
)

print("Number of antibiotics with higher/lower mean F1:")
print(counts.to_string())
print("\nSaved both comparison tables in outputs/tables.")