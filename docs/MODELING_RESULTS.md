@'
# Modeling Results: M6–M8

## Implementation

Logistic Regression, Random Forest and XGBoost were evaluated
for 16 antibiotics using the stored M5 five-fold stratified splits.

ROC-AUC remains the primary metric. Hyperparameters are fixed.
Class weighting is determined from each training fold.

## Threshold experiment

Both fixed 0.5 and selected-threshold results are retained.
Threshold selection maximizes resistance F1 using three-fold
inner cross-validation within each outer training fold.
Outer validation data are excluded from threshold selection.

## Mean results across antibiotics

| Model | Accuracy fixed/selected | Recall fixed/selected | F1 fixed/selected |
|---|---|---|---|
| Logistic Regression | 0.8524 / 0.8396 | 0.7518 / 0.8241 | 0.7850 / 0.7952 |
| Random Forest | 0.8423 / 0.8299 | 0.6792 / 0.8593 | 0.7424 / 0.7845 |
| XGBoost | 0.8146 / 0.7995 | 0.7049 / 0.8109 | 0.7305 / 0.7506 |

Selected thresholds improved mean F1 for 10/16 LR targets,
14/16 RF targets and 11/16 XGBoost targets.
Individual targets can decline; statistical significance
has not been established.

Comparison figures use selected-threshold predictions.
Both rules are saved in outputs/tables/threshold_comparison_folds.csv
and outputs/tables/threshold_comparison_by_antibiotic.csv.

Grouped MLST and K-locus validation are supplementary analyses.
They do not replace the primary stored M5 folds.

## Reproduction

Run from the project root after installing project dependencies:

```bash
python scripts/06_train_models.py
python scripts/07_evaluate_models.py
python scripts/09_compare_thresholds.py
python -m pytest tests/test_modeling.py -q