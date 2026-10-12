# Modeling Protocol

## Task definition

Train a separate binary classifier for each eligible antibiotic:

`AMR genomic features -> Resistant (1) / Susceptible (0)`

Intermediate samples are not used as binary ground truth.

## Eligibility

An antibiotic is modelled only if both classes have enough samples. The default minimum in `project.yaml` is 10 per class, but this is a minimum viability rule, not a guarantee of reliable generalization.

## Validation

Preferred initial validation: stratified 5-fold cross-validation.

Because the cohort has only 100 isolates, a single 80/20 split can be unstable. If lineage information is available, run an additional lineage-aware/grouped validation as a robustness analysis.

## Models

1. Logistic Regression — interpretable baseline.
2. Random Forest — nonlinear interactions and feature importance.
3. XGBoost — strong sparse/tabular baseline.

## Metrics

Primary: ROC-AUC.  
Also: Accuracy, Precision, Recall, F1, confusion matrix.

If class imbalance is strong, emphasize F1/Recall/Precision and consider PR-AUC as an extra metric.
