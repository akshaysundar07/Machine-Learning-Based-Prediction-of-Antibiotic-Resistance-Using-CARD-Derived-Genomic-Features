# Module Tracker

Update status only after acceptance criteria are satisfied.

| Module | Status | Acceptance criteria | Evidence / file |
|---|---|---|---|
| M0 Setup | COMPLETE | Directory exists; config loads; workspace check passes | `configs/project.yaml` |
| M1 Data Collection | COMPLETE | 100-isolate source files and provenance manifest collected | `data/dataset_manifest.csv` |
| M2 Data Audit & QC | COMPLETE | Unique isolate IDs, joins validated, missingness report generated | `data/reports/data_audit.csv` |
| M3 AST Processing | COMPLETE | 100 unique isolates; 19 primary antibiotics; S=0/R=1; I preserved/excluded from binary targets; counts, eligibility and ranking validated; CST audit-only; raw SHA-256 unchanged; deterministic rerun; 39 M3 / 42 total tests passed | `data/interim/ast/ast_clean_long.csv`; `data/processed/per_antibiotic/*_labels.csv`; `data/reports/ast_binary_summary.csv`; `data/reports/ast_antibiotic_ranking.csv`; `tests/test_ast_processing.py` |
| M4 Feature Engineering | COMPLETE | 100 unique isolates matching M3; 138 source columns -> 152 encoded features (136 genes/alleles, 16 variants); 7 all-zero features removed -> 145 modeling-ready features; no AST/metadata leakage; dictionary/provenance and byte-identical rerun verified; raw/M3 unchanged; 68 M4 / 109 total tests passed | `data/interim/amr/amr_feature_matrix_full.csv`; `data/processed/master/amr_feature_matrix.csv`; `data/reports/amr_feature_dictionary.csv`; `data/reports/amr_feature_summary.csv`; `tests/test_amr_features.py` |
| M5 Dataset Assembly | COMPLETE | 16 M3-eligible targets; all 145 M4 features/values/order preserved; exact S/R joins; configured 5-fold stratified CV (shuffle=true, seed=42); 80 fold-summary rows, both classes in every fold; provenance and byte-identical rerun verified; M3/M4 unchanged; 57 M5 / 165 total tests passed | `data/processed/ml_datasets/*_ml.csv`; `data/processed/splits/*_folds.csv`; `data/reports/ml_dataset_summary.csv`; `data/reports/cv_fold_summary.csv`; `tests/test_split.py` |
| M6 Logistic Regression | BLOCKED | CV metrics + coefficients saved | `models/metrics/logistic_regression.csv` |
| M7 RF + XGBoost | BLOCKED | CV metrics + trained models saved | `models/metrics/` |
| M8 Evaluation | BLOCKED | Comparison table, ROC-AUC, F1 and confusion matrices generated | `outputs/tables/model_comparison.csv` |
| M9 Explainability | BLOCKED | Ranked AMR determinants per antibiotic | `outputs/tables/feature_importance.csv` |
| M10 Reporting | BLOCKED | Final figures/tables exported and reproducible | `outputs/reports/` |
