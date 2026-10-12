# Directory Map

```text
KP_AMR_Prediction_Project/
├── README.md                         # Project overview and start point
├── AGENTS.md                         # Rules for future agents/coders
├── MODULE_TRACKER.md                 # Module completion status
├── DIRECTORY_MAP.md                  # This file
├── CHANGELOG.md                      # Record major project changes
├── .gitignore
├── requirements.txt
├── environment.yml
├── Makefile
│
├── configs/
│   └── project.yaml                  # Central project/model configuration
│
├── data/
│   ├── README.md
│   ├── dataset_manifest.csv          # Source and derived-dataset provenance + SHA-256
│   ├── external/                     # CARD/reference resources if later needed
│   ├── raw/                          # NEVER MODIFY
│   │   ├── metadata/                 # Original isolate metadata
│   │   ├── genomes/                  # Optional FASTA files / accession list
│   │   ├── amr/                      # Published AMR genes/alleles
│   │   ├── ast/                      # Original AST phenotype tables
│   │   └── published_study/
│   │       ├── mgen-9-967-s001.pdf
│   │       ├── mgen-9-967-s002.xlsx   # Table S1: M3/M4 source of truth
│   │       ├── mgen-9-967-s003.xlsx
│   │       └── mgen-9-967-s004.xlsx
│   ├── interim/                      # Cleaned, standardized data
│   │   ├── ast/ast_clean_long.csv     # Original/standardized AST; I/missing retained; CST audit-only
│   │   └── amr/amr_feature_matrix_full.csv # M4: all encoded genomic features, including constants
│   ├── processed/
│   │   ├── master/                   # Master isolate/AST/AMR tables
│   │   │   └── amr_feature_matrix.csv # M4: binary genomic features; zero-variance features removed
│   │   ├── per_antibiotic/           # M3 binary labels per antibiotic
│   │   │   └── *_labels.csv           # M3: 19 S/R target files; CT_labels.csv retains C/T internally
│   │   ├── ml_datasets/
│   │   │   └── *_ml.csv               # M5: 16 eligible targets; 145 locked genomic features + binary target
│   │   └── splits/                    # Saved CV/train-test fold assignments
│   │       └── *_folds.csv            # M5: five stratified validation folds (0-4), seed 42
│   └── reports/                      # Audit/missingness/class balance reports
│       ├── ast_binary_summary.csv     # Counts, eligibility thresholds and reasons; CST audit-only
│       ├── ast_antibiotic_ranking.csv # Dataset suitability; eligible targets first
│       ├── amr_feature_dictionary.csv # M4: safe IDs, exact biological names, sources and retention
│       ├── amr_feature_summary.csv    # M4: counts, encoding decisions and genomic feature QC
│       ├── ml_dataset_summary.csv    # M5: dimensions, source paths, class counts and CV settings
│       └── cv_fold_summary.csv       # M5: class balance for each antibiotic/validation fold
│
├── docs/
│   ├── MODULE_PLAN.md                # Detailed module-by-module workflow
│   ├── PROJECT_ARCHITECTURE.md
│   ├── DATA_DICTIONARY.md
│   ├── DATA_COLLECTION_CHECKLIST.md
│   ├── MODELING_PROTOCOL.md
│   └── DECISIONS.md
│
├── notebooks/
│   ├── 01_data_audit.ipynb
│   ├── 02_ast_exploration.ipynb
│   ├── 03_feature_matrix_qc.ipynb
│   ├── 04_model_baselines.ipynb
│   └── 05_results_review.ipynb
│
├── scripts/
│   ├── 00_check_workspace.py
│   ├── 01_build_manifest.py
│   ├── 02_audit_raw_data.py
│   ├── 03_clean_ast.py               # Existing M3 implementation/entry point
│   ├── 03_process_ast.py             # Preferred M3 command, delegates to 03_clean_ast.py
│   ├── 04_build_amr_features.py
│   ├── 05_build_ml_datasets.py
│   ├── 06_train_models.py
│   ├── 07_evaluate_models.py
│   └── 08_export_results.py
│
├── src/kpamr/
│   ├── __init__.py
│   ├── config.py
│   ├── data/
│   │   ├── __init__.py
│   │   ├── io.py
│   │   ├── ast.py
│   │   └── validation.py
│   ├── features/
│   │   ├── __init__.py
│   │   └── amr.py
│   ├── modeling/
│   │   ├── __init__.py
│   │   ├── split.py
│   │   ├── train.py
│   │   └── evaluate.py
│   └── utils/
│       ├── __init__.py
│       └── paths.py
│
├── tests/
│   ├── test_ast.py
│   ├── test_ast_processing.py        # M3 schema, labels, counts, ranking and reproducibility
│   ├── test_amr_features.py
│   └── test_split.py
│
├── models/
│   ├── checkpoints/
│   ├── metrics/
│   └── explainability/
│
├── outputs/
│   ├── figures/
│   ├── tables/
│   └── reports/
│
└── logs/
```
