# Project Architecture

```text
SOURCE LAYER
  PRJNA717739 + publication supplementary files
       |
       +-- isolate metadata / accessions
       +-- AST phenotypes
       +-- published AMR determinants
       |
       v
RAW IMMUTABLE LAYER
  data/raw/*
       |
       v
STANDARDIZATION LAYER
  isolate keys / antibiotic names / S-I-R / AMR names
       |
       v
FEATURE + LABEL LAYER
  X = binary AMR determinant matrix
  y = antibiotic-specific R/S phenotype
       |
       v
VALIDATION LAYER
  stratified CV folds (and optional lineage-aware robustness split)
       |
       v
MODELING LAYER
  Logistic Regression | Random Forest | XGBoost
       |
       v
EVALUATION + INTERPRETATION
  metrics | confusion matrices | ROC | feature importance
       |
       v
FINAL OUTPUTS
  tables | figures | report-ready evidence
```

## Important architecture decision

The initial project does **not** require de novo resistance-gene calling from all genomes because the selected publication already provides AMRFinderPlus/ARIBA-derived AMR determinants. This reduces unnecessary computation and keeps the first implementation aligned with the published panel.
