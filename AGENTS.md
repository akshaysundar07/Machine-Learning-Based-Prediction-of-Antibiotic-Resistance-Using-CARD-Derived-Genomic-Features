# AGENTS.md

This file tells any future coding/research agent how to work inside this project.

## Project truth

- Species: *Klebsiella pneumoniae*
- Source cohort: MRSN diversity panel, 100 isolates selected from a larger collection
- NCBI BioProject: PRJNA717739
- Primary inputs: isolate metadata, genome accessions, published AMR determinants, antibiotic susceptibility phenotypes
- Primary feature representation: AMR gene/allele presence/absence and resistance-associated variants
- Primary labels: antibiotic-specific Resistant/Susceptible
- Intermediate phenotype: preserve in source data, exclude from primary binary ML dataset
- Required models: Logistic Regression, Random Forest, XGBoost
- Optional models should not replace the required three without team approval

## Working rules

1. Never overwrite raw source data.
2. Every transformed dataset must have an upstream source documented in `data/dataset_manifest.csv`.
3. Use isolate IDs as the canonical join key whenever possible.
4. Do not convert Intermediate to Resistant or Susceptible in the primary analysis.
5. Do not report overall accuracy without class counts and at least F1/ROC-AUC.
6. Avoid leakage: feature extraction must not use the target AST label.
7. Perform preprocessing inside training folds where preprocessing is learned from data.
8. Prefer stratified cross-validation because the cohort is small.
9. If sequence type / lineage metadata is available, evaluate whether related isolates could be split across train and test; optionally add grouped validation as a robustness check.
10. Save model metrics, fold assignments and random seeds.

## File placement

- Raw downloaded files -> `data/raw/`
- External reference files -> `data/external/`
- Clean intermediate files -> `data/interim/`
- Final ML-ready tables -> `data/processed/`
- Model files -> `models/`
- Figures and final tables -> `outputs/`
- Experiment logs -> `logs/`

## Current phase

Start at M1 Data Collection unless `MODULE_TRACKER.md` marks it complete.
