# Machine Learning-Based Antibiotic Resistance Prediction in *Klebsiella pneumoniae*

**Course:** DNA Sequencing  
**Team:** 20  
**Dataset:** MRSN *Klebsiella pneumoniae* diversity panel  
**NCBI BioProject:** PRJNA717739  
**Project objective:** Predict antibiotic-specific Resistant/Susceptible phenotypes from published genomic AMR determinants using interpretable machine-learning models.

## Core project decision

The main model is **binary classification**:

- Resistant (R) = 1
- Susceptible (S) = 0
- Intermediate (I) = retained in raw/interim data but excluded from the primary binary training/evaluation dataset

The project will use the **published AMR genes/alleles already identified in the original study using AMRFinderPlus and ARIBA** as the primary genomic features. Re-running AMR callers is an optional validation extension, not required for the first implementation.

## Main pipeline

```text
Published isolate metadata + genome accessions + AST phenotypes + AMR determinants
                                 |
                                 v
                     Data audit and standardization
                                 |
                                 v
          Antibiotic-specific R/S labels + binary AMR feature matrix
                                 |
                                 v
             Stratified split / stratified cross-validation
                                 |
                                 v
              Logistic Regression / Random Forest / XGBoost
                                 |
                                 v
        Accuracy / Precision / Recall / F1 / ROC-AUC / Confusion Matrix
                                 |
                                 v
                 Feature importance / biological interpretation
```

## Module order

| Module | Name | Main output |
|---|---|---|
| M0 | Project setup & provenance | Reproducible workspace |
| M1 | Data collection | Raw source files + manifest |
| M2 | Data audit & QC | Validated isolate-level tables |
| M3 | AST phenotype processing | Clean antibiotic-specific R/S labels |
| M4 | AMR feature engineering | Binary AMR feature matrix |
| M5 | ML-ready dataset assembly | Per-antibiotic datasets and splits |
| M6 | Baseline modeling | Logistic Regression results |
| M7 | Tree models | Random Forest and XGBoost results |
| M8 | Evaluation | Metrics, confusion matrices, ROC curves |
| M9 | Explainability | Ranked resistance determinants |
| M10 | Final reporting | Figures, tables, report-ready outputs |

Read `docs/MODULE_PLAN.md` before starting. Update `MODULE_TRACKER.md` after each completed module.

## Recommended first command

```bash
python scripts/00_check_workspace.py
```

Then begin with **Module 1** only. Do not move to later modules until the preceding module passes its acceptance checks.
