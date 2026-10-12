# Project Decisions

## D001 — Species
Use *Klebsiella pneumoniae* only.

## D002 — Cohort
Use the published 100-isolate MRSN diversity panel associated with PRJNA717739.

## D003 — Features
Primary features come from published AMRFinderPlus/ARIBA AMR gene/allele calls. Genome re-analysis is optional later.

## D004 — Labels
Antibiotic-specific AST phenotypes are the ground truth.

## D005 — Intermediate phenotype
Preserve Intermediate records in the source/interim data but exclude them from primary R/S binary model training and evaluation.

## D006 — Models
Required comparison: Logistic Regression, Random Forest, XGBoost. SVM and deep learning are not part of the core plan.

## D007 — Validation
Prefer stratified cross-validation because the dataset is small. Add lineage-aware validation if lineage/ST fields permit it.
