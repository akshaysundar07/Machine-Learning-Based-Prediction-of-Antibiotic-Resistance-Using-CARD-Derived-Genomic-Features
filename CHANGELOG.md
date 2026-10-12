# Changelog

## 2026-09-29 — Module 5 dataset assembly and reproducible splits complete

- Assembled 16 datasets using M3 eligibility and S/R labels, preserving all 145 locked M4 genomic features and their order. CST and ineligible targets are excluded; phenotype text and metadata do not enter datasets.
- Saved configured five-fold StratifiedKFold assignments (shuffle=true, seed=42), dataset summaries and 80 fold-summary rows. Every usable isolate is assigned once; every fold contains both classes. Registered all 34 outputs with source SHA-256 provenance.
- Workspace check, 57 M5 tests and all 165 project tests passed; independent joins, feature values/order, exact sklearn assignments and byte-identical rerun verified. All 36 protected predecessor/configuration files and existing manifest records unchanged. No models trained; Module 6 not started.

## 2026-09-29 — Module 4 AMR feature engineering complete

- Encoded Table S1 published AMRFinderPlus/ARIBA gene calls and Kleborate variant annotations for 100 isolates matching M3: 138 source columns, 152 full features (136 genes/alleles, 16 variants), and 145 modeling-ready features after removing seven all-zero features. AST, phenotype summaries and metadata are excluded.
- Preserved the distinct cell-named `aac(6')-Ib` allele, combined the duplicate `blaSHV-11` columns by presence, and encoded exact semicolon-separated variant tokens. The feature dictionary preserves source names/columns; QC reports spelling differences and retained identical vectors. Rare nonconstant features remain; percentage variant tokens are not reinterpreted.
- Registered four outputs with SHA-256 provenance. Workspace check, 68 M4 tests and all 109 project tests passed; independent feature reconstruction and byte-identical rerun passed. Raw data, M3 outputs and configuration unchanged; Module 5 not started.

## 2026-09-24 — Module 3 AST processing complete

- Processed Table S1 directly: 100 unique isolates, 19 primary antibiotics, original phenotypes preserved, strict unknown-label rejection, and S=0/R=1 targets excluding I/missing. CST remains audit-only (93 missing).
- Generated cleaned long data, 19 label files, balance/eligibility and dataset-suitability ranking reports; registered source/output SHA-256 provenance. The configured 70 usable / 20 minority / 0.20 minority-fraction rule selects 16 eligible targets.
- Workspace check, 39 M3 tests and all 42 project tests passed; independent source-to-output checks and byte-identical rerun passed; raw workbook SHA-256 unchanged. Prior module statuses retained; Module 4 not started.

## 2026-09-24 — Project workspace initialized

- Finalized project around *Klebsiella pneumoniae* MRSN diversity panel.
- BioProject fixed to PRJNA717739.
- Primary task fixed to antibiotic-specific binary R/S prediction.
- Intermediate AST values retained in source data but excluded from primary binary modeling.
- Primary features fixed to published AMR genes/alleles and resistance-associated variants.
- Required models fixed to Logistic Regression, Random Forest and XGBoost.
- Created modular project structure from data collection through final reporting.
