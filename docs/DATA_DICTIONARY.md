# Data Dictionary

This is the target schema. Actual source column names may differ and must be mapped during M2/M3.

## Isolate metadata

| Field | Type | Meaning |
|---|---|---|
| isolate_id | string | Canonical isolate identifier |
| biosample | string | NCBI BioSample accession if available |
| assembly_accession | string | Genome assembly accession if available |
| collection_year | integer | Collection year |
| country | string | Geographic source |
| facility | string | Source facility if published |
| sequence_type | string/int | Optional lineage/ST variable |

## AST long-format table

| Field | Type | Meaning |
|---|---|---|
| isolate_id | string | Join key |
| antibiotic | string | Standardized antibiotic name |
| phenotype_raw | string | Original AST interpretation |
| phenotype_std | category | S, I or R |
| label_binary | nullable int | S=0, R=1, I=NA |

## AMR determinant long-format table

| Field | Type | Meaning |
|---|---|---|
| isolate_id | string | Join key |
| determinant | string | Gene/allele/variant name |
| caller | string | AMRFinderPlus / ARIBA / published combined result |
| present | int | 1 if detected |

## ML matrix

Rows = isolates. Columns = AMR determinants. Values = 0/1. Target column is antibiotic-specific and must not be included among predictors.
