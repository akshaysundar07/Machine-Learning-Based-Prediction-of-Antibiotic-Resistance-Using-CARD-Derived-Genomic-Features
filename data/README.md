# Data Area

## Raw data rule

`data/raw/` is immutable. Never edit source files in place.

Expected source categories:

1. **metadata/** — 100-isolate metadata and genome accessions.
2. **ast/** — published antibiotic susceptibility phenotypes.
3. **amr/** — published AMRFinderPlus/ARIBA AMR gene/allele calls.
4. **genomes/** — optional FASTA files if sequence-level validation is later added.

## Manifest

Every downloaded source file must be registered in `dataset_manifest.csv` with source URL/accession, retrieval date, description and checksum where possible.

## Intermediate phenotype policy

Do not delete Intermediate records from raw data. Preserve them, then exclude them only when generating the primary binary R/S modeling tables.
