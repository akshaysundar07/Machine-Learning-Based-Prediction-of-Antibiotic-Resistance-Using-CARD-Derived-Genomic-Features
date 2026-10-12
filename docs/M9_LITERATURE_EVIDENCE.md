# M9: Isolate case studies and literature evidence

Scope: two selected cases from the 100-isolate panel. This is a targeted literature comparison, not a systematic novelty review or independent external model validation. Local results below are taken from the project's reported outputs; the scripts and checkpoint files were not independently rerun for this document.

## Isolate identities and observed phenotypes

| Isolate | BioSample | GenBank assembly | MLST | Relevant recorded determinants | IPM / MEM / ETP |
|---|---|---|---|---|---|
| MRSN517281 | SAMN18874744 | GCA_019928005.1 | 147 | blaNDM-5, blaOXA-181; also recorded porin variants OmpK35-76% and OmpK36TD | R / R / R |
| MRSN681054 | SAMN18874767 | GCA_019927525.1 | 340 | blaNDM-7; other determinants include blaCTX-M-15 | R / R / R |

IPM = imipenem; MEM = meropenem; ETP = ertapenem. R/S labels are experimental phenotypes inherited from the source cohort, not newly measured in this project. Gene calls are recorded genomic features, not new functional assays.

## Literature evidence table

| Project case | Reference and location | Published evidence | Relationship to our case | Permitted conclusion |
|---|---|---|---|---|
| Both cases: cohort provenance | Martin et al. (2023), A panel of diverse Klebsiella pneumoniae clinical isolates for research and development. DOI: 10.1099/mgen.0.000967 | Publication describing the source 100-isolate MRSN panel and its phenotypic/genomic characterization | Source cohort used by this project | Establishes provenance; reuse of these labels is not independent validation |
| MRSN517281 | Rojas et al. (2017), NDM-5 and OXA-181 Beta-Lactamases, a Significant Threat Continues To Spread in the Americas. DOI: 10.1128/AAC.00454-17; abstract and isolate characterization | Carbapenem-resistant K. pneumoniae ST147 carrying blaNDM-5 and blaOXA-181 | Same species, sequence type and gene combination; the publication does not identify this case as MRSN517281 | Supports biological plausibility; matching MLST and genes does not establish identical isolate identity |
| MRSN681054 | Shao et al. (2020), Genotypic and Phenotypic Characterization of blaNDM-7-Harboring IncX3 Plasmid in a ST11 Klebsiella pneumoniae Isolated From a Pediatric Patient in China. DOI: 10.3389/fmicb.2020.576823; Table 1 and Results | B14 (ST11) carried blaNDM-7; imipenem, meropenem and ertapenem MICs were each >=32 micrograms/mL. Transfer of its NDM-7 plasmid produced a carbapenem-resistant transconjugant | Same species and allele, different isolate and sequence type | Supports a resistance role for NDM-7 in a separate experimental context; does not independently verify MRSN681054 |
| Both exact strain names | Ali et al., Engineered antibodies bypass bacterial immune evasion to drive complement-mediated protection against lethal infections; Cambridge repository manuscript, Table S1, PDF page 40 | MRSN517281 and MRSN681054 are listed as MDR clones among strains used for an antibody ELISA binding assay. Classification uses a Kleborate resistance score or MDR-associated lineage | Exact strain-name matches; this table gives classification rather than drug-specific MICs | Confirms use of these named strains in another study; the table cannot be used as an independent IPM/MEM/ETP susceptibility confirmation |

No independent drug-specific AST confirmation for these exact two isolates was established in this targeted search. A paper mentioning the same strain or genes is not automatically an independent phenotype measurement.

## Training coverage

| Target | Isolate | Allele | Total carriers in target dataset | Carriers in this case's outer training fold |
|---|---|---|---|---|
| IPM | MRSN517281 | blaNDM-5 | 1 | 0 |
| IPM | MRSN517281 | blaOXA-181 | 1 | 0 |
| MEM and ETP (each) | MRSN517281 | blaNDM-5 | 1 | 0 |
| MEM and ETP (each) | MRSN517281 | blaOXA-181 | 2 | 1 |
| IPM, MEM and ETP (each) | MRSN681054 | blaNDM-7 | 1 | 0 |

These allele columns provided no observed positive training examples where training coverage was zero. Other genomic features could still support predictions. Missing allele coverage does not by itself prove the cause of a prediction error.

## Out-of-fold prediction comparison

All actual labels below are resistant. Entries show fixed-0.5 prediction -> training-only selected-threshold prediction. These compare the same predicted probabilities; threshold selection did not retrain a different classifier for this comparison.

| Isolate | Target | Logistic regression | Random Forest | XGBoost |
|---|---|---|---|---|
| MRSN517281 | IPM | R -> R | S -> R | R -> R |
| MRSN517281 | MEM | R -> R | R -> R | R -> R |
| MRSN517281 | ETP | R -> R | R -> R | R -> R |
| MRSN681054 | IPM | S -> S | S -> S | S -> S |
| MRSN681054 | MEM | R -> R | S -> S | R -> R |
| MRSN681054 | ETP | S -> R | S -> R | S -> R |

For MRSN517281, threshold selection corrected the Random Forest IPM prediction. For MRSN681054, it corrected all three ETP predictions but did not correct any IPM prediction or the Random Forest MEM prediction. These cases were selected after inspection and do not estimate generalization performance on their own; overall cross-validation results remain the performance evidence.

## Interpretation and contribution

The case studies connect allele-specific training coverage, known resistance biology and out-of-fold errors. They illustrate that training-only threshold selection can correct some resistant predictions, while errors can persist in isolates carrying rare resistance alleles. Correct predictions despite absent allele-specific training examples also demonstrate the possible contribution of other features.

The project has not discovered a new carbapenemase or established a new resistance mechanism. Its current contribution is a reproducible analysis of threshold trade-offs and rare-allele coverage limitations in this cohort. Novelty relative to prior ML studies is not established by this targeted biological literature search.

A possible future experiment is to add predefined gene-family features (for example, any NDM allele) alongside allele-specific columns and compare them using the same leakage-controlled evaluation. Improvement is a hypothesis, not a demonstrated result. Because these cases have already been inspected, a confirmatory claim would benefit from an independent cohort.

## Sources

1. Martin et al. (2023): https://doi.org/10.1099/mgen.0.000967
2. Rojas et al. (2017): https://doi.org/10.1128/AAC.00454-17
3. Shao et al. (2020): https://doi.org/10.3389/fmicb.2020.576823
4. Ali et al., Cambridge repository manuscript (publication date/status not verified here): https://api.repository.cam.ac.uk/server/api/core/bitstreams/616c12f9-3cbc-4422-a45a-fea4f8bd353a/content

Search date: 11 October 2026. Targeted queries included both exact MRSN strain names and the relevant resistance alleles. This evidence table is not an exhaustive literature review.
