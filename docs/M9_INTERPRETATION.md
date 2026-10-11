# M9 — Genomic Feature Interpretation

## Methods

Feature rankings were generated for 16 antibiotics and three
models using 145 genomic features.

- Logistic Regression: absolute coefficient magnitude for ranking;
  coefficient sign retained for direction.
- Random Forest: mean decrease in impurity.
- XGBoost: average split gain.
- Feature IDs mapped to biological names through the feature dictionary.

Full-data fitted models provide descriptive rankings.
Additional models fitted on each outer training fold provide
top-10 frequency, mean rank and rank variability.

Raw importance values are not comparable across model types.
Tree importance does not indicate association direction.

## Observed patterns

- blaCTX-M-15 ranks first in LR and RF for ATM, CAZ, CRO and FEP.
- aac(3)-IIe ranks first across all three GEN models.
- sul1 and sul2 rank highly across the SXT models.
- tet(A) and tet(D) rank highly across the TET models.
- Carbapenem models include blaOXA-48, blaNDM-1, KPC genes
  and, for some targets, porin-associated variants.

These are observed predictive patterns. Biological explanations
and novelty claims require a separate literature assessment.

## IPM fold stability

blaOXA-48 and blaNDM-1 occur in the top 10 in all five training
folds for every model. Their LR coefficients are positive in
all five folds.

dfrA12 occurs in the top 10 in four LR folds, five RF folds
and four XGBoost folds.

Rare features require caution. blaKPC-3 has only one carrier
in some training folds and substantial rank variability.

## IPM co-occurrence analysis

Among 15 dfrA12-positive isolates:
- 9 are IPM-resistant and 6 are susceptible.
- 7 resistant isolates carry at least one of blaOXA-48,
  blaNDM-1, blaKPC-2 or blaKPC-3.
- The other 8 isolates comprise 2 resistant and 6 susceptible isolates.

Co-occurrence may contribute to dfrA12's predictive importance.
An independent contribution has not been established.

MRSN517281 and MRSN681054 are resistant and lack these four
checked determinants. Other genes and variants must be examined
before drawing conclusions about their resistance mechanisms.

## Interpretation limits

- Importance and coefficient associations do not establish causation.
- Negative coefficients do not demonstrate biological protection.
- Correlated features and bacterial lineages can influence rankings.
- Rare features have limited supporting observations.
- The five training folds overlap; stability is not independent replication.
- Fold stability does not establish generalization to external isolates.
- XGBoost average gain can emphasize features used in few splits.
- No new resistance mechanism has been demonstrated.

## Outputs

- outputs/tables/feature_importance.csv
- outputs/tables/top10_features.csv
- outputs/tables/feature_importance_by_fold.csv
- outputs/tables/feature_stability.csv
- outputs/tables/IPM_dfrA12_cooccurrence.csv
- outputs/figures/feature_importance/

Output checks passed for biological mappings, feature ranks,
fold summaries and the presence of all 16 plots.

## Remaining interpretation work

Inspect selected isolates, compare findings with published
evidence, and distinguish independent susceptibility testing
from reuse of the original panel's phenotype data.

## Rare-allele case analysis

MRSN517281 carries blaNDM-5 and blaOXA-181, while MRSN681054
carries blaNDM-7. Both isolates are resistant to IPM, MEM and ETP.

Within the IPM dataset, each of these three alleles occurs in
only one isolate. Consequently, when its carrier is held out,
the corresponding outer training set contains zero carriers.

All three models correctly classify MRSN517281 as IPM-resistant,
but all three misclassify MRSN681054 as susceptible.

The models cannot learn a presence-versus-absence association
for these singleton alleles in the relevant training folds.
The correct prediction for MRSN517281 therefore depends on
other features and the decision rule.

These cases illustrate a limitation of allele-specific features
for rare-variant generalization. They do not demonstrate a new
resistance mechanism or establish the sole cause of the error.

A future experiment could supplement allele-specific features
with biologically curated gene-family features, evaluated using
the same outer folds and training-only model selection.

## Cross-antibiotic case comparison

Correct predictions among three models, fixed versus selected thresholds:

| Isolate | IPM | MEM | ETP |
|---|---|---|---|
| MRSN517281 | 2/3 -> 3/3 | 3/3 -> 3/3 | 3/3 -> 3/3 |
| MRSN681054 | 0/3 -> 0/3 | 2/3 -> 2/3 | 0/3 -> 3/3 |

NDM-5 and NDM-7 have zero training carriers in the relevant
folds for all three antibiotics.

OXA-181 has zero training carriers for IPM, but one resistant
training carrier for MEM and ETP. Carrier counts therefore
depend on antibiotic-specific sample eligibility.

Selected thresholds correct the RF IPM prediction for MRSN517281
and all three ETP predictions for MRSN681054. All three models
still miss IPM resistance in MRSN681054.

These cases demonstrate antibiotic-specific prediction behavior
and the limits of threshold adjustment. They do not establish
that rare alleles alone caused the errors or demonstrate a
new resistance mechanism.

Evidence tables:
- outputs/tables/rare_allele_training_coverage.csv
- outputs/tables/rare_allele_case_predictions.csv
## Literature-supported case: MRSN681054

MRSN681054 carries blaNDM-7 and has recorded resistant phenotypes
for imipenem, meropenem and ertapenem.

Shao et al. (2020) studied a different blaNDM-7-positive
K. pneumoniae isolate, B14 (ST11), and reported MICs of
>=32 micrograms/mL for all three carbapenems.

This provides supporting biological evidence, rather than
independent phenotype confirmation for MRSN681054 (ST340).

In our cross-validation experiment, blaNDM-7 had no carriers
in this isolate's outer training fold. All three models missed
its imipenem resistance. For ertapenem, training-only threshold
selection changed all three predictions from susceptible to resistant.

These observations illustrate a potential limitation of rare
allele representation in a small dataset. They do not establish
that the unseen allele alone caused the prediction errors or
demonstrate a new resistance mechanism.

Reference:
Shao et al. (2020). Genotypic and Phenotypic Characterization
of blaNDM-7-Harboring IncX3 Plasmid in a ST11 Klebsiella
pneumoniae Isolated From a Pediatric Patient in China.
https://doi.org/10.3389/fmicb.2020.576823