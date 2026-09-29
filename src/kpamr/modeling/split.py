"""M5 assembly and CV: preserve locked genomic columns, use published binary labels."""

from dataclasses import dataclass
from pathlib import Path
import csv
import re

import pandas as pd
from sklearn.model_selection import StratifiedKFold

FEATURE_PATH = Path("data/processed/master/amr_feature_matrix.csv")
DICTIONARY_PATH = Path("data/reports/amr_feature_dictionary.csv")
ELIGIBILITY_PATH = Path("data/reports/ast_binary_summary.csv")
RANKING_PATH = Path("data/reports/ast_antibiotic_ranking.csv")
CONFIG_PATH = Path("configs/project.yaml")
FEATURE_COUNT = 145
COHORT_SIZE = 100


def read_csv_checked(path: Path) -> pd.DataFrame:
    """Reject duplicate CSV headers before pandas can silently rename them."""
    with Path(path).open(encoding="utf-8", newline="") as stream:
        header = next(csv.reader(stream), [])
    if not header or len(header) != len(set(header)):
        raise ValueError(f"Empty or duplicate CSV header: {path}")
    return pd.read_csv(path, dtype={"isolate_id": str}, keep_default_na=False)


def validate_ids(ids: pd.Series) -> None:
    if ids.empty or ids.isna().any() or ids.duplicated().any():
        raise ValueError("Missing or duplicate isolate IDs")
    if not ids.map(lambda value: isinstance(value, str) and bool(value) and value == value.strip()).all():
        raise ValueError("Isolate IDs must be nonempty canonical strings")


def validate_binary(values, description: str) -> None:
    if values.isna().to_numpy().any() or not values.isin([0, 1]).to_numpy().all():
        raise ValueError(f"{description} must contain only 0/1 with no missing values")


def validate_features(features: pd.DataFrame) -> list[str]:
    if not features.columns.is_unique or features.columns[0] != "isolate_id":
        raise ValueError("Invalid locked M4 column schema")
    columns = features.columns[1:].tolist()
    if len(features) != COHORT_SIZE or len(columns) != FEATURE_COUNT:
        raise ValueError("Locked M4 matrix must have 100 isolates and 145 features")
    if not all(re.fullmatch(r"(?:gene|variant)__[A-Za-z0-9_]+", name) for name in columns):
        raise ValueError("Non-genomic or unsafe column in locked M4 matrix")
    validate_ids(features.isolate_id)
    validate_binary(features[columns], "Genomic features")
    return columns


def validate_locked_features(features: pd.DataFrame, dictionary: pd.DataFrame) -> list[str]:
    columns = validate_features(features)
    required = {"feature_id", "feature_type", "retained_in_model_matrix"}
    if required - set(dictionary.columns) or dictionary.feature_id.duplicated().any():
        raise ValueError("Invalid M4 feature dictionary")
    if not pd.api.types.is_bool_dtype(dictionary.retained_in_model_matrix) or dictionary.retained_in_model_matrix.isna().any():
        raise ValueError("Dictionary retention flags must be boolean")
    retained = dictionary.loc[dictionary.retained_in_model_matrix]
    if retained.feature_id.tolist() != columns:
        raise ValueError("Locked M4 feature names/order disagree with dictionary")
    if not retained.feature_type.isin(["gene_or_allele", "resistance_variant"]).all():
        raise ValueError("Dictionary contains non-genomic feature types")
    return columns


def eligible_targets(summary: pd.DataFrame, ranking: pd.DataFrame) -> pd.DataFrame:
    """Use M3 eligibility flags verbatim; cross-check both reports without rethresholding."""
    fields = ["antibiotic", "primary_antibiotic", "total_isolates", "susceptible_count",
              "intermediate_count", "resistant_count", "missing_count",
              "usable_binary_count", "eligible_for_ml"]
    for report in (summary, ranking):
        if not report.columns.is_unique or set(fields) - set(report.columns):
            raise ValueError("Missing or duplicate M3 eligibility report columns")
        if report.antibiotic.isna().any() or report.antibiotic.duplicated().any():
            raise ValueError("Missing/duplicate antibiotics in M3 report")
        for field in ("primary_antibiotic", "eligible_for_ml"):
            if not pd.api.types.is_bool_dtype(report[field]) or report[field].isna().any():
                raise ValueError("M3 eligibility/primary flags must be boolean")
        if len(report) != 20 or report.primary_antibiotic.sum() != 19:
            raise ValueError("Expected 19 primary antibiotics and audit-only CST")
        if report.loc[~report.primary_antibiotic, "antibiotic"].tolist() != ["CST"]:
            raise ValueError("CST must be the audit-only antibiotic")
        if report.loc[~report.primary_antibiotic, "eligible_for_ml"].any():
            raise ValueError("CST/audit-only targets cannot be eligible")
        for field in fields[2:-1]:
            values = report[field]
            if not pd.api.types.is_numeric_dtype(values) or values.isna().any() or (values < 0).any() or (values % 1 != 0).any():
                raise ValueError(f"Invalid M3 count column: {field}")
        if not report.total_isolates.eq(COHORT_SIZE).all():
            raise ValueError("M3 totals must reconcile to 100 isolates")
        if not report[["susceptible_count", "intermediate_count", "resistant_count", "missing_count"]].sum(axis=1).eq(COHORT_SIZE).all():
            raise ValueError("M3 phenotype counts do not reconcile")
        if not (report.susceptible_count + report.resistant_count).eq(report.usable_binary_count).all():
            raise ValueError("M3 usable counts do not reconcile")
    left = summary[fields].set_index("antibiotic").sort_index()
    right = ranking[fields].set_index("antibiotic").sort_index()
    if not left.equals(right):
        raise ValueError("M3 summary and ranking disagree")
    if "initial_experiment_rank" not in ranking:
        raise ValueError("Missing M3 initial experiment ranks")
    selected = ranking.loc[ranking.eligible_for_ml].copy()
    if selected.empty:
        raise ValueError("No M3-eligible primary targets")
    ranks = pd.to_numeric(selected.initial_experiment_rank, errors="raise")
    if ranks.tolist() != list(range(1, len(selected) + 1)):
        raise ValueError("M3 eligible ranks must be consecutive and ordered")
    return selected.reset_index(drop=True)


def antibiotic_stem(antibiotic: str) -> str:
    if antibiotic == "C/T":
        return "CT"
    if antibiotic == "CST" or not isinstance(antibiotic, str) or not re.fullmatch(r"[A-Z]{3}", antibiotic):
        raise ValueError(f"Unsupported primary antibiotic: {antibiotic!r}")
    return antibiotic


def assemble_dataset(features: pd.DataFrame, labels: pd.DataFrame, stats: pd.Series) -> pd.DataFrame:
    columns = validate_features(features)
    antibiotic = stats["antibiotic"]
    antibiotic_stem(antibiotic)
    if not stats["eligible_for_ml"] or not stats["primary_antibiotic"]:
        raise ValueError("Only M3-eligible primary targets may be assembled")
    required = {"isolate_id", "phenotype", "binary_label"}
    if (not labels.columns.is_unique or required - set(labels.columns)
            or set(labels.columns) - (required | {"antibiotic"})):
        raise ValueError("Unexpected M3 label-file schema")
    validate_ids(labels.isolate_id)
    validate_binary(labels.binary_label, "Binary labels")
    if not labels.phenotype.isin(["S", "R"]).all():
        raise ValueError("M3 binary labels must contain only S/R phenotypes")
    if not labels.phenotype.map({"S": 0, "R": 1}).eq(labels.binary_label).all():
        raise ValueError("M3 phenotype and binary label disagree")
    if "antibiotic" in labels and not labels.antibiotic.eq(antibiotic).all():
        raise ValueError("Label-file antibiotic disagrees with M3 report")
    if len(labels) != stats["usable_binary_count"]:
        raise ValueError("Label row count differs from M3 usable_binary_count")
    for target, field in ((0, "susceptible_count"), (1, "resistant_count")):
        if int(labels.binary_label.eq(target).sum()) != stats[field]:
            raise ValueError(f"Label class counts disagree with M3: {field}")
    extra = set(labels.isolate_id) - set(features.isolate_id)
    if extra:
        raise ValueError(f"Label isolates missing from locked M4 matrix: {sorted(extra)}")
    # A validated one-to-one join follows M4 row order; features never depend on y.
    dataset = features.merge(labels[["isolate_id", "binary_label"]], on="isolate_id",
                             how="inner", sort=False, validate="one_to_one")
    if len(dataset) != stats["usable_binary_count"]:
        raise ValueError("Inner join lost expected usable isolates")
    dataset.insert(1, "antibiotic", antibiotic)
    dataset = dataset[["isolate_id", "antibiotic", *columns, "binary_label"]]
    validate_binary(dataset[columns], "Joined features")
    validate_binary(dataset.binary_label, "Joined target")
    return dataset


@dataclass(frozen=True)
class CVSettings:
    n_splits: int
    shuffle: bool
    seed: int

    def __post_init__(self):
        if type(self.n_splits) is not int or self.n_splits < 2:
            raise ValueError("n_splits must be an integer >= 2")
        if type(self.shuffle) is not bool:
            raise ValueError("shuffle must be boolean")
        if type(self.seed) is not int or not 0 <= self.seed < 2**32:
            raise ValueError("random_seed must be a nonnegative 32-bit integer")

    @classmethod
    def from_config(cls, config: dict):
        validation = config["modeling"]["validation"]
        if validation["preferred"] != "stratified_cross_validation":
            raise ValueError("M5 requires the approved stratified_cross_validation protocol")
        return cls(validation["n_splits"], validation["shuffle"], config["project"]["random_seed"])


def make_stratified_folds(
    ids: pd.Series, y: pd.Series, n_splits: int = 5, seed: int = 42, shuffle: bool = True,
) -> pd.DataFrame:
    """Assign each canonical ID once; sklearn sees only dummy X and binary y."""
    settings = CVSettings(n_splits, shuffle, seed)
    validate_ids(ids)
    validate_binary(y, "CV target")
    if len(ids) != len(y):
        raise ValueError("CV IDs and target lengths differ")
    if set(y) != {0, 1} or y.value_counts().min() < n_splits:
        raise ValueError("Each class needs at least n_splits isolates")
    rows = pd.DataFrame({"isolate_id": ids.to_numpy(), "binary_label": y.to_numpy()})
    canonical = rows.sort_values("isolate_id", kind="stable").reset_index(drop=True)
    canonical["fold_id"] = -1
    splitter = StratifiedKFold(n_splits=settings.n_splits, shuffle=settings.shuffle,
                              random_state=settings.seed if settings.shuffle else None)
    for fold_id, (_, validation_indices) in enumerate(splitter.split([[0]] * len(canonical), canonical.binary_label)):
        canonical.loc[validation_indices, "fold_id"] = fold_id
    # Save in dataset row order, while making membership independent of caller row order.
    return rows[["isolate_id"]].merge(canonical, on="isolate_id", how="left", validate="one_to_one", sort=False)


def validate_folds(dataset: pd.DataFrame, folds: pd.DataFrame, n_splits: int) -> pd.DataFrame:
    if folds.columns.tolist() != ["isolate_id", "antibiotic", "binary_label", "fold_id"]:
        raise ValueError("Unexpected fold-file columns")
    validate_ids(folds.isolate_id)
    validate_binary(folds.binary_label, "Fold target")
    if not pd.api.types.is_integer_dtype(folds.fold_id) or set(folds.fold_id) != set(range(n_splits)):
        raise ValueError("Invalid or incomplete fold IDs")
    expected = dataset[["isolate_id", "antibiotic", "binary_label"]].set_index("isolate_id").sort_index()
    observed = folds[expected.columns.tolist() + ["isolate_id"]].set_index("isolate_id").sort_index()
    if not expected.equals(observed):
        raise ValueError("Fold IDs/antibiotic/targets disagree with assembled dataset")
    counts = pd.crosstab(folds.fold_id, folds.binary_label).reindex(index=range(n_splits), columns=[0, 1], fill_value=0)
    if (counts <= 0).any().any():
        raise ValueError("Every validation fold must contain both classes")
    sizes = counts.sum(axis=1)
    if sizes.max() - sizes.min() > 1 or ((counts.max() - counts.min()) > 1).any():
        raise ValueError("Fold sizes or class allocation are not stratified")
    return pd.DataFrame({
        "antibiotic": dataset.antibiotic.iloc[0], "fold_id": range(n_splits),
        "fold_size": sizes.to_numpy(), "susceptible_count": counts[0].to_numpy(),
        "resistant_count": counts[1].to_numpy(),
        "resistant_fraction": (counts[1] / sizes).to_numpy(),
    })

