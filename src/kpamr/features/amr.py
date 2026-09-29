"""Encode reviewed Table S1 genomic calls; phenotype columns are never features."""

from collections import defaultdict
from pathlib import Path
import hashlib
import json
import re

import pandas as pd
from openpyxl import load_workbook
from openpyxl.utils import get_column_letter

from kpamr.data.ast import SOURCE_PATH

# Reviewed row-2 schema, including the repeated published blaSHV-11 header.
# Footnote q/r column letters are stale; named headers define the actual blocks.
NONFEATURE_HEADERS = (
    "MRSN IDa", "Regionb", "Datec", "Year", "Sample Typed",
    "MLSTe", "K_locusf", "O_locusg", "O_typeh", "VIR Scorei",
    "Yersiniabactinj", "Colibactink", "Aerobactinl", "Salmochelinm", "RmpADCn",
    "rmpA2o", "AMR Phenotypep", "AMK", "GEN", "TOB",
    "SAM", "ATM", "CAZ", "CRO", "FEP",
    "CZA", "C/T", "ETP", "IPM", "MEM",
    "CIP", "LVX", "TZP", "SXT", "TET",
    "TGC", "CST", "MinION data", "All Resistance Genesq",
)
GENE_HEADERS = (
    "aac(2')-Iia", "aac(3)-IId", "aac(3)-IIe", "aac(3)-IIg", "aac(3)-IVa",
    "aac(6')-Ib4", "aac(6')-Ib-AKT", "aac(6')-Ib-cr5", "aac(6')-Ib-G", "aac(6')-IIc",
    "aac(6')-Il", "aadA1", "aadA16", "aadA2", "ant(2'')-Ia",
    "aph(3')-Ia", "aph(3'')-Ib", "aph(3')-VI", "aph(3')-Vib", "aph(3')-XV",
    "aph(4)-Ia", "aph(6)-Id", "aphA16", "armA", "arr-2",
    "arr-269927220", "arr-3", "arr-8-like", "blaCARB-2", "blaCMY-4",
    "blaCTX-M-14", "blaCTX-M-15", "blaCTX-M-3", "blaCTX-M-55", "blaCTX-M-8",
    "blaDHA-1", "blaFOX-5", "blaGES-5", "blaIMP-15", "blaKPC-2",
    "blaKPC-3", "blaNDM-1", "blaNDM-5", "blaNDM-7", "blaOXA-1",
    "blaOXA-10", "blaOXA-9", "blaOXA-181", "blaOXA-232", "blaOXA-48",
    "blaOXA-2", "blaSHV-1", "blaSHV-11", "blaSHV-110", "blaSHV-119",
    "blaSHV-142", "blaSHV-168", "blaSHV-187", "blaSHV-202", "blaSHV-207",
    "blaSHV-26", "blaSHV-27", "blaSHV-28", "blaSHV-32", "blaSHV-33",
    "blaSHV-171", "blaSHV-108", "blaSHV-36", "blaSHV-52", "blaSHV-60",
    "blaSHV-71", "blaSHV-75", "blaSHV-76", "blaSHV-77", "blaSHV-5",
    "blaSHV-7", "blaSHV-11", "blaSHV-12", "blaTEM-1", "blaVIM-1",
    "blaVIM-4", "catA1", "catA2", "catB", "catB11",
    "catB2", "catB3", "cmlA1", "cmlA5", "cmlA6",
    "dfr7", "dfrA1", "dfrA12", "dfrA14", "dfrA18",
    "dfrA19", "dfrA26", "dfrA27", "dfrA32", "dfrA5",
    "dfrA7", "dfrA8", "erm(42)", "erm(B)", "ermC",
    "floR", "floR2", "fosA_gen", "fosA5", "fosA7",
    "fosA7.3", "mcr-1.1", "mcr-8.1", "mph(A)", "mph(B)",
    "mph(E)", "msr(E)", "qnrA1", "qnrB1", "qnrB19",
    "qnrB4", "qnrB6", "qnrB9", "qnrS1", "rmtB1",
    "rmtF1", "rmtH", "sul1", "sul2", "sul3",
    "tet(A)", "tet(B)", "tet(D)", "tet(G)", "tet(M)",
    "tet(X)",
)
VARIANT_HEADERS = ("OmpK35/OmpK36 mutations", "MgrB mutations")
GENE_COLUMNS = tuple(get_column_letter(i) for i in range(40, 40 + len(GENE_HEADERS)))
VARIANT_COLUMNS = ("FT", "FU")
SOURCE_COLUMNS = (*GENE_COLUMNS, *VARIANT_COLUMNS)
AST_PATH = Path("data/interim/ast/ast_clean_long.csv")

# One reviewed cell names a distinct allele under the Ib-AKT header.
# Preserve the named allele, never infer Ib-AKT presence from a nonempty cell.
CELL_ALLELE_EXCEPTIONS = {"AT": ("aac(6')-Ib",)}
VARIANT_PATTERNS = {
    "FT": re.compile(r"(?:OmpK35-(?:[0-9]|[1-9][0-9])%|OmpK36(?:GD|TD))"),
    "FU": re.compile(r"MgrB-(?:[0-9]|[1-9][0-9])%"),
}


def validate_isolate_ids(ids: pd.Series, expected_count: int = 100) -> None:
    if (len(ids) != expected_count or ids.isna().any() or ids.duplicated().any()
            or not ids.map(lambda value: isinstance(value, str) and bool(re.fullmatch(r"[0-9]+", value))).all()):
        raise ValueError(f"Expected {expected_count} unique numeric-string isolate IDs")


def read_amr_table(path: str | Path) -> pd.DataFrame:
    """Validate the named published schema and project only IDs and genomic calls."""
    workbook = load_workbook(path, read_only=True, data_only=False)
    try:
        if "Table S1" not in workbook.sheetnames:
            raise ValueError("Missing worksheet 'Table S1'")
        rows = list(workbook["Table S1"].iter_rows(values_only=True))
    finally:
        workbook.close()
    expected = (*NONFEATURE_HEADERS, *GENE_HEADERS, *VARIANT_HEADERS)
    if len(rows) < 3 or tuple(rows[1]) != expected:
        raise ValueError("Table S1 schema mismatch: reviewed metadata/AST/gene/variant headers required")
    if rows[0][39] != "aac(2')" or rows[0][175] != "Mutationsr":
        raise ValueError("Table S1 AMR family/mutation header semantics changed")
    boundaries = [i for i, row in enumerate(rows[2:], start=2)
                  if isinstance(row[0], str) and row[0].startswith("Abbreviations -")]
    if len(boundaries) != 1:
        raise ValueError("Expected exactly one Table S1 footnote boundary")
    boundary = boundaries[0]
    footer = rows[boundary:]
    if any(any(value is not None for value in row[1:]) for row in footer):
        raise ValueError("Unexpected data after Table S1 footnote boundary")
    if any(isinstance(row[0], (int, float)) or
           (isinstance(row[0], str) and row[0].strip().isdigit()) for row in footer):
        raise ValueError("Unexpected isolate after Table S1 footnote boundary")
    notes = [row[0] for row in footer if isinstance(row[0], str)]
    gene_note = next((value for value in notes if value.startswith("qAll resistance genes")), "")
    variant_note = next((value for value in notes if value.startswith("rMutations")), "")
    if not all(value in gene_note for value in ("AMRFinderPlus", "ARIBA", "gene name in the row indicates that it is present")):
        raise ValueError("Missing/changed published gene presence semantics")
    if not all(value in variant_note for value in ("ompK35/ompK36", "mgrB", "Kleborate")):
        raise ValueError("Missing/changed published variant provenance")
    records = [[str(row[0]).strip() if row[0] is not None else "", *row[39:177]]
               for row in rows[2:boundary]]
    source = pd.DataFrame(records, columns=["isolate_id", *SOURCE_COLUMNS])
    validate_isolate_ids(source["isolate_id"])
    return source.sort_values("isolate_id", kind="stable").reset_index(drop=True)


def is_absent(value) -> bool:
    """An empty published call is absence of a reported call, not inferred wild type."""
    return bool(pd.isna(value)) or (isinstance(value, str) and not value.strip())


def gene_call(column: str, header: str, value) -> str | None:
    if is_absent(value):
        return None
    if isinstance(value, str):
        value = value.strip()
        if value.casefold() == header.casefold():
            return header
        if value in CELL_ALLELE_EXCEPTIONS.get(column, ()):
            return value
    raise ValueError(f"Unexpected determinant in {column} ({header}): {value!r}")


def variant_tokens(column: str, value) -> tuple[str, ...]:
    """Split the observed semicolon syntax; keep exact published variant tokens.

    Percentage tokens are source annotations, not inferred effect sizes. No
    explicit wild-type tokens occur in this workbook; unreviewed tokens fail.
    """
    if column not in VARIANT_PATTERNS:
        raise ValueError(f"Unreviewed variant column: {column}")
    if is_absent(value):
        return ()
    if not isinstance(value, str):
        raise ValueError(f"Unexpected variant value in {column}: {value!r}")
    tokens = tuple(token.strip() for token in value.split(";"))
    if any(not VARIANT_PATTERNS[column].fullmatch(token) for token in tokens):
        raise ValueError(f"Unexpected variant token in {column}: {value!r}")
    return tuple(sorted(set(tokens)))


def machine_feature_id(feature_type: str, display_name: str) -> str:
    """Readable slug plus exact-name digest avoids punctuation-normalization collisions."""
    prefix = {"gene_or_allele": "gene", "resistance_variant": "variant"}[feature_type]
    slug = re.sub(r"[^A-Za-z0-9]+", "_", display_name).strip("_")
    if not slug:
        raise ValueError("Empty biological feature name")
    digest = hashlib.sha256(display_name.encode("utf-8")).hexdigest()[:12]
    return f"{prefix}__{slug}__{digest}"


def validate_binary_matrix(matrix: pd.DataFrame) -> None:
    if not matrix.columns.is_unique or matrix.columns[0] != "isolate_id":
        raise ValueError("Invalid matrix column schema")
    validate_isolate_ids(matrix["isolate_id"])
    features = matrix.drop(columns="isolate_id")
    if features.isna().any().any() or not features.isin([0, 1]).all().all():
        raise ValueError("AMR features must be binary with no missing values")
    if not all(re.fullmatch(r"(?:gene|variant)__[A-Za-z0-9_]+", col) for col in features):
        raise ValueError("Unsafe or non-genomic matrix column")


def remove_invariant_features(df: pd.DataFrame, id_col: str = "isolate_id") -> tuple[pd.DataFrame, list[str]]:
    validate_binary_matrix(df)
    if id_col != "isolate_id":
        raise ValueError("The canonical isolate key is isolate_id")
    invariant = [column for column in df.columns if column != id_col and df[column].nunique() == 1]
    return df.drop(columns=invariant), invariant


def build_binary_feature_matrix(source: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return full matrix and dictionary from validated, projected genomic columns.

    Identically named genes across duplicate source columns are OR-combined.
    Biologically different genes are never merged based on matching vectors.
    """
    if tuple(source.columns) != ("isolate_id", *SOURCE_COLUMNS):
        raise ValueError("Only the reviewed genomic source columns may enter encoding")
    validate_isolate_ids(source["isolate_id"])
    source = source.sort_values("isolate_id", kind="stable").reset_index(drop=True)
    vectors = {}
    definitions = {}
    aliases = defaultdict(set)
    raw_values = defaultdict(set)

    def register(feature_type, name, column, header):
        key = (feature_type, name)
        if key not in vectors:
            vectors[key] = [0] * len(source)
            definitions[key] = []
        pair = (column, header)
        if pair not in definitions[key]:
            definitions[key].append(pair)
        return key

    for column, header in zip(GENE_COLUMNS, GENE_HEADERS):
        register("gene_or_allele", header, column, header)  # Retain all-zero candidates.
        for index, value in enumerate(source[column]):
            try:
                name = gene_call(column, header, value)
            except ValueError as error:
                raise ValueError(f"Isolate {source.at[index, 'isolate_id']}: {error}") from error
            if name is not None:
                key = register("gene_or_allele", name, column, header)
                vectors[key][index] = 1
                aliases[key].add(value)
                raw_values[key].add(value)

    for column, header in zip(VARIANT_COLUMNS, VARIANT_HEADERS):
        calls = []
        for index, value in enumerate(source[column]):
            try:
                calls.append(variant_tokens(column, value))
            except ValueError as error:
                raise ValueError(f"Isolate {source.at[index, 'isolate_id']}: {error}") from error
        for token in sorted({token for call in calls for token in call}):
            key = register("resistance_variant", token, column, header)
            for index, call in enumerate(calls):
                if token in call:
                    vectors[key][index] = 1
                    aliases[key].add(token)
                    raw_values[key].add(source.at[index, column])

    encoded = {"isolate_id": source["isolate_id"]}
    dictionary = []
    for (feature_type, name), vector in vectors.items():
        key = (feature_type, name)
        feature_id = machine_feature_id(feature_type, name)
        if feature_id in encoded:
            raise ValueError(f"Feature ID collision: {feature_id}")
        encoded[feature_id] = pd.Series(vector, dtype="int8")
        count = sum(vector)
        retained = 0 < count < len(source)
        is_gene = feature_type == "gene_or_allele"
        dictionary.append({
            "feature_id": feature_id, "display_name": name,
            "source_column": json.dumps([header for _, header in definitions[key]]),
            "source_excel_columns": json.dumps([column for column, _ in definitions[key]]),
            "feature_type": feature_type,
            "encoding": ("named_cell_presence_case_normalized; OR_identical_determinant_columns; blank=0"
                         if is_gene else "exact_semicolon_token_presence; blank=no_published_variant_call"),
            "published_names": json.dumps(sorted(aliases[key] | {name})),
            "published_cell_values": json.dumps(sorted(raw_values[key])),
            "published_callers": "AMRFinderPlus;ARIBA" if is_gene else "Kleborate",
            "prevalence_count": count, "prevalence_fraction": count / len(source),
            "retained_in_model_matrix": retained,
            "drop_reason": "" if retained else ("zero_variance_absent" if count == 0 else "zero_variance_present"),
        })
    full = pd.DataFrame(encoded)
    validate_binary_matrix(full)
    return full, pd.DataFrame(dictionary)


def validate_m3_isolate_set(matrix: pd.DataFrame, ast_path: str | Path) -> None:
    """Read only the canonical ID column from M3; never import AST labels."""
    validate_isolate_ids(matrix["isolate_id"])
    ast_ids = pd.read_csv(ast_path, usecols=["isolate_id"], dtype=str, keep_default_na=False)["isolate_id"]
    unique_ids = ast_ids.drop_duplicates().reset_index(drop=True)
    validate_isolate_ids(unique_ids)
    if set(matrix["isolate_id"]) != set(unique_ids):
        missing = sorted(set(unique_ids) - set(matrix["isolate_id"]))
        extra = sorted(set(matrix["isolate_id"]) - set(unique_ids))
        raise ValueError(f"M3 isolate set mismatch: missing={missing}; extra={extra}")


def summarize_features(source, full, model, dictionary) -> pd.DataFrame:
    validate_binary_matrix(full)
    validate_binary_matrix(model)
    if dictionary["feature_id"].tolist() != full.columns[1:].tolist():
        raise ValueError("Feature dictionary does not map every full-matrix column")
    counts = full.iloc[:, 1:].sum()
    if dictionary["prevalence_count"].tolist() != counts.tolist():
        raise ValueError("Feature dictionary prevalence mismatch")
    if dictionary.loc[dictionary.retained_in_model_matrix, "feature_id"].tolist() != model.columns[1:].tolist():
        raise ValueError("Feature dictionary retention mismatch")
    identical_groups = defaultdict(list)
    for feature in model.columns[1:]:
        identical_groups[tuple(model[feature])].append(feature)
    groups = [features for features in identical_groups.values() if len(features) > 1]
    case_cells = 0
    distinct_cells = 0
    for column, header in zip(GENE_COLUMNS, GENE_HEADERS):
        for value in source[column]:
            if not is_absent(value):
                case_cells += int(value.strip() != header and value.strip().casefold() == header.casefold())
                distinct_cells += int(value.strip().casefold() != header.casefold())
    return pd.DataFrame([{
        "total_isolates": len(full),
        "candidate_source_columns": len(SOURCE_COLUMNS),
        "gene_allele_source_columns": len(GENE_COLUMNS),
        "variant_source_columns": len(VARIANT_COLUMNS),
        "full_encoded_features": full.shape[1] - 1,
        "gene_allele_features": int(dictionary.feature_type.eq("gene_or_allele").sum()),
        "variant_features": int(dictionary.feature_type.eq("resistance_variant").sum()),
        "zero_variance_absent_features": int(counts.eq(0).sum()),
        "zero_variance_present_features": int(counts.eq(len(full)).sum()),
        "modeling_ready_features": model.shape[1] - 1,
        "missing_encoded_values": int(full.iloc[:, 1:].isna().sum().sum()),
        "duplicate_isolate_ids": int(full.isolate_id.duplicated().sum()),
        "duplicate_gene_header_columns_merged": len(GENE_HEADERS) - len(set(GENE_HEADERS)),
        "case_normalized_gene_cells": case_cells,
        "distinct_cell_alleles_under_other_header": distinct_cells,
        "identical_vector_groups_retained": len(groups),
        "identical_vector_feature_ids": json.dumps(groups),
        "m3_isolate_set_match": True,
        "source_sheet": "Table S1",
        "schema_note": "Row-2 named headers override stale column-letter references in footnotes q/r",
        "gene_policy": "Blank=absence of published call; named cells=presence; exact same determinant columns OR-combined",
        "variant_policy": "Exact published semicolon tokens; blanks=no published call, not inferred wild type; percentage tokens uninterpreted",
        "filter_policy": "Remove only zero variance; preserve rare features and distinct determinants with identical vectors",
    }])

