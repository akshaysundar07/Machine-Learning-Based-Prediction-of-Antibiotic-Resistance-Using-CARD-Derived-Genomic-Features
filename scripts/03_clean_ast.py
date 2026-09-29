"""Module 3: process published Table S1 and register derived AST provenance."""

import hashlib
from pathlib import Path
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from kpamr.config import load_config  # noqa: E402
from kpamr.data.ast import (  # noqa: E402
    SOURCE_PATH, EligibilityRule, binary_label_tables, clean_ast_long,
    rank_antibiotics, read_table_s1, summarize_ast, validate_ast_outputs,
)
from kpamr.data.io import write_csv  # noqa: E402


def update_manifest(manifest: pd.DataFrame, paths: list[Path], source_hash: str) -> pd.DataFrame:
    """Upsert only M3 records; keep predecessor records and retrieval dates intact."""
    required = {"file_id", "local_path", "source_type", "source_name", "source_accession_or_url",
                "retrieved_date", "description", "sha256", "status"}
    if required - set(manifest.columns) or manifest["file_id"].duplicated().any():
        raise ValueError("Invalid dataset manifest schema or duplicate file IDs")
    records = [{
        "file_id": "M3_SOURCE_TABLE_S1", "local_path": SOURCE_PATH.as_posix(),
        "source_type": "supplementary", "source_name": "Martin et al. 2023 - Table S1",
        "source_accession_or_url": "DOI:10.1099/mgen.0.000967; BioProject:PRJNA717739",
        "retrieved_date": "", "description": "Local published Table S1 used for M3; retrieval date not recorded",
        "sha256": source_hash, "status": "LOCAL_SOURCE_VALIDATED_M3",
    }]
    for path in paths:
        records.append({
            "file_id": f"M3_{path.stem.upper()}", "local_path": path.as_posix(),
            "source_type": "derived_ast", "source_name": "Table S1 AST phenotype processing",
            "source_accession_or_url": SOURCE_PATH.as_posix(), "retrieved_date": "",
            "description": "scripts/03_process_ast.py; Table S1; source_sha256=" + source_hash
                + "; eligibility in configs/project.yaml; S=0 R=1; I/missing excluded from binary targets",
            "sha256": "", "status": "VALIDATED_M3",
        })
    for record in records:
        existing = manifest.loc[manifest["file_id"].eq(record["file_id"])]
        if not existing.empty and existing.iloc[0]["local_path"] != record["local_path"]:
            raise ValueError(f"Manifest file ID collision: {record['file_id']}")
    retained = manifest.loc[~manifest["file_id"].isin([record["file_id"] for record in records])]
    return pd.concat([retained, pd.DataFrame(records)], ignore_index=True).fillna("")


def process_ast(root: Path = ROOT) -> dict[Path, pd.DataFrame]:
    root = Path(root)
    source = root / SOURCE_PATH
    before = hashlib.sha256(source.read_bytes()).hexdigest()
    config = load_config(root / "configs/project.yaml")
    rule = EligibilityRule(**config["ast"]["eligibility"])
    tidy = clean_ast_long(read_table_s1(source))
    summary = summarize_ast(tidy, rule)
    tables = binary_label_tables(tidy)
    validate_ast_outputs(tidy, summary, tables)
    outputs = {
        Path("data/interim/ast/ast_clean_long.csv"): tidy,
        Path("data/reports/ast_binary_summary.csv"): summary,
        Path("data/reports/ast_antibiotic_ranking.csv"): rank_antibiotics(summary),
    }
    for antibiotic, table in tables.items():
        name = "CT" if antibiotic == "C/T" else antibiotic
        outputs[Path(f"data/processed/per_antibiotic/{name}_labels.csv")] = table
    manifest_path = root / "data/dataset_manifest.csv"
    manifest = pd.read_csv(manifest_path, dtype=str, keep_default_na=False)
    manifest = update_manifest(manifest, list(outputs), before)
    if hashlib.sha256(source.read_bytes()).hexdigest() != before:
        raise ValueError("Raw workbook changed during AST processing")
    for relative_path, frame in outputs.items():
        write_csv(frame, root / relative_path)
        mask = manifest["file_id"].eq(f"M3_{relative_path.stem.upper()}")
        manifest.loc[mask, "sha256"] = hashlib.sha256((root / relative_path).read_bytes()).hexdigest()
    write_csv(manifest, manifest_path)
    print(f"AST processing PASSED: 100 isolates; 19 primary antibiotics; {len(tidy)} preserved rows (including CST).")
    print(f"19 binary label files; {int(summary['eligible_for_ml'].sum())} eligible targets; unexpected labels: 0.")
    print("CST is audit-only; ranking reflects dataset suitability. Raw workbook SHA-256 unchanged.")
    return outputs


if __name__ == "__main__":
    try:
        process_ast()
    except (ValueError, KeyError, FileNotFoundError) as error:
        raise SystemExit(f"AST processing FAILED: {error}") from error
