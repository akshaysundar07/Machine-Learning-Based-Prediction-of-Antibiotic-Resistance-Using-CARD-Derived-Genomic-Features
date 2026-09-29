"""M4 entry point: published genomic AMR features without AST targets or metadata."""

import hashlib
from pathlib import Path
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from kpamr.config import load_config  # noqa: E402
from kpamr.data.io import write_csv  # noqa: E402
from kpamr.features.amr import (  # noqa: E402
    AST_PATH, SOURCE_PATH, build_binary_feature_matrix, read_amr_table,
    remove_invariant_features, summarize_features, validate_m3_isolate_set,
)

OUTPUT_PATHS = (
    Path("data/interim/amr/amr_feature_matrix_full.csv"),
    Path("data/processed/master/amr_feature_matrix.csv"),
    Path("data/reports/amr_feature_dictionary.csv"),
    Path("data/reports/amr_feature_summary.csv"),
)


def prepare_manifest(manifest, paths, source_hash, isolate_hash):
    required = {"file_id", "local_path", "source_type", "source_name", "source_accession_or_url",
                "retrieved_date", "description", "sha256", "status"}
    if required - set(manifest.columns) or manifest.file_id.duplicated().any():
        raise ValueError("Invalid dataset manifest schema or duplicate IDs")
    source = manifest.loc[manifest.local_path.eq(SOURCE_PATH.as_posix())]
    if len(source) != 1 or source.sha256.item().lower() != source_hash:
        raise ValueError("Table S1 source must have one matching SHA-256 provenance record")
    records = []
    for path in paths:
        file_id = "M4_" + path.stem.upper()
        existing = manifest.loc[manifest.file_id.eq(file_id)]
        if not existing.empty and existing.local_path.item() != path.as_posix():
            raise ValueError(f"Manifest ID collision: {file_id}")
        same_path = manifest.loc[manifest.local_path.eq(path.as_posix())]
        if not same_path.empty and (len(same_path) != 1 or same_path.file_id.item() != file_id):
            raise ValueError(f"Manifest path collision: {path}")
        records.append({
            "file_id": file_id, "local_path": path.as_posix(), "source_type": "derived_amr",
            "source_name": "Table S1 published genomic AMR determinants and variants",
            "source_accession_or_url": SOURCE_PATH.as_posix(), "retrieved_date": "",
            "description": "scripts/04_build_amr_features.py; source_sha256=" + source_hash
                + "; M3 isolate-ID-only validation: " + AST_PATH.as_posix()
                + "; isolate_set_sha256=" + isolate_hash
                + "; genes: AMRFinderPlus/ARIBA; variants: Kleborate; no caller reruns; no AST features",
            "sha256": "", "status": "VALIDATED_M4",
        })
    retained = manifest.loc[~manifest.file_id.isin([row["file_id"] for row in records])]
    return pd.concat([retained, pd.DataFrame(records)], ignore_index=True).fillna("")


def process_amr(root: Path = ROOT):
    root = Path(root)
    config = load_config(root / "configs/project.yaml")
    required = {
        "primary_source": "published_AMRFinderPlus_and_ARIBA_determinants",
        "binary_presence_absence": True,
        "include_resistance_variants": True,
        "rerun_amr_callers": False,
    }
    if any(config["features"].get(key) != value for key, value in required.items()):
        raise ValueError("M4 requires the configured published binary gene/variant workflow")
    source_path, ast_path = root / SOURCE_PATH, root / AST_PATH
    hashes = {path: hashlib.sha256(path.read_bytes()).hexdigest() for path in (source_path, ast_path)}
    source = read_amr_table(source_path)
    validate_m3_isolate_set(source, ast_path)
    full, dictionary = build_binary_feature_matrix(source)
    model, removed = remove_invariant_features(full)
    summary = summarize_features(source, full, model, dictionary)
    source_hash = hashes[source_path]
    isolate_hash = hashlib.sha256(("\n".join(sorted(source.isolate_id)) + "\n").encode("utf-8")).hexdigest()
    summary["source_sha256"] = source_hash
    summary["isolate_set_sha256"] = isolate_hash
    outputs = dict(zip(OUTPUT_PATHS, (full, model, dictionary, summary)))
    manifest_path = root / "data/dataset_manifest.csv"
    manifest = prepare_manifest(
        pd.read_csv(manifest_path, dtype=str, keep_default_na=False),
        list(outputs), source_hash, isolate_hash,
    )
    if any(hashlib.sha256(path.read_bytes()).hexdigest() != digest for path, digest in hashes.items()):
        raise ValueError("Protected source or M3 AST file changed during M4 processing")
    for path, frame in outputs.items():
        write_csv(frame, root / path)
        mask = manifest.file_id.eq("M4_" + path.stem.upper())
        manifest.loc[mask, "sha256"] = hashlib.sha256((root / path).read_bytes()).hexdigest()
    write_csv(manifest, manifest_path)
    print(f"M4 PASSED: {len(full)} isolates matching M3; {full.shape[1] - 1} full features; "
          f"{len(removed)} zero-variance features removed; {model.shape[1] - 1} modeling-ready features.")
    print("Binary genomic features only; raw workbook and M3 data unchanged; provenance registered.")
    return outputs


if __name__ == "__main__":
    try:
        process_amr()
    except (ValueError, KeyError, FileNotFoundError) as error:
        raise SystemExit(f"AMR processing FAILED: {error}") from error

