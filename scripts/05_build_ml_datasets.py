"""M5: join locked M3/M4 artifacts and save configured stratified CV assignments."""

import hashlib
import json
from pathlib import Path
import sys

import pandas as pd
import sklearn

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from kpamr.config import load_config  # noqa: E402
from kpamr.data.io import write_csv  # noqa: E402
from kpamr.modeling.split import (  # noqa: E402
    CONFIG_PATH, DICTIONARY_PATH, ELIGIBILITY_PATH, FEATURE_PATH, RANKING_PATH,
    CVSettings, antibiotic_stem, assemble_dataset, eligible_targets,
    make_stratified_folds, read_csv_checked, validate_folds, validate_locked_features,
)

SUMMARY_PATH = Path("data/reports/ml_dataset_summary.csv")
FOLD_SUMMARY_PATH = Path("data/reports/cv_fold_summary.csv")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_input_provenance(root, manifest, inputs):
    required = {"file_id", "local_path", "source_type", "source_name", "source_accession_or_url",
                "retrieved_date", "description", "sha256", "status"}
    if required - set(manifest.columns) or manifest.file_id.duplicated().any():
        raise ValueError("Invalid dataset manifest schema or duplicate IDs")
    hashes = {}
    for path in inputs:
        row = manifest.loc[manifest.local_path.eq(path.as_posix())]
        digest = sha256(root / path)
        if len(row) != 1 or row.sha256.item() != digest:
            raise ValueError(f"Locked input provenance mismatch: {path}")
        hashes[path] = digest
    return hashes


def prepare_manifest(manifest, upstream, hashes, settings):
    records = []
    for path, sources in upstream.items():
        file_id = "M5_" + path.stem.upper()
        by_id = manifest.loc[manifest.file_id.eq(file_id)]
        by_path = manifest.loc[manifest.local_path.eq(path.as_posix())]
        if (not by_id.empty and by_id.local_path.item() != path.as_posix()
                or not by_path.empty and (len(by_path) != 1 or by_path.file_id.item() != file_id)):
            raise ValueError(f"M5 manifest ID/path collision: {path}")
        records.append({
            "file_id": file_id, "local_path": path.as_posix(), "source_type": "derived_ml_dataset_or_cv",
            "source_name": "Locked M4 genomic features plus M3 binary labels",
            "source_accession_or_url": json.dumps([source.as_posix() for source in sources]),
            "retrieved_date": "",
            "description": "scripts/05_build_ml_datasets.py; input_sha256="
                + json.dumps({source.as_posix(): hashes[source] for source in sources}, sort_keys=True)
                + f"; stratified CV: n_splits={settings.n_splits}, shuffle={settings.shuffle}, seed={settings.seed}",
            "sha256": "", "status": "VALIDATED_M5",
        })
    retained = manifest.loc[~manifest.file_id.isin([row["file_id"] for row in records])]
    return pd.concat([retained, pd.DataFrame(records)], ignore_index=True).fillna("")


def process_ml_datasets(root: Path = ROOT):
    root = Path(root)
    config = load_config(root / CONFIG_PATH)
    settings = CVSettings.from_config(config)
    features = read_csv_checked(root / FEATURE_PATH)
    columns = validate_locked_features(features, read_csv_checked(root / DICTIONARY_PATH))
    selected = eligible_targets(read_csv_checked(root / ELIGIBILITY_PATH), read_csv_checked(root / RANKING_PATH))
    label_paths = {
        antibiotic: Path(f"data/processed/per_antibiotic/{antibiotic_stem(antibiotic)}_labels.csv")
        for antibiotic in selected.antibiotic
    }
    if len(set(label_paths.values())) != len(label_paths):
        raise ValueError("Antibiotic filenames collide")
    shared = [FEATURE_PATH, DICTIONARY_PATH, ELIGIBILITY_PATH, RANKING_PATH]
    manifest_path = root / "data/dataset_manifest.csv"
    manifest = read_csv_checked(manifest_path)
    hashes = verify_input_provenance(root, manifest, [*shared, *label_paths.values()])
    hashes[CONFIG_PATH] = sha256(root / CONFIG_PATH)
    outputs, upstream, summaries, fold_summaries = {}, {}, [], []
    feature_order_hash = hashlib.sha256(("\n".join(columns) + "\n").encode("utf-8")).hexdigest()
    for _, stats in selected.iterrows():
        antibiotic = stats.antibiotic
        label_path = label_paths[antibiotic]
        dataset = assemble_dataset(features, read_csv_checked(root / label_path), stats)
        folds = make_stratified_folds(
            dataset.isolate_id, dataset.binary_label, settings.n_splits, settings.seed, settings.shuffle,
        )
        folds.insert(1, "antibiotic", antibiotic)
        fold_summary = validate_folds(dataset, folds, settings.n_splits)
        stem = antibiotic_stem(antibiotic)
        dataset_path = Path(f"data/processed/ml_datasets/{stem}_ml.csv")
        fold_path = Path(f"data/processed/splits/{stem}_folds.csv")
        outputs[dataset_path], outputs[fold_path] = dataset, folds
        for path in (dataset_path, fold_path):
            upstream[path] = [*shared, label_path, CONFIG_PATH]
        fold_summaries.append(fold_summary)
        susceptible = int(dataset.binary_label.eq(0).sum())
        resistant = int(dataset.binary_label.eq(1).sum())
        summaries.append({
            "antibiotic": antibiotic, "filename": dataset_path.name, "total_samples": len(dataset),
            "susceptible_count": susceptible, "resistant_count": resistant, "feature_count": len(columns),
            "positive_fraction": resistant / len(dataset), "negative_fraction": susceptible / len(dataset),
            "n_splits": settings.n_splits,
            "minimum_fold_size": int(fold_summary.fold_size.min()),
            "maximum_fold_size": int(fold_summary.fold_size.max()),
            "minimum_resistant_per_fold": int(fold_summary.resistant_count.min()),
            "minimum_susceptible_per_fold": int(fold_summary.susceptible_count.min()),
            "class_ratio": max(susceptible, resistant) / min(susceptible, resistant),
            "source_feature_matrix": FEATURE_PATH.as_posix(), "source_label_file": label_path.as_posix(),
            "eligible_for_ml": True, "random_seed": settings.seed, "shuffle": settings.shuffle,
            "validation_strategy": "stratified_cross_validation",
            "class_ratio_definition": "majority_count/minority_count",
            "feature_order_sha256": feature_order_hash, "sklearn_version": sklearn.__version__,
            "fold_interpretation": "fold_id == k: validation; all other folds: training",
        })
    outputs[SUMMARY_PATH] = pd.DataFrame(summaries)
    outputs[FOLD_SUMMARY_PATH] = pd.concat(fold_summaries, ignore_index=True)
    for path in (SUMMARY_PATH, FOLD_SUMMARY_PATH):
        upstream[path] = [*shared, *label_paths.values(), CONFIG_PATH]
    # An eligibility change must not leave stale targets mixed into the output set.
    for directory, pattern in (("data/processed/ml_datasets", "*_ml.csv"), ("data/processed/splits", "*_folds.csv")):
        stale = [path.relative_to(root).as_posix() for path in (root / directory).glob(pattern)
                 if path.relative_to(root) not in outputs]
        if stale:
            raise ValueError(f"Unexpected/stale M5 outputs; review before regeneration: {stale}")
    manifest = prepare_manifest(manifest, upstream, hashes, settings)
    if any(sha256(root / path) != digest for path, digest in hashes.items()):
        raise ValueError("A locked input changed during M5 assembly")
    # All targets and fold checks finish before any persistent output is replaced.
    for path, frame in outputs.items():
        write_csv(frame, root / path)
        manifest.loc[manifest.file_id.eq("M5_" + path.stem.upper()), "sha256"] = sha256(root / path)
    write_csv(manifest, manifest_path)
    print(f"M5 PASSED: {len(selected)} eligible antibiotics; {len(columns)} unchanged genomic features each; "
          f"{settings.n_splits} folds; {len(outputs[FOLD_SUMMARY_PATH])} fold-summary rows.")
    print("Every validation fold contains both classes; source inputs unchanged; no models trained.")
    return outputs


if __name__ == "__main__":
    try:
        process_ml_datasets()
    except (ValueError, KeyError, FileNotFoundError) as error:
        raise SystemExit(f"M5 assembly FAILED: {error}") from error

