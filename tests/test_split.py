"""M5 regression checks for locked features, exact joins, and deterministic CV."""

import hashlib
from pathlib import Path
import runpy
import sys

import pandas as pd
import pytest
from sklearn.model_selection import StratifiedKFold

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from kpamr.modeling.split import (  # noqa: E402
    CONFIG_PATH, DICTIONARY_PATH, ELIGIBILITY_PATH, FEATURE_PATH, RANKING_PATH,
    CVSettings, antibiotic_stem, assemble_dataset, eligible_targets, make_stratified_folds,
    read_csv_checked, validate_folds, validate_locked_features,
)


@pytest.fixture
def inputs():
    columns = [f"gene__determinant_{i:03d}" for i in range(144, -1, -1)]
    features = pd.DataFrame({
        "isolate_id": [str(1000 + i) for i in range(100)],
        **{name: [int((row + i) % (i % 11 + 2) == 0) for row in range(100)]
           for i, name in enumerate(columns)},
    })
    features[columns[0]] = [0] * 99 + [1]  # Constant after the target subset; must remain.
    labels = pd.DataFrame({
        "isolate_id": features.isolate_id.iloc[:80],
        "antibiotic": "CIP", "phenotype": ["S"] * 40 + ["R"] * 40,
        "binary_label": [0] * 40 + [1] * 40,
    })
    stats = pd.Series({"antibiotic": "CIP", "eligible_for_ml": True, "primary_antibiotic": True,
                       "usable_binary_count": 80, "susceptible_count": 40, "resistant_count": 40})
    dictionary = pd.DataFrame({"feature_id": columns, "feature_type": "gene_or_allele",
                               "retained_in_model_matrix": True})
    return features, labels, stats, dictionary


def make_reports():
    antibiotics = ["AMK", "GEN", "TOB", "SAM", "ATM", "CAZ", "CRO", "FEP", "CZA", "C/T",
                   "ETP", "IPM", "MEM", "CIP", "LVX", "TZP", "SXT", "TET", "TGC", "CST"]
    rows = []
    for antibiotic in antibiotics:
        rows.append({
            "antibiotic": antibiotic, "primary_antibiotic": antibiotic != "CST",
            "total_isolates": 100, "susceptible_count": 70, "intermediate_count": 0,
            "resistant_count": 30, "missing_count": 0, "usable_binary_count": 100,
            "eligible_for_ml": antibiotic not in {"AMK", "CZA", "TGC", "CST"},
        })
    summary = pd.DataFrame(rows)
    ranking = summary.sort_values(["eligible_for_ml", "antibiotic"], ascending=[False, True]).reset_index(drop=True)
    ranking.insert(0, "initial_experiment_rank", [
        i + 1 if eligible else "" for i, eligible in enumerate(ranking.eligible_for_ml)
    ])
    return summary, ranking


def test_exact_join_preserves_145_features_order_and_inputs(inputs):
    features, labels, stats, dictionary = inputs
    before_features, before_labels = features.copy(deep=True), labels.copy(deep=True)
    columns = validate_locked_features(features, dictionary)
    dataset = assemble_dataset(features, labels.sample(frac=1, random_state=3), stats)
    assert dataset.columns.tolist() == ["isolate_id", "antibiotic", *columns, "binary_label"]
    assert len(columns) == 145 and dataset.shape == (80, 148)
    assert dataset.isolate_id.tolist() == features.isolate_id.iloc[:80].tolist()
    assert dataset[columns[0]].eq(0).all()  # No per-target constant/rare-feature removal.
    pd.testing.assert_frame_equal(dataset[columns], features.iloc[:80][columns])
    assert dataset.binary_label.tolist() == labels.binary_label.tolist()
    assert not dataset.isna().any().any()
    assert dataset[columns].isin([0, 1]).all().all()
    assert "phenotype" not in dataset
    pd.testing.assert_frame_equal(features, before_features)
    pd.testing.assert_frame_equal(labels, before_labels)


def test_target_changes_do_not_transform_or_select_features(inputs):
    features, labels, stats, _ = inputs
    original = assemble_dataset(features, labels, stats)
    labels["binary_label"] = 1 - labels.binary_label
    labels["phenotype"] = labels.binary_label.map({0: "S", 1: "R"})
    changed = assemble_dataset(features, labels, stats)
    pd.testing.assert_frame_equal(original.drop(columns="binary_label"), changed.drop(columns="binary_label"))
    assert not original.binary_label.equals(changed.binary_label)


@pytest.mark.parametrize("mutation", [
    "duplicate_feature_ids", "duplicate_label_ids", "missing_id", "missing_label", "missing_label_column",
    "nonbinary_label", "unknown_label_isolate", "wrong_count", "wrong_class_counts",
    "intermediate", "label_phenotype_mismatch", "wrong_antibiotic", "phenotype_metadata",
    "missing_feature", "extra_feature", "feature_nan", "nonbinary_feature", "metadata_feature",
])
def test_invalid_joins_and_input_schemas_fail_loudly(inputs, mutation):
    features, labels, stats, _ = inputs
    if mutation == "duplicate_feature_ids":
        features.loc[1, "isolate_id"] = features.loc[0, "isolate_id"]
    elif mutation == "duplicate_label_ids":
        labels.loc[1, "isolate_id"] = labels.loc[0, "isolate_id"]
    elif mutation == "missing_id":
        labels.loc[0, "isolate_id"] = ""
    elif mutation == "missing_label":
        labels.loc[0, "binary_label"] = None
    elif mutation == "missing_label_column":
        labels = labels.drop(columns="binary_label")
    elif mutation == "nonbinary_label":
        labels.loc[0, "binary_label"] = 2
    elif mutation == "unknown_label_isolate":
        labels.loc[0, "isolate_id"] = "absent-from-M4"
    elif mutation == "wrong_count":
        labels = labels.iloc[:-1]
    elif mutation == "wrong_class_counts":
        stats["susceptible_count"] = 39
    elif mutation == "intermediate":
        labels.loc[0, "phenotype"] = "I"
    elif mutation == "label_phenotype_mismatch":
        labels.loc[0, "binary_label"] = 1
    elif mutation == "wrong_antibiotic":
        labels.loc[0, "antibiotic"] = "TOB"
    elif mutation == "phenotype_metadata":
        labels["MDR"] = "yes"
    elif mutation == "missing_feature":
        features = features.drop(columns=features.columns[1])
    elif mutation == "extra_feature":
        features["gene__extra"] = 0
    elif mutation == "feature_nan":
        features.loc[0, features.columns[1]] = None
    elif mutation == "nonbinary_feature":
        features.loc[0, features.columns[1]] = 2
    elif mutation == "metadata_feature":
        features = features.rename(columns={features.columns[1]: "country"})
    with pytest.raises(ValueError):
        assemble_dataset(features, labels, stats)


@pytest.mark.parametrize("mutation", ["reorder", "rename", "duplicate_dictionary", "invalid_type", "string_flag"])
def test_locked_dictionary_rejects_feature_schema_changes(inputs, mutation):
    features, _, _, dictionary = inputs
    if mutation == "reorder":
        features = features[["isolate_id", *features.columns[1:][::-1]]]
    elif mutation == "rename":
        features = features.rename(columns={features.columns[1]: "gene__unreviewed"})
    elif mutation == "duplicate_dictionary":
        dictionary.loc[1, "feature_id"] = dictionary.loc[0, "feature_id"]
    elif mutation == "invalid_type":
        dictionary.loc[0, "feature_type"] = "phenotype"
    else:
        dictionary["retained_in_model_matrix"] = "False"
    with pytest.raises(ValueError):
        validate_locked_features(features, dictionary)


def test_duplicate_csv_headers_rejected_before_pandas_renaming(tmp_path):
    path = tmp_path / "duplicate.csv"
    path.write_text("isolate_id,gene__x,gene__x\n1000,0,1\n", encoding="utf-8")
    with pytest.raises(ValueError, match="duplicate CSV header"):
        read_csv_checked(path)


def test_m3_flags_drive_eligibility_without_new_thresholds():
    summary, ranking = make_reports()
    selected = eligible_targets(summary, ranking)
    assert len(selected) == 16
    assert not set(selected.antibiotic) & {"AMK", "CZA", "TGC", "CST"}
    # Ineligible rows are balanced too; do not decide eligibility from class counts.
    for report in (summary, ranking):
        mask = report.antibiotic.eq("CIP")
        report.loc[mask, ["susceptible_count", "resistant_count"]] = [90, 10]
    assert "CIP" in eligible_targets(summary, ranking).antibiotic.tolist()


@pytest.mark.parametrize("mutation", ["cst_eligible", "disagreement", "duplicate", "string_flag", "bad_rank", "bad_count"])
def test_invalid_eligibility_reports_fail_loudly(mutation):
    summary, ranking = make_reports()
    if mutation == "cst_eligible":
        for report in (summary, ranking):
            report.loc[report.antibiotic.eq("CST"), "eligible_for_ml"] = True
    elif mutation == "disagreement":
        ranking.loc[ranking.antibiotic.eq("CIP"), "eligible_for_ml"] = False
    elif mutation == "duplicate":
        summary.loc[1, "antibiotic"] = summary.loc[0, "antibiotic"]
    elif mutation == "string_flag":
        summary["eligible_for_ml"] = "False"
    elif mutation == "bad_rank":
        ranking.loc[0, "initial_experiment_rank"] = 2
    elif mutation == "bad_count":
        summary.loc[0, "usable_binary_count"] -= 1
    with pytest.raises(ValueError):
        eligible_targets(summary, ranking)


def test_ineligible_and_cst_cannot_be_assembled(inputs):
    features, labels, stats, _ = inputs
    stats["eligible_for_ml"] = False
    with pytest.raises(ValueError, match="Only M3-eligible"):
        assemble_dataset(features, labels, stats)
    with pytest.raises(ValueError, match="Unsupported primary antibiotic"):
        antibiotic_stem("CST")


def test_ct_filename_keeps_actual_antibiotic(inputs):
    features, labels, stats, _ = inputs
    stats["antibiotic"], labels["antibiotic"] = "C/T", "C/T"
    result = assemble_dataset(features, labels, stats)
    assert antibiotic_stem("C/T") == "CT"
    assert result.antibiotic.eq("C/T").all()


@pytest.mark.parametrize("n_splits,shuffle,seed", [(5, True, 42), (4, True, 73), (5, False, 42)])
def test_cv_matches_configured_sklearn_and_class_totals(inputs, n_splits, shuffle, seed):
    features, labels, stats, _ = inputs
    dataset = assemble_dataset(features, labels, stats)
    folds = make_stratified_folds(dataset.isolate_id, dataset.binary_label, n_splits, seed, shuffle)
    folds.insert(1, "antibiotic", "CIP")
    summary = validate_folds(dataset, folds, n_splits)
    assert set(folds.fold_id) == set(range(n_splits))
    assert folds.isolate_id.is_unique and set(folds.isolate_id) == set(dataset.isolate_id)
    assert len(folds) == len(dataset) and summary.fold_size.sum() == len(dataset)
    assert summary.susceptible_count.sum() == summary.resistant_count.sum() == 40
    assert summary[["susceptible_count", "resistant_count"]].gt(0).all().all()
    splitter = StratifiedKFold(n_splits=n_splits, shuffle=shuffle, random_state=seed if shuffle else None)
    expected = pd.Series(-1, index=dataset.index)
    for fold_id, (train, validation) in enumerate(splitter.split([[0]] * len(dataset), dataset.binary_label)):
        assert not set(train) & set(validation)
        assert set(train) | set(validation) == set(dataset.index)
        expected.iloc[validation] = fold_id
    assert folds.fold_id.tolist() == expected.tolist()
    second = make_stratified_folds(dataset.isolate_id, dataset.binary_label, n_splits, seed, shuffle)
    assert folds.fold_id.tolist() == second.fold_id.tolist()


def test_cv_is_order_stable_and_shuffle_false_ignores_seed(inputs):
    _, labels, _, _ = inputs
    first = make_stratified_folds(labels.isolate_id, labels.binary_label)
    shuffled = labels.sample(frac=1, random_state=17)
    repeated = make_stratified_folds(shuffled.isolate_id, shuffled.binary_label)
    pd.testing.assert_frame_equal(first.set_index("isolate_id").sort_index(), repeated.set_index("isolate_id").sort_index())
    no_shuffle = make_stratified_folds(labels.isolate_id, labels.binary_label, seed=42, shuffle=False)
    different_seed = make_stratified_folds(labels.isolate_id, labels.binary_label, seed=73, shuffle=False)
    pd.testing.assert_frame_equal(no_shuffle, different_seed)


@pytest.mark.parametrize("n_splits,shuffle,seed", [(1, True, 42), (5.0, True, 42), (5, "true", 42), (5, True, -1)])
def test_invalid_cv_settings_rejected(n_splits, shuffle, seed):
    with pytest.raises(ValueError):
        CVSettings(n_splits, shuffle, seed)


def test_config_protocol_is_used_without_defaults():
    config = {"project": {"random_seed": 73},
              "modeling": {"validation": {"preferred": "stratified_cross_validation", "n_splits": 4, "shuffle": False}}}
    assert CVSettings.from_config(config) == CVSettings(4, False, 73)
    config["modeling"]["validation"]["preferred"] = "phylogenetic"
    with pytest.raises(ValueError, match="approved"):
        CVSettings.from_config(config)


@pytest.mark.parametrize("targets", [[0] * 10, [0] * 10 + [1] * 4, [0] * 10 + [2] * 10])
def test_cv_rejects_absent_insufficient_or_invalid_classes(targets):
    with pytest.raises(ValueError):
        make_stratified_folds(pd.Series([str(i) for i in range(len(targets))]), pd.Series(targets))


@pytest.mark.parametrize("mutation", ["duplicate", "missing", "wrong_target", "wrong_antibiotic", "bad_fold", "single_class", "imbalanced"])
def test_fold_validator_rejects_incomplete_or_leaking_assignments(inputs, mutation):
    features, labels, stats, _ = inputs
    dataset = assemble_dataset(features, labels, stats)
    folds = make_stratified_folds(dataset.isolate_id, dataset.binary_label)
    folds.insert(1, "antibiotic", "CIP")
    if mutation == "duplicate":
        folds.loc[1, "isolate_id"] = folds.loc[0, "isolate_id"]
    elif mutation == "missing":
        folds = folds.iloc[:-1]
    elif mutation == "wrong_target":
        folds.loc[0, "binary_label"] = 1
    elif mutation == "wrong_antibiotic":
        folds.loc[0, "antibiotic"] = "TOB"
    elif mutation == "bad_fold":
        folds.loc[0, "fold_id"] = 7
    elif mutation == "single_class":
        folds.loc[folds.binary_label.eq(1) & folds.fold_id.eq(0), "fold_id"] = 1
    elif mutation == "imbalanced":
        moved = folds.index[folds.fold_id.eq(0) & folds.binary_label.eq(0)][:2]
        folds.loc[moved, "fold_id"] = 1
    with pytest.raises(ValueError):
        validate_folds(dataset, folds, 5)


@pytest.fixture
def pipeline_workspace(tmp_path):
    # Integration uses only the supplied locked CSVs and configuration, never raw data.
    paths = [FEATURE_PATH, DICTIONARY_PATH, ELIGIBILITY_PATH, RANKING_PATH, CONFIG_PATH]
    paths.extend(path.relative_to(ROOT) for path in sorted((ROOT / "data/processed/per_antibiotic").glob("*_labels.csv")))
    source_rows = []
    columns = ["file_id", "local_path", "source_type", "source_name", "source_accession_or_url",
               "retrieved_date", "description", "sha256", "status"]
    for index, path in enumerate(paths):
        destination = tmp_path / path
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes((ROOT / path).read_bytes())
        if path != CONFIG_PATH:
            source_rows.append([f"LOCKED_{index}", path.as_posix(), "locked", "fixture source", "upstream",
                                "2026-01-01", "preserve existing record",
                                hashlib.sha256(destination.read_bytes()).hexdigest(), "VALIDATED"])
    manifest_path = tmp_path / "data/dataset_manifest.csv"
    pd.DataFrame(source_rows, columns=columns).to_csv(manifest_path, index=False)
    process = runpy.run_path(str(ROOT / "scripts/05_build_ml_datasets.py"))["process_ml_datasets"]
    return tmp_path, paths, manifest_path, process


def test_pipeline_exclusions_provenance_protected_inputs_and_deterministic_rerun(pipeline_workspace):
    root, inputs, manifest_path, process = pipeline_workspace
    protected = {path: (root / path).read_bytes() for path in inputs}
    original_manifest = pd.read_csv(manifest_path, keep_default_na=False)
    outputs = process(root)
    summary = pd.read_csv(root / "data/reports/ml_dataset_summary.csv")
    folds = pd.read_csv(root / "data/reports/cv_fold_summary.csv")
    eligible = pd.read_csv(root / ELIGIBILITY_PATH)
    expected = set(eligible.loc[eligible.eligible_for_ml, "antibiotic"])
    assert set(summary.antibiotic) == expected
    assert len(outputs) == 2 * len(expected) + 2
    assert len(folds) == len(expected) * 5
    assert not set(summary.antibiotic) & {"AMK", "CZA", "TGC", "CST"}
    assert summary.feature_count.eq(145).all() and summary.random_seed.eq(42).all()
    assert folds[["susceptible_count", "resistant_count"]].gt(0).all().all()
    assert {path.name for path in (root / "data/processed/ml_datasets").glob("*.csv")} == {antibiotic_stem(name) + "_ml.csv" for name in expected}
    ct = pd.read_csv(root / "data/processed/ml_datasets/CT_ml.csv")
    assert ct.antibiotic.eq("C/T").all()
    snapshot = {path: (root / path).read_bytes() for path in outputs}
    snapshot[Path("data/dataset_manifest.csv")] = manifest_path.read_bytes()
    process(root)
    assert all((root / path).read_bytes() == before for path, before in snapshot.items())
    assert all((root / path).read_bytes() == before for path, before in protected.items())
    updated = pd.read_csv(manifest_path, keep_default_na=False)
    assert updated.file_id.is_unique
    pd.testing.assert_frame_equal(updated.iloc[:len(original_manifest)].reset_index(drop=True), original_manifest)
    for row in updated.iloc[len(original_manifest):].itertuples():
        assert row.sha256 == hashlib.sha256((root / row.local_path).read_bytes()).hexdigest()
        assert FEATURE_PATH.as_posix() in row.source_accession_or_url
        assert "_labels.csv" in row.source_accession_or_url
        assert "scripts/05_build_ml_datasets.py" in row.description


def test_locked_source_change_fails_before_any_output_write(pipeline_workspace):
    root, _, manifest_path, process = pipeline_workspace
    original_manifest = manifest_path.read_bytes()
    label_path = root / "data/processed/per_antibiotic/CIP_labels.csv"
    with label_path.open("a", encoding="utf-8") as stream:
        stream.write("\n")
    with pytest.raises(ValueError, match="Locked input provenance mismatch"):
        process(root)
    assert not (root / "data/processed/ml_datasets").exists()
    assert manifest_path.read_bytes() == original_manifest


def test_stale_ineligible_output_is_not_silently_left_in_output_set(pipeline_workspace):
    root, _, manifest_path, process = pipeline_workspace
    directory = root / "data/processed/ml_datasets"
    directory.mkdir()
    stale = directory / "CST_ml.csv"
    stale.write_text("unrelated old file", encoding="utf-8")
    before = manifest_path.read_bytes()
    with pytest.raises(ValueError, match="Unexpected/stale M5 outputs"):
        process(root)
    assert stale.read_text(encoding="utf-8") == "unrelated old file"
    assert manifest_path.read_bytes() == before
    assert len(list(directory.iterdir())) == 1

