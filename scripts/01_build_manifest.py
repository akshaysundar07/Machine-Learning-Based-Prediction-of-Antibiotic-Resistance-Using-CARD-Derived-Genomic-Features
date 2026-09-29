"""Register downloaded source files in data/dataset_manifest.csv.

This script intentionally does not scrape publication supplements automatically.
Module 1 requires verified source-file selection and provenance first.
"""
from pathlib import Path
import hashlib
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
manifest_path = ROOT / "data" / "dataset_manifest.csv"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

manifest = pd.read_csv(manifest_path)
print(manifest.to_string(index=False))
print("\nPlace verified source files under data/raw/* and update this manifest.")
