"""Module 2 entry point.

Populate column mappings after the actual supplementary files are collected.
The script refuses to invent source schemas.
"""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
raw = ROOT / "data" / "raw"
files = [p for p in raw.rglob("*") if p.is_file() and p.name != ".gitkeep"]
if not files:
    raise SystemExit("No raw source files found. Complete Module 1 first.")
print("Raw files found:")
for p in files:
    print(" -", p.relative_to(ROOT))
print("\nNext: map actual source columns before running phenotype/feature processing.")
