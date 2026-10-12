import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from kpamr.data.ast import standardize_phenotype


def test_standardize_phenotype():
    assert standardize_phenotype("Resistant") == "R"
    assert standardize_phenotype("S") == "S"
    assert standardize_phenotype("intermediate") == "I"
