from pathlib import Path
import sys
import yaml

ROOT = Path(__file__).resolve().parents[1]
required = [
    ROOT / "configs/project.yaml",
    ROOT / "data/dataset_manifest.csv",
    ROOT / "docs/MODULE_PLAN.md",
    ROOT / "src/kpamr",
]
missing = [str(p) for p in required if not p.exists()]
if missing:
    print("Workspace check FAILED. Missing:")
    print("\n".join(missing))
    sys.exit(1)

with (ROOT / "configs/project.yaml").open("r", encoding="utf-8") as f:
    cfg = yaml.safe_load(f)
assert cfg["project"]["bioproject"] == "PRJNA717739"
assert cfg["labels"]["intermediate_policy"] == "exclude_from_primary_binary_model"
print("Workspace check PASSED")
print(f"Project: {cfg['project']['name']}")
print(f"Species: {cfg['project']['species']}")
print(f"BioProject: {cfg['project']['bioproject']}")
