from pathlib import Path
import yaml

ROOT = Path(__file__).resolve().parents[2]


def load_config(path: str | Path | None = None) -> dict:
    cfg_path = Path(path) if path else ROOT / "configs" / "project.yaml"
    with cfg_path.open("r", encoding="utf-8") as f:
        return yaml.safe_load(f)
