"""Preferred M3 command; delegates to the existing AST entry point."""

from pathlib import Path
import runpy

if __name__ == "__main__":
    runpy.run_path(str(Path(__file__).with_name("03_clean_ast.py")), run_name="__main__")
