.PHONY: check test audit clean-ast features datasets train evaluate

check:
	python scripts/00_check_workspace.py

test:
	pytest -q

audit:
	python scripts/02_audit_raw_data.py

clean-ast:
	python scripts/03_clean_ast.py

features:
	python scripts/04_build_amr_features.py

datasets:
	python scripts/05_build_ml_datasets.py

train:
	python scripts/06_train_models.py

evaluate:
	python scripts/07_evaluate_models.py
