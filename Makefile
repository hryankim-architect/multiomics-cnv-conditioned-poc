.PHONY: canary test lint english-only check

canary:        ## deterministic <1s smoke
	PYTHONPATH=src python -m mocnv.canary

test:          ## unit tests (synthetic fixtures; no torch/data needed for v0.0)
	PYTHONPATH=src python -m pytest -q

lint:          ## ruff
	ruff check .

english-only:  ## fail on any CJK in tracked text artifacts
	PYTHONPATH=src python scripts/check_english_only.py

check: lint english-only test  ## the full pre-push gate
