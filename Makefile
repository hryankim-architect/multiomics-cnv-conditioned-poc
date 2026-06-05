.PHONY: canary test lint check

canary:        ## deterministic <1s smoke
	PYTHONPATH=src python -m mocnv.canary

test:          ## unit tests (synthetic fixtures; no torch/data needed for v0.0)
	PYTHONPATH=src python -m pytest -q

lint:          ## ruff
	ruff check .

check: lint test  ## the full pre-push gate
