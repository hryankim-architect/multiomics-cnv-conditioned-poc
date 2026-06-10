.PHONY: canary test lint check compare-cnv

canary:        ## deterministic <1s smoke
	PYTHONPATH=src python -m mocnv.canary

compare-cnv:   ## CNV amplicon prior vs unsupervised baselines (needs sibling dmoi labels)
	PYTHONPATH=src python scripts/compare_cnv_prior.py

test:          ## unit tests (synthetic fixtures; no torch/data needed for v0.0)
	PYTHONPATH=src python -m pytest -q

lint:          ## ruff
	ruff check .

check: lint test  ## the full pre-push gate
