PYTHON ?= .venv/bin/python
PIP    ?= .venv/bin/pip
JUPYTER ?= .venv/bin/jupyter

.PHONY: setup test scripts notebook reproduce clean

## Create the environment and install pinned dependencies
setup:
	python3 -m venv .venv
	$(PIP) install -r requirements.txt

## Check the biological expectations and dropout maths
test:
	$(PYTHON) -m pytest -q tests

## Rebuild the provenance record, figures and tables
scripts:
	$(PYTHON) scripts/00_dataset_provenance.py
	$(PYTHON) scripts/01_build_atlas.py
	$(PYTHON) scripts/02_cd20_audit.py

## Re-execute the notebook in place so saved outputs match the scripts
notebook:
	$(JUPYTER) nbconvert --execute --to notebook --inplace notebooks/cd20_analysis.ipynb

## Full reproducibility: environment -> tests -> scripts -> notebook
reproduce: setup test scripts notebook

clean:
	rm -rf .venv
