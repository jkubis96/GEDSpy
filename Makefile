.PHONY: format lint all

format:
	isort GEDSpy
	black GEDSpy
	isort tests
	black tests


lint:
	pylint --exit-zero GEDSpy/DataPrepare.py
	pylint --exit-zero GEDSpy/Enrichment.py



all: format lint
