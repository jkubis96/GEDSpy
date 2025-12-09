.PHONY: format lint all

format:
	isort gedspy
	black gedspy
	isort tests
	black tests


lint:
	pylint --exit-zero gedspy/DataPrepare.py
	pylint --exit-zero gedspy/Enrichment.py



all: format lint
