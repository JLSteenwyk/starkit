install:
	python3 setup.py install

develop:
	python3 setup.py develop

test: test.unit test.integration

test.unit:
	python3 -m pytest -m "not integration"

test.integration:
	rm -rf output/
	mkdir output/
	python3 -m pytest --basetemp=output -m "integration"

test.fast:
	python3 -m pytest -m "not (integration or slow)"

test.coverage: coverage.unit coverage.integration

coverage.unit:
	python3 -m pytest --cov=./ -m "not integration" --cov-report=xml:unit.coverage.xml

coverage.integration:
	rm -rf output/
	mkdir output/
	python3 -m pytest --basetemp=output --cov=./ -m "integration" --cov-report=xml:integration.coverage.xml
