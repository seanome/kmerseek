PYTHON ?= python
IMAGE ?= kmerseek
VERSION ?= $(shell grep '^version' Cargo.toml | head -1 | sed 's/.*"\(.*\)".*/\1/')

all: build

.PHONY:

build: .PHONY
	$(PYTHON) -m pip install -e .

clean:
	rm -rf tests/testdata/fasta/*sig.zip*
	rm -rf tests/testdata/fasta/*manysketch.csv

test: clean
	pytest

docker-build:
	docker build -t $(IMAGE):$(VERSION) -t $(IMAGE):latest .

docker-push:
	docker push $(IMAGE):$(VERSION)
	docker push $(IMAGE):latest