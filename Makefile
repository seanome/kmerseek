PYTHON ?= python

all: build

.PHONY:

build: .PHONY
	$(PYTHON) -m pip install -e .

clean:
	rm -rf tests/testdata/fasta/*sig.zip*

test: clean
	pytest