IMAGE ?= kmerseek
VERSION ?= $(shell grep '^version' Cargo.toml | head -1 | sed 's/.*"\(.*\)".*/\1/')

all: build

.PHONY:

build: .PHONY
	cargo build --release --no-default-features

clean:
	rm -rf tests/testdata/fasta/*sig.zip*
	rm -rf tests/testdata/fasta/*manysketch.csv

test: clean
	cargo test --no-default-features --lib -- --test-threads=2

install: build
	cargo install --path . --no-default-features

docker-build:
	docker build -t $(IMAGE):$(VERSION) -t $(IMAGE):latest .

docker-push:
	docker push $(IMAGE):$(VERSION)
	docker push $(IMAGE):latest
