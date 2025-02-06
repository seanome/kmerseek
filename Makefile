clean:
	rm -rf tests/testdata/fasta/*sig.zip*

test: clean
	pytest