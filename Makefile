obj = GENBANK_FASTA

$(obj) : test.o
	cc -o GENBANK_FASTA test.o -L. $(pwd)  -lfasta -lgenbank -I.

test.o : test.c
	cc -c test.c

clean : 
	rm GENBANK_FASTA test.o
