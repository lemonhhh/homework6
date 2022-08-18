obj = GENBANK_FASTA

$(obj) : test.o
	cc -o GENBANK_FASTA test.o -L./lib -I./src/include -lfasta -lgenbank

test.o : test.c
	cc -c test.c

clean : 
	rm GENBANK_FASTA test.o
