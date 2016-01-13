# wgMLST
tools developed to perform a bacterial wgMLST

Dependencies:
* drmma	http://drmaa-python.github.io/
* biopython http://biopython.org/wiki/Main_Page
* HTSeq http://www-huber.embl.de/users/anders/HTSeq/doc/overview.html
* BLAST
* Prodigal http://prodigal.ornl.gov/
* CommonFastaFunctions.py
* Create_Genome_Blastdb.py


suggested folder structure:

1. main folder - scripts and txt files
 1. sub folder - genomes - all genomes fasta files
 2. sub folder - genes - all genes fasta files

**lists of files MUST contain FULL PATH!**

=============
#CreateSchema.py

dependencies:
* CommonFastaFunctions.py
* biopython
* HTSeq

Given a concatenated ffn file, removes genes that are substring of bigger genes and genes smaller than choosen in the -g parameter. Blasts all the genes against each other and saves the bigger genes, removing the smaller genes with a 0.6>BSR

	% CreateSchema.py -i allffnfile.fasta -g 200
	
Output:

* proteins.fasta containing the transaltion of all the genes from the given ffn file, without substring genes
* x.fasta large set of .fasta files, 1 per gene



=============
##alleleCalling_probe_based_main2.py

Given a list of genomes and a list of alleles, the program will perform an allele call using the defined alleles as probes

	% alleleCalling_probe_based_main2.py -i listGenomes.txt -g listGenes.txt -o outputFileName.txt -p True
	
`-i` path to the list of genomes file

`-g` path to the list of alleles file

`-o` output file name

`-p` (optional and recomended) parameter to return a phyloviz output file type and a statistics.txt file - Default = True

short example statistics file:

* EXC - allele has exact match (100% identity)
* NA - new allele found
* undefined - allele found is contained in an already defined allele but match size is more than 2 bases different from the defined allele
* LNF - locus not found
* LOT - locus is on the tip of the contig of the genome
* PLOT - locus is possibly on the tip of the contig of the genome
* incomplete - match size is less than 80% of the allele size and identity% is smaller than 50%
* small - match size is less than 50% of the allele size
```
Stats:	EXC	NA	undefined	LNF	LOT	PLOT	incomplete	small
NC_017162.fna	1026	4	0	0	0	0	0	0	
NC_009085.fna	248	1	0	0	0	0	0	0	

```
short example phyloviz file output:
```
FILE	Unique_Acinetobacter_baumannii_1656-2.1.peg.gi_384132717_ref_YP_005515329.1_.fasta	Unique_Acinetobacter_baumannii_1656-2.1.peg.gi_384133246_ref_YP_005515858.1_.fasta
10_S10_L001.fasta	7	1
2_S12_L001.fasta	7	1
7_S7_L001.fasta	NA6:-8	NA5:-7
```

=============
##alleleCalling_ORFbased_protein_main2.py

Given a list of genomes and a list of alleles, the program will perform an allele call using the defined alleles translated aminoacid sequences and the ORF translated sequences of a genome obtained using Prodigal

	% alleleCalling_ORFbased_protein_main2.py -i listgenomes.txt -g listgenes.txt -o output_file_name.txt -p True
	
`-i` path to the list of genomes file

`-g` path to the list of alleles file

`-o` output file name

`-p` (optional and recomended) parameter to return a phyloviz output file type and a statistics.txt file - Default = True


short example statistics file:

* EXC - allele has exact match (100% identity)
* INF - infered allele with prodigal
* LNF - locus not found
* LOT - locus on the tip of the contig
* PLOT - locus possibly on the tip of the contig
* NIPL - Non informative paralog locus - two or more good blast matches for the protein
* ALM - allele larger than the gene match - match > allele match size + allele match size * 0.2
* ASM - allele smaller than the gene match - match < allele match size - allele match size * 0.2

```
Stats:	EXC	INF	LNF	LOT	PLOT	NIPL	ALM	ASM
NC_017162.fna	892	2319	1909	0	0	104	5	37	
NC_011586.fna	1563	1697	1809	0	0	116	6	75	
```

short example phyloviz file output:

```
FILE	gi_126640115_ref_NC_009085.1_:1032446-1033294.fasta	gi_126640115_ref_NC_009085.1_:103903-104649.fasta	gi_126640115_ref_NC_009085.1_:1056402-1057004.fasta	gi_12664011510_S10_L001.fasta
NC_017162.fna	INF-2	LNF
NC_011586.fna	INF-3	LNF
NC_011595.fna	3	LNF
```

=============
## GetCleanLoci4Phyloviz.py

Dependencies:
* matplotlib
* numpy

Clean a raw output file from an allele calling to a phyloviz readable file. Keep the locus with only Exact matches or new alleles found for all genomes.

Basic usage:

	% GetCleanLoci4Phyloviz.py -i rawDataToClean.txt -g cleanedOutput.txt
	
`-i` raw output file from an allele calling

`-g` name of the clean file

`-o` information file, tab separated

`-p` which property # to group by

`-s` which group value to select by

Using an info file for a specific genome selection:

	% GetCleanLoci4Phyloviz.py -i rawDataToClean.txt -g cleanedOutput.txt -o info.tsv -p 4 -s 3
	
In this example we will be using the property #4 (serotype) and all genomes that have a serotype "3". Other genome profiles will be ignored.
	
short example info file:

```	
file	ST	CC	serotype	study
NC_003028.fna	205	34	4	completegenome
NC_003098.fna	595	378	2	completegenome
...
```	
