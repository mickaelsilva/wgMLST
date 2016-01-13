#!/usr/bin/python
import sys
import os

def main(questionDB,directory,genomeFile,nucleotide):
	

	base = os.path.basename(questionDB)
	dirname = os.path.dirname( questionDB )

	
	if len(dirname)==0:
		dirname='.'
	basename = os.path.splitext(base)[0]

	name = directory+"/"+genomeFile+"_db"
	
	print name
	
	if nucleotide:
		os.system( "makeblastdb -in " + questionDB + " -out " + name + " -dbtype nucl -logfile " + name + "_blast.log" )
	else:
		os.system( "makeblastdb -in " + questionDB + " -out " + name + " -dbtype prot -logfile " + name + "_blast.log" )
	
	
	
	
	return True
	
if __name__ == "__main__":
    main()
