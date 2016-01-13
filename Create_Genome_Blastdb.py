#!/usr/bin/python
import sys
import os
import os.path






def main():
	
	
	try:
		questionDB = sys.argv[1]
		directory = sys.argv[2]
		genomeFile = sys.argv[3]
	except IndexError:
		print "usage: Sum Ting Wong"
	
	
	
	base = os.path.basename(questionDB)
	dirname = os.path.dirname( questionDB )

	
	if len(dirname)==0:
		dirname='.'
	basename = os.path.splitext(base)[0]

	name = directory+"/"+genomeFile+"_db"
	

	
	try:
		nucleotide=sys.argv[4]
		os.system( "makeblastdb -in " + questionDB + " -out " + name + " -dbtype nucl -logfile " + name + "_blast.log" )
	except:
		os.system( "makeblastdb -in " + questionDB + " -out " + name + " -dbtype prot -logfile " + name + "_blast.log" )
	
	
	
	
	return True
	
if __name__ == "__main__":
    main()
