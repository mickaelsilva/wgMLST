#!/usr/bin/python
import sys
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastpCommandline

import os
import argparse
import time
import drmaa
import pickle
import shutil

# ================================================ MAIN ================================================ #

def main():
	
	#time ./alleleCalling_ORFbased_main.py -i asd.txt -g allffn.txt -o out_all_fnn_spades.txt -p True
	
	parser = argparse.ArgumentParser(description="This program screens a set of genes in a fasta file.")
	parser.add_argument('-i', nargs='?', type=str, help='List of genome files (list of fasta files)', required=True)
	parser.add_argument('-g', nargs='?', type=str, help='List of genes (fasta)', required=True)
	parser.add_argument('-o', nargs='?', type=str, help="Name of the output files", required=True)

	args = parser.parse_args()
	
	genomeFiles = args.i
	genes = args.g
	phylovizinput=True	
	
	print ("Starting Script at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))

	listOfCDSDicts = []
	listOfGenomes = []
	listOfGenomesDict = []

	fp = open(genomeFiles, 'r')

	for genomeFile in fp:

		genomeFile = genomeFile.rstrip('\n')
		genomeFile = genomeFile.rstrip('\r')
		listOfGenomes.append( genomeFile )
		genomeDict = {}

	fp.close()
	
	
	
	
	
	gene_fp = open( genes, 'r')

	genepath=''
	basepath=''
	lGenesFiles = []
	argumentsList = []
	for gene in gene_fp:
		gene = gene.rstrip('\n')
		lGenesFiles.append( gene )
		genepath=os.path.dirname(gene)
		basepath=os.path.join(genepath, "temp")
		if not os.path.exists(basepath):
			os.makedirs(basepath)
		filepath=os.path.join(basepath,str(os.path.basename(gene))+"_argList.txt")
		with open(filepath, 'wb') as f:
			var = [gene, listOfGenomes]
			pickle.dump(var, f)
		argumentsList.append(filepath)
		

	gene_fp.close()
	
	
	print ("Starting Genome Blast Db creation at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))
	joblist =[]
	with drmaa.Session() as s:
		for genomeFile in listOfGenomes:
			print genomeFile
			filepath=genomeFile
			os.makedirs(os.path.join(basepath,str(os.path.basename(genomeFile)) ))
			
			#Create_Blastdb2( filepath,os.path.join(basepath,str(os.path.basename(genomeFile)) ),str(os.path.basename(genomeFile)) )

			jt = s.createJobTemplate()
			jt.remoteCommand = os.path.join(os.getcwd(), 'Create_Genome_Blastdb.py')
			jt.args = [filepath,os.path.join(basepath,str(os.path.basename(genomeFile)) ),str(os.path.basename(genomeFile)),"nucl" ]
			jt.joinFiles=True
			jt.nativeSpecification='-V'
			jobid = s.runJob(jt)
			joblist.append(jobid)
			with open("jobsid.txt","a") as f:
				f.write(str(genomeFile)+"\n"+str(jobid))
			print('Your job has been submitted with ID %s' % jobid)

			s.deleteJobTemplate(jt)
		s.synchronize(joblist, drmaa.Session.TIMEOUT_WAIT_FOREVER, True)
		

	print
	print ("Starting Allele Calling at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))
	

	totloci= len(argumentsList)
	joblist =[]
	with drmaa.Session() as s:
		for argList in argumentsList:
			jt = s.createJobTemplate()
			jt.remoteCommand = os.path.join(os.getcwd(), 'callAlleles_probe_based2.py')
			jt.args = [str(argList),basepath]
			jt.joinFiles=True
			jt.nativeSpecification='-V'
			jobid = s.runJob(jt)
			joblist.append(jobid)
			with open("jobsid.txt","a") as f:
				f.write(str(argList)+"\n"+str(jobid))
			print('Your job has been submitted with ID %s' % jobid)

			s.deleteJobTemplate(jt)
		for curjob in joblist:
			print 'Collecting job ' + curjob
			retval = s.wait(curjob, drmaa.Session.TIMEOUT_WAIT_FOREVER)
			print 'Job: ' + str(retval.jobId) + ' finished with status ' + str(retval.hasExited)
	

	print ("Finished Allele Calling at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))
	
	output=[]
	for gene in lGenesFiles:
		filepath=os.path.join(basepath, os.path.basename(gene)+"_result.txt")
		with open(filepath,'rb') as f:
			var = pickle.load(f)
			output.append(var)
	
	output2=[]
	for gene in lGenesFiles:
		filepath2=os.path.join(basepath, os.path.basename(gene)+"_result2.txt")
		with open(filepath2,'rb') as f:
			var = pickle.load(f)
			output2.append(var)
	
			
	shutil.rmtree(basepath)
	
	print "##################################################\n %s genomes used for %s loci" % (len(output[0][0]),len(output) )
	numberexactmatches=0
	for gene in output:
		for gAllele in gene[0]:
			if("EXC:" in gAllele):
				numberexactmatches+=1
					
				
	print "\n %s exact matches found out of %s" % (numberexactmatches,(len(output[0][0])*len(output)) )	
	print "\n %s percent of exact matches \n##################################################" % (float((numberexactmatches*100)/(len(output[0][0])*len(output))) )	
		
	print "\nWriting output files\n"
	args.o = '/' + args.o
	

	try:
		phylovout=[]
		phylovout2=[]
		genesnames=[]
		statistics=[]
		
		for gene in lGenesFiles:
			
			genename=gene.split("/")
			#genename=genename[len(genename)-1].split(".")
			genename=genename[len(genename)-1]
			genesnames.append(genename)
		for geneOut in output:
			gene=0
			alleleschema=[]
			while gene<len(output[0][0]): 

				genename=(geneOut[1][gene]).split("_")

				if(len(genename)!=1):
					alleleschema.append(genename[1])
				else:
					alleleschema.append(genename[0])
				gene+=1
			phylovout.append(alleleschema)
		
		for geneOut in output2:
			#print str(geneOut)
			gene=0
			alleleschema=[]
			while gene<len(output2[0]): 

				genename=(geneOut[gene])
				#print genename
				#if(len(genename)!=1):
				#	alleleschema.append(genename[1])
				#else:
				
				alleleschema.append(genename)
				gene+=1
			phylovout2.append(alleleschema)

		
		genome=0
		finalphylovinput= "FILE"+ "\t" 
		finalphylovinput2= "FILE"+ "\t" 
		for geneid in genesnames:
			finalphylovinput+= str(geneid)+ "\t"
			finalphylovinput2+= str(geneid)+ "\t"
			
		
		while genome<len(listOfGenomes):
			currentGenome = os.path.basename(listOfGenomes[genome])
			statsaux=[0]*8 # EXC NA undef LNF LOT incomplete SAC 
			finalphylovinput+= "\n" + currentGenome + "\t"
			for gene in phylovout:
				
				val= str(gene[genome])
				finalphylovinput+= val + "\t"
				if "NA" in val:
					statsaux[1]+=1
				elif "undefined" in val:
					statsaux[2]+=1
				elif "LNF" in val:
					statsaux[3]+=1
				elif "PLOT" in val:
					statsaux[4]+=1
				elif "LOT" in val:
					statsaux[5]+=1
				elif "incomplete" in val:
					statsaux[6]+=1 
				elif "small" in val:
					statsaux[7]+=1
				else:
					statsaux[0]+=1
				
			genome+=1
			statistics.append(statsaux)
			
		genome=0	
		while genome<len(listOfGenomes):
			currentGenome = os.path.basename(listOfGenomes[genome])
			finalphylovinput2+= "\n" + currentGenome + "\t"
			for gene in phylovout2:
				
				val= str(gene[genome])
				finalphylovinput2+= val + "\t"
			
			genome+=1
				
			
		gOutFile = os.path.dirname( "./")
		gOutFile2 = os.path.dirname( "./")
		gOutFile += args.o
		gOutFile2 += "contigsInfo.txt"
		statswrite='Stats:\tEXC\tNA\tundefined\tLNF\tPLOT\tLOT\tincomplete\tsmall'
		i=0
		genome=0
		while genome<len(listOfGenomes):
			currentGenome = os.path.basename(listOfGenomes[genome])
			statsaux=[0]*8 # EXC NA undef LNF LOT incomplete SAC
			statswrite+= "\n" + currentGenome + "\t"
			for k in statistics[i]:
				statswrite+= str(k) + "\t"
			i+=1	
			genome+=1
				
		print statswrite
		with open(gOutFile, 'w') as f:
			f.write(finalphylovinput)
		statoutfile=os.path.dirname( "./")
		with open("stastics.txt", 'w') as f:
			f.write(str(statswrite))
		
		with open(gOutFile2, 'a') as f:
			f.write(str(finalphylovinput2))	
		
	except Exception as e:
		print e
		
		exc_type, exc_obj, tb = sys.exc_info()
		f = tb.tb_frame
		lineno = tb.tb_lineno
		print lineno
			
		

	print ("Finished Script at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))

if __name__ == "__main__":
    main()
	
