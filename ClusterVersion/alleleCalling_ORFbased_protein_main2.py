#!/usr/bin/python
import HTSeq
import sys
from Bio.Seq import Seq
import os
import argparse
from CommonFastaFunctions import runBlastParser
import time
import drmaa
import pickle
import shutil


def reverseComplement(strDNA):

	basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'R':'R', 'N':'N', 'K':'K'}
        strDNArevC = ''
        for l in strDNA:

        	strDNArevC += basecomplement[l]

        return strDNArevC[::-1]

def translateSeq(DNASeq):
	seq=DNASeq
	try:
		myseq= Seq(seq)
		#print myseq
		protseq=Seq.translate(myseq, table=11,cds=True)
	except:
		try:
			seq=reverseComplement(seq)
			myseq= Seq(seq)
			#print myseq
			protseq=Seq.translate(myseq, table=11,cds=True)
						
		except:
			try:
				seq=seq[::-1]
				myseq= Seq(seq)
				#print myseq
				protseq=Seq.translate(myseq, table=11,cds=True)
			except:
				try:
					seq=seq[::-1]							
					seq=reverseComplement(seq)
					myseq= Seq(seq)
					#print myseq
					protseq=Seq.translate(myseq, table=11,cds=True)
				except Exception as e:
					print "translated error"
					print e
					raise
	return protseq
	
def loci_translation (genesList,listOfGenomes2):
	
	gene_fp = open( genesList, 'r')

	genepath=''
	basepath=''
	lGenesFiles = []
	argumentsList = []
	
	for gene in gene_fp:
		k=0
		gene = gene.rstrip('\n')
		multiple=True
		gene_fp2 = HTSeq.FastaReader(gene)
		for allele in gene_fp2: 
			k+=1
			if (len(allele.seq) % 3 != 0):
				multiple=False
				print "allele "+str(k)+" is not multiple of 3: "+str(gene)+" this gene is to be removed ffs"
				break
			else:
				try:
					protseq=translateSeq(allele.seq)
				except:
					multiple=False
					print "allele "+str(k)+" is not translating : "+ str(gene)+" this gene is to be removed ffs"
					break
		if multiple:
			lGenesFiles.append( gene )
			genepath=os.path.dirname(gene)
			basepath=os.path.join(genepath, "temp")
			if not os.path.exists(basepath):
				os.makedirs(basepath)
			filepath=os.path.join(basepath,str(os.path.basename(gene))+"_argList.txt")
			with open(filepath, 'wb') as f:
				var = [gene, listOfGenomes2]
				pickle.dump(var, f)
			argumentsList.append(filepath)

	

	gene_fp.close()	
	return	genepath,basepath,lGenesFiles,argumentsList

# ================================================ MAIN ================================================ #

def main():
	
	print(os.path.dirname(os.path.abspath(__file__)))
	
	parser = argparse.ArgumentParser(description="This program screens a set of genes in a fasta file.")
	parser.add_argument('-i', nargs='?', type=str, help='List of genome files (list of fasta files)', required=True)
	parser.add_argument('-g', nargs='?', type=str, help='List of genes (fasta)', required=True)
	parser.add_argument('-o', nargs='?', type=str, help="Name of the output files", required=True)
	parser.add_argument('-p', nargs='?', type=str, help="Path to prodigal exec file", required=True)

	args = parser.parse_args()
	
	genomeFiles = args.i
	genes = args.g
	prodigalPath = args.p

		
	
	
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
	
	genepath,basepath,lGenesFiles,argumentsList=loci_translation (genes,listOfGenomes)

	if len(argumentsList) ==0:
		print "ERROR! AT LEAST ONE GENE FILE MUST CONTAIN AN ALLELE REPRESENTING A CDS"
		raise ValueError('ERROR! AT LEAST ONE GENE FILE MUST CONTAIN AN ALLELE REPRESENTING A CDS')
		
	# ------------------------------------------------- #
	#           RUN PRODIGAL OVER ALL GENOMES           #
	# ------------------------------------------------- #




	print ("Starting Prodigal at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))

	totgenomes= len(listOfGenomes)
	

	joblist =[]
	with drmaa.Session() as s:
		for genome in listOfGenomes:
			jt = s.createJobTemplate()
			jt.remoteCommand = os.path.join(os.getcwd(), 'runProdigal.py')
			jt.args = [str(genome),basepath]
			jt.joinFiles=True
			jt.nativeSpecification='-V'
			jobid = s.runJob(jt)
			joblist.append(jobid)
			with open("jobsid.txt","a") as f:
				f.write(str(genome)+"\n"+str(jobid))
			print('Your job has been submitted with ID %s' % jobid)

			s.deleteJobTemplate(jt)
		s.synchronize(joblist, drmaa.Session.TIMEOUT_WAIT_FOREVER, True)

	
	print ("Finishing Prodigal at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))
	
	#---CDS to protein---#
	listOFProt=[]
	listAllCDS=[]
	
	i = 0
	j=0
	
	
	#translate the genome CDSs, load them into dictionaries and fasta files to be used further ahead
	for genomeFile in listOfGenomes:
		listOfCDS={}
		genomeProts=""
		currentCDSDict = {}
		currentGenomeDict = {}
		filepath=os.path.join(basepath,str(os.path.basename(genomeFile))+"_ORF.txt")
		with open(filepath,'rb') as f:
			currentCDSDict = pickle.load(f)
			
		g_fp = HTSeq.FastaReader( genomeFile )
		for contig in g_fp:
			sequence=str(contig.seq)
			genomeDict[ contig.name ] = sequence
		
		currentGenomeDict = genomeDict
		
		i+=1
		for contigTag,value in currentCDSDict.iteritems():
			for protein in value:
				try:
					seq= currentGenomeDict[ contigTag ][ protein[0]:protein[1] ].upper()
					protseq=translateSeq(seq)
					idstr=">"+contigTag+"&protein"+str(j)+"&"+str(protein[0])+"-"+str(protein[1])
					genomeProts+=idstr+"\n"
					listOfCDS[idstr]=seq
					genomeProts+=str(protseq)+"\n"
					
				except Exception as e:
					print str(e)+" "+str(genomeFile)
					pass
				j+=1
		filepath=os.path.join(basepath,str(os.path.basename(genomeFile))+"_ORF_Protein.txt")
		with open(filepath, 'wb') as f:
			var = listOfCDS
			pickle.dump(var, f)

		filepath=os.path.join(basepath,str(os.path.basename(genomeFile))+"_Protein.fasta")
		with open(filepath, 'wb') as f:
			f.write(genomeProts)
			
	
	print ("Starting Genome Blast Db creation at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))
	
	#creation of the Databases for each genome, one genome per core using n-2 cores (n number of cores)

	
	with drmaa.Session() as s:
		for genomeFile in listOfGenomes:
			
			filepath=os.path.join(basepath,str(os.path.basename(genomeFile))+"_Protein.fasta")
			os.makedirs(os.path.join(basepath,str(os.path.basename(genomeFile)) ))
			
			jt = s.createJobTemplate()
			jt.remoteCommand = os.path.join(os.getcwd(), 'Create_Genome_Blastdb.py')
			jt.args = [filepath,os.path.join(basepath,str(os.path.basename(genomeFile)) ),str(os.path.basename(genomeFile)) ]
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
	
	# Run the allele call, one gene per core using n-2 cores (n number of cores)


	totloci= len(argumentsList)
	joblist =[]
	with drmaa.Session() as s:
		for argList in argumentsList:
			jt = s.createJobTemplate()
			jt.remoteCommand = os.path.join(os.getcwd(), 'callAlleles_protein2.py')
			jt.args = [str(argList),basepath]
			jt.joinFiles=True
			jt.nativeSpecification='-V'
			jobid = s.runJob(jt)
			joblist.append(jobid)
			with open("jobsid.txt","a") as f:
				f.write(str(argList)+"\n"+str(jobid))
			print('Your job has been submitted with ID %s' % jobid)

			s.deleteJobTemplate(jt)
		s.synchronize(joblist, drmaa.Session.TIMEOUT_WAIT_FOREVER, True)

	
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


	print ("Finished Allele Calling at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))
	
	#delete all temp files
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
	
	#wrapping up the results into a matrix and save in a file
	try:
		phylovout=[]
		phylovout2=[]
		genesnames=[]
		statistics=[]
		
		
		for gene in lGenesFiles:
			
			genename=gene.split("/")
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
			gene=0
			alleleschema=[]
			while gene<len(output2[0]): 

				genename=(geneOut[gene])
				
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
			statsaux=[0]*8 # EXC INF LNF LOT incomplete SAC
			finalphylovinput+= "\n" + currentGenome + "\t"
			for gene in phylovout:
				
				val= str(gene[genome])
				finalphylovinput+= val + "\t"
				if "INF" in val:
					statsaux[1]+=1
				elif "LNF" in val:
					statsaux[2]+=1
				elif "PLOT" in val:
					statsaux[4]+=1
				elif "LOT" in val:
					statsaux[3]+=1 
				elif "NIPL" in val:
					statsaux[5]+=1
				elif "ALM" in val:
					statsaux[6]+=1
				elif "ASM" in val:
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
		gOutFile += args.o
		statswrite='Stats:\tEXC\tINF\tLNF\tLOT\tPLOT\tNIPL\tALM\tASM'
		i=0
		genome=0
		while genome<len(listOfGenomes):
			currentGenome = os.path.basename(listOfGenomes[genome])
			statsaux=[0]*8 # EXC NA INF LNF LOT incomplete SAC
			statswrite+= "\n" + currentGenome + "\t"
			for k in statistics[i]:
				statswrite+= str(k) + "\t"
			i+=1	
			genome+=1
		
		
		print statswrite
		with open(gOutFile, 'a') as f:
			f.write(finalphylovinput)
		statoutfile=os.path.dirname( "./")
		with open("stastics.txt", 'a') as f:
			f.write(str(statswrite))
		
		with open("contigsInfo.txt", 'a') as f:
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
