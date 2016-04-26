#!/usr/bin/python
import HTSeq
import sys
from Bio.Seq import Seq
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast.Applications import NcbiblastpCommandline
from Counter import Counter 
import os
import argparse
from CommonFastaFunctions import Create_Blastdb
from CommonFastaFunctions import runBlastParser
import time
import pickle
import shutil

def getBlastScoreRatios(genefile,basepath,doAll):
	
	gene_fp = HTSeq.FastaReader(genefile)
	alleleI=0
	allelescores=[]
	alleleProt=''
	alleleAllProt=''
	alleleList=[]
	for allele in gene_fp: #new db for each allele to blast it against himself
		alleleI+=1
		genome=-1
		alleleList.append(allele.seq)
		translatedSequence,x,y=translateSeq(allele.seq)
		
		if translatedSequence =='':
			pass
			
		else:	
			alleleProt=">"+str(alleleI)+"\n"+str(translatedSequence+"\n")
			alleleAllProt+=">"+str(alleleI)+"\n"+str(translatedSequence+"\n")
			proteinfastaPath=os.path.join(basepath,str(os.path.basename(genefile)+'_protein2.fasta'))
			
			with open(proteinfastaPath, "wb") as f:
				f.write(alleleProt)
			Gene_Blast_DB_name = Create_Blastdb( proteinfastaPath, 1, True )
			if doAll:
				
				blast_out_file = os.path.join(basepath,'blastdbs/temp.xml')
				print ("Starting Blast alleles at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))
				

				# --- get BLAST score ratio --- #
				cline = NcbiblastpCommandline(query=proteinfastaPath, db=Gene_Blast_DB_name, evalue=0.001, out=blast_out_file, outfmt=5)
				allelescore=0
			
				blast_records = runBlastParser(cline,blast_out_file, alleleProt)
			
				print ("Blasted alleles at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))
			
				for blast_record in blast_records:

					for alignment in blast_record.alignments:

						for match in alignment.hsps:
								
							allelescores.append(int(match.score))
							
				geneScorePickle=os.path.abspath(genefile)+'_bsr.txt'
				print "________"
				var=[alleleI,allelescores]
				with open(geneScorePickle,'wb') as f:
					pickle.dump(var, f)			
			
			else:
				geneScorePickle=os.path.abspath(genefile)+'_bsr.txt'
				with open(geneScorePickle,'rb') as f:
					var = pickle.load(f)
					allelescores=var[1]
				
	proteinfastaPath=os.path.join(basepath,str(os.path.basename(genefile)+'_protein.fasta'))
	with open(proteinfastaPath, "wb") as f:
			f.write(alleleAllProt)
			
			
	return int(alleleI),allelescores,alleleList
	
def reDogetBlastScoreRatios(genefile,basepath,alleleI,allelescores2,newGene_Blast_DB_name,alleleList2,picklepath):
	
	gene_fp = HTSeq.FastaReader(genefile)

	alleleProt=''
	
	alleleI+=1
		
	proteinfastaPath=genefile
	
	print ("Re-starting Blast alleles at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))
	
	blast_out_file2 = os.path.join(basepath,'blastdbs/temp.xml')

	cline = NcbiblastpCommandline(query=proteinfastaPath, db=newGene_Blast_DB_name, evalue=0.001, out=blast_out_file2, outfmt=5)
	allelescore=0
	blast_records = runBlastParser(cline,blast_out_file2, proteinfastaPath)
	
	print ("Blasted alleles at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))
	
	found =False
	for blast_record in blast_records:
		
		for alignment in blast_record.alignments:
			
			
			for match in alignment.hsps:
				allelescores2.append(int(match.score))
				

	var=[alleleI,allelescores2]
	with open(picklepath,'wb') as f:
		currentCDSDict = pickle.dump(var, f)
	
	return int(alleleI),allelescores2,alleleList2

def reverseComplement(strDNA):

	basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        strDNArevC = ''
        for l in strDNA:

        	strDNArevC += basecomplement[l]

        return strDNArevC[::-1]

def translateSeq(DNASeq):
	seq=DNASeq
	reversedSeq=False
	try:
		myseq= Seq(seq)
		protseq=Seq.translate(myseq, table=11,cds=True)
	except:
		reversedSeq=True
		try:
			seq=reverseComplement(seq)
			myseq= Seq(seq)
			protseq=Seq.translate(myseq, table=11,cds=True)
						
		except:
			try:
				seq=seq[::-1]
				myseq= Seq(seq)
				protseq=Seq.translate(myseq, table=11,cds=True)
			except:
				reversedSeq=False
				try:
					seq=seq[::-1]							
					seq=reverseComplement(seq)
					myseq= Seq(seq)
					protseq=Seq.translate(myseq, table=11,cds=True)
				except Exception as e:
					print "translated error"
					print e
					protseq=""
	return protseq,seq,reversedSeq


# ======================================================== #
#            Allele calling and classification             #
# ======================================================== #
def main():
	print ("Starting script at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))
	try:
		input_file = sys.argv[1]
		temppath = sys.argv[2]
	except IndexError:
		print "usage: list_pickle_obj"

	argumentList=[]
	with open(input_file,'rb') as f:
		argumentList = pickle.load(f)
	
	geneFile = argumentList[0]
	genomesList = argumentList[1]
	
	basepath=os.path.join(temppath,os.path.splitext(geneFile)[0])
	if not os.path.exists(basepath):
			os.makedirs(basepath)

	gene_fp = HTSeq.FastaReader(geneFile)
	alleleI = 0

	resultsList = []
	i = 0
	perfectMatchIdAllele=[]
	perfectMatchIdAllele2=[]
	allelescores=[]
	
	print ("Getting BSR at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))

	geneScorePickle=os.path.abspath(geneFile)+'_bsr.txt'
	
	#check if bsr as arealdy been calculated and recalculate it

	if os.path.isfile(geneScorePickle) :
		
		alleleI,allelescores,alleleList=getBlastScoreRatios(geneFile,basepath,False)
		
	else:	
		alleleI,allelescores,alleleList=getBlastScoreRatios(geneFile,basepath,True)
		
			
			
	print ("Finished BSR at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))
	genome=-1	
	
	genomeDict = {}
	print ("starting allele call at: "+time.strftime("%H:%M:%S-%d/%m/%Y"))
	for genomeFile in genomesList:
		print genomeFile
		bestmatch=[0,0,False,'',0] #score, score ratio, perfectmatch, key name of the DNA sequence string, allele ID
		currentGenomeDict={}
		currentCDSDict={}
		
		# load the translated CDS from the genome to a dictionary
		filepath=os.path.join(temppath,str(os.path.basename(genomeFile))+"_ORF_Protein.txt")
		with open(filepath,'rb') as f:
			currentCDSDict = pickle.load(f)
		
		#load the contig info of the genome to a dictionary
		g_fp = HTSeq.FastaReader( genomeFile )
		for contig in g_fp:
			sequence=str(contig.seq)
			genomeDict[ contig.name ] = sequence
		
		currentGenomeDict = genomeDict

		genome+=1
		listOfCDS=currentCDSDict
		genomeProteinfastaPath=os.path.join(temppath,str(os.path.basename(genomeFile)+'_Protein.fasta'))
		
		print ("Blasting alleles on genome at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))
		
		blast_out_file = os.path.join(basepath,"blastdbs/"+os.path.basename(geneFile)+ '_List.xml')

		Gene_Blast_DB_name = os.path.join(temppath,str(os.path.basename(genomeFile))+"/"+str(os.path.basename(genomeFile))+"_db")

		proteinfastaPath=os.path.join(basepath,str(os.path.basename(geneFile)+'_protein.fasta'))
		
		
		#blast the genome CDS against the translated locus
		cline = NcbiblastpCommandline(query=proteinfastaPath, db=Gene_Blast_DB_name, evalue=0.001, out=blast_out_file, outfmt=5)
			
		blast_records = runBlastParser(cline, blast_out_file, proteinfastaPath)
		print ("Blasted alleles on genome at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))
		
		alleleSizes=[]
		for allele in alleleList:
			alleleSizes.append(len(allele))
		
		biggestSizeAllele=0
		
		moda=max(set(alleleSizes), key=alleleSizes.count)
		contador= Counter(alleleSizes).most_common()
		
		if (contador[0])[1] ==1:
			moda= alleleSizes[0]

		try:
			
			# iterate through the blast results
			for blast_record in blast_records:
					
				locationcontigs=[]
				
				for alignment in blast_record.alignments:
					
					# select the best match
					for match in alignment.hsps:
						
						alleleMatchid=str(blast_record.query_id).split("_")[1]
						
						scoreRatio=float(match.score)/float(allelescores[int(alleleMatchid)-1])

						cdsStrName=((alignment.title).split(" "))[1]
						
						DNAstr=listOfCDS[">"+cdsStrName]

						AlleleDNAstr=alleleList[int(alleleMatchid)-1]
						if len(AlleleDNAstr)>biggestSizeAllele:
							biggestSizeAllele=len(AlleleDNAstr)
							
						compare=False
						
						#compare the DNA match and the allele DNA sequence (protein sequences may be equal and DNA different)
						if DNAstr==AlleleDNAstr is False:
							try:
								DNAstr=reverseComplement(DNAstr)
								if DNAstr==AlleleDNAstr is False:
									pass
								else:
									compare=True
							except:
								pass
						else:
							compare=True
						
						if scoreRatio>0.6:
							locationcontigs.append(cdsStrName)
							
						if "N" in DNAstr or "K" in DNAstr or "R" in DNAstr:
							pass
							
						elif(scoreRatio == 1 and bestmatch[2] is False and compare is True):
							bestmatch=[match.score,scoreRatio,True,cdsStrName,int(alleleMatchid),match,len(AlleleDNAstr)]

						elif(scoreRatio == 1 and match.score>bestmatch[0] and compare is True):
							bestmatch=[match.score,scoreRatio,True,cdsStrName,int(alleleMatchid),match,len(AlleleDNAstr)]

						elif(scoreRatio == 1 and bestmatch[2] is False and compare is False):
							bestmatch=[match.score,scoreRatio,False,cdsStrName,int(alleleMatchid),match,len(AlleleDNAstr)]
						
						elif(scoreRatio == 1 and match.score>bestmatch[0] and compare is False):
							bestmatch=[match.score,scoreRatio,False,cdsStrName,int(alleleMatchid),match,len(AlleleDNAstr)]

						elif(match.score>bestmatch[0] and scoreRatio>0.6 and scoreRatio>bestmatch[1] and bestmatch[2] is False):
							bestmatch=[match.score,scoreRatio,False,cdsStrName,int(alleleMatchid),match,len(AlleleDNAstr)]
							
										
			print ("Classifying the match at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))		
			
			#if no best match was found it's a Locus Not Found
			if bestmatch[0]==0 or "N" in AlleleDNAstr or "K" in AlleleDNAstr or "R" in AlleleDNAstr :
						
						###################
						# LOCUS NOT FOUND #
						###################
				if 	bestmatch[0]==0:		
					resultsList.append('LNF3:-1')
					perfectMatchIdAllele.append('LNF')
					perfectMatchIdAllele2.append('LNF')
					print "Locus not found, no matches \n"
				else:
					resultsList.append('LNFN:-1')
					perfectMatchIdAllele.append('LNF')
					perfectMatchIdAllele2.append('LNF')
					print "Locus has strange base (N, K or R) \n"
			
			#if more than one BSR >0.6 in two different CDSs it's a Non Paralog Locus
			elif len(list(set(locationcontigs)))>1:
				resultsList.append('NIPL')            
				perfectMatchIdAllele.append('NIPL')
				perfectMatchIdAllele2.append('NIPL')
				for elem in locationcontigs:
					print elem
				
			
			#in case the DNA match sequence equal to the DNA sequence of the comparing allele
			elif bestmatch[2] is True:
				contigname=bestmatch[3]	
				
				contigname=contigname.split("&")
				matchLocation=contigname[2]	
				contigname=contigname[0]	
				print contigname
				alleleStr=listOfCDS[">"+bestmatch[3]]
				protSeq,alleleStr,Reversed=translateSeq(alleleStr)
				

				#check for possible locus on tip
				match=bestmatch[5]
				matchLocation2=matchLocation.split("-")			
				seq=currentGenomeDict[ contigname ]
				bestMatchContigLen=len(seq)
				
				rightmatchContig=bestMatchContigLen-int(matchLocation2[1])	
				leftmatchContig=int(matchLocation2[0])
				
				if Reversed:
					aux=rightmatchContig
					rightmatchContig=leftmatchContig
					leftmatchContig=aux
				
				
				
				
				
				# get extra space to the right and left between the allele and match
				
				possibleExtra=int(moda)-((int(match.query_end)*3)-(int(match.query_start)*3))
				
				if possibleExtra<0:
					perfectMatchIdAllele.append(str(bestmatch[4]))
					if not Reversed:
						perfectMatchIdAllele2.append(str(contigname)+"&"+str(matchLocation)+"&"+"+")
					else:
						perfectMatchIdAllele2.append(str(contigname)+"&"+str(matchLocation)+"&"+"-")
					resultsList.append('EXC:' + str(bestmatch[4]) )
				
				else:	
					rightmatchAllele=possibleExtra
					leftmatchAllele=possibleExtra
					
					if leftmatchContig<leftmatchAllele and 	rightmatchContig < rightmatchAllele:
				
						resultsList.append('PLOTSC:-1')
						perfectMatchIdAllele.append('PLOTSC')
						perfectMatchIdAllele2.append('PLOTSC')

						print match
						print "contig extras (l,r)"
						print leftmatchContig,rightmatchContig
						print "allele extras (l,r)"
						print leftmatchAllele,rightmatchAllele
						
						print "Locus is possibly bigger than the contig \n"
					
					elif leftmatchContig<leftmatchAllele:
						
						
						resultsList.append('PLOT3:-1')
						perfectMatchIdAllele.append('PLOT3')
						perfectMatchIdAllele2.append('PLOT3')
						
						print match
						print "contig extras (l,r)"
						print leftmatchContig,rightmatchContig
						print "allele extras (l,r)"
						print leftmatchAllele,rightmatchAllele
						
						print "Locus is possibly on the 3' tip of the contig \n"
					
					
					elif 	rightmatchContig < rightmatchAllele:
						
						resultsList.append('PLOT5:-1')
						perfectMatchIdAllele.append('PLOT5')
						perfectMatchIdAllele2.append('PLOT5')
						
						print match
						print "contig extras (l,r)"
						print leftmatchContig,rightmatchContig
						print "allele extras (l,r)"
						print leftmatchAllele,rightmatchAllele

						print "Locus is possibly on the 5' tip of the contig \n"
				
					else:
						#if a perfect match was found
								
						################################################
						# EXACT MATCH --- MATCH == GENE --- GENE FOUND #
						################################################
								
						perfectMatchIdAllele.append(str(bestmatch[4]))
						if not Reversed:
							perfectMatchIdAllele2.append(str(contigname)+"&"+str(matchLocation)+"&"+"+")
						else:
							perfectMatchIdAllele2.append(str(contigname)+"&"+str(matchLocation)+"&"+"-")
						resultsList.append('EXC:' + str(bestmatch[4]) )

			
			# if match with BSR >0.6 and not equal DNA sequences
			else:
				
				match=bestmatch[5]
				geneLen=bestmatch[6]

				contigname=bestmatch[3]	
				
				contigname=contigname.split("&")
				matchLocation=contigname[2]	
				matchLocation=matchLocation.split("-")
				contigname=contigname[0]
				
				seq=currentGenomeDict[ contigname ]
				bestMatchContigLen=len(seq)
				
				alleleStr=listOfCDS[">"+bestmatch[3]]
				protSeq,alleleStr,Reversed=translateSeq(alleleStr)
				
				
				rightmatchContig=bestMatchContigLen-int(matchLocation[1])	
				leftmatchContig=int(matchLocation[0])
				
				if Reversed:
					aux=rightmatchContig
					rightmatchContig=leftmatchContig
					leftmatchContig=aux
				
				
				print rightmatchContig,leftmatchContig
				
				
				# get extra space to the right and left between the allele and match and check if it's still inside the contig
				
				rightmatchAllele=geneLen-((int(match.query_end)+1)*3)	
				leftmatchAllele=((int(match.query_start)-1)*3)
				

						###########################
						# LOCUS ON THE CONTIG TIP #
						###########################
				
				
				
				if leftmatchContig<leftmatchAllele and 	rightmatchContig < rightmatchAllele:
				
					resultsList.append('LOTSC:-1')
					perfectMatchIdAllele.append('LOTSC')
					perfectMatchIdAllele2.append('LOTSC')
					print match
					print contigname
					print geneFile
					print leftmatchAllele,rightmatchAllele
					print "Locus is bigger than the contig \n"
				
				elif leftmatchContig<leftmatchAllele:
					
					
					resultsList.append('LOT3:-1')
					perfectMatchIdAllele.append('LOT3')
					perfectMatchIdAllele2.append('LOT3')
					print match
					print contigname
					print geneFile
					print leftmatchAllele,rightmatchAllele
					print "Locus is on the 3' tip of the contig \n"
				
				
				elif 	rightmatchContig < rightmatchAllele:
					
					resultsList.append('LOT5:-1')
					perfectMatchIdAllele.append('LOT5')
					perfectMatchIdAllele2.append('LOT5')
					print match
					print contigname
					print geneFile
					print leftmatchAllele,rightmatchAllele
					print "Locus is on the 5' tip of the contig \n"
				
				
							
				elif len(alleleStr) > moda+(moda*0.2) :
					
					print moda
					print alleleStr
					resultsList.append('ALM')
					perfectMatchIdAllele.append('ALM')
					perfectMatchIdAllele2.append('ALM')
				
				elif len(alleleStr) < moda-(moda*0.2):
					
					print moda
					print alleleStr
					resultsList.append('ASM')
					perfectMatchIdAllele.append('ASM')
					perfectMatchIdAllele2.append('ASM')
			
					
				else:
							#######################
							# ADD INFERRED ALLELE #		# a new allele 
							#######################
							
													
					tagAux='INF'
					perfectMatchIdAllele.append( tagAux +"-"+str(alleleI+1))
					
					if not Reversed:
						perfectMatchIdAllele2.append(str(contigname)+"&"+str(matchLocation[0])+"-"+str(matchLocation[1])+"&"+"+")
					else:
						perfectMatchIdAllele2.append(str(contigname)+"&"+str(matchLocation[0])+"-"+str(matchLocation[1])+"&"+"-")
					
					
					print "New allele! Adding allele "+ tagAux + str(alleleI+1) +" to the database\n"
																						
					resultsList.append( tagAux + str(alleleI+1) )

												# --- add the new allele to the gene fasta --- #
					
					
					appendAllele='>allele_' + str(alleleI+1) + '_' + tagAux[:-1] +"_" + str(os.path.basename(genomesList[genome])) + '\n'
					fG = open( geneFile, 'a' )
					fG.write(appendAllele)
						
					fG.write( alleleStr + '\n')
					fG.close()
					
					fG = open( os.path.join(basepath,str(os.path.basename(geneFile)+'_protein2.fasta')), 'w' )
					fG.write('>'+str(alleleI+1)+'\n'+str(protSeq) + '\n')
					fG.close()
					fG = open( os.path.join(basepath,str(os.path.basename(geneFile)+'_protein.fasta')), 'a' )
					fG.write('>'+str(alleleI+1)+'\n'+str(protSeq) + '\n')
					fG.close()	
					
					match=bestmatch[5]
					
					# --- remake blast DB and recalculate the BSR for the locus --- #
					alleleList.append(alleleStr)
					print os.path.join(basepath,str(os.path.basename(geneFile)+'_protein.fasta'))
					genefile2= os.path.join(basepath,str(os.path.basename(geneFile)+'_protein2.fasta'))
					Gene_Blast_DB_name2 = Create_Blastdb( genefile2, 1, True )
					print ("Re-calculating BSR at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))
					alleleI,allelescores,alleleList=reDogetBlastScoreRatios(genefile2,basepath,alleleI,allelescores,Gene_Blast_DB_name2,alleleList,geneScorePickle)
					print "allele id " + str(alleleI)
					print ("Done Re-calculating BSR at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))
		
		except Exception as e:
			print "some error occurred"
			print e
			print 'Error on line {}'.format(sys.exc_info()[-1].tb_lineno)
			perfectMatchIdAllele2.append("ERROR")
			perfectMatchIdAllele.append("ERROR")
			resultsList.append('ERROR')  
		
	
	final =	(resultsList,perfectMatchIdAllele)	
	print ("Finished allele calling at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))
	filepath=os.path.join(temppath , os.path.basename(geneFile)+"_result.txt")
	filepath2=os.path.join(temppath , os.path.basename(geneFile)+"_result2.txt")
	with open(filepath, 'wb') as f:
		pickle.dump(final, f)
	with open(filepath2, 'wb') as f:
		pickle.dump(perfectMatchIdAllele2, f)
	shutil.rmtree(basepath)
	return True
	
if __name__ == "__main__":
    main()
