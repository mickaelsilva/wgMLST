#!/usr/bin/python
import HTSeq
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
import os
from CommonFastaFunctions import runBlastParser
import time
import pickle

def reverseComplement(strDNA):

	basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        strDNArevC = ''
        for l in strDNA:

        	strDNArevC += basecomplement[l]

        return strDNArevC[::-1]



def printinfo(genome, gene):
	
	print "Genome : "+ str(os.path.basename(genome)) 
	print "Locus : "+ str(os.path.basename(gene)) 

# ======================================================== #
#            Allele calling and classification             #
# ======================================================== #
def main():
	
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

	basepath=temppath+"/"+os.path.basename(geneFile)

	if not os.path.exists(basepath+"/blastdbs/"):
		os.makedirs(basepath+"/blastdbs/")
	
	
	gene_fp = HTSeq.FastaReader(geneFile)
	geneDict = {}
	alleleI = 1
	inverted=False
	orderedAlleleNames=[]
	biggestAllelelen=0
	smallestAllelelen=999999
	for allele in gene_fp:
		if allele.seq in geneDict:
			print "\nWARNING: this file contains a repeated allele, it should be checked. Ignoring it now!\n", geneFile
		else:
			if len(allele.seq)>biggestAllelelen:
				biggestAllelelen=len(allele.seq)
			if len(allele.seq)<smallestAllelelen:
				smallestAllelelen=len(allele.seq)
			orderedAlleleNames.append(str(alleleI))
			geneDict[ allele.seq ] = alleleI
			alleleI += 1

	# --- make 1st blast DB --- #

	geneF = os.path.basename(geneFile)
	blast_out_file = os.path.dirname(geneFile)+"/blastdbs/"+geneF + '.xml'

	# list of results - the output of the function
	i = 0
	perfectMatchIdAllele=[]
	perfectMatchIdAllele2=[]
	genomeDict = {}
	genome=-1	
	resultsList=[]
	print genomesList
	for genomeFile in genomesList:
		print "_______________________________________________________"
		print perfectMatchIdAllele
		printinfo(genomeFile,geneFile)
		#currentCDSDict = listOfCDSDicts[i]
		
		g_fp = HTSeq.FastaReader( genomeFile )
		for contig in g_fp:
			sequence=str(contig.seq)
			genomeDict[ contig.name ] = sequence
		
		currentGenomeDict = genomeDict
		
		genome+=1
		
		print ("Blasting alleles on genome at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))
		
		blast_out_file = os.path.join(basepath,"blastdbs/"+os.path.basename(geneFile)+ '_List.xml')


		Gene_Blast_DB_name = os.path.join(temppath,str(os.path.basename(genomeFile))+"/"+str(os.path.basename(genomeFile))+"_db")

		
		cline = NcbiblastnCommandline(query=geneFile, db=Gene_Blast_DB_name, evalue=0.001, out=blast_out_file, outfmt=5)

			
		blast_records = runBlastParser(cline, blast_out_file, geneFile)
		print ("Blasted alleles on genome at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))

		# ------ DETERMINING BEST MATCH ------ #

		bestMatch = ['','', 0]
		bestMatchContig=''
		bestMatchContigLen=''
		bestalignlen=0
		perfectMatch=False
		bmAlleleLen2=0
		bmAllele=''
		for blast_record in blast_records:

			if perfectMatch==True:
				break
			try:
				hspC = blast_record.alignments[0]
				
				if bestMatch[0] == '' and bestMatch[1] == '':
					bestMatch[0] = blast_record.query
					bestMatch[1] = hspC
			except IndexError:
				continue


			# --- the contig tag is used in the progigal function --- #

			contigTag = blast_record.query
			
			
			# --- brute force parsing of the contig tag - better solution is advisable --- #			

			j=0
			for l in contigTag:
				if l == ' ':
					break
				j+=1

			contigTag = contigTag[:j]

			contigLen = blast_record.query_letters
			
			# --- iterating over all the results to determine the best match --- #
			for alignment in blast_record.alignments:
				contigTag = alignment.hit_def
				contigTag=(contigTag.split(" "))[0]

				index=orderedAlleleNames.index(str(blast_record.query_id).split("_")[1])
				
				for k, v in geneDict.iteritems():
					if v == index+1:
						bmAlleleLen2= len(k)
					
				if perfectMatch:
					break
				for match in alignment.hsps:

					scoreRatio = float(match.score) / float(bmAlleleLen2)
					

					#if #identities is the same as the length of the allele and it has no gaps or N's
					if (int(match.identities)==int(bmAlleleLen2) and int(match.identities)==int(len(match.query)) and "N" not in match.sbjct and "K" not in match.sbjct and "Y" not in match.sbjct and "R" not in match.sbjct ): 
						
						index=orderedAlleleNames.index(str(blast_record.query_id).split("_")[1])
						for seq, alleleid in geneDict.iteritems():
							if alleleid == index+1:
								bmAllele=seq
								break
						bmAlleleLen= len(bmAllele)
						
						lenratio=float(len(match.sbjct))/float(bmAlleleLen)
						bestMatch = [blast_record.query, match, scoreRatio, blast_record.query_id,lenratio,bmAlleleLen]
						bestMatchContig=contigTag
						perfectMatch=True
						index=orderedAlleleNames.index(str(blast_record.query_id).split("_")[1])
						bmAlleleLen= len(geneDict.keys()[index])
						break
						
					#chose the match with the best score ratio (score/length of allele)
					elif scoreRatio > bestMatch[2]:
						index=orderedAlleleNames.index(str(blast_record.query_id).split("_")[1])
						for seq, alleleid in geneDict.iteritems():
							if alleleid == index+1:
								bmAllele=seq
								break
						bmAlleleLen= len(bmAllele)
						lenratio=float(len(match.sbjct))/float(bmAlleleLen)
						bestMatch = [blast_record.query, match, scoreRatio, blast_record.query_id,lenratio,bmAlleleLen]
						bestMatchContig=contigTag
						bestMatchContigLen=len(currentGenomeDict[contigTag])
						print contigTag
						bestalignlen=alignment.length						
						
					
					if perfectMatch==True:
						break


		# ---------- ALLELE CALLING AFTER DETERMINING BEST MATCH ---------- #
		print ("Finished choosing best match at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))
		
		
		try:
			match = bestMatch[1]
			bestMatchStart=match.sbjct_start
			bestMatchEnd=match.sbjct_end
			if match.query_start > match.query_end:
				bestMatchEnd=match.sbjct_start
				bestMatchStart=match.sbjct_end
			
			print match

			geneLen = bestMatch[5]
			alleleStr = match.sbjct
			nIdentities = match.identities
			idPercent = float(nIdentities) / float(geneLen)
			scoreRatio = bestMatch[2]
			lenRatio = bestMatch[4]
		
		except:
			#if no best match was found
			
			###################
			# LOCUS NOT FOUND #
			###################
			
			perfectMatchIdAllele.append('LNF')
			perfectMatchIdAllele2.append('LNF')
			
			print "Locus not found, no matches \n"
			continue
		
		print "is perfect match true?" +str(perfectMatch)
		if perfectMatch is True:
			
			#if a perfect match was found
			

			try:
				alleleNumber = geneDict[ alleleStr ]
			except:
				alleleStr=reverseComplement(alleleStr)
				alleleNumber = geneDict[ alleleStr ]
			
			################################################
			# EXACT MATCH --- MATCH == GENE --- GENE FOUND #
			################################################
			if "_" in bestMatch[3]:
				a=bestMatch[3].split("_")
				perfectMatchIdAllele.append(a[1])
				perfectMatchIdAllele2.append(str(bestMatchContig)+"&"+str(bestMatchStart)+"-"+str(bestMatchEnd)+"&"+"+")
			else:
				perfectMatchIdAllele.append(bestMatch[3])
				perfectMatchIdAllele2.append(str(bestMatchContig)+"&"+str(bestMatchStart)+"-"+str(bestMatchEnd)+"&"+"+")
			printinfo(genomeFile,geneFile)
			print "Exact match \n"
			continue
						
			

		else:
						
			#if a best match was found but it's not an exact match	

					###########################
					# LOCUS ON THE CONTIG TIP #
					###########################
			print geneLen
			if bestMatchContigLen <= geneLen:
								
				perfectMatchIdAllele.append('LOTSC')
				perfectMatchIdAllele2.append(str(bestMatchContig)+"&"+str(bestMatchStart)+"-"+str(bestMatchEnd)+"&"+"+")
				printinfo(genomeFile,geneFile)
				print "Locus is bigger than the contig \n"
			
			
			elif (match.sbjct_start ==1 and len(match.query) < geneLen) or (match.sbjct_start == bestMatchContigLen and len(match.query) < bestMatchContigLen and match.sbjct_start > match.sbjct_end):
			
				perfectMatchIdAllele.append('LOT5')
				perfectMatchIdAllele2.append(str(bestMatchContig)+"&"+str(bestMatchStart)+"-"+str(bestMatchEnd)+"&"+"+")
				printinfo(genomeFile,geneFile)
				
				print "Locus is on the 5' tip of the contig \n"
			
			elif (match.sbjct_end ==1 and len(match.query) < geneLen and match.sbjct_start > match.sbjct_end) or (match.sbjct_end == bestMatchContigLen and len(match.query) < bestMatchContigLen):
			
				perfectMatchIdAllele.append('LOT3')
				perfectMatchIdAllele2.append(str(bestMatchContig)+"&"+str(bestMatchStart)+"-"+str(bestMatchEnd)+"&"+"+")
				printinfo(genomeFile,geneFile)
				
				print "Locus is on the 3' tip of the contig \n"
			
	

			elif 'N' in alleleStr or "K" in alleleStr or "R" in alleleStr or "Y" in alleleStr:

					#####################
					# ALLELE NOT FOUND  #		# N base found!
					#####################
				
				geneFile2= os.path.splitext(geneFile)[0] + "LNFN.fasta"
				with open(geneFile2, 'a') as f:
					f.write(">"+ (str(os.path.basename(genomeFile)))+"|"+(str(os.path.basename(geneFile)))+"\n")
					f.write((alleleStr) +"\n")
				perfectMatchIdAllele.append('LNFN')
				perfectMatchIdAllele2.append('LNFN')
				printinfo(genomeFile,geneFile) 
				print "LNFN, contains strange (N,K,R) bases! \n"
			
			
			
			else:
				
				print "new allele?"
				#removing gaps
					
				alleleStr = alleleStr.replace('-', '')
				lenExtraThresh=int(biggestAllelelen*0.2)
			

				#else: #check if best match without gaps are contained inside an already defined allele

				isContainedDefinedAllele = False	
				definedAllele=''
				definedAlleleName=''

				for k in geneDict.keys():
					if alleleStr in k:
						definedAllele=k
						isContainedDefinedAllele = True
						definedAlleleName=geneDict.get(k)
						break
				print "is contained? " + str(isContainedDefinedAllele)
				print idPercent
				print geneLen
				print lenExtraThresh
				print lenRatio
				
				if isContainedDefinedAllele  and int(len(match.sbjct))<=int(len(definedAllele))+lenExtraThresh and int(len(match.sbjct))>=int(len(definedAllele))-lenExtraThresh :
					#allele without gaps is contained in a defined allele
					#best match with gaps has same size +1/-1 base as the defined allele
					
					isnewallele=False
						
					if int(len(alleleStr))==int(len(definedAllele)): # if match without gaps has same size as the defined allele 
						tagAux = 'NA?:'
						printinfo(genomeFile,geneFile) 
						perfectMatchIdAllele.append("NA?-"+str(alleleI))
						perfectMatchIdAllele2.append(str(bestMatchContig)+"&"+str(bestMatchStart)+"-"+str(bestMatchEnd)+"&"+"+")
						isnewallele=True
						
					elif int(len(alleleStr))==int(len(definedAllele))-1 : # if match without gaps has minus one base than the defined allele
						
						tagAux = 'NA2:'
						printinfo(genomeFile,geneFile) 
						perfectMatchIdAllele.append("NA2-"+str(alleleI))
						perfectMatchIdAllele2.append(str(bestMatchContig)+"&"+str(bestMatchStart)+"-"+str(bestMatchEnd)+"&"+"+")
						isnewallele=True

						
					else:
							extraleft=0
							extraright=0
							tS=0
							tE=0

							handle = open(genomeFile, "rU")
							record_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
							handle.close()
							record= record_dict[bestMatchContig]
							
							# if match without gaps has more than one base missing comparing to the defined allele 
							if (1<int(match.query_start) and 1<int(match.query_end)):
								
								if match.query_start > match.query_end:
									extraleft=match.query_end-1
									
								else:
									extraleft=match.query_start-1
							
							print 	extraleft, 	extraright		
									
							if (int(geneLen)>int(match.query_start) and int(geneLen)>int(match.query_end) ): # if 3' tip bases of the allele are missing on the match
								
								
								if match.query_start > match.query_end:
									extraright=geneLen-match.query_start
									
								else:
									extraright=geneLen-match.query_end
									
							print 	extraleft, 	extraright
							
							
							if match.sbjct_start > match.sbjct_end:
								tE=match.sbjct_start+extraleft
								tS=match.sbjct_end-extraright-1
								alleleStr=str(record.seq[tS:tE])
								alleleStr = reverseComplement(alleleStr)
							else:
								tS=match.sbjct_start-extraleft-1
								tE=match.sbjct_end+extraright
								alleleStr=str(record.seq[tS:tE])
							
							print tS
							print tE
							print "allele is:"
							print alleleStr
							
							if tE> bestMatchContigLen:
								perfectMatchIdAllele.append('LOT3B')
								perfectMatchIdAllele2.append(str(bestMatchContig)+"&"+str(tS)+"-"+str(bestMatchContigLen)+"&"+"+")
								printinfo(genomeFile,geneFile)
								
								print "Locus is on the 3B' tip of the contig \n"
							
							elif tS<0:
								perfectMatchIdAllele.append('LOT5B')
								perfectMatchIdAllele2.append(str(bestMatchContig)+"&"+str(0)+"-"+str(tE)+"&"+"+")
								printinfo(genomeFile,geneFile)
								
								print "Locus is on the 5B' tip of the contig \n"
						
						
							else:
						
								tagAux = 'NA2:'
								printinfo(genomeFile,geneFile) 
								perfectMatchIdAllele.append("NA2-"+str(alleleI))
								perfectMatchIdAllele2.append(str(bestMatchContig)+"&"+str(tS)+"-"+str(tE)+"&"+"+")
								isnewallele=True
						
					
					if isnewallele:
						print "New allele found! Adding allele "+ tagAux + str(alleleI) +" to the database"
						geneDict[alleleStr] = alleleI
							
							
						orderedAlleleNames.append(str(alleleI))						
						# --- add the new allele to the gene fasta --- #
							
						fG = open( geneFile, 'a' )
						fG.write('>allele_' + str(alleleI) + '_' + tagAux[:-1] +'_' + str(os.path.basename(genomeFile)) + '\n')
						fG.write( alleleStr + '\n')
						fG.close()
						alleleI += 1
					
				#if best match is not contained in an already defined allele, check if it has similar size with the match allele and has 0.8 similarity
				
				elif not isContainedDefinedAllele and idPercent >= 0.8 and int(len(match.sbjct))<=int(geneLen)+lenExtraThresh and int(len(match.sbjct))>=int(geneLen)-lenExtraThresh :
					#best match with gaps has 80% identity
					#best match with gaps is the same size or +1/-1 as the defined allele
					
					ratio=float(len(alleleStr)) / float(geneLen)
					
					if ratio>=0.8 and ratio<=1.2: # if match without gaps has same size as the best match allele and 80%similarity
						
						tagAux = ''
						extraleft=0
						extraright=0
						tS=0
						tE=0

						handle = open(genomeFile, "rU")
						record_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
						handle.close()
						record= record_dict[bestMatchContig]

						#if len(match.sbjct)<geneLen and "-" not in match.sbjct:  #if the allele is not fully covered against the match, compensate the tips
						try:
							print match
							if (1<int(match.query_start) and 1<int(match.query_end)):
								
								if match.query_start > match.query_end:
									extraleft=match.query_end-1
									
								else:
									extraleft=match.query_start-1
							
							print 	extraleft, 	extraright		
									
							if (int(geneLen)>int(match.query_start) and int(geneLen)>int(match.query_end) ): # if 3' tip bases of the allele are missing on the match
								
								
								if match.query_start > match.query_end:
									extraright=geneLen-match.query_start
									
								else:
									extraright=geneLen-match.query_end
									
							print 	extraleft, 	extraright
							
							
							if match.sbjct_start > match.sbjct_end:
								tE=match.sbjct_start+extraleft
								tS=match.sbjct_end-extraright-1
								alleleStr=str(record.seq[tS:tE])
								alleleStr = reverseComplement(alleleStr)
							else:
								tS=match.sbjct_start-extraleft-1
								tE=match.sbjct_end+extraright
								alleleStr=str(record.seq[tS:tE])
							
							print tS
							print tE
							print "allele is:"
							print alleleStr
							
							if tE> bestMatchContigLen:
								perfectMatchIdAllele.append('LOT3C')
								perfectMatchIdAllele2.append(str(bestMatchContig)+"&"+str(tS)+"-"+str(bestMatchContigLen)+"&"+"+")
								printinfo(genomeFile,geneFile)
								
								print "Locus is on the 3C' tip of the contig \n"
							
							elif tS<0:
								perfectMatchIdAllele.append('LOT5C')
								perfectMatchIdAllele2.append(str(bestMatchContig)+"&"+str(0)+"-"+str(tE)+"&"+"+")
								printinfo(genomeFile,geneFile)
								
								print "Locus is on the 5C' tip of the contig \n"
							
							else:
								tagAux = 'NA3:'
								printinfo(genomeFile,geneFile) 
								perfectMatchIdAllele.append("NA3-"+str(alleleI))
								perfectMatchIdAllele2.append(str(bestMatchContig)+"&"+str(tS)+"-"+str(tE)+"&"+"+")
									
								print "New allele found! Adding allele "+ tagAux + str(alleleI) +" to the database"
								geneDict[alleleStr] = alleleI
									
									
								orderedAlleleNames.append(str(alleleI))
								# --- add the new allele to the gene fasta --- #
									
								fG = open( geneFile, 'a' )
								fG.write('>allele_' + str(alleleI) + '_' + tagAux[:-1] +'_' + str(os.path.basename(genomeFile)) + '\n')
								fG.write( alleleStr + '\n')
								fG.close()
								alleleI += 1

						except Exception as e:
							##################
							#       LNF      #
							##################
							print e
							geneFile2= os.path.splitext(geneFile)[0] + "LNF3.fasta"
							print geneFile2
							with open(geneFile2, 'a') as f:
								f.write(">"+ (str(os.path.basename(genomeFile)))+"|"+(str(os.path.basename(geneFile)))+" | "+str(bestMatchContig)+"\n")
								f.write((alleleStr) +"\n")
								f.write(">Allele\n")
								f.write((bmAllele)+"\n")
							printinfo(genomeFile,geneFile) 
							perfectMatchIdAllele.append("LNF3")
							perfectMatchIdAllele2.append("LNF3")
							print "No allele found"
					else:
						##################
						#       LNF      #
						##################
						geneFile2= os.path.splitext(geneFile)[0] + "LNF4.fasta"
						print geneFile2
						with open(geneFile2, 'a') as f:
							f.write(">"+ (str(os.path.basename(genomeFile)))+"|"+(str(os.path.basename(geneFile)))+" | "+str(bestMatchContig)+"\n")
							f.write((alleleStr) +"\n")
							f.write(">Allele\n")
							f.write((bmAllele)+"\n")
						printinfo(genomeFile,geneFile) 
						perfectMatchIdAllele.append("LNF4")
						perfectMatchIdAllele2.append("LNF4")
						print "No allele found"
					

						
				elif isContainedDefinedAllele:
					####################
					# UNDEFINED ALLELE #		# it is contained in another allele
					####################
						
					alleleStr=match.query

					perfectMatchIdAllele.append("undefined allele")
					perfectMatchIdAllele2.append("undefined allele")
					printinfo(genomeFile,geneFile) 
					print "Undefined allele \n"
					
					geneFile2= os.path.splitext(geneFile)[0] + "undefined.fasta"
					print geneFile2
					

				elif lenRatio < 0.5:
						
					###############
					# SMALL MATCH #
					###############
								
					perfectMatchIdAllele.append('small match')
					perfectMatchIdAllele2.append('small match')
					printinfo(genomeFile,geneFile) 
					print "lower than 50% match \n"	
							
				elif lenRatio < 0.8 and idPercent < 0.5:
					#####################
					# INCOMPLETE ALLELE #		# it was not possible to extend it to at least 80% of the length of the gene
					#####################
					perfectMatchIdAllele.append('allele incomplete')
					perfectMatchIdAllele2.append('allele incomplete')
					printinfo(genomeFile,geneFile)
					print "Incomplete allele\n"
						
				else:	
					##################
					#       LNF      #
					##################
					
					printinfo(genomeFile,geneFile) 
					perfectMatchIdAllele.append("LNF5")
					perfectMatchIdAllele2.append("LNF5")
					print "Locus not found"
						
							
	final =	(resultsList,perfectMatchIdAllele)	
	print ("Finished allele calling at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))
	filepath=os.path.join(temppath , os.path.basename(geneFile)+"_result.txt")
	filepath2=os.path.join(temppath , os.path.basename(geneFile)+"_result2.txt")
	with open(filepath, 'wb') as f:
		pickle.dump(final, f)
	with open(filepath2, 'wb') as f:
		pickle.dump(perfectMatchIdAllele2, f)
	return True
	
if __name__ == "__main__":
    main()
