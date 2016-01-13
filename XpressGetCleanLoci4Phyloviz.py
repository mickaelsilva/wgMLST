import sys
import csv
import numpy as np
from numpy import array
import argparse
import collections
from collections import OrderedDict
#import HTSeq
import matplotlib.pyplot as plt



def clean (inputfile,outputfile,totaldeletedgenes,rangeFloat,toremovegenes):
	
	#open the raw file to be clean
	
	with open(inputfile) as f:
		reader = csv.reader(f, delimiter="\t")
		d = list(reader)
	
	d2 = array(d)
	
	
	genomeslist= d2[1:,:1]
	
	d2=d2.T
	rowid=1
	deleted=0
	badAlleles=False
	lnfdel=0
	balldel=0
	
	
	
	numbergenomes=len(genomeslist)
	pontuationmatrix=[0]*numbergenomes	
	lostgenesList=[]
	
	while rowid< d2.shape[0]:
		columnid=0
		genomeindex=0
		
		
		#print badAlleles
		while columnid< d2.shape[1]:
			#print d2[rowid][0]
			if d2[rowid][0] in toremovegenes:
				
				d2=np.delete(d2, rowid, 0)
				totaldeletedgenes+=1
				deleted+=1
				rowid-=1
				columnid=1
				lnfdel+=1
				break
			
			elif ( "PLOT" in d2[rowid][columnid] or "NIPL" in d2[rowid][columnid] or "LNF" in d2[rowid][columnid] or "LOT" in d2[rowid][columnid]  or "ALM" in d2[rowid][columnid]  or "ASM" in d2[rowid][columnid] or "ERROR" in d2[rowid][columnid]  or "ABM" in d2[rowid][columnid]) or "undefined" in d2[rowid][columnid] or "small match" in d2[rowid][columnid] or "allele incomplete" in d2[rowid][columnid]:	
				
				d2=np.delete(d2, rowid, 0)
				totaldeletedgenes+=1
				deleted+=1
				rowid-=1
				columnid=1
				lnfdel+=1
				break
			
			if badAlleles :
				allele=d2[rowid][columnid]
				allele=allele.replace("INF-","",)
				allele=allele.replace('NA1-', '')
				allele=allele.replace('NA2-', '')
				allele=allele.replace('NA3-', '')
				allele=allele.replace('NA4-', '')
				allele=allele.replace('NA5-', '')
				try:
					if int(allele) in badAlleles:

						d2=np.delete(d2, rowid, 0)
						totaldeletedgenes+=1
						deleted+=1
						rowid-=1
						columnid=1
						balldel+=1
						break
				except Exception as e:
					print e
					pass
			
			columnid+=1
			genomeindex+=1
			#pontuationmatrix=[x + y for x, y in zip(pontuationmatrix, aux)]
		rowid+=1
	
	
	
	d2=d2.T
	d2=d2.tolist()
	
	
	with open(outputfile, "wb") as f:
		writer = csv.writer(f,delimiter='	')
		writer.writerows(d2)
	
	

	file = open(outputfile)
	contents = file.read()
	contents = contents.replace('INF-', '')
	contents = contents.replace('INF1:-', '')
	contents = contents.replace('INF2:-', '')
	contents = contents.replace('INF3:-', '')
	contents = contents.replace('INF4:-', '')
	contents = contents.replace('INF5:-', '')
	contents = contents.replace('INF6:-', '')
	contents = contents.replace('INF7:-', '')
	contents = contents.replace('NA1-', '')
	contents = contents.replace('NA2-', '')
	contents = contents.replace('NA3-', '')
	contents = contents.replace('NA4-', '')
	contents = contents.replace('NA5-', '')
	
	with open(outputfile, 'w') as f:
			f.write(contents)
	
	#if deleted==0:
	print str(lnfdel),str(balldel)
	print "deleted : %s genes" % totaldeletedgenes
	print "total genes remaining : "+ str(rowid)
	#else:
		#clean(outputfile,outputfile,totaldeletedgenes,shortGenomeList,genesDirect,rangeFloat)
	

def main():

	parser = argparse.ArgumentParser(description="This program cleans an output file for phyloviz")
	parser.add_argument('-i', nargs='?', type=str, help='output to clean', required=True)
	parser.add_argument('-g', nargs='?', type=str, help='name of the clean file', required=True)
	parser.add_argument('-r', nargs='?', type=str, help='listgenes to remove', required=False)
	
	args = parser.parse_args()

	
	pathOutputfile = args.i
	newfile = args.g
	
	
	
	
	genesToRemove=[]
	
	try:
		genesToRemoveFile = (args.r)
		fp = open(genesToRemoveFile, 'r')
		
		for geneFile in fp:

			geneFile = geneFile.rstrip('\n')
			geneFile = geneFile.rstrip('\r')
			
			genesToRemove.append( geneFile )
			
	except:		
		pass
	
	#print matrix
	
	clean(pathOutputfile,newfile,0,0.2,genesToRemove)

	

	
	
	
	
		
if __name__ == "__main__":
	main()
