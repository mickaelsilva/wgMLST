import sys
import csv
import numpy as np
from numpy import array
import argparse
import collections
from collections import OrderedDict
#import HTSeq
import matplotlib.pyplot as plt


#build presence/abscense matrix per locus considering it abscent if LNF
def presence (d2):	

	d2c=np.copy(d2)
	
	#print "Warning: all profiles must have same length"


	genomeslist= d2c[1:,:1]
	
	d2c = d2c[1:,1:-1]
	
	row=0
	while row<d2c.shape[0]:
		column=0
		while column<d2c.shape[1]:
			if d2c[row,column] == "LNF":
				d2c[row,column]=0
			else:
				d2c[row,column]=1
			column+=1
		row+=1
	
	d2c=d2c.T

	genomeslist=(genomeslist.tolist())

	return d2c,genomeslist


# function to report the most problematic locus/genomes and report the number of good locus by % of genomes
# a list of genomes with a sum of problems higher than the given ythreshold will be returned
def presence3(d2,ythreshold,vector):
	d2d=np.copy(d2)
	

	genomeslist= d2d[1:,:1]
	geneslist = d2d[0]
	
	d2d = d2d[1:,1:-1]
	
	column=0
	totals=[]
	allbadgenomes=[]
	reallybadgenomes=[]
	allverygood=[]
	allverybad=[]
	plus95=0
	plus99=0
	plus995=0
	
	# if a locus has one problematic call, check how many genomes have problem in that call
	
	while column<d2d.shape[1]:
		row=0
		notfound=0
		badgenomes=[]
		while row<d2d.shape[0]:
			if d2d[row,column] == "LNF" or d2d[row,column] =="NIPL" or "LOT" in d2d[row][column]  or "ALM" in d2d[row][column]  or "ASM" in d2d[row][column] or "ABM" in d2d[row][column] or "ERROR" in d2d[row][column] or "undefined" in d2d[row][column] or "small match" in d2d[row][column] or "allele incomplete" in d2d[row][column]:	
			#if d2d[row,column] == "LNF":
				d2d[row,column]=0
				notfound+=1
				badgenomes.append(row)
				
			else:
				d2d[row,column]=1
				pass
			row+=1
		
		value=float(d2d.shape[0]-notfound)/float(d2d.shape[0])
		if len(genomeslist)>100:
			xthreshold=0.99
		else:
			xthreshold=0.95

		if(value>xthreshold and value<1):
			for badgenome in badgenomes:
				allbadgenomes.append((genomeslist[int(badgenome)])[0])
		elif value==1:
			allverygood.append(geneslist[column+1])
		elif value==0:
			allverybad.append(geneslist[column+1])
		if value>0.95:
			plus95+=1
		if value>0.99:
			plus99+=1
		if value>0.995:
			plus995+=1

		totals.append(float(float(d2d.shape[0]-notfound)/float(d2d.shape[0])))
		column+=1
	
	counter=collections.Counter(allbadgenomes)
	
	#print "Ordered count of genomes with bad calls, considering only a bad call when the locus is present in more than "+str(xthreshold*100)+"% of the genomes"
	
	for elem in counter.most_common():
		if int(elem[1])>ythreshold:
				reallybadgenomes.append((elem)[0])
	#print counter
	print
	print "number genes in 100% genomes: " +str(len(allverygood))
	print "number genes above 99%: " +str(plus99)
	print "number genes above 95%: " +str(plus95)
	print "number genes above 99.5%: " +str(plus995)
	print "number genes in 0% genomes: " +str(len(allverybad))
	print len (totals)
	
	
	
	d2d=d2d.T
	
	
	#number of used genomes
	vector[0].append(len(genomeslist))
	#print vector
	
	#number of loci at 95%
	vector[1].append(plus95)
	
	#number of loci at 99%
	vector[2].append(plus99)
	
	#number of loci at 100%
	vector[3].append(len(allverygood))
	
	#number of loci at 0%
	vector[4].append(len(allverybad))
	
	#number of to be removed genomes%
	try:
		vector[5].append(len(reallybadgenomes)+((vector[5])[-1]))
	except:
		vector[5].append(len(reallybadgenomes))
	
		#number of loci at 99.5%
	vector[6].append(plus995)
	
	return d2d,reallybadgenomes,vector
	

def lostGenesPontuation(matrix,genomeslist):
	
	numbergenomes=len(genomeslist)
	pontuationmatrix=[0]*numbergenomes
	
	for gene in matrix:
		aux=[0]*numbergenomes
		presentgenes=0
		genomeindex=0
		notpresentgenes=0
		for value in gene:
			if int(value)==0:
				aux[genomeindex]=aux[genomeindex]-1
				notpresentgenes+=1
				
			else:
				presentgenes+=1
			genomeindex+=1

		if notpresentgenes==1:		
			pontuationmatrix=[x + y for x, y in zip(pontuationmatrix, aux)]
				
			
	
	cleangenomelist=[]
	for genome in genomeslist:
		cleangenomelist.append(genome[0])
		
	ordered=OrderedDict(zip(cleangenomelist,pontuationmatrix ))
	ordered = OrderedDict(sorted(ordered.items(), key=lambda(k,v):(v,k)))
	ordered=OrderedDict(ordered.items()[::-1])
	#print "List of negative pontuation by locus loss responsability. -1 point per genome when that genome is the sole responsible for that locus loss"
	#print ordered
	

	#for genome in ordered.items():
		#print genome
		
	return False
	
def removegenomes(d2a,shortGenomeList):
	
	rowid=1
	deleted=0
	
	while rowid< d2a.shape[0]:
		inside=False
		for genome in shortGenomeList:
			
			if genome == d2a[rowid][0] :	
				inside=True
				break

		if not inside:
			d2a=np.delete(d2a, rowid, axis=0)
			deleted+=1
		rowid+=1

	
	if deleted>0:
		return removegenomes(d2a,shortGenomeList)
	else:
		return d2a
		
def removegenes(d2a,genesToRemove):
	d2a=d2a.T

	linenmbr=0
	for line in d2a:

		if line[0] in genesToRemove:
			
			d2a=np.delete(d2a,linenmbr,axis=0)
			linenmbr-=1
			
		linenmbr+=1
	return d2a.T


def clean (inputfile,totaldeletedgenes,shortGenomeList,genesDirect,rangeFloat,toremovegenes,iterations,ythreshold):
	
	#open the raw file to be clean
	
	with open(inputfile) as f:
		reader = csv.reader(f, delimiter="\t")
		d = list(reader)
	
	d2 = array(d)
	shortGenomeList2=[]
	if shortGenomeList:
		d2=removegenomes(d2,shortGenomeList)
		shortGenomeList= d2[1:,:1]
		shortGenomeList=shortGenomeList.tolist()
		
		for elem in shortGenomeList:
			shortGenomeList2.append(elem[0])
	else:
		shortGenomeList= d2[1:,:1]
		shortGenomeList=shortGenomeList.tolist()
		
		for elem in shortGenomeList:
			shortGenomeList2.append(elem[0])
	
	if len(toremovegenes)>1:
		d2=removegenes(d2,toremovegenes)
		
	
	i=0
	removedlistgenomes=[]
	toremovegenomes=[]
	
	statsvector=[[] for x in range(7)]
	
	#run a function that gives information about the usable locus and the possible bad genomes
	lastremovedgenomesCount=0
	iterationStabilizedat=None
	isStable=False
	while i<=iterations:
		
		print "\n########## ITERATION NUMBER %s  ##########  \n" % str(i)

		for elem in toremovegenomes:
			shortGenomeList2.remove(elem)
			removedlistgenomes.append(elem)

		if len(removedlistgenomes)>lastremovedgenomesCount and i >0:
			lastremovedgenomesCount=len(removedlistgenomes)
		elif iterationStabilizedat is None and i >0:
			iterationStabilizedat=i
			print "stabilized at "+str(i)
			isStable=True
		if not isStable:
			print "total removed genomes :" +str(len(removedlistgenomes))
			d2=removegenomes(d2,shortGenomeList2)
			matrix3,toremovegenomes,statsvector=presence3(d2,ythreshold,statsvector)
		
		else:
			for vector in statsvector:
				vector.append(vector[-1])
		
		i+=1
		

	matrix2,genomeslist2=presence (d2)
	
	

	#lostGenesPontuation(matrix3,genomeslist2)
	
	#run a function that returns a matrix with 0 and 1 depending on wheter the call is a Locus Not Found or not
	
	with open("removedGenomes2.txt", "a") as f:
		f.write("using a threshold of "+ str(ythreshold)+" at iteration number " +str(i)+"\n")
		
		for x in removedlistgenomes:
			f.write(x+"\n")
	
	
	
	return statsvector,iterationStabilizedat

def main():

	parser = argparse.ArgumentParser(description="This program cleans an output file for phyloviz")
	parser.add_argument('-i', nargs='?', type=str, help='output to clean', required=True)
	parser.add_argument('-o', nargs='?', type=str, help='info file', required=False)
	parser.add_argument('-p', nargs='?', type=int, help='group by property', required=False)
	parser.add_argument('-s', nargs='?', type=str, help='specify group', required=False)
	parser.add_argument('-d', nargs='?', type=str, help='genes directory', required=False)
	parser.add_argument('-r', nargs='?', type=str, help='listgenes to remove', required=False)
	parser.add_argument('-n', nargs='?', type=int, help='number of iterations', required=True)
	parser.add_argument('-t', nargs='?', type=int, help='threshold number of bad calls above 99%', required=True)
	
	args = parser.parse_args()

	
	pathOutputfile = args.i
	iterationNumber=int(args.n)
	thresholdBadCalls=int(args.t)
	
	
	groupselected=False
	
	#check if a specific property was selected
	
	try:
		info = args.o
		property = (args.p)-1
		groupselected = args.s
		groupselected=groupselected.split(',')
	except:
		pass
	genesdirect=False
	try:
		genesdirect = (args.d)
	except:
		pass
	print groupselected
	
	
	#if a specific property was selected, get the list of genomes that meet that selected property
	
	if groupselected:
		with open(info) as f:
			reader = csv.reader(f, delimiter="\t")
			oldlist = list(reader)
			
		print "sorted by : " +str((oldlist[0])[property])
		oldlist.remove(oldlist[0])

		
		values = set(map(lambda x:x[property], oldlist))
		newlist = [[y[0] for y in oldlist if y[property]==x] for x in values]
	
		shortGenomeList=[]
		values=list(values)
		for value in groupselected:
			index=values.index(value)
			listvalue=newlist[index]
			for genome in listvalue:
				shortGenomeList.append(genome)
	else:
		shortGenomeList=False
		
	print "list of selected genomes :"
	print shortGenomeList
	
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
	
	
	allresults=[]
	threshold=5
	thresholdlist=[]
	listStableIter=[]
	
	#for each threshold run a clean function on the dataset, using the previous run output (to be removed genomes) as input for the new one
	
	while threshold<=thresholdBadCalls:
		
		thresholdlist.append(threshold)
		print "########## USING A THRESHOLD AT " +str(threshold) + " ########"
		result,stabilizedIter=clean(pathOutputfile,0,shortGenomeList,genesdirect,0.2,genesToRemove,iterationNumber,threshold)
		listStableIter.append(stabilizedIter)
		allresults.append(result)
		
		threshold+=5
		
	
	x=list(range(0,len(allresults)))
	

	i=0
	
	
	labels=["number of genomes",
	"Number of Loci present in 95% genomes",
	"Number of Loci present in 99% genomes",
	"Number of Loci present in 100% genomes",
	#"Number of Loci present in 0% genomes",
	#"Selected genomes to be removed",
	"Number of Loci present in 99.5% genomes"]
		
	while i<=iterationNumber:
		fig, ax1 = plt.subplots(figsize=(20.5,10.0))
		threshindex=0
		
		aux2=[]
		for resultPerThresh in allresults:
			aux=[]
			
			for result in resultPerThresh:
				#resultPerThresh[i]
				aux.append(result[i])
			aux2.append(aux)
			
		d2=np.asarray(aux2)
			
		d2=d2.T
		linenmbr=0

		d2 = np.delete (d2,(4), axis=0)
		d2 = np.delete (d2,(4), axis=0)
		
		ax2 = ax1.twinx()
		
		
		for line in d2:
			
			if(linenmbr==0):
				ax2.plot(thresholdlist,line,linestyle='-',color='#663300', marker='v')
			else:	
				ax1.plot(thresholdlist,line, label = labels[linenmbr], linestyle='--', linewidth=1, marker='o')

			linenmbr+=1
		
		
		for tl in ax2.get_yticklabels():
			tl.set_color('#663300')
			
		plt.xticks(thresholdlist)
		
		ax1.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05))
		ax1.set_xlabel('Threshold')

		ax2.set_ylabel('Number of genomes in use', color='#663300')
		plt.title('Iteration number '+str(i),loc='left')

		
		
		i+=1
	
	i=0
	for stableiter in listStableIter:
		
		print "At threshold "+str(thresholdlist[i]) + " it stabilized at the iteration number "+str(stableiter)
		i+=1
	plt.show()
	plt.close()		
	
	
	
		
if __name__ == "__main__":
	main()
