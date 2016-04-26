import csv
import numpy as np
from numpy import array
import argparse
import collections
from collections import OrderedDict
import matplotlib.pyplot as plt
import Counter


# function to report the most problematic locus/genomes and report the number of good locus by % of genomes
# a list of genomes with a sum of problems higher than the given ythreshold will be returned
def presence3(d2,ythreshold,vector,abscenceMatrix):
	print "calculating presence abscence ..."
	d2d=d2
	
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
			
			if not abscenceMatrix:
				
				if  "LNF" in d2d[row,column] or d2d[row,column] =="NIPL" or "LOT" in d2d[row][column]  or "ALM" in d2d[row][column]  or "ASM" in d2d[row][column] or "ABM" in d2d[row][column] or "ERROR" in d2d[row][column] or "undefined" in d2d[row][column] or "small match" in d2d[row][column] or "allele incomplete" in d2d[row][column]:	

				#if d2d[row,column] == "LNF" :
					d2d[row,column]=0
					notfound+=1
					badgenomes.append(row)
					
				else:
					d2d[row,column]=1
					pass
			else:
				if int(d2d[row,column]) == 0:
					notfound+=1
					badgenomes.append(row)
					
				else:
					d2d[row,column]=1
					pass
				
				
			row+=1
		
		value=float(d2d.shape[0]-notfound)/float(d2d.shape[0])
		if len(genomeslist)>500:
			xthreshold=0.99
		elif len(genomeslist)>200:
			xthreshold=0.97
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
	
	counter=Counter.Counter(allbadgenomes)
	
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
	
	
	#print "this genomes have more than " +str(ythreshold)+" bad calls in locus present in more than "+str(xthreshold*100)+"% of the genomes :"
	#for x in reallybadgenomes:
	#	print x
	

	
	d2d=d2d.T
	
	
	#number of used genomes
	vector[0].append(len(genomeslist))
	
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
	
	print "presence abscence calculated"
	
	return d2d,reallybadgenomes,vector,True
	

	
def removegenomes(d2a,bagenomeslist):
	
	rowid=1
	deleted=0
	genomesList= (d2a[1:,:1]).tolist()
	
	print len(genomesList)
	print len(bagenomeslist)
	print "removing genomes..."
	badgenomeslist2=[]
	for elem in bagenomeslist:
		badgenomeslist2.append(elem)
	
	if len(bagenomeslist)>0:
		i=0
		for genome in genomesList:
			i+=1
			if genome[0] in badgenomeslist2:
				d2a=np.delete(d2a, i, axis=0)
				deleted+=1
				i-=1

	print "genomes removed"
	print d2a.shape
	return d2a


def clean (inputfile,iterations,ythreshold):
	
	#open the raw file to be clean
	
	"""with open(inputfile) as f:
		reader = csv.reader(f, delimiter="\t")
		d = list(reader)"""
	
	print "will try to open file..."
	#d2=np.loadtxt(inputfile, delimiter='\t')
	d2 = np.genfromtxt(inputfile, delimiter='\t', dtype=None)
	print "file was read"
	#d2 = array(d)
	
	

		
	
	i=0
	removedlistgenomes=[]
	toremovegenomes=[]
	
	statsvector=[[] for x in range(7)]
	
	#run a function that gives information about the usable locus and the possible bad genomes
	lastremovedgenomesCount=0
	iterationStabilizedat=None
	isStable=False
	abscencematrix=False
	while i<=iterations:
		
		print "\n########## ITERATION NUMBER %s  ##########  \n" % str(i)
		#print "genomes to be removed:"
		#print toremovegenomes
		
		for elem in toremovegenomes:
			removedlistgenomes.append(elem)

		if len(removedlistgenomes)>lastremovedgenomesCount and i >0:
			lastremovedgenomesCount=len(removedlistgenomes)
		elif iterationStabilizedat is None and i >0:
			iterationStabilizedat=i
			print "stabilized at "+str(i)
			isStable=True
		if not isStable:
			print "total removed genomes :" +str(len(removedlistgenomes))
			d2=removegenomes(d2,toremovegenomes)
			matrix3,toremovegenomes,statsvector,abscencematrix=presence3(d2,ythreshold,statsvector,abscencematrix)
		
		else:
			for vector in statsvector:
				vector.append(vector[-1])
		
		i+=1
		

	
	#run a function that returns a matrix with 0 and 1 depending on wheter the call is a Locus Not Found or not
	
	with open("removedGenomes.txt", "a") as f:
		f.write("using a threshold of "+ str(ythreshold)+" at iteration number " +str(i)+"\n")
		
		for x in removedlistgenomes:
			f.write(x+"\n")
	
	
	
	return statsvector,iterationStabilizedat

def main():

	parser = argparse.ArgumentParser(description="This program cleans an output file for phyloviz")
	parser.add_argument('-i', nargs='?', type=str, help='output to clean', required=True)
	parser.add_argument('-n', nargs='?', type=int, help='number of iterations', required=True)
	parser.add_argument('-t', nargs='?', type=int, help='threshold number of bad calls above 99%', required=True)
	parser.add_argument('-s', nargs='?', type=int, help='step', required=True)
	
	args = parser.parse_args()

	
	pathOutputfile = args.i
	iterationNumber=int(args.n)
	thresholdBadCalls=int(args.t)
	step=int(args.s)
	
	
	#print matrix
	
	allresults=[]
	threshold=5
	thresholdlist=[]
	listStableIter=[]
	
	#for each threshold run a clean function on the dataset, using the previous run output (to be removed genomes) as input for the new one
	
	while threshold<thresholdBadCalls:
		
		thresholdlist.append(threshold)
		print "########## USING A THRESHOLD AT " +str(threshold) + " ########"
		result,stabilizedIter=clean(pathOutputfile,iterationNumber,threshold)
		listStableIter.append(stabilizedIter)
		allresults.append(result)
		
		#step=20
		threshold+=step
		
	
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
		fig, ax1 = plt.subplots(figsize=(30.5,20.0))
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

		
		plt.savefig("Iteration_"+str(i)+'.png')
		i+=1
		
	
	i=0
	for stableiter in listStableIter:
		
		print "At threshold "+str(thresholdlist[i]) + " it stabilized at the iteration number "+str(stableiter)
		i+=1
	
	
	#plt.show()
	plt.close()		
	
	
	
		
if __name__ == "__main__":
	main()
