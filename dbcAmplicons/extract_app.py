#extract_app

import sys, os, time, gzip
import time
start_time = time.time()
from dbcAmplicons import OneReadIlluminaRun
#from multiprocessing import Pool
import multiprocessing
#from itertools import product
from functools import partial



def extractFastq(vGoodReads, batch, name, batchName):
	vContinue=False					
	first_variable=""
	print "tset2"
	check_first_varable=True
	write_list=[]
	x=0
	with open(batchName, "r") as file:
		for read in (file): #loops over the file. 

		 	 	######## @M may not be a good variable as it may change with speceis. save first two characters of line 1 and sort by that?
			if check_first_varable==True: #this should fix problem above, but only if header names are relatively consistant throughout
		 	 	first_variable=read[:1]
		 		check_first_varable=False

		 	if first_variable in read:  #While '@' is used as a quality score, '@M' has so far only been found in header
		 	 	if read.split(" ",1)[0][1:] in vGoodReads: #checks to see if the header is in the good reads
		 	 		vContinue=True
					vGoodReads.remove(read.split(" ",1)[0][1:]) #if header is in good reads, it will be removed from good reads to save memory

		 	if vContinue==True: #If this is turned on, the next four lines (which is the length of the sample), will be recorded
		 	 	write_list.append(read)
		 	 	x+=1
			 	if x ==4:
					x=0
					vContinue=False

			if len(vGoodReads)==0 and vContinue==False:
				break
	file.close()
	fFileWrite('extOut_'+name+"_"+str(batch),write_list)
	print("test")

def fFileLength(file):
	#Check that files are the same length
	with open(file,"r") as file:
		for i, l in enumerate(file):
			pass
		return i + 1

def fRemoveFiles(vFastqLoops, vFixrankLoops, batch):
	i=0
	while i<=vFastqLoops:
		os.remove('extFastqTemp_'+str(batch)+"_"+str(i))
		i+=1
	i=0
	while i<=vFixrankLoops:
		os.remove('extFixrankTemp_'+str(batch)+"_"+str(i))
		i+=1

def fFileWrite(vFileName,vToWrite, vSolution=False, vBatchNumber="", overwrite=False):
	#Writes vDict to vFileName file.
	if os.path.isfile(vFileName+"_"+vBatchNumber) ==False: #checks to see if file exists and if it should be overwriten or appended
		vTempFile=open(vFileName+"_"+vBatchNumber,'w')
	elif overwrite==True:
		vTempFile=open(vFileName+"_"+vBatchNumber,'w')
	else:
		vTempFile=open(vFileName+"_"+vBatchNumber,'a')
	if vSolution==True:			#writes a file containing a multiple array 
		for i in range(len(vToWrite)):
			for y in range(4):
				vTempFile.write(vToWrite[i][y])
		vTempFile.close()
	else:
		for i in vToWrite:
			vTempFile.write(i)	
		vTempFile.close()

class extractApp():

	def __init__(self):
		self.verbose = False

	def start(self, species, fixrank, fastq, threshold, output, batchsize, order ="", family="", genus="", procs=1):
		starttime = time.time()

		#Format input files:
		print("---\nSpecies target is: " + species)
		vFixrankFile = open(fixrank, "r")
		if fastq[-3:] == '.gz':
			vFastqFile = gzip.open(fastq, "r")
		else:
			vFastqFile = open(fastq, "r")

		loops_fixrank=0
		loops_fastq=0

		#if __name__ == '__main__'

		#Run through the batches
		print "Begining to break fixed rank into batches"
		vFixBatchList = []
		vBatch_fixrank,n = 0,0
		loops_fixrank=0
		continue_loop=True
		while continue_loop==True: #breaks the fix rank file into a number of temporary files with lengths equal to the batch size
			with vFixrankFile as file: 
				for classification in file.read().split('\n'):	#writes a file containing as many lines of the fix rank file as the batch demands
					if n < batchsize and classification !="":
						if classification!="":
							vFixBatchList.append(classification+'\n')
						n+=1						
					if n == batchsize: #If a batch size is reached, the batch will be written to file, the batch will be cleared, 
										#and the counter for the number of loops for the fixed rank will increase
						fFileWrite('extFixrankTemp_'+str(batchsize), vFixBatchList, False, str(loops_fixrank), True)
						vFixBatchList = []
						loops_fixrank+=1
						n=0
						print "Fix rank batch #"+str(loops_fixrank)+" written"
					
				if len(vFixBatchList) >0:					
					vBatch_fixrank =n
					fFileWrite('extFixrankTemp_'+str(batchsize), vFixBatchList, False, str(loops_fixrank), True)
					vFixBatchList=[] #clear batch to save memory
					print "Fix rank batch #"+str(loops_fixrank+1)+" written"
					continue_loop=False
				else:
					continue_loop=False
			file.close()
			#Generate fastq batches
		n = 0
		k = 0
		j = 0
		print "Begining to break Fastq into batches"
		vList = []
		vTotal=[]
		continue_loop=True
		i=0
		vBatch_fastq=0
		loops_fastq=0
		while continue_loop ==True: #breaks the fastq file into a number of temporary files with lengths equal to the batch size
			with vFastqFile as file:     
				for line in file.readlines():
					if n < batchsize:
						if line!="":
							vList.append(line)
						k+=1
						if k % 4 == 0 and k!=0:  ## stores each of the four lines of the sample into one array 
							vTotal.append(vList)
							vList = []
							j+=1
							n+=1
					if n == batchsize:	#If a batch size is reached, the batch will be written to file, the batch will be cleared, 
											#and the counter for the number of loops for the fastq rank will increase
						fFileWrite('extFastqTemp_'+str(batchsize), vTotal, True, str(loops_fastq), True)
						n = 0
						j = 0
						vList = []
						vTotal = []
						loops_fastq+=1
						print "Fastq batch #"+str(loops_fastq)+" written"
				if len(vTotal)>0: #if the file ends before the batch limit is reached, will recorded
					fFileWrite('extFastqTemp_'+str(batchsize), vTotal, True, str(loops_fastq), True)
					vList = []
					vTotal = []
					vBatch_fastq=n
					print "Fastq batch #"+str(loops_fastq+1)+" written"
					continue_loop=False
				else:
					vList = []
					vTotal = []
					continue_loop=False
		file.close()	
				




		# Extract
		vList = []
		vFilePlace = 0		
		vCounter = loops_fixrank

		print "Finding valid samples"
		#global vGoodReads
		vGoodReads=[]
		while vCounter>=0: #Loops through all of the split fix rank files and records all good reads
			with open('extFixrankTemp_'+str(batchsize)+"_"+str(vCounter)) as vFile:
				for line in vFile.readlines():
					if (genus == "" or genus in line) and (family =="" or family in line) and (order =="" or order in line):
			 	 		if species in line and float(line.split("\t",26)[25]) >= float(threshold): #checks threshold and if the sample is above or at the score
			 	 				vGoodReads.append(line.split("|",1)[0]) 		   #it will add the sample to the good reads file	
		 		vCounter-=1
		 	vFile.close()
		
		vFile.close()
		vContinue=False					
		vCounter=loops_fastq
		first_variable=""
		check_first_varable=True

		print str(len(vGoodReads))+" valid samples found"


		print "Extracting the samples. This may take a while"

		vList=[]
		temp=""
		while vCounter>=0:	#Loops through all of the fastq files and writes all of the good reads
		 	temp='extFastqTemp_'+str(batchsize)+"_"+str(vCounter)
		 	vList.append(temp)
			vCounter-=1
		
		pool=multiprocessing.Pool(procs)
		func =partial(extractFastq, vGoodReads, batchsize, species)
		pool.map(func, vList)
		

		print "Samples extracted"
		print "Cleaning up"
		fRemoveFiles(loops_fastq, loops_fixrank, batchsize)
		print "Done"
		print("--- %s seconds ---" % (time.time() - start_time))

