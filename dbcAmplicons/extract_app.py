#extract_app

import sys, os, time, gzip
from dbcAmplicons import OneReadIlluminaRun

def fFileLength(file):
	#Check that files are the same length
	with open(file,"r") as file:
		for i, l in enumerate(file):
			pass
		return i + 1

def fFileWrite(vFileName,vToWrite):
	#Writes vDict to vFileName file.
	if os.path.isfile(vFileName):
		open(vFileName,'w')
	vTempFile = open(vFileName,'a')
	try:
		#Test if dictionary
		if vToWrite.keys() == True:
			for i in vToWrite:
				for j in vToWrite[i]:
					vTempFile.write(j)
	except:
		#if fails, assume its a list
		for i in vToWrite:
			vTempFile.write(i)	

class extractApp():

	def __init__(self):
		self.verbose = False

	def start(self, taxon, fixrank, fastq, threshold, output, batchsize):

		#Format input files:
		print("---\nTaxon target is: " + taxon)
		vFixrankFile = open(fixrank, "r")
		if fastq[-3:] == '.gz':
			vFastqFile = gzip.open(fastq, "r")
		else:
			vFastqFile = open(fastq, "r")

		#Run through the batches
		block = False
		if block == False:

			vFixBatchList = []
			vBatch,n = 0,0
			with vFixrankFile as file:
				for classification in file.read().split('\n'):

					if n < batchsize:
						vFixBatchList.append(classification+'\n')
						n+=1
					
					if n == batchsize:
						fFileWrite('extFixrankTemp_'+str(vBatch), vFixBatchList)
						vBatch+=1
						vFixBatchList = []
						n = 0

				if n != batchsize:
					if len(vFixBatchList) != 0:
						fFileWrite('extFixrankTemp_'+str(vBatch), vFixBatchList)
						vBatch +=1

		#Generate fastq batches
		block = True
		if block == False:

			n = 0
			k = 0
			j = 0
			vBatch = 0
			vList = []
			vDict = {}
			with vFastqFile as file:
				for line in file.readlines():
					if n < batchsize:
						vList.append(line)
						k+=1
						if k % 4 == 0:
							vDict[j] = vList
							vList = []
							j+=1
							n+=1
					if n == batchsize:
						fFileWrite('extFastqTemp_'+str(vBatch), vDict)
						print(str(vBatch) + ' complete.')
						n = 0
						j = 0
						vBatch +=1
						vList = []
						vDict = {}
				if len(vDict) != 0:
					fFileWrite('extFastqTemp_'+str(vBatch), vDict)

		#Compares file lengths
		block = True
		if block == False:

			for batch in range(vBatch):
				with open('extFixrankTemp__'+str(batch),"r") as file:
					vFixrankLength = (len(file.readlines()))
					file.close()
				with open('extFastqTemp_'+str(vBatch),"r") as file:
					vFastqLength = (len(file.readlines())/4)
					file.close()
				if vFixrankLength == vFastqLength:
					print(str(batch)+' '+str(vFixrankLength)+':'+str(vFastqLength)+' Files align. Proceeding with extraction.')
				else:
					print(str(batch)+' '+str(vFixrankLength)+':'+str(vFastqLength)+' [ERROR]: Files do not align.')
		
		# Extract
		vList = []
		vFilePlace = 0
		vCounter = 0
		for batch in range(vBatch):
			print('Batch '+str(batch)+' begun')
			with open('extFixrankTemp_'+str(batch)) as vFile:
			 	for line in vFile.readlines():
			 		vFastqPlace = (vFilePlace*4)
	 	 			if line.isspace() == True:
	 	 				continue
	 	 			elif line[0:6] == '[FAIL]':
	 	 				vFilePlace +=1
	 	 			else:
	 	 				if line.split()[22] == taxon:
	 	 					if float(line.split("\t",26)[25]) >= float(threshold):
	 	 						#Extract Read from fastq
	 	 						loop = True
	 	 						with open(fastq, "r") as file:
	 	 							for i, read in enumerate(file):
	 	 								if vCounter > 0:
	 	 									vList.append(read)
	 	 									vCounter -= 1
	 	 									loop = False

	 	 								if loop == True:
		 	 								if i == vFastqPlace:
		 	 									if read.split()[0].split('|')[0] == ('@'+line.split()[0].split('|')[0]):
		 	 										vList.append(read)
		 	 										vCounter = 3
	 	 							vFilePlace+=1

	 	 			 	else:
			 		 		vFilePlace+=1
			if len(vList) > 0:
				fFileWrite('extOut_'+str(batch),vList)
				vList = []
			print('Batch '+str(batch)+' complete')
