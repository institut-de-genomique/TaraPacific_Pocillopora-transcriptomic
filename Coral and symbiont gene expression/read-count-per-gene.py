#!/usr/bin/python
# -*- coding: iso-8859-1 -*-
########################################################################################################################
## To have mean coverage, mean identity and count reads, that mapped on transcriptome reference, order by transcripts ##
## From: Sam-Carac results; 					   	                                              ##
## Author: Julie Lê-Hoang		 		                                                              ##
## Contact : lehoang.julie@hotmail.fr                                                                                 ##
## Date: 20/03/2018			                                                             		      ##	
########################################################################################################################
import os
import re
import sys

###############
## Variables ##
###############
# General
myDict = dict()
split = list()
tempID = ''
l = -1.0
tempStat = ''
tempID = ''
myFile = ''
myIsland = ''
#myDictStat = dict()

# Coverage
refStart = 0
refEnd = 0
length = 0
mySeq = list()
myTempSeq = ''
previousID = ''
count1 = 0
count0 = 0
cover = 0.0

# Identity
identity = list()
matches = 0.0
totErrors = 0.0
calc = 0.0
iden = 0.0

# CountID
countID = 1

###############
## Read file ##
###############
myFile = sys.argv[1]
mySize = sys.argv[2]
mySort = sys.argv[3]
split = re.split("_",myFile)
myIsland = split[0]

##############################
## Selection ID transcripts ##
##############################
handlerTranscripts = open(mySize,"rU")
myTrans = handlerTranscripts.readlines()

for trans in myTrans:
	split = re.split("\t", trans)
	myDict[split[0]] = str(split[1].replace("\n",""))

##############
## Out file ##
##############
sortFile = open(mySort,"w")
sortFile.write("QuerySeqID\tCoverage\tIdentity\tMapCount\tLength\tOrigin\n")

#######################
## Coverage Identity ##
## and Count reads   ##
#######################
handlerStat = open(myFile, "rU")
myStat = handlerStat.readlines()

for stat in myStat:
	split = re.split("\t",stat)
	tempID = split[1]
	
	if tempID != "chr_name":
		if previousID == '':
			# Coverage
			length = myDict[tempID]
			refStart = split[3]
			refEnd = split[4]
			myTempSeq = '0' * (int(refStart)) +  '1' * (int(refEnd)-int(refStart)) + '0' * (int(length)-int(refEnd))
			mySeq = list(myTempSeq)
			previousID = tempID
			
			# Identity
			matches = split[10]
			totErrors = split[14]
			calc = (float(matches) / (float(matches) + float(totErrors))) * 100
			identity.append(calc)
			
		else:
			if tempID == previousID:
				# Coverage
				length = myDict[tempID]
				refStart = int(split[3])
				refEnd = int(split[4])

				for i in range(int(refStart-1), int(refEnd-1), 1):
					mySeq[i] = '1'

				previousID = tempID
				myNewSeq = ''

				# Identity
				matches = split[10]
				totErrors = split[14]
				calc = (float(matches) / (float(matches) + float(totErrors))) * 100
				identity.append(calc)

				# CountID
				countID += 1

			else:
				#Coverage
				for i in range(0, int(length), 1):
					if mySeq[i] == '1':
						count1 = count1 + 1.0
					else:
						count0 = count0 + 1.0
				cover = (count1 / int(length)) * 100
				
				# Identity
				for i in range(0, len(identity), 1):
					iden = iden + identity[i]
				calc = float(iden) / float(len(identity))

				sortFile.write(previousID+"\t"+str(cover)+"\t"+str(calc)+"\t"+str(countID)+"\t"+str(length)+"\t"+myIsland+"\n")				
	
				#myDictStat[previousID] = 0
				# Pass
				# Coverage
				count1 = 0.0
				count0 = 0.0

				previousID = tempID

				refStart = split[3]
				refEnd = split[4]
				readAlign = split[8]
				length = myDict[tempID]
				

				myTempSeq = '0' * (int(refStart)) +  '1' * (int(refEnd)-int(refStart)) + '0' * (int(length)-int(refEnd))
				mySeq = list(myTempSeq)
				# Identity
				iden = 0.0
				identity = list()
				matches = split[10]
				totErrors = split[14]
				calc = (float(matches) / (float(matches) + float(totErrors))) * 100
				identity.append(calc)

				# CountID
				countID = 1

# For the last data
# Coverage 
for i in range(0, int(length), 1):
	if mySeq[i] == '1':
		count1 = count1 + 1.0
	else:
		count0 = count0 + 1.0
cover = (count1 / int(length)) * 100

# Identity
for i in range(0, len(identity), 1):
	iden = iden + identity[i]
calc = float(iden) / float(len(identity))
sortFile.write(previousID+"\t"+str(cover)+"\t"+str(calc)+"\t"+str(countID)+"\t"+str(length)+"\t"+myIsland+"\n")
#myDictStat[previousID] = 0

######################
## Transcripts that ##
## are note present ##
######################
#keysDict = myDict.keys()
#Search = ''

#for i in range(0, len(keysDict), 1):
#	Search = myDictStat.has_key(keysDict[i])
	#tempID = myDict[keysDict[i]]
#	length = myDict[keysDict[i]]
#	if Search == False:
#		sortFile.write(keysDict[i]+"\t0.0\t0.0\t0.0\t"+str(length)+"\t"+myIsland+"\n")

#################
## Close files ##
#################
sortFile.close()
handlerStat.close()
handlerTranscripts.close()


