#!/usr/bin/python
# -*- coding: iso-8859-1 -*-
############################################################
## To have SNP filtered on the coverage and Biallelic SNP ##
## Author: Julie Lê-Hoang																	##
## Contact : lehoang.julie@hotmail.fr                     ##
## Date: 11/12/2020																		    ##	
############################################################
import os
import re
import sys

###############
## Read file ##
###############
myFrequency = sys.argv[1]
myFrequencySort = sys.argv[2]

# Variables
nbNA = 0
sortFrequencyNAFiltered = open(myFrequencySort,"w")

## Program 2 ## Removed SNP with NA in each samples
handlerFreq = open(myFrequency,"U")
myFreq = handlerFreq.readlines()

for freq in myFreq:

	if "CHROM" in freq:
		sortFrequencyNAFiltered.write(freq)

	else:
		split = re.split("\t", freq.replace("\n",""))

		for i in range(9, len(split), 1):
			if "NA" in split[i]:
				nbNA = nbNA + 1

		if nbNA == 0:
			sortFrequencyNAFiltered.write(freq)
		nbNA = 0
			
#################
## Close files ##
#################
sortFrequencyNAFiltered.close()

