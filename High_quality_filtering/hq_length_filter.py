##!/usr/bin/python
#10 November 2016
#this code for filter out sequence that have more than 50 b length
 
import sys
import gzip

argvs = sys.argv 
title = str(argvs[1].split(".trim.fq")[0])
in1	  = argvs[1]

if in1.endswith('.gz'): 
	outSuffix='.fq.gz'
else:
	outSuffix='.fq'

#print out the ussage command if the user type the wrong input
if len(argvs) != 2:
	print ("ussage : python hq_length.filter.py sequence_after_cutadapt")
	sys.exit(0)

#defining functions
#open file function
def myopen(infile, mode='r'):
	if infile.endswith('.gz'):
		return gzip.open(infile,mode=mode)
	else:
		return open(infile,mode=mode)

with myopen(title+'.length.filter'+outSuffix,'w') as f2:
	with myopen(in1) as f1:

		while True:
			
			line1=f1.readline()
			if not line1: break
			line2=f1.readline()
			line3=f1.readline()	
			line4=f1.readline()	
			
				
			#FOR WRITING UP THE FILE
			#extract the line that contain more than 50b sequence
			if len(line2)>50 :
				f2.write(line1)
				f2.write(line2)
				f2.write(line3)
				f2.write(line4)

		



