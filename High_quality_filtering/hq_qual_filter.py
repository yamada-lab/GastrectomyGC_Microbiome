##!/usr/bin/python
#10 November 2016
#this code sequence that have more than 25 average quality score

import sys
import gzip

argvs = sys.argv 
title = str(argvs[1].split(".length.filter")[0])
in1 = argvs[1]

if in1.endswith('.gz'): 
	outSuffix='.fq.gz'
else:
	outSuffix='.fq'

#print out the ussage command if the user type the wrong input
if len(argvs) != 3:
	print ("ussage : python hq_qual_filter.py file-after-length_filter ascii-table")
	sys.exit(0)

ascii_dict={}
asciifile = open(sys.argv[2],"r")
for line in asciifile:
    ascii_dict[line.split("\t")[0]] = int(line.split("\t")[1]) - 33

#defining functions
#open file function
def myopen(infile, mode='r'):
	if infile.endswith('.gz'):
		return gzip.open(infile,mode=mode)
	else:
		return open(infile,mode=mode)

with myopen(title+'.qual.filter'+outSuffix,'w') as f2:
	with myopen(in1) as f1:

		while True:
			
			line1=f1.readline()
			if not line1: break
			line2=f1.readline()
			line3=f1.readline()	
			line4=f1.readline()	
			qc_sum=float(0)
			qc_num=float(0)
			qc_val=float(0)
			for Q in list(line4[:-1]):
			#	if Q in ascii_dict:
					#qc_sum += float(ascii_dict[Q])
				qc_sum += float(ord(Q)-33)
				qc_num += float(1)
				qc_val =float(qc_sum/qc_num)
				
			#FOR WRITING UP THE FILE
			#extract the line that contain more than or equal 25 quality 
			if qc_val>=25.00 :
				f2.write(line1)
				f2.write(line2)
				f2.write(line3)
				f2.write(line4)

#last edited 11 May 2016 --> using ASCII table rather than ord function


