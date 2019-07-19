##!/usr/bin/python
#9 November 2016
#this code for N-contained sequence filtering
#high quality sequence : not containt N reads


import sys
import gzip

argvs = sys.argv 
title = str(argvs[1].split(".fq")[0])
in1 = argvs[1]

if in1.endswith('.gz'): 
	outSuffix='.fq.gz'
else:
	outSuffix='.fq'

#print out the ussage command if the user type the wrong input
if len(argvs) != 2:
	print ("ussage : python hq_N.filter.py fastq_raw_reads_file")
	sys.exit(0)	

#defining functions
#open file function
def myopen(infile, mode='r'):
	if infile.endswith('.gz'):
		return gzip.open(infile,mode=mode)
	else:
		return open(infile,mode=mode)



with myopen(title+'.N.filter'+outSuffix,'w') as f2:
	with myopen(in1) as f1:
		while True:
			
			name=f1.readline()
			if not name: break
			seq=f1.readline()
			name2=f1.readline()	
			qual=f1.readline()	

			if "N" not in seq:
				f2.write(name)
				f2.write(seq)
				f2.write(name2)
				f2.write(qual)



#
		




