#!/usr/bin/python
#24 May 2018
#perform matrix generations to combine results from each sample
#usage : python CombineMatrix.py path_to_the_results file_reference_convert_seq_to_stage
#e.g python CombineMatrix.py ~/gut_microbiome_analysis/stomach_removed/results/ ~/gut_microbiome_analysis/stomach_removed_analysis/list_data/160706_seq-state_convertion.txt pattern


import sys
import pandas as pd 
import glob
import os
from datetime import date

argvs 		= sys.argv
path 		=argvs[1]
dirpath 	= os.getcwd()
fileref 	=open(argvs[2])
file 		= str(argvs[3])
pattern 	='*/'+file
refname 	={}
today  		= str(date.today())
if '*' in file:
	file1=file.split('*')[-1]
elif '/' in file:
	file1 = file.split('/')
else:
	file1=file
title_list 	= os.path.join(dirpath,today+'_'+file1+'_analysis.txt')
title_df 	= os.path.join(dirpath,today+'_'+file1+'.tsv')


list_file=glob.glob(os.path.join(path,pattern)) #KO-USCG.txt
with open(title_list,'w') as writeup:
	for x in list_file:
		writeup.write(x+"\n")

for tab in fileref:
	tab1	= tab.split("\t")[0]
	tab2	= (tab.split("\t")[3]).strip()
	#change the ID into IDnumber.Stage
	value	=tab1+'.'+tab.split("\t")[1]
	
	if tab1 not in refname:
		refname[tab2]=value


name_g=""
df_list=[]
for file1 in list_file:
	name_g3=file1.split("/")[-2]
	'''
	#remove the suffix (April, 25)
	if "." in name_g31:
		name_g3= name_g31.split('.')[0]
	else:
		name_g3=name_g31
	'''
	if name_g3 in refname:
		name_g=refname[name_g3]
	if file1.endswith('.gz'):
		df1 =pd.read_table(file1, compression='gzip', index_col=0, header=None, names=["",name_g], usecols=["",name_g])
	else:
		df1 =pd.read_table(file1,index_col=0, header=None, names=["",name_g], usecols=["",name_g])
	df2=df1.loc[(df1!=name_g3).all(1)]
	df_list.append(df2)
	big_df=df_list[0]


for n in range(len(df_list)):
	if n!=(len(df_list)-1):
		big_df=big_df.join(df_list[n+1], how='outer')
del df_list

big_df_f=big_df.fillna(0)
del big_df

#remove Nan columns
row = []
[row.append(a) for a,b in (big_df_f.astype(bool).sum(axis=1)).iteritems() if b<=1]
#remove row with zero value
df_table_drop=big_df_f.drop(row)
df_table_drop1=df_table_drop.loc[(df_table_drop!='0').any(1)]
df_table_drop2=df_table_drop1.loc[(df_table_drop1!='motus.processing.1').all(1)]
df_table_drop2.to_csv((title_df), header=True, index=True, sep='\t')




####example###
#time python ~/gut_microbiome_analysis/script/matrix_generator/CombineMatrix.py ~/gut_microbiome_analysis/gastrectomy_early/Seq_data/ ../../list_data/2018-07-18_sampleID_seqID.txt *.diamond.aln.id40.sc70.cov80.parsed.ko.relative.tsv
