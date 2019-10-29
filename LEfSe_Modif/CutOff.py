#!/usr/bin/python
#May 26 2018
#this code is created to remove the rows that the have mean less than thresshold
#generating the file for LEfSe
#ussage : time python cutt_off.py rel_abundance_file cut_off_appearance cut_off_means <str to remove from columns>

import sys
import pandas as pd 
import numpy as np

argvs = sys.argv
path=argvs[1]
c=int(argvs[2])
d=float(argvs[3])
RefFile=open(argvs[4])
columns =int(argvs[5])


def Refname_Gen(RefFile,columns):
	refname = {}
	for tab in RefFile:
		tab1	= tab.split("\t")[0]
		#tab2	= (tab.split("\t")[3]).strip()
		##new data
		key		= tab1+'.'+tab.split("\t")[1]
		##old data
		#key 	= tab.split("\t")[1]+'.'+tab1
		value	= tab.split('\t')[columns]
		if key not in refname:
			refname[key]=value
	return refname
#pandas dataframe 

refname=Refname_Gen(RefFile,columns)

df_table1 = pd.read_table(path,index_col=0)

if len(argvs)== 5:
	drop = (argvs[4])
	df_table = df_table1.drop([drop])
else:
	df_table=df_table1

del df_table1
list_samples = df_table.columns.values.tolist()

#remove by appearance
c1 = int(len(list_samples)*c//100)
row = []
[row.append(a) for a,b in (df_table.astype(bool).sum(axis=1)).iteritems() if b<=c1]
df_table_drop=df_table.drop(row)

#add the mean and standard deviation
mean_list =[]
std_list =[]
list_KO=df_table_drop.index.values
for KO in list_KO:
	mean_KO=(df_table_drop.loc[KO].values).mean()
	std_KO=(df_table_drop.loc[KO].values).std()
	mean_list.append(mean_KO)
	std_list.append(std_KO)

df=df_table_drop.describe()
se1=pd.Series(mean_list)
se2=pd.Series(std_list)
df_table_drop['mean']=se1.values
df_table_drop['std']=se2.values

#remove the row that have the mean value less than define in c
df_remove1=df_table_drop[df_table_drop['mean'] >= d]
df_remove=df_remove1.drop(['mean','std'],axis=1)
del df_remove1

columns_name 	=(list(df_remove.columns.values))
add_row 		= []

for columns in columns_name:
	add 		= refname[columns]
	add_row.append(add.strip())
try:
	del df_remove.index.name
except:
	print "No index name"

title =str(argvs[1].split(".tsv")[0])
outtitle=title+'.'+str(c)+'.'+str(d)+".up.tsv"
outtitle2=title+'.'+str(c)+'.'+str(d)+".upLEfSe.txt"

df_remove.to_csv(outtitle,header=True, index=True, sep='\t')
df_remove.columns = pd.MultiIndex.from_tuples(zip(add_row, columns_name),names=['Class','Subject_ID'])
df_remove.to_csv(outtitle2,header=True, index=True, sep='\t')

