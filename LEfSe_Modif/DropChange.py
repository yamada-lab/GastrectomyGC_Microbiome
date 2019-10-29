#!/usr/bin/python
#Erawijantari-Yamada Laboratory-Titech-Jul 25, 2018
#create file for pairwise and change the status
#generating the file for LEfSe; exclude; or selecting based on the list

import sys
import pandas as pd 
import numpy as np
import os
import argparse

#Parameters to input
def read_params(args):
	parser = argparse.ArgumentParser(description='DropID-DF')
	parser.add_argument('input_file', metavar='INPUT_FILE', type=str, help="the input file")
	parser.add_argument('out_file', metavar='OUTPUT_FILE', type=str, help="the input file")
	parser.add_argument('RefFile', metavar='ReferenceFile', type=str, help="the reference file")
	#parser.add_argument('-r',dest="ReferenceFile", metavar='str', type=str, default='',
	#help="prepare file with ID to retrive or the stage to drop)")
	parser.add_argument('-d',dest="DropList", metavar='str', type=str, default='',
	help="prepare file with ID to retrive or the stage to drop)")
	parser.add_argument('-g',dest="Get list", metavar='str', type=str, default='',
	help="prepare file with ID to retrive or the stage to drop)")
	parser.add_argument('-c',dest="ColumnsOnRefFile", metavar='str', type=int, default=1,
	help="define the drop ID)")
	args = parser.parse_args() 
	params = vars(args)
	return params 

def read_dataframe(path_matrix):
	matrix_raw 			= path_matrix
	matrix_bef 			= pd.read_table(matrix_raw,index_col=0)
	feature_list		= matrix_bef.index.values

	if 'C00158_Citrate' in feature_list:
		matrix1 		= matrix_bef.drop(['C00158_Citrate'],axis=0, inplace=False)
	else:
		matrix1 		= matrix_bef

	if '-1' in feature_list:
		matrix1 		= matrix1.drop(['-1'],axis=0, inplace=False)
	else:
		matrix1 		= matrix_bef

	#removing if the 'mean, std' were in dataframe
	try:
		matrix2  		= matrix1.drop(['mean','std'],axis=1,inplace=False)
	except:
		matrix2 		= matrix1
		matrix 	 		= matrix2.fillna(0)
		#matrix 			= matrix2.apply(lambda x:x-x.mean())
	return matrix

def IDdropStage(df,dictOfElements,valueToFind):
	listOfKeys = list()
	listOfItems = dictOfElements.items()
	df_columns =list(df.columns.values) 
	for item  in listOfItems:
		if item[0] in df_columns:  
			if item[1]== valueToFind:
				listOfKeys.append(item[0])
	df_remove = df.drop(listOfKeys,axis=1)
	return df_remove

def IDdropList(df,catdropL):
	listID		=[x.strip() for x in catdropL]
	df_remove = df.drop(listID,axis=1)
	return df_remove

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

def IDgetList(df,dictOfElements,catdropL):
	listID		=[x.strip() for x in catdropL]
	listget 	=[]
	for a,b in dictOfElements.items():
		##new data
		tab1 	= (a.split(".")[0])
		##old data
		#tab1 = (a.split(".")[1]).strip()
		if tab1 in listID:
			key		= a
			listget.append(key)

	df_remove = df[listget]
	return df_remove

def RelAbCalc(df):
	#calculate the sum abundance for each sample
 	columnsDf = 	(list(df_remove.columns.values)) 

	#apply relative abudance calc


	#return pandas dataframe

if __name__ == '__main__':
	params = read_params(sys.argv)
	df_table= read_dataframe(params['input_file'])
	RefFile = open(params['RefFile'])
	drop = params['DropList']
	title 	= str(params['out_file'])
	columns = params['ColumnsOnRefFile']	
	refname=Refname_Gen(RefFile,columns)
	
	if drop=='No':
		df_remove 	= df_table
		outtitle 	= title+'.'+"tsv"
		outtitle2 	= title+'.'+"upLEfSe.txt"
	else:
		try :
			dropfile = open(drop)
			df_remove= IDgetList(df_table,refname,dropfile)
			outtitle 	= title+'.'+'Specify.'+"tsv"
			outtitle2 	= title+'.'+'Specify.'+"upLEfSe.txt"
			
		except:
			df_remove=IDdropStage(df_table,refname,drop)
			outtitle 	= title+'.'+'state.'+str(columns)+'drop.'+drop+".tsv"
			outtitle2 	= title+'.'+'state.'+str(columns)+'drop.'+drop+"upLEfSe.txt"

	
	

	columns_name 	=(list(df_remove.columns.values))
	add_row 		= []

	for columns in columns_name:
		add 		= refname[columns]
		add_row.append(add.strip())
	try:
		del df_remove.index.name
	except:
		print "No index name"
	df_remove.to_csv(outtitle,header=True, index=True, sep='\t')
	df_remove.columns = pd.MultiIndex.from_tuples(zip(add_row, columns_name),names=['Class','Subject_ID'])
	print outtitle2
	df_remove.to_csv(outtitle2,header=True, index=True, sep='\t')


###example###
#python ~/gut_microbiome_analysis/script/matrix_generator/DropChange.py 2018-07-25_annotated.mOTU.abundances.7.0.0001.up.tsv ../../list_data/2018-07-18_sampleID_seqID.txt 1 EMR_ESD	





