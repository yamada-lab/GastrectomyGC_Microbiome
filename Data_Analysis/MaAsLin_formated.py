#!/usr/bin/python
#Erawijantari-Yamada Laboratory-Titech-Jul 25, 2018
#create file for pairwise and change the status
#generating the file MaAsLin as a input for merge and transpose file

import sys
import pandas as pd 
import numpy as np
import os
import argparse


#Parameters to input
def read_params(args):
	parser = argparse.ArgumentParser(description='MaAsLinfile_Preparation')
	parser.add_argument('input_file', metavar='INPUT_FILE',type=str, help="the input file")
	parser.add_argument('metadata_file', metavar='METADATA_FILE',type=str, help="the input file")
	parser.add_argument('out_file', metavar='OUTPUT_FILE',type=str, help="the input file")
	parser.add_argument('RefFile', metavar='ReferenceFile', type=str, help="the reference file")
	parser.add_argument('-m',dest="MetadataList", metavar='METADATA_LIST,str', type=str, default='all',
	help="the list of metadata that want to be retrieved)")
	parser.add_argument('-f',dest="Feature", metavar='str', type=str, default='Feature',
	help="Determine the microbial feature to test)")
	args = parser.parse_args() 
	params = vars(args)
	return params 

#generating the formated features for MaAsLin
def df_parser(df):
	listID 	= (list(df.columns.values))
	Feature =[1]*len(listID)
	df_add 	= pd.DataFrame({str(params['Feature']): Feature}, index=listID)
	df_append = (df_add.T).append(df)
	return df_append
	

#generate dictionary
def Refname_Gen(RefFile):
	refname = {}
	for tab in RefFile:
		key	= tab.split("\t")[0]
		value= key+'.'+tab.split("\t")[1]
		if key not in refname:
			refname[key]=value
	return refname

#generate the metadata dataframe
def metadata_parser(metadata_df,ListMeta,Ref_dict):
	try:
		retrive_meta 	= open(ListMeta)
		metadata_l		= [x.strip() for x in list_meta]
		df_retrieve 	= metadata_df[metadata_l]
		df_T 			= df_retrieve
	except:
		df_T 			= metadata_df

	list_ID 	= (list(df_T.index.values))
	RefIDnew 	= []

	for x in list_ID:
		if str(x)in Ref_dict:
			RefIDnew.append(Ref_dict[str(x)])
	df_T.index 		= RefIDnew
	return df_T

def ReadConfGen(metadata_df,feat_df,outfile):
	metadata 	= (list(metadata_df.columns.values))
	metadata.sort()
	features	= (list(feat_df.index.values))
	features.sort()
	last_elementM= metadata[-1]
	first_elementF = features[0]
	with open(outfile,'w') as f1:
		f1.write('Matrix: Metadata'+'\n')
		f1.write('Read_PCL_Rows: -'+str(last_elementM)+'\n')
		f1.write(''+'\n')
		f1.write('Matrix: Abundance'+'\n')
		f1.write('Read_PCL_Rows: '+str(first_elementF)+'-')



if __name__ == '__main__':
	params 			= read_params(sys.argv)
	
	#metadata formated generation
	metadata_df 	= pd.read_table(open(params['metadata_file']),index_col=0)
	feature_df 		= pd.read_table(open(params['input_file']),index_col=0)
	retrive_meta 	= params['MetadataList']
	Ref_File 		= open (params['RefFile'])
	Ref_Meta 		= Refname_Gen(Ref_File)
	Metadata_New 	= metadata_parser(metadata_df,retrive_meta,Ref_Meta)
	outMeta 		= str(params['out_file'])+'.Metadata.tsv'
	outConfig 		= str(params['out_file'])+'.read.config'
	ReadConfGen(metadata_df, feature_df,outConfig)
	Metadata_New = Metadata_New.fillna('NA')
	Metadata_New.to_csv(outMeta,header=True, index=True, sep='\t')
	
	#example
	#time python ~/gut_microbiome_analysis/script/ReRunOrganized/StastViz/MaAsLIn_formated.py ../trial/181031_RemovalNot/Removed/181031_NCBI.species.abundances.tsv.Specifiy.tsv ../Metadata/181102_Metadata_selectedExcluded.txt Test_metadata_run.txt ../list_data/2018-10-21_sampleID_seqIDNoEMR
