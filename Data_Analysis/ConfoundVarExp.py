#!/usr/bin/python
#Erawijantari-Yamada Laboratory-Titech-Jan,7 2019
#Plot the variance explain

import sys
import matplotlib as plotlt
plotlt.use("Agg")
import pylab
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import argparse
import re

#Parameters to input
def read_params(args):
	parser = argparse.ArgumentParser(description='DF-modifier')
	parser.add_argument('input_file1', metavar='INPUT_FILE1', type=str, help="ExplainVar by status")
	parser.add_argument('input_file2', metavar='INPUT_FILE1', type=str, help="ExplainVar by metadata")
	parser.add_argument('out_file', metavar='OUTPUT_FILE', type=str, help="the output file")
	parser.add_argument('-r',dest="Rel_abundance", metavar='str', type=str, default='',
	help="LEfSe with LDA score")
	parser.add_argument('-s',dest="Stats_results", metavar='str', type=str, default='',
	help="Kruskall wallis results")
	parser.add_argument('-a',dest="Absolute", metavar='str', type=str, default='',
	help="absolut value")
	args = parser.parse_args() 
	params = vars(args)
	return params 


#string matching in list ignoring the special characters
def string_com(char,list1,list_N):
	#remove the character from list
	char_N = re.sub('[^a-zA-Z0-9]+', '', char)
	if char_N in list_N:
		ind=list_N.index(char_N)
		element_in = list1[ind]
	else:
		element_in = 'None'
	return element_in
	 

#Explain variance data reader
def exp_var(inpt1, inpt2,absolute=None):
	ExpVar_dict1 	= {}
	ExpVar_dict2 	= {}
	inpt1_df 		= pd.read_table(open(inpt1),index_col=0)
	inpt2_df 		= pd.read_table(open(inpt2),index_col=0)
	index 			= inpt1_df["Feature"].tolist()
	absolute 		= str(params['Absolute'])
	for x in index:
		stats_var	= inpt1_df.loc[inpt1_df['Feature'] == x, 'Coefficient'].item()
		meta_var 	= inpt2_df.loc[inpt2_df['Feature'] == x, 'Coefficient'].item()
		#ExpVar_dict[x]=[stats_var,meta_var]
		if absolute =='Yes':
			ExpVar_dict1[x]=abs(stats_var)
			ExpVar_dict2[x]=abs(meta_var)
		else:
			ExpVar_dict1[x]=(stats_var)
			ExpVar_dict2[x]=(meta_var)
	return ExpVar_dict1,ExpVar_dict2,index

#get the p-value
def pVal_get(stats, index_check,list_N):
	pval_dict = {}
	pval_df 		= pd.read_table(open(stats),index_col=0)
	index_l 		= list(pval_df.index.values)
	for x1 in index_l:
		if "md_M" in x1:
			x = x1[3:]
		else:
			x=x1
		new_x=string_com(x,index_check,list_N)
		if new_x != 'None':
			FDR 		= pval_df.loc[x1,"FDR"].item()
			if FDR < 0.1:
				colour = 'q value <0.1'
			else:
				colour = 'q value >0.1'
			pval_dict[new_x]= colour
	return pval_dict

#scale relab
def relab_scale(relab,index_check,list_N):
	#get the relab average and turn into list
	relab_df 		= pd.read_table(open(relab),index_col=0)
	index_l 		= list(relab_df.index.values)
	relab_dict 		={}
	for i in index_l:
		new_x=string_com(i,index_check,list_N)
		if new_x != 'None':
			mean_rel	= np.mean(relab_df.loc[i].values)
			if mean_rel < 0.000001:
				scale = 10
			elif 0.000001 <= mean_rel < 0.00001:
				scale = 50
			elif 0.00001 <= mean_rel < 0.001:
				scale = 100
			elif mean_rel> 0.001:
				scale = 200
			relab_dict[new_x]	= scale
	#df_relabdict 	= pd.DataFrame.from_dict(relab_dict, orient='index',columns=['mean','scale'])	
	return relab_dict


#combine the dictionary into one dataframe
def combine_df(ExpVar_dict1,ExpVar_dict2,pval_dict,scale_dict):
	#combine multiple dict into pandas dataframe
	df_com = pd.DataFrame([ExpVar_dict1,ExpVar_dict2,pval_dict,scale_dict]).T
	#df_com.columns = ['d{}'.format(i) for i, col in enumerate(df_com, 1), column=['Expvar_stats', 'Expvar_metadata', 'pval','scale']]
	df_com.columns = ['Status', 'Metadata', 'pval','scale']
	return df_com


#plot the scatter plot
def scatter_plot(df,outtitle,absolute=None):
	plotlt.rcParams['pdf.fonttype'] = 42
	plotlt.rcParams['ps.fonttype'] = 42
	#hue --> give colour for FDR<0.1
	#size scale based on relab
	colour =['r','black']
	sns.set(style="white")
	sns.plotting_context(font_scale=0.5)
	size_col = df['scale'].tolist()
	ax=sns.scatterplot(x='Status', y='Metadata',data=df,hue='pval',palette=colour)#, sizes='scale', size='scale')#, data=df
	ax.legend(loc='center left', bbox_to_anchor=(1, 1),fontsize=7,markerscale=0.3)
	ax.plot([0, 1], [0, 1], transform=ax.transAxes, c=".3", linewidth=0.2,zorder=1)
	# control x and y limits
	#axes = ax.axes
	absolute 		= str(params['Absolute'])
	concat_list= (df['Status'].tolist())+(df['Metadata'].tolist())
	max_val = max(concat_list)+ 0.001
	min_val = min(concat_list) - 0.001
	ax.set_xlim(min_val,max_val)
	ax.set_ylim(min_val,max_val)
	if absolute == 'Yes':
		#ax.set_xlim(-0.0001, 0.01)
		#ax.set_ylim(-0.0001, 0.01)
		ax.set_xlabel('Status',fontsize=8)
		ax.set_ylabel('Metadata',fontsize=8)
	else:
		#ax.set_xlim(-0.1, 0.05)
		#ax.set_ylim(-0.1, 0.05)
		#ax.set_xlim(min_val,max_val)
		#ax.set_ylim(min_val,max_val)
		ax.set_xlabel('Crude coefficient',fontsize=8)
		ax.set_ylabel('Adjusted coefficient',fontsize=8)
	
	plt.savefig(outtitle, dpi=400)

def boxplot_stats(dfin,title):
	plotlt.rcParams['pdf.fonttype'] = 42
	plotlt.rcParams['ps.fonttype'] = 42
	#hue --> give colour for FDR<0.1
	#size scale based on relab
	colour =['r','black']
	sns.set(style="white")
	sns.plotting_context(font_scale=0.5)
	


if __name__ == '__main__':
	params 	= read_params(sys.argv)
	input1 	= params['input_file1']
	input2 	= params['input_file2']
	title	= str(params['out_file'])
	Explain_variance1,Explain_variance2,index_list = exp_var(input1,input2)
	index_mod = [re.sub('[^a-zA-Z0-9]+', '', _) for _ in index_list]
	Stats 	= params['Stats_results']
	pVal 	= pVal_get(Stats, index_list,index_mod)
	reldf 	= params['Rel_abundance']
	scale	= relab_scale(reldf,index_list,index_mod)
	df_combine = combine_df(Explain_variance1, Explain_variance2, pVal, scale)
	scatter_plot(df_combine,title)
	


	


