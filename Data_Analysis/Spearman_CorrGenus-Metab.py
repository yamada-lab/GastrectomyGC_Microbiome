#!/usr/bin/python
#11 April 2017
# This code was created to perform the correlation calculation among features
# After the correlation --> create the heatmap and get the p_value
# ussage:
#python spearman_rank_corr.py the_rel_abundance_table1 the_rel_abundance_table2 the_significance-hypothesis_file1 the_significance-hypothesis_file1

import matplotlib as plotlt
plotlt.use("Agg")
import pylab
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from string import letters
from scipy.stats import spearmanr
from statsmodels.sandbox.stats.multicomp import multipletests
import math


#Parameters to input
def read_params(args):
	parser = argparse.ArgumentParser(description='Spearman-correlation analysis')
	parser.add_argument('input_file1', metavar='INPUT_FILE1', type=str, help="the input file")
	parser.add_argument('input_file2', metavar='INPUT_FILE2', type=str, help="the input file")
	parser.add_argument('out_file', metavar='OUTPUT_FILE', type=str, help="the input file")
	parser.add_argument('RefFile1', metavar='ReferenceFile1', type=str, help="the reference file")
	parser.add_argument('RefFile2', metavar='ReferenceFile2', type=str, help="the reference file")
	parser.add_argument('-h1',dest="Hyphothesis1", metavar='str', type=str, default='',
	help="prepare file with ID to retrive or the stage to drop)")
	parser.add_argument('-h2',dest="Hypothesis2", metavar='str', type=str, default='',
	help="prepare file with ID to retrive or the stage to drop)")
	args = parser.parse_args() 
	params = vars(args)
	return params 


def read_dataframe(path_matrix):
	matrix_raw 			= path_matrix
	matrix_bef 			= pd.read_table(matrix_raw,index_col=0)
	feature_list		= matrix_bef.index.values

	if 'C00158_Citrate' in feature_list:
		matrix1 		= matrix_bef.drop(['C00158_Citrate'],axis=0, inplace=False)
		
		try:
			matrix2  		= matrix1.drop(['mean','std'],axis=1,inplace=False)
		except:
			matrix2 		= matrix1
		#matrix 			= matrix2.apply(lambda x:x-x.mean())
		matrix 	 		= matrix2.fillna(0)
	else:
		matrix1 		= matrix_bef
		try:
			matrix2  		= matrix1.drop(['mean','std'],axis=1,inplace=False)
		except:
			matrix2 		= matrix1
		matrix	 		= matrix2.fillna(0)
		#matrix 			= matrix2.apply(lambda x:x-x.mean())
	return matrix



def sig_list(dataframe=None,list_inp=None):
	list_sig 	= []
	list_lable 	= []
	if list_inp==None:
		feature_list = dataframe.index.values
		for i in feature_list:
			p_value 	= (dataframe.loc[i,'corrected_p_value'])
			if p_value < 0.05:
				list_sig.append(i)
			if " " not in i:
				if '_' in i:
					i_list=i.split('_')
					if len(i_list)==2:
						i2 = i_list[1]
					elif  len(i_list)==3:
						i2 = i_list[1]+'_'+i_list[2]
					list_lable.append(i2)
				elif '[' in i:
					a  = i.split()
					i1 = (a[0].split('[')[1]).split(']')[0]
					i2 = i1+' '+a[1]
					list_lable.append(i2)
			else :
				list_lable.append(i)
	elif dataframe==None:
		list_sig = list_inp
		for i in list_sig:
			if " " not in i:
				if '_' in i:
					i_list=i.split('_')
					if len(i_list)==2:
						i2 = i_list[1]
					elif  len(i_list)==3:
						i2 = i_list[1]+'_'+i_list[2]
					list_lable.append(i2)
				elif '[' in i:
					a  = i.split()
					i1 = (a[0].split('[')[1]).split(']')[0]
					i2 = i1+' '+a[1]
					list_lable.append(i2)
			else :
				list_lable.append(i)
	return list_sig,list_lable


def concatinate_df(df1_all,df2_all,stage=None):
	#list_ID2 = list(df2_all.columns.values)
	if stage==None:
		#df1_all =df1_all1[list_ID2]
		df1 	= df1_all
		df2 	= df2_all
	else:
		#df1_all =df1_all1[list_ID2]
		df1=(df1_all[df1_all.columns[df1_all.columns.to_series().str.contains(stage)]])
		df2=(df2_all[df2_all.columns[df2_all.columns.to_series().str.contains(stage)]])
	df_con=pd.concat([df1, df2], axis=0, join='inner')
	df_con.fillna(0.0, inplace=True)
	return df_con


def spearman_rank(dataframe_all,list_cat1,list_cat2,stage=None):

	if stage == None:
		dataframe = dataframe_all
		title 	= 'all'
	else:
		dataframe=(dataframe_all[dataframe_all.columns[dataframe_all.columns.to_series().str.contains(stage)]])
		title 	= stage
	df_T  		= dataframe.transpose()
	feature 	= dataframe.index.values #78

	del dataframe
	list_app 	= [] #78
	list_app.extend(list_cat1)
	list_app.extend(list_cat2)

	matrix_df1a = df_T[list_app]

	#print matrix_df1a.columns #80

	#matrix_df1b = matrix_df1a.loc[list_cat1]
	matrix_df 	= matrix_df1a.as_matrix()
	#del matrix_df1a
	rho,p_value = spearmanr(matrix_df)

	numcols 	= len(p_value[0])
	p_flat 		= p_value.flatten()
	mask 		= np.isfinite(p_flat)
	pval_corrected = np.empty(p_flat.shape)
	pval_corrected.fill(np.nan) 
	pval_corrected[mask] = multipletests(p_flat[mask],method='fdr_bh')[1]

	#p_corrected1= (multipletests(p_flat,method='fdr_bh'))[1]
	p_corrected2= np.reshape(pval_corrected,(-1,numcols))
	df_rho		= pd.DataFrame(rho,index=list(matrix_df1a.columns.values), columns=list(matrix_df1a.columns.values))
	df_p1a 		= pd.DataFrame(p_corrected2,index=list(matrix_df1a.columns.values), columns=list(matrix_df1a.columns.values))
	df_p1b 		= df_p1a.loc[list_cat2]
	df_p1 		= df_p1b[list_cat1]
	del df_p1a
	del df_p1b
	df_p_sig1 	= (df_p1[(df_p1 <=0.1)])#.dropna(axis=0,how='all'))
	df_p_sig 	= (df_p_sig1[(df_p_sig1 <=0.1)])#.dropna(axis=0,how='all'))
	#df_p_sig 	= df_p_sig2.fillna(value=np.nan)
	sig_index 	= df_p_sig.index.values
	sig_column 	= df_p_sig.columns.values
	df_rho1 	= df_rho.loc[sig_index]
	df_rho_sig 	= df_rho1[sig_column]
	#df_rho_sig 	= (df_rho_sig1[(df_rho_sig1 <=(0.6))].dropna(axis=1,how='all'))
	df_rho2 	= df_rho.loc[list_cat2]
	df_rho_non 	= df_rho2[list_cat1]
	return df_rho_sig, df_rho_non, df_p_sig, df_p1, sig_index,sig_column,title

def plot_heatmap(coor_matrix,list_lablex,list_labley,out_heat_sig,df_mask=None):
	
	plotlt.rcParams['pdf.fonttype'] = 42
	plotlt.rcParams['ps.fonttype'] = 42
	sns.set(style="white")
	sns.plotting_context(font_scale=0.6)
	f, ax= plt.subplots(figsize=(5,10))
	#cmap=sns.palplot(sns.diverging_palette(220, 20, n=7))
	#cmap 	= sns.diverging_palette(220, 10, as_cmap=True)
	#cmap 	= sns.diverging_palette("RdBu_r", as_cmap=True)

	if df_mask is not None:
		mask_val = df_mask.isnull()
		g=sns.clustermap(coor_matrix,
   	 	xticklabels=list_lablex, yticklabels=list_labley,
    	method='average',metric='euclidean',mask=mask_val,cmap="RdBu_r")
		plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0)
		hm = g.ax_heatmap.get_position()
		g.ax_heatmap.set_position([hm.x0, hm.y0, hm.width*0.50, hm.height*0.50])
		col = g.ax_col_dendrogram.get_position()
		g.ax_col_dendrogram.set_position([col.x0, col.y0, col.width*0.25, col.height*0.75])
		g.ax_row_dendrogram.set_visible(False)
		g.ax_col_dendrogram.set_visible(True)

		#linewidths=.1,
		
	else:
		#sns.heatmap(coor_matrix, vmax=1.0,
		#square=True, xticklabels=list_lablex, yticklabels=list_labley,
		#x=ax,linewidths=.3, cbar_kws={"shrink": 0.3})
		g=sns.clustermap(coor_matrix,
   	 	xticklabels=list_lablex, yticklabels=list_labley,
    	method='average',metric='euclidean',figsize=(10,10),cmap="RdBu_r",vmin=-0.5, vmax=0.5)
    	#rotate the y lables

    	#linewidths=.3,
		plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0)

		#remove the dendogram on the row site
		g.ax_row_dendrogram.set_visible(True) 
		g.ax_col_dendrogram.set_visible(True)

	plt.savefig(out_heat_sig,bbox_inches='tight', dpi=300) #format='eps'


def g(text):
    replacements = {
        ' ' :'',
        '.':'_',
        '-':'_',
        '/':'_',
        ':':'_',
        '+':'_',
        '(':'_',
        ')':'_',
        ',':'_',
        "'":''

    }
    text = "".join([replacements.get(c, c) for c in text])
    return text

#generate  a custom diverging colormap
#read dataframe
argvs 			= sys.argv
category1_raw 	= read_dataframe(argvs[1])
category2_raw1 	= read_dataframe(argvs[2])

#change the index in dataframe 2
category2_raw_index = (list(category2_raw1.index.values))
index_new = {}
for y in category2_raw_index:
	y_new = g(y)
	index_new[y]=y_new
category2_raw = category2_raw1.rename(index=index_new)
del category2_raw1

columns_name 	=(list(category2_raw.columns.values))
add_row 		= []
for columns in columns_name:
	stage = columns.split('.')[1].strip()
	add_row.append(stage)	

#category2_raw.columns = pd.MultiIndex.from_tuples(zip(add_row, columns_name),names=['Class','Subject_ID'])

#category2_raw.to_csv("181109_Metabolomic_5.up.1e-06.up.tsv",header=True, index=True, sep='\t')

#df_1 			= read_dataframe(argvs[3])
#df_2 			= read_dataframe(argvs[4])

'''
list_cat1=["Eggerthella",
"Veillonella",
"Lactobacillus",
"Coprobacillus",
"Prevotella",
"Streptococcus",
"Alistipes",
"Roseburia",
"Bifidobacterium",
"Klebsiella",
"Akkermansia",
"Odoribacter",
"Anaerotruncus",
"Atopobium",
"Coprococcus",
"Blautia"]
'''

#create color based on range of value
def color_genVal(a,b):

	return dict_col

def dict_feat_colour(input1,colour_dict):
	 return dict_feat

list_cat1 =[]
list_cat1File = open((argvs[3]))
for x in list_cat1File:
	list_cat1.append(x.strip())


list_cat2 =[]
list_cat2File = open((argvs[4]))
for x in list_cat2File:
	list_cat2.append(x.strip())


#category1_sig,category1_lable 	= sig_list(dataframe=df_1)
category1_sig,category1_lable	= sig_list(list_inp=list_cat1)
#category2_sig,category2_lable 	= sig_list(dataframe=df_2)
category2_sig,category2_lable 	= sig_list(list_inp=list_cat2)

category1_df 	= category1_raw.loc[category1_sig]
category2_df 	= category2_raw.loc[category2_sig]
correlation_df 	= concatinate_df(category1_df,category2_df)

correlation_sig, correlation_non, dfp_value_sig, df_p_value,index_sig,column_sig,title3 = spearman_rank(correlation_df,category1_sig,category2_sig)#,stage='Gastrectomy')

#correlation_sig, correlation_non, dfp_value_sig, df_p_value,index_sig,column_sig,title3 = spearman_rank(correlation_df,category1_sig,category2_sig,stage='Healthy')
#title0 			= '1807masked'

del category1_raw
del category2_raw
title0 			= '190823'
title1 			= 'Metadata'
#title1 			= ((argvs[1]).split('_')[1]).split('.')[0]
title2 			= 'SpeciesmOTU'
out_rho 	 	= title0+'_'+title1+'_'+title2+'.'+title3+'.rhodf_LDA3.relab.tsv'
outitle_sig_df 	= title0+'_'+title1+'_'+title2+'.'+title3+'.sigdf_LDA3.relab.tsv'
outfig_df_sig 	= title0+'_'+title1+'_'+title2+'.'+title3+'.corrdf_mask.pdf'
#outfig_df_sig 	= title0+'_'+title1+'_'+title2+'.'+title3+'.corrdf_mask01LDA30.Relab.pdf'
dfp_value_sig.to_csv(outitle_sig_df,header=True, index=True, sep='\t')
correlation_sig.to_csv(out_rho,header=True, index=True, sep='\t')
#cat1_sig,cat1_lable 	= sig_list(list_inp=column_sig)
#cat2_sig,cat2_lable 	= sig_list(list_inp=index_sig)

#plot_heatmap(correlation_sig,category1_lable)
#plot_heatmap(correlation_non,category1_lable,category2_lable,outfig_df)
#unmasked heatmap
correlation_sig.dropna(axis=0, how='any', thresh=None, subset=None, inplace=True)
index_same = correlation_sig.index.values
dfp_value_sig = dfp_value_sig.loc[index_same]
#plot_heatmap(correlation_sig,column_sig,index_sig,outfig_df_sig)#df_mask=dfp_value_sig
#masked heatmap
plot_heatmap(correlation_sig,column_sig,index_sig,outfig_df_sig,df_mask=dfp_value_sig)