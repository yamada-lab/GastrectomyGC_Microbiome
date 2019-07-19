#!/usr/bin/python
#Erawijantari-Titech; April, 2018
#Code for statistical analysis on metadata MWU for numerical and Fisher for categorical 
#The boxplot will be generated for numerical metadata
# usage:
#python Metadata_stats.py metadata_table reference columns_to_test(e.g BMI, age,etc)

'''
example of reference
ID  BMI age 
1.Gastrectomy   20  52
2.Healthy   18  67

'''

import matplotlib as plotlt
plotlt.use("Agg")
import pylab
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.sandbox.stats.multicomp import multipletests
import math
from scipy.stats import mannwhitneyu
#from scipy.stats import fisher_exact
from math import log
from math import sqrt
from math import isnan
from collections import Counter
from FisherExact import fisher_exact

argvs 			= sys.argv
file_conv		= {}
#generate the name converter
fileref			= open(argvs[2])
columns 		= argvs[3]

#for numerical data : Mann Whitney U test
def stat_diff_num(dictionary,Healthy,Gastrectomy):
	calc_dict={}
	for x,y in dictionary.iteritems():
		if 'Healthy' in x:
			Healthy.append(y)
		elif 'Gastrectomy'in x:
			Gastrectomy.append(y)
	hyp_test=mannwhitneyu(np.array(Healthy),np.array(Gastrectomy),alternative='two-sided')
	p_value = hyp_test[1]
	calc_dict['p_value']=p_value
	mean_Healthy=np.array(Healthy).mean()
	std_Healthy=np.array(Healthy).std()
	calc_dict['mean Healthy']=mean_Healthy
	calc_dict['std Healthy']=std_Healthy
	mean_Gastrectomy=np.array(Gastrectomy).mean()
	std_Gastrectomy=np.array(Gastrectomy).std()
	calc_dict['mean Gastrectomy']= mean_Gastrectomy
	calc_dict['std Gastrectomy']=std_Gastrectomy
	var1=np.array(Healthy).var()
	var2=np.array(Gastrectomy).var()
	n1,n2=len(Healthy), len(Gastrectomy)
	diff=float(mean_Gastrectomy - mean_Healthy)
	pooled_var=((n2*var2)+(n1*var1))/(n1+n2)
	odds_score= float(0)
	effect_size=diff/sqrt(pooled_var)
	calc_dict['effect size']=effect_size
	print len(Healthy), len(Gastrectomy)
	return calc_dict

#for categorical data: Fisher exact test
def stat_diff_cat(dictionary,Healthy,Gastrectomy):
	calc_dict={}
	for x,y in dictionary.iteritems():
		if 'Healthy' in x:
			Healthy.append(y)
		elif 'Gastrectomy'in x:
			Gastrectomy.append(y)
	healthy_c 		= Counter(Healthy)
	gastrectomy_c 	= Counter(Gastrectomy)
	array_dict		={}
	array_dict['Healthy'] = healthy_c
	array_dict['Gastrectomy']= gastrectomy_c
	df1 				= (pd.DataFrame.from_dict(array_dict))
	df 					= df1.fillna(0)
	print df
	#array 			= df.as_matrix(df)
	hyp_test=fisher_exact(df)
	
	calc_dict['p_value']=hyp_test
	#calc_dict['oddsratio']=hyp_test[0]

	return calc_dict
#old
#for tab in fileref:
#	tab2		= tab.split("\t")[1][:-1]
#	pat_ID 		= tab2.split('.')[1].split('-')[0]
#	if pat_ID not in file_conv:
#		file_conv[pat_ID]=tab2

#new
patID_l1 =[]


for tab in fileref:
	pat_ID 		= (tab.split("\t")[0]).strip()
	patID_l1.append(pat_ID)
	stage 		= (tab.split("\t")[1]).strip()
	newID 		= pat_ID+'_'+stage
	if pat_ID not in file_conv:
		file_conv[pat_ID]=newID
patID_l = [int(y) for y in patID_l1]
print patID_l
path=open(argvs[1])
df_table1=((pd.read_table(path,index_col=0,error_bad_lines=False))).loc[patID_l]
df_table = df_table1.transpose()
#change the name
row=(df_table.columns.values)
row2 = [str(x) for x in row]

for y in row2:
	if y in file_conv:
		df_table.columns = df_table.columns.astype(str)
		df_table.columns=df_table.columns.str.replace(y,file_conv[y])
	else:
		print "None"

#remove all of the row that contain zero value
df_table_Tpose = df_table.transpose()

list_sample = df_table_Tpose.index.values

list_stage_n 	= []
for x in list_sample:
	stage 		= str(x).split('_')[1]
	list_stage_n.append(stage)
df_table_Tpose['category'] = list_stage_n

#create the boxplot
#list of feature
if argvs[3]			== "all":
	list_feature 		= df_table_Tpose.columns.values
else:
	list_feature 		= [argvs[3]]
for cat in list_feature:
	outfig 				= cat+"metadata.eps"
	dict_df1 			= df_table_Tpose[cat].to_dict()
	Healthy 			= []
	Gastrectomy 		= []
	x 					= (dict_df1.values()[1])
	
	if isinstance(x, str)==True:
		calcres 		=stat_diff_cat(dict_df1,Healthy,Gastrectomy)
	else:
		dict_df			= {k: dict_df1[k] for k in dict_df1 if not isnan(dict_df1[k])}
		calcres			= stat_diff_num(dict_df,Healthy,Gastrectomy)
		
		plt.figure()
		fig,ax 				= plt.subplots()
		pal 				=['#0072BC','#F7941D']
		sns.set(style='ticks',font_scale=0.6)
		val 				= str(cat)
		df_table_ext 		= df_table_Tpose[[cat,'category']]
		df_long 			= pd.melt(df_table_ext, 'category', var_name="feature", value_name=val)
		sns.factorplot(y=val, hue='category', x="feature", data=df_long, kind="box",palette=pal,linewidth=1.00,width=0.5,legend=False,showfliers=False)
		sns.despine(offset=10, trim=True)
		plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
		plt.savefig(outfig, bbox_inches='tight',format='eps', dpi=1000)
		plt.close()
	
	out_calc 			= str(cat)+'.statscalc'+'.tsv'
	df_calc 			= pd.DataFrame.from_dict(calcres,orient='index')	
	df_calc.to_csv(out_calc,header=True, index=True, sep='\t')
