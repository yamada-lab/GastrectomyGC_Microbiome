#!/usr/bin/python
#13 February 2017
#This code was perform to calculate the alpha diversity followed by statistcal test
#only for interested phylum
#ussage python alpha_phylum.py matrix_table reference_table

from scipy.stats import mannwhitneyu
from math import log
from math import sqrt
import sys
import pandas as pd
import numpy as np
from skbio.diversity import alpha_diversity
from skbio.diversity import alpha
import matplotlib as plotlt
plotlt.use("Agg")
import seaborn as sns
import matplotlib.pyplot as plt
from statsmodels.sandbox.stats.multicomp import multipletests

#get the essentials information of the inputed data
argvs 			= sys.argv
matrix_raw 		= argvs[1]
reference 	 	= open(argvs[2])
method_select 	= str(argvs[3])
matrix_bef 		= pd.read_table(matrix_raw,index_col=0)
title 			= str(argvs[1].split(".csv")[0])
outtitle1		= title+'.phylum_alpha_'+method_select+'.1.csv'
outtitle2 		= title+'.phylum_stat.alpha_'+method_select+'.1.csv'
outfig	 		= title+'phylum_boxplot.'+method_select+'.1.eps'
#writer 		= pd.ExcelWriter(title+"variance_scatter.xlsx", engine='xlsxwriter')
list_feature 	= matrix_bef.index.values
stage 			= matrix_bef.columns.values

#dataframe
try:
	matrix 			= matrix_bef.drop(['mean','std'],axis=1)
except:
	matrix 			= matrix_bef
matrix_Tpose 	= matrix.transpose()
#sample_matrix 	= matrix_Tpose.as_matrix()
list_sample		= matrix_Tpose.index.values
normal_1		= (matrix[matrix.columns[matrix.columns.to_series().str.contains('Healthy')]])
surgery_1 		= (matrix[matrix.columns[matrix.columns.to_series().str.contains('Gastrectomy')]])

#cretate the dictionary with key : Family; value: list of specides
#fam_list		= ['Bacteroidetes','Firmicutes','Proteobacteria','Actinobacteria','Euryarchaeota','Fusobacteria','Verrucomicrobia']
fam_list		= ['Bacteroidetes','Firmicutes','Proteobacteria','Actinobacteria','Fusobacteria']
fam_dict 		= {}
all_fam 		= {}
for column in reference:
	spec_code 	= (column[:-1].split("\t"))[-1]
	family 	 	= ((column[:-1].split("\t"))[2])

	if family in fam_list and spec_code in list_feature:
		if family not in fam_dict:
			fam_dict[family] = []
			if ("NA" in spec_code)!=True and spec_code not in fam_dict[family]:
				fam_dict[family].append(spec_code)
		else:
			if ("NA" in spec_code)!=True and spec_code not in fam_dict[family]:
				fam_dict[family].append(spec_code)


#function
#calculate the alpha diversity for each family
#dataframe is the after the selected feature only
def alpha_div(dataframe,list_s,fam_code,methode):
	df_matrix 	= dataframe.as_matrix()
	results 	= ((alpha_diversity(methode,df_matrix,ids=list_s,validate=True)))
	df_results	= pd.DataFrame(results,columns=[fam_code])
	df_results=df_results.fillna(0)
	return df_results, results

#the statistic differeces
def stat_diff(dictionary,normal,surgery,family,method):
	stat_list = []
	for x,y in dictionary.iteritems():
		if 'Healthy' in x:
			normal.append(y)
		else:
			surgery.append(y)
	hyp_test 			= mannwhitneyu(np.array(normal),np.array(surgery))
	p_value 			= hyp_test[1]
	stat_list.append(p_value)
	#print p_value
	mean_normal 		= np.array(normal).mean()
	mean_surgery 		= np.array(surgery).mean()
	var1				= np.array(normal).var()
	var2				= np.array(surgery).var()
	n1,n2				= len(normal), len(surgery)
	diff 				= float(mean_surgery - mean_normal)
	pooled_var 			= ((n2*var2)+(n1*var1))/(n1+n2)
	odds_score 			= float(0)
	effect_size 		= diff/sqrt(pooled_var)
	stat_list.append(effect_size)
	try:
		odds_score=np.log2((mean_surgery/mean_normal))
		stat_list.append(odds_score)
		#diff=float((mean_normal - mean_surgery)/std_dev)
	except:
		stat_list.append(0.0)
	stat_df  			= pd.DataFrame(stat_list,index=['p_value', 'effect size','log_odds_score'],columns=[family])
	return stat_df,p_value

#initialize the dataframe for the alpha diversity
df_alpha 				= pd.DataFrame([])
df_stat 				= pd.DataFrame([])

#initialize the dataframe for the statistical analysis
#list of p_value
p_val_list  			= []
for key in fam_list:
	each 				= fam_dict[key]
	df_select			= matrix_Tpose[each]
	df_fam,result_fam 	= alpha_div(df_select,list_sample,key,method_select)
	df_alpha 			= df_alpha.join(df_fam,how='outer')
	dict_result 		= dict(result_fam)
	normal_fam 			= []
	surgery_fam 		= []
	stat_fam,p_val_fam 	= stat_diff(dict_result,normal_fam,surgery_fam,key,method_select)
	df_stat 			= df_stat.join(stat_fam,how='outer')
	p_val_list.append(p_val_fam)

#get the hue and write to file
list_stage_n 			= []
for x in list_sample:
	stage 		= x.split('.')[1]
	list_stage_n.append(stage)
df_alpha['category'] 	= list_stage_n
df_alpha.to_csv(outtitle1,header=True, index=True, sep='\t')

#get the multiple test p_value correction
correction_p 			= (multipletests(p_val_list,method='fdr_bh'))[1]
stats_column 			= df_stat.columns.values
df_correct  			= (pd.DataFrame(correction_p, index=stats_column,columns=['corrected_p_value'])).transpose()
df_stat 	 			= df_stat.append(df_correct)
df_stat.to_csv(outtitle2,header=True, index=True, sep='\t')
writer 					= pd.ExcelWriter(title+'.phylum_alpha_'+method_select+".xlsx", engine='xlsxwriter')
df_alpha.to_excel(writer,'alpha_diversity')
df_stat.to_excel(writer,'statistical_test')

#create the boxplot
plt.figure()
fig,ax 				= plt.subplots()
pal 				=['#0072BC','#F7941D']
sns.set(style='ticks',font_scale=0.6)
df_long 			= pd.melt(df_alpha, 'category', var_name="phylum", value_name="alpha_diversity")
sns.factorplot("phylum", hue='category', y="alpha_diversity", data=df_long, kind="box",palette=pal,linewidth=1.00,width=0.5,legend=False)
sns.despine(offset=10, trim=True)
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.savefig(outfig, bbox_inches='tight',format='eps', dpi=1000)