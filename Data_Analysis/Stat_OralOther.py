#!/usr/bin/python
# 6 January 2017
#This code was created to calculate the ratio of oralGutratio, and create the boxplot 
#ussage python Stat_OralOther.py <species_table> <oral_spec_table>

from scipy.stats import mannwhitneyu
from math import log
from math import sqrt
import sys
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

argvs = sys.argv
path = (argvs[1])
oral_microbes = (argvs[2])
title=str(argvs[1].split(".csv")[0])

#create the list of the oral microbes
oral_microbes_df = pd.read_excel(oral_microbes,index_col=0)
df_table = pd.read_table(path,index_col=0)
list_numb = oral_microbes_df.index.values
species_oral =[]
for i in list_numb:
	species = (oral_microbes_df.loc[i,'Species'])
	genus 	= (oral_microbes_df.loc[i,'Genus'])
	#name 	= genus + ' ' + species
	#metaphlan2
	name 	= genus + '_' + species
	species_oral.append(name)

#list up the sample name
df_table_Tpose = df_table.transpose()
list_sample=df_table_Tpose.index.values

# function to get the oral microbes and other proportion and the ratio of it (oral/other)
def tendency(dictionary_tend,series_cat,list_loop,list_ref,list_oral):
	for i in list_loop:
		if i in list_ref:
			rel_oral = series_cat.loc[i]
			dictionary_tend['oral']+=rel_oral
			if i not in list_oral and rel_oral!=0.0:
				list_oral.append(i)
		else:
			rel_other = series_cat.loc[i]
			dictionary_tend['other']+=rel_other
	dictionary_tend['ratio']=dictionary_tend['oral']/dictionary_tend['other']
	return dictionary_tend

#function to get the statistical difference
def stat_diff(dictionary_normal,dictionary_surgery,normal,surgery,or_normal,or_surgery):
	calc_dict={}
	for x,y in dictionary_normal.iteritems():
		ratio=y['ratio']
		normal.append(ratio)
		oral =y['oral']
		or_normal.append(oral)
	for x,y in dictionary_surgery.iteritems():
		ratio=y['ratio']
		surgery.append(ratio)
		oral = y['oral']
		or_surgery.append(oral)
	#print normal
	#print surgery
	hyp_test=mannwhitneyu(np.array(normal),np.array(surgery),alternative='two-sided')
	p_value = hyp_test[1]
	hyp_test_oral=mannwhitneyu(np.array(or_normal),np.array(or_surgery))
	p_value_oral = hyp_test_oral[1]
	print p_value
	print p_value_oral
	mean_normal=np.array(normal).mean()
	mean_surgery=np.array(surgery).mean()
	var1=np.array(normal).var()
	var2=np.array(surgery).var()
	n1,n2=len(normal), len(surgery)
	diff=float(mean_surgery - mean_normal)
	pooled_var=((n2*var2)+(n1*var1))/(n1+n2)
	odds_score= float(0)
	effect_size=diff/sqrt(pooled_var)
	calc_dict['p_value_ratio']=p_value
	calc_dict['p_value_oral']=p_value_oral
	calc_dict['mean normal']=mean_normal
	calc_dict['std normal']=var1
	calc_dict['mean surgery']=mean_surgery
	calc_dict['std surgery']=var2
	calc_dict['effect size']=effect_size
	return calc_dict
			

#get the oral and other composition
tendency_dict_normal={}
tendency_dict_surgery={}
each_dict_normal={'oral':0.0,'other':0.0}
each_dict_surgery={'oral':0.0,'other':0.0}
normal_oral = []
surgery_oral = []
for x in list_sample:
	#if x.startswith('Normal'):
	if 'Healthy' in x:
		each_sample 	= df_table[x]
		list_species	= each_sample.index.values
		tendency_dict_normal[x] = tendency(each_dict_normal,each_sample,list_species,species_oral,normal_oral)
		each_dict_normal={'oral':0.0,'other':0.0}
for x in list_sample:
	if 'Gastrectomy' in x:
	#if x.startswith('Stomach'):
		each_sample 	= df_table[x]
		list_species	= each_sample.index.values
		tendency_dict_surgery[x] = tendency(each_dict_surgery,each_sample,list_species,species_oral,surgery_oral)
		each_dict_surgery={'oral':0.0,'other':0.0}

normal		= []
surgery 	= []
or_normal 	= []
or_surgery 	= []
out_calc 	= title+'.calc' +'OralOther'+'.tsv'
calcres=stat_diff(tendency_dict_normal,tendency_dict_surgery, normal,surgery,or_normal,or_surgery)
df_calc 		= pd.DataFrame.from_dict(calcres,orient='index')	
df_calc.to_csv(out_calc,header=True, index=True, sep='\t')
#print normal_oral
#print surgery_oral

tendency_dict = dict(tendency_dict_normal, **tendency_dict_surgery)
df_new = pd.DataFrame(tendency_dict)
df_new_Tpose = df_new.transpose()
#df_new_Tpose.index.name = 'Sample'
#write to file
outitle = title+'.OralOtherRatio'+".tsv"
outfig_oral = title+'.Oral'+".pdf"
outfig_rat 	= title+'.OralOtherRatio'+".pdf"
df_new_Tpose.to_csv(outitle,header=True,index=True,sep='\t')
writer = pd.ExcelWriter(title+'.OralOtherRatio'+".xlsx", engine='xlsxwriter')
df_new_Tpose.to_excel(writer)
#create the boxplot
def boxplot(dataframe,list_stage,purpose,outfig):
	list_stage_pur	= dataframe.index.values
	for x in list_stage_pur:
		stage 		= x.split('.')[1]
		list_stage.append(stage)
	dataframe['stage']=list_stage
	pal =['#0072BC','#F7941D']	
	sns.set(style='ticks')
	sns.boxplot(x="stage", y=purpose, data=dataframe,linewidth=1.00,width=0.3,palette=pal)
	sns.despine(offset=10, trim=True)
	plt.savefig(outfig, dpi=500)

#create the boxplot of oral bacteria and ratio
#oral
list_stage_oral = []
purpose_oral 	= 'oral'
boxplot(df_new_Tpose,list_stage_oral,purpose_oral,outfig_oral)
# clear the plot
plt.clf()
#ratio
list_stage_ratio = []
purpose_ratio 	= 'ratio'
boxplot(df_new_Tpose,list_stage_ratio,purpose_ratio,outfig_rat)



