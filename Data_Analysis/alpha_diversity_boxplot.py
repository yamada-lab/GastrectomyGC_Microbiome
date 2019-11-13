#!/usr/bin/python
#16 April 2018
#This code was perform to calculate the alpha diversity followed by statistcal test
#ussage python alpha_diversity_boxplot.py matrix_table metric

from scipy.stats import mannwhitneyu
from math import log
from math import sqrt
import sys
import pandas as pd
import numpy as np
from skbio.diversity import alpha_diversity
from skbio.diversity import alpha
import seaborn as sns
import matplotlib.pyplot as plt


argvs 		= sys.argv
path 		= argvs[1]
index 		= argvs[2]
title 		= str(argvs[1].split(".csv")[0])
outtitle1	= title+'.alpha_diversity_' + str(index) +'.tsv'
outfig_1 	= title+'_'+str(index)+'.boxplot.eps'
out_calc 	= title+'.calc' + str(index) +'.tsv'

#writer = pd.ExcelWriter(title+"alpha_diver.xlsx", engine='xlsxwriter')
df_table = pd.read_table(path,index_col=0)
try:
	matrix2  		= df_table.drop(['mean','std'],axis=1,inplace=False)
except:
	matrix2 		= df_table
df_table_Tpose = matrix2.transpose()
list_sample=df_table_Tpose.index.values
sample_matrix = df_table_Tpose.as_matrix()


def stat_diff(dictionary,normal,surgery):
	calc_dict={}
	for x,y in dictionary.iteritems():
		if 'Healthy' in x:
			normal.append(y)
		else:
			surgery.append(y)
	hyp_test=mannwhitneyu(np.array(normal),np.array(surgery),alternative='two-sided')
	p_value = hyp_test[1]
	calc_dict['p_value']=p_value
	mean_normal=np.array(normal).mean()
	std_normal=np.array(normal).std()
	calc_dict['mean normal']=mean_normal
	calc_dict['std normal']=std_normal
	mean_surgery=np.array(surgery).mean()
	std_surgery=np.array(surgery).std()
	calc_dict['mean surgery']= mean_surgery
	calc_dict['std surgery']=std_surgery
	var1=np.array(normal).var()
	var2=np.array(surgery).var()
	n1,n2=len(normal), len(surgery)
	diff=float(mean_surgery - mean_normal)
	pooled_var=((n2*var2)+(n1*var1))/(n1+n2)
	odds_score= float(0)
	effect_size=diff/sqrt(pooled_var)
	calc_dict['effect size']=effect_size
	return calc_dict
normal = []
surgery = []
#p_value = 0.0
results_index 	= dict(alpha_diversity(str(index),sample_matrix,ids=list_sample,validate=True))
calcres			=stat_diff(results_index,normal,surgery)
y 				= str(index)+'_index'
column_name 	= ['Sample', y]

df_final1 		= pd.DataFrame(results_index.items(),columns=column_name)
df_final1.to_csv(outtitle1,header=True, index=True, sep='\t')
list_stage 		=[]
list_stage_alpha = df_final1.index.values
df_calc 		= pd.DataFrame.from_dict(calcres,orient='index')	
df_calc.to_csv(out_calc,header=True, index=True, sep='\t')

for x in list_stage_alpha:
	stage_x 	= (df_final1.loc[x,'Sample'])
	stage 		=stage_x.split('.')[1]
	list_stage.append(stage)
df_final1['stage']=list_stage
pal =['#0072BC','#F7941D']
sns.set(style='ticks')
sns.boxplot(x="stage", y=y, data=df_final1,linewidth=1.00,width=0.3,palette=pal)
sns.despine(offset=10, trim=True)
plt.savefig(outfig_1, format='eps', dpi=1000)
#df_final2.to_excel(writer,'shannon')
#iterate over column to get the alpha diversity


#http://scikit-bio.org/docs/0.1.3/math.diversity.alpha.html
