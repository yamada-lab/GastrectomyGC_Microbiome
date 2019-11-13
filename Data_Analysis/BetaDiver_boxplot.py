#!/usr/bin/python
#This code was perform to calculate the beta diversity and plot the boxplot
#ussage python PCoAboxplot_beta_diversity.py matrix_table

from scipy.stats import mannwhitneyu
from math import log
from math import sqrt
import sys
import pandas as pd
import numpy as np
from skbio.diversity import beta_diversity
from skbio.diversity import beta
from skbio.stats.distance import mantel
from skbio.stats.ordination import pcoa
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy.spatial.distance as distance
import seaborn as sns 



argvs 		= sys.argv
path 		= argvs[1]
title 		= str(argvs[1].split(".csv")[0])
outtitle 	= title+'.beta_diversity_bray_curtis.csv'
outtitle2 	= title+'.alpha_diversity_unweight.csv'
writer 		= pd.ExcelWriter(title+".beta_diver.xlsx", engine='xlsxwriter')
df_table1 	= pd.read_table(path,index_col=0)
outfig 		= title+'.beta_diversity_bray_curtisbox.eps'

#drop unecessary features
try:
	df_table = df_table1.drop(['mean','std'],axis=1,inplace=False)
except:
	df_table = df_table1

df_table_Tpose	= df_table.transpose()
list_sample		= df_table_Tpose.index.values
list_otu 		= df_table.index.values
sample_matrix 	= df_table_Tpose.as_matrix()

#dataframe partition (depend on group)
normal 		=(df_table[df_table.columns[df_table.columns.to_series().str.contains('Healthy')]])
surgery 	=(df_table[df_table.columns[df_table.columns.to_series().str.contains('Gastrectomy')]])

#bray curtis pairwise distance matrix (scikit-bio)
pw=beta_diversity('braycurtis', sample_matrix, ids=list_sample, validate=True, pairwise_func=None) #return class

#convert matrix into square dataframe
df=pd.DataFrame(data=pw[0:,0:],index=list_sample,columns=list_sample)

#remove duplicate values
df_values = list(distance.squareform(df))

#make the group matrix
gdf_values = list()
for i in range(len(list_sample)):
	for j in range(i, len(list_sample)):
		g1 = list_sample[i].split('.')[1]
		if g1 == 'Healthy':
			g1= 'Healthy'
		else:
			g1= 'Gastrectomy'
		g2 = list_sample[j].split('.')[1]
		if g2 == 'Healthy':
			g2= 'Healthy'
		else:
			g2= 'Gastrectomy'
		g_g = '-'.join(sorted([g1, g2]))

		gdf_values.append(g_g)

#boxplot
plt.figure()
fig,ax 				= plt.subplots()
pal 				=['#0072BC','#F0E442','#F7941D']
sns.set(style='ticks',font_scale=0.6)

df_for_box1 		= pd.DataFrame([df_values, gdf_values], index=['Bray curtis distance', 'group']).T
#df_for_box1=df_for_box.fillna(0)
#df_for_box = df_for_box1[(df_for_box1.T != 0).any()]
df_for_box 			= df_for_box1.dropna(axis=0, how='any')

df_for_box.to_csv('boxplot_distance.tsv',header=True, index=True, sep='\t')
sns.factorplot(y="Bray curtis distance", x="group", data=df_for_box, kind="box", palette=pal, linewidth=1.00, width=0.5, legend=False, showfliers=True)
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.savefig(outfig, bbox_inches='tight',format='eps', dpi=1000)
df.to_excel(writer,'Bray_curtis')
df.to_csv(outtitle,header=True, index=True, sep='\t')

#statistic analysis
def stat_diff(dataframe,var1_group,var2_group):
	var1 = (dataframe.loc[dataframe['group'] == var1_group, 'Bray curtis distance']).tolist()
	var2 = (dataframe.loc[dataframe['group'] == var2_group, 'Bray curtis distance']).tolist()
	hyp_test=mannwhitneyu(np.array(var1),np.array(var2),alternative='two-sided')
	p_value = hyp_test[1]
	mean_var1=np.array(var1).mean()
	std_var1=np.array(var1).std()
	mean_var2=np.array(var2).mean()
	std_var2=np.array(var2).std()
	var1=np.array(var1).var()
	var2=np.array(var2).var()
	value = list()
	value.extend([p_value,mean_var1,std_var1,mean_var2,std_var2])
	return value

HH_GH 	=stat_diff(df_for_box,"Healthy-Healthy", "Gastrectomy-Healthy")
GH_GG 	=list(stat_diff(df_for_box,"Gastrectomy-Healthy", "Gastrectomy-Gastrectomy"))
HH_GG 	=list(stat_diff(df_for_box,"Healthy-Healthy", "Gastrectomy-Gastrectomy"))

df_stats 	= pd.DataFrame([HH_GH,GH_GG,HH_GG], columns=['p_value','mean_var1','std_var1','mean_var2','std_var2'], index=['HH_GH','GH_GG','HH_GG'])
df_stats.to_csv('Stats_bray_curtis.tsv',header=True, index=True, sep='\t')