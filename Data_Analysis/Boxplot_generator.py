#!/usr/bin/python
#This code was perform to create boxplot for intended list
#ussage python boxplot_list_verAll.py list_of_features_to_plot

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


#get the essentials information of the inputed data
argvs 			= sys.argv
matrix_raw 		= argvs[1]
#feature_inp		= open(argvs[2])
matrix_bef 		= pd.read_table(matrix_raw,index_col=0)


#dataframe and elimination of citrate
feature_list	= matrix_bef.index.values
if 'C00158__Citrate' in feature_list:
	matrix 	= matrix_bef.drop(['C00158__Citrate'],axis=0, inplace=False)
	
else:
	matrix  	= matrix_bef
df_table_Tpose 	= matrix.transpose()
list_sample 	= df_table_Tpose.index.values
#list of feature
list_feature 	= []
feature_inp = df_table_Tpose.columns.values

list_stage_n 	= []
for x in list_sample:
	stage 		= x.split('.')[1]
	list_stage_n.append(stage)


for feature in feature_inp:
	#list_feature.append(x.split("\n")[0])
	outfig 			= feature+'boxplot.eps'
	df_feature 		= pd.DataFrame(df_table_Tpose[feature])
	
	df_feature['category'] = list_stage_n
	#create the boxplot
	import seaborn as sns 
	plt.figure()
	fig,ax 				= plt.subplots()
	pal 				=['#0072BC',"#F7941D"]
	sns.set(style='ticks',font_scale=0.6)
	sns.boxplot(x='category', y=feature, data=df_feature,linewidth=1.00,width=0.1,palette=pal,showfliers=True)
	sns.despine(offset=10, trim=True)
	#plt.xticks(rotation=90)
#plt.xlim(0.0,4000)
	plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
	plt.savefig(outfig, bbox_inches='tight',format='eps', dpi=310)


#modify the x and y axis while playing with orient