#!/usr/bin/python
#Erawijantari-Yamada Laboratory-Titech-Aug,2 2019
#Distribution plot samplewise

from __future__ import print_function
from __future__ import unicode_literals
# encoding=utf8  
import sys  
reload(sys)  
sys.setdefaultencoding('utf8')
import pandas as pd 
import numpy as np
import os
import argparse
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as plotlt



def read_params(args):
	parser = argparse.ArgumentParser(description='distribution_plot')
	parser.add_argument('input_file', metavar='INPUT_FILE', type=str, help="relative abundance file")
	parser.add_argument('metadata', metavar='METADATA_TABLE', type=str, help="the metadata file")
	parser.add_argument('out_file', metavar='OUTPUT_FILE', type=str, help="the output file")
	parser.add_argument('-f',dest="feature", metavar='str', type=str, default='',
	help="features selected for the visualization")
	parser.add_argument('-o',dest="norm", metavar='str', type=bool, default=False,
	help="normalization into log scale or not")
	parser.add_argument('-s',dest="sort", metavar='str', type=str, default='',
	help="whether to short the samples based on relab value or not")
	parser.add_argument('-c',dest="colormap", metavar='str', type=str, default='',
	help="whether to short the samples based on relab value or not")
	parser.add_argument('-r',dest="orientation", metavar='str', type=str, default='v',
	help="orientation of the plot")
	parser.add_argument('-t',dest="type", metavar='str', type=str, default='barplot',
	help="boxplot or barplot")
	args = parser.parse_args() 
	params = vars(args)
	return params 

#generate the reference for the file (contain the type of surgery and recontr)
def Refname_Gen(RefFile):
	refname = {}
	for tab in RefFile:
		tab1	= tab.split("\t")[0]
		key		= tab1+'.'+tab.split("\t")[1]
		value1	= (tab.split('\t')[2]).strip() #type of surgery
		value2	= (tab.split('\t')[3]).strip()#reconstruction
		if key not in refname:
			refname[key]=[value1,value2]
	return refname

def log_transform(df):
	#add pseudocount 
	pseudo = float(0.000000010)
	#add all values with that
	df_return1 = df+pseudo
	#log
	df_return= df_return1.apply(np.log10)
	return df_return

def df_read(df,feat,refname,log):
	#read dataframe convert to long form with surgery type and reconstruction
	feature_list	= df.index.values
	if 'C00158__Citrate' in feature_list:
		matrix1 	= df.drop(['C00158__Citrate'],axis=0, inplace=False)	
	else:
		matrix1 = df
	if log == True:
		print("Yes")
		matrix=log_transform(matrix1)
	else:
		print('No')
		matrix=matrix1
	list_sample 	= (list(matrix.columns.values))
	type_surgery 	= []
	reconstruction 	= []
	for x in list_sample:
		t_surgery = refname[x][0]
		r_surgery = refname[x][1]
		type_surgery.append(t_surgery)
		reconstruction.append(r_surgery)
	df_feature 		= matrix.loc[feat]
	df_feature.columns = pd.MultiIndex.from_tuples(zip(list_sample, type_surgery, reconstruction),names=['Subject','Type_Surgery','Reconstruction'])
	#remove the unknown

	return df_feature


def boxplot(df_inp,outitle,orient):
	#create the boxplot
	plotlt.rcParams['pdf.fonttype'] = 42
	plotlt.rcParams['ps.fonttype'] = 42
	plt.figure()
	#pal 				=['#0072BC',"#F7941D"]
	sns.set(style='ticks',font_scale=0.6)
	index=list(df_inp.index.values)

	for ind in index:
		outfig=outitle+ind+'boxplot'+'.pdf'
		df_select 		= df_inp.loc[[ind]]
		df=pd.melt(df_select,value_name=ind)
		#df.sort_values(by=[ind],inplace=True)
		df.sort_values(by=['Type_Surgery'],inplace=True)
		flierprops = dict(markerfacecolor='0.75', markersize=5,linestyle='none')
		#pal=sns.color_palette(['#0072BC','#F7941D','#c57617'])
		pal=sns.color_palette(['#0072BC','#d4dc38','#dc3882','#b679e7','#dc9238','#82dc38','#7b4a0e'])
		if orient=='v':
			sns.catplot(x='Type_Surgery', y=ind,hue='Reconstruction',kind='box',data=df,palette=pal,linewidth=0.5,legend=False,orient='v',flierprops={"marker": "o","markersize":0.5})
		else:
			sns.catplot(x='ind', y='Type_Surgery',hue='Reconstruction',kind='box',data=df,palette=pal,linewidth=0.5,legend=False,orient='h')#flierprops={"marker": "o","markersize":0.5}
		sns.despine(offset=10, trim=True)
		plt.xticks(rotation=90)
		#plt.xlim(0.0,4000)
		plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
		plt.savefig(outfig, bbox_inches='tight',format='pdf', dpi=300)

def boxplot_gas(df_inp,outitle,orient):
	#create the boxplot
	plotlt.rcParams['pdf.fonttype'] = 42
	plotlt.rcParams['ps.fonttype'] = 42
	plt.figure()
	#pal 				=['#0072BC',"#F7941D"]
	sns.set(style='ticks',font_scale=0.6)
	index=list(df_inp.index.values)

	for ind in index:
		outfig=outitle+ind+'boxplot.gastrectomy'+'.pdf'
		df_select 		= df_inp.loc[[ind]]
		df=pd.melt(df_select,value_name=ind)
		#df.sort_values(by=[ind],inplace=True)
		df.sort_values(by=['Type_Surgery'],inplace=True)
		flierprops = dict(markerfacecolor='0.75', markersize=5,linestyle='none')
		#pal=sns.color_palette(['#0072BC','#F7941D','#c57617'])
		pal=sns.color_palette(['#0072BC','#d4dc38','#dc3882','#b679e7','#dc9238','#82dc38','#7b4a0e'])
		if orient=='v':
			sns.catplot(x='Type_Surgery', y=ind,kind='box',data=df,palette=pal,linewidth=0.5,legend=False,orient='v',flierprops={"marker": "o","markersize":0.5})
		else:
			sns.catplot(x='ind', y='Type_Surgery',kind='box',data=df,palette=pal,linewidth=0.5,legend=False,orient='h')#flierprops={"marker": "o","markersize":0.5}
		sns.despine(offset=10, trim=True)
		plt.xticks(rotation=90)
		#plt.xlim(0.0,4000)
		plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
		plt.savefig(outfig, bbox_inches='tight',format='pdf', dpi=300)

	
def viz(df_multi,outitle):
	plotlt.rcParams['pdf.fonttype'] = 42
	plotlt.rcParams['ps.fonttype'] = 42
	plt.figure()
	sns.set(style='ticks',font_scale=0.6)
	sns.despine(offset=10, trim=True)
	index=list(df_multi.index.values)

	for ind in index:
		#pal=sns.color_palette(['#0072BC','#F7941D','#f9b460','#fbc98e','#c57617','#945811','#7b4a0e'])
		pal=sns.color_palette(['#0072BC','#d4dc38','#dc3882','#b679e7','#dc9238','#82dc38','#7b4a0e'])	
		outfig=outitle+ind+'barplot'+'.pdf'
		df_select 		= df_multi.loc[[ind]]
		df1=pd.melt(df_select,value_name=ind)
		
		df2=df1.sort_values(by=[ind])
		df=df2.sort_values(by=['Type_Surgery'])
		ax = sns.catplot(x="Subject", y=ind, hue="Reconstruction",
                     data=df,kind="bar",dodge=False) #col='Type_Surgery'

		#g=sns.FacetGrid(df,col='Type_Surgery',sharex=False,row_order=df["Subject"])
		#g=g.map(sns.barplot,'Subject',ind,'Reconstruction',dodge=False,palette=pal)#hue_order=df["Reconstruction"])#hue_order=np.unique(df["Reconstruction"])
		#g.add_legend()
		sns.despine(offset=10, trim=True)
		plt.xticks(rotation=90)
		plt.tick_params(labelsize=5)
		#plt.xlim(0.0,4000)
		plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
		plt.savefig(outfig, dpi=300)
	#loop based on the list of columns


def main( ):
	params = read_params(sys.argv)
	df_input = pd.read_table(open(params['input_file']),index_col=0)
	metadata = open(params['metadata'])
	featureF = open(params['feature'])
	feature=[]
	for i in featureF:
		feature.append(i.strip())
	out= str(params['out_file'])
	orientation=str(params['orientation'])
	normalize = bool(params['norm'])
	typeplot = str(params['type'])

	#generate the reference
	Reffile=Refname_Gen(metadata)
	#prepare the datafame
	df_ready=df_read(df_input,feature,Reffile,normalize)
	#plot either barplot or boxplot
	if typeplot=='barplot':
		viz(df_ready,out)
	else:
		boxplot(df_ready,out,orientation)
		boxplot_gas(df_ready,out,orientation)

if __name__ == '__main__':
	main( )