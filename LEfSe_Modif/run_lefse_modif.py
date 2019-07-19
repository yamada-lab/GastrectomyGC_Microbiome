#!/usr/bin/env python

import os,sys,math,pickle
from lefse_modif import *

####Notes###
#Modification by Erawijantari-Titech
#Script were modified from the original LEfSe script to enable the FDR calculations and writing to the output
#All of the modifications were marked with "###"

def read_params(args):
        parser = argparse.ArgumentParser(description='LEfSe 1.0')
        parser.add_argument('input_file', metavar='INPUT_FILE', type=str, help="the input file")
        parser.add_argument('output_file', metavar='OUTPUT_FILE', type=str,
                help="the output file containing the data for the visualization module")
        parser.add_argument('-o',dest="out_text_file", metavar='str', type=str, default="",
                help="set the file for exporting the result (only concise textual form)")
        parser.add_argument('-a',dest="anova_alpha", metavar='float', type=float, default=0.05,
                help="set the alpha value for the Anova test (default 0.05)")
        parser.add_argument('-w',dest="wilcoxon_alpha", metavar='float', type=float, default=0.05,
                help="set the alpha value for the Wilcoxon test (default 0.05)")
        parser.add_argument('-l',dest="lda_abs_th", metavar='float', type=float, default=2.0,
                help="set the threshold on the absolute value of the logarithmic LDA score (default 2.0)")
        parser.add_argument('--nlogs',dest="nlogs", metavar='int', type=int, default=3,
		help="max log ingluence of LDA coeff")
        parser.add_argument('--verbose',dest="verbose", metavar='int', choices=[0,1], type=int, default=0,
		help="verbose execution (default 0)")
        parser.add_argument('--wilc',dest="wilc", metavar='int', choices=[0,1], type=int, default=1,
		help="wheter to perform the Wicoxon step (default 1)")
	parser.add_argument('-r',dest="rank_tec", metavar='str', choices=['lda','svm'], type=str, default='lda',
		help="select LDA or SVM for effect size (default LDA)")
	parser.add_argument('--svm_norm',dest="svm_norm", metavar='int', choices=[0,1], type=int, default=1,
		help="whether to normalize the data in [0,1] for SVM feature waiting (default 1 strongly suggested)")
        parser.add_argument('-b',dest="n_boots", metavar='int', type=int, default=30,
                help="set the number of bootstrap iteration for LDA (default 30)")
        parser.add_argument('-e',dest="only_same_subcl", metavar='int', type=int, default=0,
                help="set whether perform the wilcoxon test only among the subclasses with the same name (default 0)")
        parser.add_argument('-c',dest="curv", metavar='int', type=int, default=0,
                help="set whether perform the wilcoxon test ing the Curtis's approach [BETA VERSION] (default 0)")
        parser.add_argument('-f',dest="f_boots", metavar='float', type=float, default=0.67,
                help="set the subsampling fraction value for each bootstrap iteration (default 0.66666)")
        parser.add_argument('-s',dest="strict", choices=[0,1,2], type=int, default=0,
                help="set the multiple testing correction options. 0 no correction (more strict, default), 1 correction for independent comparisons, 2 correction for independent comparison")
#       parser.add_argument('-m',dest="m_boots", type=int, default=5, 
#               help="minimum cardinality of classes in each bootstrapping iteration (default 5)")
        parser.add_argument('--min_c',dest="min_c", metavar='int', type=int, default=10,
                help="minimum number of samples per subclass for performing wilcoxon test (default 10)")
        parser.add_argument('-t',dest="title", metavar='str', type=str, default="",
                help="set the title of the analysis (default input file without extension)")
        parser.add_argument('-y',dest="multiclass_strat", choices=[0,1], type=int, default=0,
                help="(for multiclass tasks) set whether the test is performed in a one-against-one ( 1 - more strict!) or in a one-against-all setting ( 0 - less strict) (default 0)")
        parser.add_argument('-q',dest="FDR_alpha", metavar='float', type=float, default=0.1,
                help="set the alpha value for the Anova test (default 0.1)")
        args = parser.parse_args()
          
        params = vars(args)
        if params['title'] == "": params['title'] = params['input_file'].split("/")[-1].split('.')[0]
        return params 


if __name__ == '__main__':
	init()
	params = read_params(sys.argv)
	import pandas as pd
	#cls is the class or dictionary that save the variables
	feats,cls,class_sl,subclass_sl,class_hierarchy = load_data(params['input_file'])
	###feats is the dictionary of features and relative abundance
	kord,cls_means = get_class_means(class_sl,feats)
	wilcoxon_res = {}
	wilcoxon_resval = {}
	wilcoxon_resadj = {}
	kw_n_ok = 0
	nf = 0
	###record in the dictionary
	Feat_pval ={} ###key is the features; value is the kw status and p-value
	for feat_name,feat_values in feats.items():
		###feat name is the features names to test
		if params['verbose']:
			print "Testing feature",str(nf),": ",feat_name,
			nf += 1
		###p-value is pv here
		kw_ok,pv = test_kw_r(cls,feat_values,params['anova_alpha'],sorted(cls.keys()))
		#print kw_ok,pv #kw_ok is the T or F; pv is the pv
		###try to write the feat_name, pv to file 
		###save the p-value to the dict

		###for the False kruskall wallis (results does not violate the H0) then
		###the wilcoxon test are not performed, therefore the pval reported as "-"
####BIG MOCIFICATIONS SHOULD PERFORMED HERE TO SELECT THE FEATURES PASSED TO WILCOXON AND LDA
		if feat_name not in Feat_pval:
			Feat_pval[feat_name]=[kw_ok,pv]
	FDRdict 	= FDRcorrection(Feat_pval,params['anova_alpha'],params['FDR_alpha'])
	outtitle1 	= params['output_file']+'Stats'
	FDR_df 		= (pd.DataFrame(FDRdict, index=['State_kw','p_value_kw','State_FDR','p_value_FDR'])).T
	FDR_df.to_csv(outtitle1,header=True, index=True, sep='\t')
	for a,FDR_ok in FDRdict.items(): 
		if (not FDR_ok[2] ):  #failed detect here #and not FDR_ok[0]
			if params['verbose']: print "\tkw ko" 
			del feats[a]
			wilcoxon_res[a] = "-" ###this value should be remained until the end; equal
			continue
		if params['verbose']: print "\tkw ok\t",
		kw_n_ok += 1 ###number of features which had p-value <= params set
		if not params['wilc']: continue
		
		res_wilcoxon_rep = test_rep_wilcoxon_r(subclass_sl,class_hierarchy,feats[a],params['wilcoxon_alpha'],params['multiclass_strat'],params['strict'],a,params['min_c'],params['only_same_subcl'],params['curv'])
		wilcoxon_resval[a] = str(FDR_ok[1])if res_wilcoxon_rep else "-"
		wilcoxon_resadj[a] = str(FDR_ok[3]) if res_wilcoxon_rep else "-"
		####the problem is here; every wilcoxon return TRUE
		if not res_wilcoxon_rep:
			if params['verbose']: print "wilc ko" 
			del feats[a]
		elif params['verbose']: print "wilc ok\t"

	if len(feats) > 0:
		print "Number of significantly discriminative features:", len(feats), "(", kw_n_ok, ") before internal wilcoxon"
		if params['lda_abs_th'] < 0.0:
			lda_res,lda_res_th = dict([(k,0.0) for k,v in feats.items()]), dict([(k,v) for k,v in feats.items()])
		else:
			if params['rank_tec'] == 'lda': lda_res,lda_res_th = test_lda_r(cls,feats,class_sl,params['n_boots'],params['f_boots'],params['lda_abs_th'],0.0000000001,params['nlogs'])
			elif params['rank_tec'] == 'svm': lda_res,lda_res_th = test_svm(cls,feats,class_sl,params['n_boots'],params['f_boots'],params['lda_abs_th'],0.0,params['svm_norm'])	
			else: lda_res,lda_res_th = dict([(k,0.0) for k,v in feats.items()]), dict([(k,v) for k,v in feats.items()])
	else: 
		print "Number of significantly discriminative features:", len(feats), "(", kw_n_ok, ") before internal wilcoxon"
		print "No features with significant differences between the two classes"
		lda_res,lda_res_th = {},{}
	###return p_adjust
	
	outres = {}
	outres['lda_res_th'] = lda_res_th
	outres['lda_res'] = lda_res
	outres['cls_means'] = cls_means
	outres['cls_means_kord'] = kord
	outres['wilcoxon_res'] = wilcoxon_res
	outres['wilcox_resval'] = wilcoxon_resval
	outres['wilcox_resadj'] = wilcoxon_resadj
	print "Number of discriminative features with abs LDA score >",params['lda_abs_th'],":",len(lda_res_th) 
	
	save_res_modif(outres,params["output_file"])

	#####
	###Additional corrected the p-value and write to file


	###



