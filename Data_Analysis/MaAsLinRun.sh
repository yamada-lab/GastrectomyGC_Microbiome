#!/bin/bash
#$ -S /bin/bash
#$ -o ./log/
#$ -e ./log/
#$ -r no


#=====MaAsLin run for several posible combination====
#ERAWIJANTARI-YAMADA LABORATORY-2019


relab=$1 #relative abundance table
metadata=$2 #metadata table
title=$3 #for title
referenceID=$4
cwd=$(pwd)
statscombine=$5
adjustfile=$6

#dir="MaAslin_Adjusted"
#mkdir $dir
#cd $dir

step1(){
	#Data preparations
time python ~/gut_microbiome_analysis/script/ReRunOrganized/StastViz/MaAsLIn_formated.py $relab $metadata $title $referenceID
time python $MaAsLin/exec/merge_metadata.py $title.Metadata.tsv < $relab > $relab.MaAsLin.pcl
time python $MaAsin/exec/transpose.py < $relab.MaAsLin.pcl >  $relab.MaAsLin.tsv
}

step2(){
#adjustment for all factors
Rscript $Maaslin/R/Maaslin.R -i $title.read.config -r 0.0 -p 0.05 $relab.MaAsLin.tsv ${title}_All
#adjustment for selected factors in manuscripts
Rscript $Maaslin/R/Maaslin.R -i ${cwd}/$title.read.config -F "BMI,Total_Cholesterol,DiabetesMed,Gastric_acid_medication,Age,Gender" -r 0.0 -p 0.05 ${cwd}/$relab.MaAsLin.tsv ${cwd}/${title}_adjAll
}

step3(){
#ploting adjusted and crude coeffient
#adjusted
mkdir MaAsLin_Adjusted/confound_factor
while read p; do
		time python ~/gut_microbiome_analysis/script/ReRunOrganized/StastViz/ConfoundVarExp.py ${dir}/${title}_All/${relab}-Status.txt ${dir}/${title}_adj.${p}/${relab}-Status.txt confound_factor/deconfound.${p}.pdf -r $relab -s $statscombine
		time python ~/gut_microbiome_analysis/script/ReRunOrganized/StastViz/ConfoundVarExp.py ${dir}/${title}_All/${relab}-Status.txt ${dir}/${title}_All/${relab}-${BMI}.txt confound_factor/${p}.pdf -r $relab -s $statscombine -a "Yes"
	done < $adjustfile
}


step1
step2
step3


#example
#MetaPhlan2
#bash ~/gut_microbiome_analysis/script/ReRunOrganized/StastViz/MaAsLinRun.sh 2018-12-07_Metaphlan2_Species.5.0.001.up.tsv ~/gut_microbiome_analysis/gastrectomy_early/181103_2Group/181109_Metadata/MetaDisBlood/190822_MaAsLinMetadata.tsv 190822MetaPhlan2 ~/gut_microbiome_analysis/gastrectomy_early/list_data/2018-10-21_sampleID_seqIDNoEMR ../190111_Stats_Combine_SpeciesMetaphlan2.tsv 

#mOTU
#bash ~/gut_microbiome_analysis/script/ReRunOrganized/StastViz/MaAsLinRun.sh 2018-11-03_NCBI.species.abundances.5.1e-05.up.tsv ~/gut_microbiome_analysis/gastrectomy_early/181103_2Group/181109_Metadata/MetaDisBlood/190822_MaAsLinMetadata.tsv 190822mOTU ~/gut_microbiome_analysis/gastrectomy_early/list_data/2018-10-21_sampleID_seqIDNoEMR ../190111_Stats_Combine_SpeciesNCBI.tsv

#Humann2
#bash ~/gut_microbiome_analysis/script/ReRunOrganized/StastViz/MaAsLinRun.sh 181218_Modules.Specifiy.5.1e-06.up.tsv ~/gut_microbiome_analysis/gastrectomy_early/181103_2Group/181109_Metadata/MetaDisBlood/190822_MaAsLinMetadata.tsv 190822Humann2 ~/gut_microbiome_analysis/gastrectomy_early/list_data/2018-10-21_sampleID_seqIDNoEMR ../190111_Stats_Combine_ModuleUniref.tsv


#KEGG
#bash ~/gut_microbiome_analysis/script/ReRunOrganized/StastViz/MaAsLinRun.sh 190117_KEGG_modules.5.1e-06.up.tsv ~/gut_microbiome_analysis/gastrectomy_early/181103_2Group/181109_Metadata/MetaDisBlood/190822_MaAsLinMetadata.tsv 190822KEGG ~/gut_microbiome_analysis/gastrectomy_early/list_data/2018-10-21_sampleID_seqIDNoEMR ../190111_Stats_Combine_ModuleKEGG.tsv

#Metabolome
#bash ~/gut_microbiome_analysis/script/ReRunOrganized/StastViz/MaAsLinRun.sh 190402_Metabolome.4.1e-06.upEdited.tsv ~/gut_microbiome_analysis/gastrectomy_early/181103_2Group/181109_Metadata/MetaDisBlood/190822_MaAsLinMetadata.tsv 190822Metabolome ~/gut_microbiome_analysis/gastrectomy_early/list_data/2018-10-21_sampleID_seqIDNoEMR ../190402_StatsCombineMetabolome.tsv