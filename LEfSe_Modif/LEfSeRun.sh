#!/bin/bash
#$ -S /bin/bash
#$ -o ./log/
#$ -e ./log/
#$ -r no


#=====LEFSE run for several posible combination====
#ERAWIJANTARI-YAMADA LABORATORY-2019


columns=$1 #columns in the reference for annotation
relab=$2 #relab table
test=$3 #for title, what to test?e.g Healthy vs R-Y
title=$4 #

echo ${columns}_${relab}_${test}_${title}

mkdir LEFSE_${test}
cd LEFSE_${test}
#prepare dataset
step1(){
echo "extracting table"
#relable
time python DropChange.py $relab ${title}_${test} ~/gut_microbiome_analysis/gastrectomy_early/list_data/190807_sampleIDseqIDPolyps -d ~/gut_microbiome_analysis/gastrectomy_early/list_data/190327_IDnoEMRnoNA -c $columns
}

step_metabolite(){
echo "extracting table"
#prepare for the reconstruction
time python DropChange.py $relab ${title}_${test} ~/gut_microbiome_analysis/gastrectomy_early/list_data/190807_sampleIDseqIDPolyps -d ~/gut_microbiome_analysis/gastrectomy_early/list_data/190402_IDmetabolomeGastrectomynoNA -c $columns
}

#LEFSE
step2(){
#run LEfSe
time python format_input.py ${title}_${test}.Specify.upLEfSe.txt ${title}_${test}.Specify.upLEfSe.in -c 1 -u 2 -o 1000000
time python run_lefse_modif.py ${title}_${test}.Specify.upLEfSe.in ${title}_${test}.Specify.upLEfSe.res
time python StatsCombine.py ${title}_${test}.Specify.tsv ${title}_${test}_StatsCombine -l ${title}_${test}.Specify.upLEfSe.res -s ${title}_${test}.Specify.upLEfSe.resStats

}


step1
#only activate for metabolite data
#step_metabolite
step2