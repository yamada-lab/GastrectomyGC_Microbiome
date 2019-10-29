#!/bin/bash
#$ -S /bin/bash
#$ -o ./log/
#$ -e ./log/
#$ -r no


#=====LEFSE run for several posible combination====
#ERAWIJANTARI-YAMADA LABORATORY-2019


columns=$1 #columns in the reference for annotation
relab=$2 #relab table
ref=$3 #reference file
test=$4 #for title, what to test?e.g Healthy vs R-Y
title=$5 #
ID=$6
cut=$7
off=$8

echo ${columns}_${relab}_${test}_${title}

mkdir LEFSE_${test}
#prepare dataset --> make another step?
step1(){
echo "extracting table"
	#relable and cut off
time python ~/gut_microbiome_analysis/script/ReRunOrganized/CoppyMergeModif/DropChange.py $relab ${title}_${test} $ref -d $ID -c $columns
time python ~/gut_microbiome_analysis/script/ReRunOrganized/CoppyMergeModif/CutOff.py ${title}_${test}.Specify.tsv $cut $off $ref $columns

}

#LEFSE
step2(){
#run LEfSe
time python ~/Downloaded_Program/LEfse/format_input.py ${title}_${test}.Specify.${cut}.${off}.upLEfSe.txt ${title}_${test}.Specify.${cut}.${off}.upLEfSe.in -c 1 -u 2 -o 1000000
time python ~/gut_microbiome_analysis/script/ReRunOrganized/LEfSeModif/run_lefse_modif.py ${title}_${test}.Specify.${cut}.${off}.upLEfSe.in ${title}_${test}.Specify.${cut}.${off}.upLEfSe.res
time python ~/gut_microbiome_analysis/script/ReRunOrganized/StastViz/StatsCombine.py ${title}_${test}.Specify.${cut}.${off}.up.tsv ${title}_${test}_StatsCombine -l ${title}_${test}.Specify.${cut}.${off}.upLEfSe.res -s ${title}_${test}.Specify.${cut}.${off}.upLEfSe.resStats
time python ~/gut_microbiome_analysis/script/ReRunOrganized/LEfSeModif/plot_res.modif.py ${title}_${test}.Specify.${cut}.${off}.upLEfSe.res ${title}_${test}.Specify.${cut}.${off}.upLEfSe.svg --dpi 300 --format svg
mv *LEfSe* LEFSE_${test}
}

#plot cladogram
#LEFSE
step3(){
#run LEfSe
cd LEFSE_${test}
time python ~/gut_microbiome_analysis/script/ReRunOrganized/LEfSeModif/run_lefse_modif.py ${title}_${test}.Specify.${cut}.${off}.upLEfSe.in ${title}_${test}.Specify.${cut}.${off}.up3LEfSe.res -l 3
time python ~/gut_microbiome_analysis/script/ReRunOrganized/LEfSeModif/plot_cladogram.modif.py ${title}_${test}.Specify.${cut}.${off}.upLEfSe.res ${title}_${test}.Specify.${cut}.${off}.up3LEfSeClad.svg --dpi 300 --format svg
}

step1spec(){
echo "extracting table"
	#relable and cut off
time python ~/gut_microbiome_analysis/script/ReRunOrganized/CoppyMergeModif/DropChange.py $relab ${title}_${test} ~/gut_microbiome_analysis/gastrectomy_early/list_data/2018-10-12_sampleID_seqID.noexclude -d $ID -c $columns

}

#LEFSE
step2spec(){
#run LEfSe
time python ~/Downloaded_Program/LEfse/format_input.py ${title}_${test}.Specify.upLEfSe.txt ${title}_${test}.Specify.upLEfSe.in -c 1 -u 2 -o 1000000
time python ~/gut_microbiome_analysis/script/ReRunOrganized/LEfSeModif/run_lefse_modif.py ${title}_${test}.Specify.upLEfSe.in ${title}_${test}.Specify.upLEfSe.res
time python ~/gut_microbiome_analysis/script/ReRunOrganized/StastViz/StatsCombine.py ${title}_${test}.Specify.up.tsv ${title}_${test}_StatsCombine -l ${title}_${test}.Specify.upLEfSe.res -s ${title}_${test}.Specify.upLEfSe.resStats
time python ~/gut_microbiome_analysis/script/ReRunOrganized/LEfSeModif/plot_res.modif.py ${title}_${test}.Specify.upLEfSe.res ${title}_${test}.Specify.upLEfSe.svg --dpi 300 --format svg
mv *LEfSe* LEFSE_${test}
}


#plot cladogram
#LEFSE
step3spec(){
#run LEfSe
cd LEFSE_${test}
time python ~/gut_microbiome_analysis/script/ReRunOrganized/LEfSeModif/run_lefse_modif.py ${title}_${test}.Specify.upLEfSe.in ${title}_${test}.Specify.up3LEfSe.res -l 3
time python ~/gut_microbiome_analysis/script/ReRunOrganized/LEfSeModif/plot_cladogram.modif.py ${title}_${test}.Specify.upLEfSe.res ${title}_${test}.Specify.up3LEfSeClad.svg --dpi 300 --format svg
}



if [ $# -lt 8 ];
then
	echo "without CutOff"
	step1spec
	step2spec
	step3spec
else
	echo "With cutoff"
	step1
	step2
	step3
fi