#!/bin/bash
#$ -S /bin/bash
#$ -o ./log/
#$ -e ./log/
#$ -r no


#=====LEFSE run for several posible combination====
#ERAWIJANTARI-YAMADA LABORATORY-2019


columns=$1 #columns in the reference for annotation; e. g. columns 1 in metadata contains informations about surgery status
relab=$2 # raw relative abundance table
ref=$3 #reference file 
test=$4 #for title, what to test?e.g Healthy vs R-Y
title=$5 #
ID=$6
cut=$7
off=$8

echo ${columns}_${relab}_${test}_${title}

mkdir LEFSE_${test}
step1(){
echo "extracting table"
	#relable and cut off
time python ./LEfSe_Modif/DropChange.py $relab ${title}_${test} $ref -d $ID -c $columns
time python ./LEfSe_Modif/CutOff.py ${title}_${test}.Specify.tsv $cut $off $ref $columns

}

#LEFSE
step2(){
#run LEfSe
time python ./LEfSe_Modif/format_input.py ${title}_${test}.Specify.${cut}.${off}.upLEfSe.txt ${title}_${test}.Specify.${cut}.${off}.upLEfSe.in -c 1 -u 2 -o 1000000
time python ./LEfSe_Modif/run_lefse_modif.py ${title}_${test}.Specify.${cut}.${off}.upLEfSe.in ${title}_${test}.Specify.${cut}.${off}.upLEfSe.res
time python ./LEfSe_Modif/StatsCombine.py ${title}_${test}.Specify.${cut}.${off}.up.tsv ${title}_${test}_StatsCombine -l ${title}_${test}.Specify.${cut}.${off}.upLEfSe.res -s ${title}_${test}.Specify.${cut}.${off}.upLEfSe.resStats
time python ./LEfSe_Modif/plot_res.modif.py ${title}_${test}.Specify.${cut}.${off}.upLEfSe.res ${title}_${test}.Specify.${cut}.${off}.upLEfSe.svg --dpi 300 --format svg
mv *LEfSe* LEFSE_${test}
}

#plot cladogram --> optional,this step will throw error if the taxon level are not fully provided
#LEFSE
step3(){
#run LEfSe
cd LEFSE_${test}
time python ./LEfSe_Modif/run_lefse_modif.py ${title}_${test}.Specify.${cut}.${off}.upLEfSe.in ${title}_${test}.Specify.${cut}.${off}.up3LEfSe.res -l 3
time python ./LEfSe_Modif/plot_cladogram.modif.py ${title}_${test}.Specify.${cut}.${off}.upLEfSe.res ${title}_${test}.Specify.${cut}.${off}.up3LEfSeClad.svg --dpi 300 --format svg
}

step1spec(){
echo "extracting table"
	#relable and cut off
time python ./LEfSe_Modif/DropChange.py $relab ${title}_${test} ~/gut_microbiome_analysis/gastrectomy_early/list_data/2018-10-12_sampleID_seqID.noexclude -d $ID -c $columns
}

#LEFSE
step2spec(){
#run LEfSe
time python ~/Downloaded_Program/LEfse/format_input.py ${title}_${test}.Specify.upLEfSe.txt ${title}_${test}.Specify.upLEfSe.in -c 1 -u 2 -o 1000000
time python ./LEfSe_Modif/run_lefse_modif.py ${title}_${test}.Specify.upLEfSe.in ${title}_${test}.Specify.upLEfSe.res
time python ./LEfSe_Modif/StatsCombine.py ${title}_${test}.Specify.up.tsv ${title}_${test}_StatsCombine -l ${title}_${test}.Specify.upLEfSe.res -s ${title}_${test}.Specify.upLEfSe.resStats
time python ./LEfSe_Modif/plot_res.modif.py ${title}_${test}.Specify.upLEfSe.res ${title}_${test}.Specify.upLEfSe.svg --dpi 300 --format svg
mv *LEfSe* LEFSE_${test}
}


#plot cladogram --> optional, this step will throw error if the taxon level are not fully provided
#LEFSE
step3spec(){
#run LEfSe
cd LEFSE_${test}
time python ./LEfSe_Modif/run_lefse_modif.py ${title}_${test}.Specify.upLEfSe.in ${title}_${test}.Specify.up3LEfSe.res -l 3
time python ./LEfSe_Modif/plot_cladogram.modif.py ${title}_${test}.Specify.upLEfSe.res ${title}_${test}.Specify.up3LEfSeClad.svg --dpi 300 --format svg
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