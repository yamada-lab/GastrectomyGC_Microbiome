#!/bin/bash
#$ -S /bin/bash
#$ -o ./log/
#$ -e ./log/
#$ -r no

source ~/.bashrc
source ~/.bash_profile


relab_input=$1
analysis=$2 #example: SpeciesRelabmOTU; GenusCountMetaPhlan2
list_sig=$3
date=`date +%Y-%m-%d`



step1(){
echo "parse the dataframe based on the grouping"
time python ~/gut_microbiome_analysis/script/ReRunOrganized/StastViz/Parse_dfGroupList.py $relab_input Healthy/${date}_${analysis} -g Healthy -d $list_sig
time python ~/gut_microbiome_analysis/script/ReRunOrganized/StastViz/Parse_dfGroupList.py $relab_input Gastrectomy/${date}_${analysis} -g Gastrectomy -d $list_sig
}

step2(){
echo "SparCC Healthy"
source activate sparcc
time python ~/Downloaded_Program/SparCC/SparCC.py Healthy/${date}_${analysis}.Healthy.tsv -c Healthy/${date}_${analysis}.Healthy.sparcc.tsv -v Healthy/${date}_${analysis}.Healthy.sparcc.cov.tsv
python ~/Downloaded_Program/SparCC/MakeBootstraps.py -n 5000 Healthy/${date}_${analysis}.Healthy.tsv -p Healthy/${analysis}_p/ -t permuted_#
mkdir Healthy/${analysis}_boot_corr
mkdir Healthy/${analysis}_boot_cov
for i in `seq 0 4999`; do python ~/Downloaded_Program/SparCC/SparCC.py Healthy/${analysis}_p/permuted_$i -c Healthy/${analysis}_boot_corr/simulated_sparcc_$i.txt -v Healthy/${analysis}_boot_cov/simulated_sparcc_$i.txt >> Healthy/${analysis}_boot_sparcc.log; done
mkdir Healthy/${analysis}_pvals
python ~/Downloaded_Program/SparCC/PseudoPvals.py Healthy/${date}_${analysis}.Healthy.sparcc.tsv Healthy/${analysis}_boot_corr/simulated_sparcc_#.txt 5000 -o Healthy/${analysis}_pvals/two_sided5000.txt -t two_sided
}


step3(){
echo "SparCC gastrectomy"
source activate sparcc
time python ~/Downloaded_Program/SparCC/SparCC.py Gastrectomy/${date}_${analysis}.Gastrectomy.tsv -c Gastrectomy/${date}_${analysis}.Gastrectomy.sparcc.tsv -v Gastrectomy/${date}_${analysis}.Gastrectomy.sparcc.cov.tsv
python ~/Downloaded_Program/SparCC/MakeBootstraps.py -n 5000 Gastrectomy/${date}_${analysis}.Gastrectomy.tsv -p Gastrectomy/${analysis}_p/ -t permuted_#
mkdir Gastrectomy/${analysis}_boot_corr
mkdir Gastrectomy/${analysis}_boot_cov
for i in `seq 0 4999`; do python ~/Downloaded_Program/SparCC/SparCC.py Gastrectomy/${analysis}_p/permuted_$i -c Gastrectomy/${analysis}_boot_corr/simulated_sparcc_$i.txt -v Gastrectomy/${analysis}_boot_cov/simulated_sparcc_$i.txt >> Gastrectomy/${analysis}_boot_sparcc.log; done
mkdir Gastrectomy/${analysis}_pvals
python ~/Downloaded_Program/SparCC/PseudoPvals.py Gastrectomy/${date}_${analysis}.Gastrectomy.sparcc.tsv Gastrectomy/Gastrectomy/${analysis}_boot_corr/simulated_sparcc_#.txt 5000 -o Gastrectomy/${analysis}_pvals/two_sided5000.txt -t two_sided
}

step4(){
echo "calculate the pvals"
#mkdir Healthy/${analysis}_pvals
python ~/Downloaded_Program/SparCC/PseudoPvals.py Healthy/2019-06-28_${analysis}.Healthy.sparcc.tsv Healthy/${analysis}_boot_corr/simulated_sparcc_#.txt 5000 -o Healthy/${analysis}_pvals/two_sided5000.txt -t two_sided
#mkdir Gastrectomy/${analysis}_pvals
python ~/Downloaded_Program/SparCC/PseudoPvals.py Gastrectomy/2019-06-28_${analysis}.Gastrectomy.sparcc.tsv Gastrectomy/${analysis}_boot_corr/simulated_sparcc_#.txt 5000 -o Gastrectomy/${analysis}_pvals/two_sided5000.txt -t two_sided
}


#call
step1
step2
step3
step4