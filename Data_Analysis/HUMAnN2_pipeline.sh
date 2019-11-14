#!/bin/bash
#$ -S /bin/bash
#$ -o ./log/
#$ -e ./log/
#$ -r no

#created by Pande Putu Erawijantari - September, 25 2018
#Yamada Laboratory, School of Life Science and Technology
#Tokyo Institute of Technology


#STEP
#1. run the humann2 command
#2. get the relative abundance
#3. map to KO
#4. remove the intermediate files

echo "$(date "+%m%d%Y %T") : Starting analysis"

source ~/.bashrc
source ~/.bash_profile
PATH=/share1/watanabe/tool/diamond-v0.8.22-linux64/:$PATH  #diamond needed
export PATH


directory_in=$3 	#directory for the high quality reads; example: /data/project/gut_cancer/runs/processed
directory_out=$2 	#directory for the output; example: /home/pande/gut_microbiome_analysis/gastrectomy_early/Seq_data
sequence_name=$1    #sequence header/ID
run=$(pwd)

#1. run the humann2 command
step1(){
echo "$(date "+%m%d%Y %T") :run the humann2 command"
humann2 --input $directory_in/$sequence_name/QC/$sequence_name.pairend.fastq.gz --output $directory_out/$sequence_name/$sequence_name.humann2 --search-mode uniref90
}

#2. get the relative abundance
step2(){
echo "$(date "+%m%d%Y %T") :get the relative abundance"
humann2_renorm_table --input $directory_out/$sequence_name/$sequence_name.humann2/*genefamilies.tsv --output $directory_out/$sequence_name/$sequence_name.humann2/$sequence_name.pairend_genefamilies_relab.tsv --units relab
humann2_genefamilies_genus_level --input $directory_out/$sequence_name/$sequence_name.humann2/$sequence_name.pairend_genefamilies_relab.tsv --output $directory_out/$sequence_name/$sequence_name.humann2/$sequence_name.pairend_genefamilies_relab_genus_level.tsv
}

#3. map to KO
step3(){
echo "$(date "+%m%d%Y %T") :map to KO"
humann2_regroup_table --input $directory_out/$sequence_name/$sequence_name.humann2/$sequence_name.pairend_genefamilies.tsv --groups uniref90_ko --output $directory_out/$sequence_name/$sequence_name.humann2/$sequence_name.pairend_genefamilies_ko.tsv 
humann2_renorm_table --input $directory_out/$sequence_name/$sequence_name.humann2/$sequence_name.pairend_genefamilies_ko.tsv  --output $directory_out/$sequence_name/$sequence_name.humann2/$sequence_name.pairend_genefamilies_relab_ko.tsv --units relab
humann2_genefamilies_genus_level --input $directory_out/$sequence_name/$sequence_name.humann2/$sequence_name.pairend_genefamilies_relab_ko.tsv --output $directory_out/$sequence_name/$sequence_name.humann2/$sequence_name.pairend_genefamilies_relab_ko_genus_level.tsv
}


#4. remove the intermediate files
step4(){
echo "$(date "+%m%d%Y %T") :removes the intermediate files"
rm -r $directory_out/$sequence_name/$sequence_name.humann2/$sequence_name.pairend_humann2_temp
}


#modif (for modifications/customize step, could be ommited)
step5(){
#cp $directory_out/$sequence_name/$sequence_name.humann2/$sequence_name.pairend_genefamilies_relab_ko.tsv $directory_out/$sequence_name/$sequence_name.pairend_genefamilies_relab_ko.tsv
#cp $directory_out/$sequence_name/$sequence_name.humann2/$sequence_name.pairend_genefamilies_relab_ko_genus_level.tsv $directory_out/$sequence_name/$sequence_name.pairend_genefamilies_relab_ko_genus_level.tsv
#cp $directory_out/$sequence_name/$sequence_name.humann2/$sequence_name.pairend_genefamilies_ko.tsv $directory_out/$sequence_name/$sequence_name.pairend_genefamilies_ko.tsv
#cp $directory_out/$sequence_name/$sequence_name.humann2/$sequence_name.pairend_humann2_temp/$sequence_name.pairend_metaphlan_bugs_list.tsv $directory_out/$sequence_name/$sequence_name.pairend_metaphlan_bugs_list.tsv
#cp $directory_out/$sequence_name/$sequence_name.humann2/$sequence_name.pairend_genefamilies_relab.tsv $directory_out/$sequence_name/$sequence_name.pairend_genefamilies_relab.tsv
}


step1
step2
step3 #this step also can be performed after combining each sample results into one table
step4
#step5



#example usage in SGE systme: qsub -l mem=4G -pe smp 8 humman2_pipeline.sh <sequence_id> <directory_out> <directory_in>


#aditional notes
##spec only 
#grep -E "(s__)|(^ID)" merged_abundance_table.txt | grep -v "t__" | sed 's/^.*s__//g' > merged_abundance_table_species.txt
##genus only
#grep -E "(g__)|(^ID)" 2018-12-07_pairend_metaphlan_bugs_list.tsv.tsv | grep -v "s__" | sed 's/^.*g__//g' > 2018-12-07_Metaphlan2_Genus.tsv
#grep -E for the level to retrieved; grep -v to ignore printing the level after; sed writing from that specific level
##KO only
#humann2_split_stratified_table --input $directory_out/$sequence_name/$sequence_name.humann2/$sequence_name.pairend_genefamilies_relab_ko.tsv --output $directory_out/$sequence_name.humann2/

#merge output
#merge_metaphlan_tables.py *.siraeum.txt > siraeum_tracker.txt


