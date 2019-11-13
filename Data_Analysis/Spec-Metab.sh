#!/bin/bash
#$ -S /bin/bash
#$ -o ./log/
#$ -e ./log/
#$ -r no

#created by Pande Putu Erawijantari - September, 25 2018
#Yamada Laboratory, School of Life Science and Technology
#Tokyo Institute of Technology
#erawijantari.p.aa@m.titech.ac.jp


source ~/.bashrc
source ~/.bash_profile
PATH=/home/pande/Downloaded_Program/MIMOSA/MIMOSA/:$PATH
export PATH


##select the intersecting ID
time python ~/gut_microbiome_analysis/script/ReRunOrganized/CoppyMergeModif/DropChange.py 190404_HUMAnN2_stratified.tsv 190404_StraHumann2 ~/gut_microbiome_analysis/gastrectomy_early/list_data/2018-10-12_sampleID_seqID.noexclude -d ~/gut_microbiome_analysis/gastrectomy_early/list_data/190402_IDmetabolomeGastrectomy

##convert the humann2 to PICRUST
Rscript ./convert_HUMANN2-PICRUST.R


#run Mimosa
Rscript runMimosa1.R --genefile="../input/190404_GenesKEGGHumann2.tsv" --metfile="../input/190404_MimosaMetabolome.tsv" --contribs_file="../input/190404_Humann2PicrustGene.tsv" --mapformula_file="/data/db/KEGG/CURRENT/ligand/reaction/reaction_mapformula.lst" --file_prefix="190407_"  --ko_rxn_file="/data/db/KEGG/CURRENT/genes/ko/ko_reaction.list" --rxn_annots_file="/data/db/KEGG/CURRENT/ligand/reaction/reaction" --metadata_file="../input/190404_IDMimosa" --metadata_var="Gastrectomy" -y --nonzero_filt=4 --taxonomy_file="../input/taxonomy_info2"

#qsub -cwd -l mem=20G Spec-Metab.sh