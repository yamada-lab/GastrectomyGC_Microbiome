#!/bin/bash
#$ -S /bin/bash
#$ -o ./log/
#$ -e ./log/
#$ -r no


#=====THIS CODE PERFORMED THE READS FILTERING FOR THE NEXT ANALYSIS====
# input can be compressed or uncompressed file
#NOVEMBER 2016
#1. N filtering
#2. Phix filtering
#3. cutadapt
#4.	length 	filtering (pick sequence which has more than 50bp length(>50)),\
#5.	filtering (pick sequence which has the more than or equal 25 of average quality)
#6. Human genome elimination
#7. combine file and non pair end elimination
#8. Converting fastq file into fasta file

#USER INPUT NEEDED
echo "$(date "+%m%d%Y %T") : Starting analysis" 

if [ $# -lt 4 ];
then
	echo "This program default set is to analyze all of data listing in the list file"
	echo "the first column in list should be the name of sequence"
	echo "ussage hq_reads_filtering.sh keyword file_output_directory raw_data_directory yes/no"
fi

directory_in=$3
directory_out=$2
sequence_name=$1
run=$(pwd)
remove=$4

# copy the sequence to the results folder for analysis 
step0(){
save=$directory_out/$sequence_name	
mkdir -p $save
cp $directory_in/$sequence_name/*1.fastq $save/$sequence_name.1.fq || cp $directory_in/$sequence_name/*1.fastq.gz $save/$sequence_name.1.fq.gz
cp $directory_in/$sequence_name/*2.fastq $save/$sequence_name.2.fq || cp $directory_in/$sequence_name/*2.fastq.gz $save/$sequence_name.2.fq.gz
}

#N-contained sequence filter
##{sequence_name}*{1/2}.fq /fq.gz--> ${sequence_name}*{1/2}.N.filter.fq/fq.gz
step1(){
echo "analysis for the $sequence_name" 
echo "check the number of sequence before analysis forward and reverse, respectivelly"
wc -l ${directory_out:2}/$sequence_name/*1.f*
wc -l ${directory_out:2}/$sequence_name/*2.f*
cd $run
time python /home/pande/gut_microbiome_analysis/script/hq_reads/hq_N.filter.py ${directory_out:2}/$sequence_name/*1.f* &
time python /home/pande/gut_microbiome_analysis/script/hq_reads/hq_N.filter.py ${directory_out:2}/$sequence_name/*2.f* &
wait
echo "A. after N filtering forward and reverse, respectivelly"
wc -l ${directory_out:2}/$sequence_name/*1.N.filter*
wc -l ${directory_out:2}/$sequence_name/*2.N.filter*
}

#Phix contaminantion filter out
##{sequence_name}.1.N.filter.fq* --> ${sequence_name}.{1/2}.phixelim.fq
step2a(){
echo "Phix filtering"
/home/pande/gut_microbiome_analysis/program/bowtie2-2.2.9/bowtie2 --no-hd --no-sq --fast-local -x /home/pande/gut_microbiome_analysis/DB/phix174/phiX174.fasta.index -U $directory_out/$sequence_name/*1.N.filter.fq -S $directory_out/$sequence_name/1.sam --un $directory_out/$sequence_name/${sequence_name}.1.phixelim.fq &
/home/pande/gut_microbiome_analysis/program/bowtie2-2.2.9/bowtie2 --no-hd --no-sq --fast-local -x /home/pande/gut_microbiome_analysis/DB/phix174/phiX174.fasta.index -U $directory_out/$sequence_name/*2.N.filter.fq -S $directory_out/$sequence_name/2.sam --un $directory_out/$sequence_name/${sequence_name}.2.phixelim.fq &
wait
echo "B. after PhiX filtering"
wc -l $directory_out/$sequence_name/${sequence_name}.1.phixelim.fq
wc -l $directory_out/$sequence_name/${sequence_name}.2.phixelim.fq
}

step2b(){
echo "Phix filtering"
gunzip -c ${directory_out:2}/$sequence_name/*1.N.filter.fq.gz > ${directory_out:2}/$sequence_name/*1.N.filter.fq
gunzip -c ${directory_out:2}/$sequence_name/*2.N.filter.fq.gz > ${directory_out:2}/$sequence_name/*2.N.filter.fq
/home/pande/gut_microbiome_analysis/program/bowtie2-2.2.9/bowtie2 --no-hd --no-sq --fast-local -x /home/pande/gut_microbiome_analysis/DB/phix174/phiX174.fasta.index -U $directory_out/$sequence_name/*1.N.filter.fq -S $directory_out/$sequence_name/1.sam --un $directory_out/$sequence_name/${sequence_name}.1.phixelim.fq &
/home/pande/gut_microbiome_analysis/program/bowtie2-2.2.9/bowtie2 --no-hd --no-sq --fast-local -x //home/pande/gut_microbiome_analysis/DB/phix174/phiX174.fasta.index -U $directory_out/$sequence_name/*2.N.filter.fq -S $directory_out/$sequence_name/2.sam --un $directory_out/$sequence_name/${sequence_name}.2.phixelim.fq &
wait
echo "B. after PhiX filtering"
wc -l $directory_out/$sequence_name/${sequence_name}.1.phixelim.fq
wc -l $directory_out/$sequence_name/${sequence_name}.2.phixelim.fq
}

trim_result1="$directory_out/$sequence_name/${sequence_name}.1.trim.fq"
trim_result2="$directory_out/$sequence_name/${sequence_name}.2.trim.fq"

#cutadapt for Adapter and 3’ end low quality reads trimming
#instal python setup.py install
##{sequence_name}.{1/2}.phixelim.fq --> ${sequence_name}.{1/2}.trim.fq
step3(){
echo "cutadapt for sequence trimming"
/home/pande/.local/bin/cutadapt -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -O 33 -o $trim_result1 -q 17 -f fastq "$directory_out/$sequence_name/${sequence_name}.1.phixelim.fq" &
/home/pande/.local/bin/cutadapt -a CTGTCTCTTATACACATCTGACGCTGCCGACGA -O 32 -o $trim_result2 -q 17 -f fastq "$directory_out/$sequence_name/${sequence_name}.2.phixelim.fq" &
wait
echo "C. the number of sequence after cutadapt forward and reverse" 
wc -l $trim_result1 
wc -l $trim_result2
}

#reads filtering based on the length that should be more than 50 base
##{sequence_name}.1/2}.trim.fq --> ${sequence_name}.{1/2}.length.filter.fq
step4(){
echo "filtering length sequence" 
time python /home/pande/gut_microbiome_analysis/script/hq_reads/hq_length_filter.py ${trim_result1:2} & 
time python /home/pande/gut_microbiome_analysis/script/hq_reads/hq_length_filter.py ${trim_result2:2} & 
wait
echo "D. the number of sequence after length filtering forward and reverse"  
wc -l $directory_out/$sequence_name/${sequence_name}.1.length.filter.*
wc -l $directory_out/$sequence_name/${sequence_name}.2.length.filter.*
}

#reads filtering based on the length that should be more than 50 base
##{sequence_name}.{1/2}.length.filter.fq --> ${sequence_name}.{1/2}.qual.filter.fq
step5a(){
echo "filtering high quality sequence" 
time python /home/pande/gut_microbiome_analysis/script/hq_reads/hq_qual_filter.py ${directory_out:2}/$sequence_name/${sequence_name}.1.length.filter* /home/pande/gut_microbiome_analysis/DB/ascii & 
time python /home/pande/gut_microbiome_analysis/script/hq_reads/hq_qual_filter.py ${directory_out:2}/$sequence_name/${sequence_name}.2.length.filter* /home/pande/gut_microbiome_analysis/DB/ascii & 
wait
echo "D. the number of sequence after quality filtering forward and reverse" 
wc -l $directory_out/$sequence_name/${sequence_name}.1.qual.filter*
wc -l $directory_out/$sequence_name/${sequence_name}.2.qual.filter*
}
step5b(){
echo "filtering high quality sequence" 
time python /home/pande/gut_microbiome_analysis/script/hq_reads/hq_qual_filter.py ${directory_out:2}/$sequence_name/${sequence_name}.1.length.filter* /home/pande/gut_microbiome_analysis/DB/ascii & 
time python /home/pande/gut_microbiome_analysis/script/hq_reads/hq_qual_filter.py ${directory_out:2}/$sequence_name/${sequence_name}.2.length.filter* /home/pande/gut_microbiome_analysis/DB/ascii & 
wait
echo "D. the number of sequence after quality filtering forward and reverse" 
wc -l $directory_out/$sequence_name/${sequence_name}.1.qual.filter*
wc -l $directory_out/$sequence_name/${sequence_name}.2.qual.filter*
gzip $directory_out/$sequence_name/${sequence_name}.1.qual.filter.fq
gzip $directory_out/$sequence_name/${sequence_name}.2.qual.filter.fq
}

#human genome elemination
##{sequence_name}.{1/2}.qual.filter* --> ${sequence_name}.{1/2}.humgenElim.fq
step6a(){
echo "human genome elimination"
/home/pande/gut_microbiome_analysis/program/bowtie2-2.2.9/bowtie2 --no-hd --no-sq --no-unal --fast-local -x /home/pande/gut_microbiome_analysis/DB/Human_Genome/human_genome.fa.index -U $directory_out/$sequence_name/${sequence_name}.1.qual.filter.fq -S $directory_out/$sequence_name/1human.sam --un $directory_out/$sequence_name/${sequence_name}.1.humgenElim.fq &
/home/pande/gut_microbiome_analysis/program/bowtie2-2.2.9/bowtie2 --no-hd --no-sq --no-unal --fast-local -x /home/pande/gut_microbiome_analysis/DB/Human_Genome/human_genome.fa.index -U $directory_out/$sequence_name/${sequence_name}.2.qual.filter.fq -S $directory_out/$sequence_name/2human.sam --un $directory_out/$sequence_name/${sequence_name}.2.humgenElim.fq &
wait
echo "E. after human genome elimination"
wc -l $directory_out/$sequence_name/${sequence_name}.1.humgenElim.fq
wc -l $directory_out/$sequence_name/${sequence_name}.2.humgenElim.fq
}

step6b(){
echo "human genome elimination"
gunzip -c ${directory_out:2}/$sequence_name/*1.qual.filter.fq.gz > ${directory_out:2}/$sequence_name/*1.qual.filter.fq
gunzip -c ${directory_out:2}/$sequence_name/*2.qual.filter.fq.gz > ${directory_out:2}/$sequence_name/*2.qual.filter.fq
/home/pande/gut_microbiome_analysis/program/bowtie2-2.2.9/bowtie2 --no-hd --no-sq --no-unal --fast-local -x /home/pande/gut_microbiome_analysis/DB/Human_Genome/human_genome.fa.index -U $directory_out/$sequence_name/${sequence_name}.1.qual.filter.fq -S $directory_out/$sequence_name/1human.sam --un $directory_out/$sequence_name/${sequence_name}.1.humgenElim.fq &
/home/pande/gut_microbiome_analysis/program/bowtie2-2.2.9/bowtie2 --no-hd --no-sq --no-unal --fast-local -x /home/pande/gut_microbiome_analysis/DB/Human_Genome/human_genome.fa.index -U $directory_out/$sequence_name/${sequence_name}.2.qual.filter.fq -S $directory_out/$sequence_name/2human.sam --un $directory_out/$sequence_name/${sequence_name}.2.humgenElim.fq &
wait
echo "E. after human genome elimination"
wc -l $directory_out/$sequence_name/${sequence_name}.1.humgenElim.fq
wc -l $directory_out/$sequence_name/${sequence_name}.2.humgenElim.fq
}

#combining paired end reads
##{sequence_name}.{1/2}.humgenElim.fq --> ${sequence_name}.pairend.fastq
step7(){
echo "combining the paired end reads" 
time python /home/pande/gut_microbiome_analysis/script/hq_reads/hq_fastqCombPairEnd.py "${directory_out:2}/$sequence_name/${sequence_name}.1.humgenElim.fq" "${directory_out:2}/$sequence_name/${sequence_name}.2.humgenElim.fq"
echo "F. number of paired end reads" 
wc -l $directory_out/$sequence_name/${sequence_name}.pairend.fq
#coppy this file to enable the process in MOCAT
gzip -c $directory_out/$sequence_name/${sequence_name}.pairend.fq > $directory_out/$sequence_name/${sequence_name}.fq.gz || echo 'no file to compress'
}

#converting fastq file to fasta file
##{sequence_name}.pairend.fq --> $sequence_name/$sequence_name.hiseq.final.fasta
step8(){
echo "generating fasta file" 
awk 'NR%4 == 1 {print ">" substr($0, 2)} NR%4 == 2 {print}' "$directory_out/$sequence_name/${sequence_name}.pairend.fq" > "$directory_out/$sequence_name/${sequence_name}.hiseq.final.fasta"
echo "G. final number of sequence" 
wc -l $directory_out/$sequence_name/$sequence_name.hiseq.final.fasta 
}

step9(){
echo "Compressing the file"
rm $directory_out/$sequence_name/$sequence_name.1.fq 
rm $directory_out/$sequence_name/$sequence_name.2.fq
rm $directory_out/$sequence_name/*sam 
rm $directory_out/$sequence_name/$sequence_name.pairend.fq
gzip $directory_out/$sequence_name/$sequence_name.hiseq.final.fasta
gzip ${directory_out:2}/$sequence_name/*1.N.filter*
gzip ${directory_out:2}/$sequence_name/*2.N.filter*
gzip $directory_out/$sequence_name/${sequence_name}.1.phixelim.fq
gzip $directory_out/$sequence_name/${sequence_name}.2.phixelim.fq
gzip $trim_result1 
gzip $trim_result2
gzip $directory_out/$sequence_name/${sequence_name}.1.length.filter.*
gzip $directory_out/$sequence_name/${sequence_name}.2.length.filter.*
gzip $directory_out/$sequence_name/${sequence_name}.1.qual.filter*
gzip $directory_out/$sequence_name/${sequence_name}.2.qual.filter*
gzip $directory_out/$sequence_name/${sequence_name}.1.humgenElim.fq
gzip $directory_out/$sequence_name/${sequence_name}.2.humgenElim.fq
}

#call step
step0
processed_file=$directory_out/$sequence_name/*
if [[ $processed_file =~ \.gz ]];
	then
		echo "compressed file processing"
		#step0
		step1
		step2b
		step3
		step4
		step5b
		step6b
		step7
		step8
else
		echo "fastq file processing"
		#step0
		step1
		step2a
		step3
		step4
		step5a
		step6a
		step7
		step8
fi

if [[ $remove =~ "yes" ]];
	then
	echo "removing the intermediate file"
	#removing all of the intermediate file
	rm -f $directory_out/$sequence_name/${sequence_name}*.qual.filter* $directory_out/$sequence_name/${sequence_name}*.length.filter*
	rm -f $directory_out/$sequence_name/${sequence_name}*.trim.* $directory_out/$sequence_name/${sequence_name}*.phixelim.* $directory_out/$sequence_name/${sequence_name}*.N.filter*
	echo "$(date "+%m%d%Y %T") : Finished $sequence_name" 
else
	step9
	echo "$(date "+%m%d%Y %T") : Finished $sequence_name"
fi


##Feature need to be added
# configuration file contain the set up for
##1 installed program location
##2 script location (or symlink generator for script)
##3 parameter setting for each process (e.g the minimum length, quality, bowtie2 parameter, idba ud parameter, cutadapt,etc)

