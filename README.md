# Gastrectomy Project Data Analysis

## Taxonomic annotations 
Taxonomic annotations were carried out using two different tools named [mOTU](https://motu-tool.org/) and [MetaPhlan2](https://bitbucket.org/biobakery/metaphlan2/src/default/)

### mOTU
The high-quality reads were subjected to mapping against the mOTU.v1.padded database using sequence identity and alignment cutoffs of 97% and 45 bp, respectively.

### MetaPhlan2
High-quality reads were mapped to unique clade-specific marker genes identified from ~17,000 reference genomes (~13,500 bacterial and archaeal, ~3,500 viral, and ~110 eukaryotic). 
For annotations using MetaPhlan2, please see the HUMAnN2 part in the Functional annotations

## Functional annotations
Functional annotations were carried out in our in house pipeline (KEGG-based) and [HUMAnN2](https://bitbucket.org/biobakery/humann2/wiki/Home)

### KEGG-based


### HUMAnN2

Custom scripts for runing in SGE system:  `HUMAnN2_pipeline.sh`.
For runing the scripts you need to specify
The format of the high quality reads should be `Sample_ID.samenames.fastq` or `Sample_ID.samenames.fastq.gz`
* `$1` : Sample ID
* `$2` : directory for the output
* `$3` : directory for the high quality reads

Example usage

```bash
qsub -l mem=4G -pe smp 8 HUMAnN2_pipeline.sh <sequence_id> <directory_out> <directory_in>
```


## Other annotations
### Oral microbes

### Phenotype annotations (BugBase)


### Microbe-Metabolite (MIMOSA)



## Data analysis custom scripts

### Statistics
* **Univariate statistic for differential relative abundance analysis**

We run [LEfSe](https://bitbucket.org/biobakery/biobakery/wiki/lefse) and MaAsLin [MaAsLin R package](https://bitbucket.org/biobakery/maaslin/src/default/) for univariate statistis as explained below. 
    * **LEfSe**

    The differences in microbiome features (taxonomy, KO modules) relative abundance were calculated by linear discriminant analysis (LDA) effect size [LEfSe](https://bitbucket.org/biobakery/biobakery/wiki/lefse). LEfSe first identifies features that are statistically different between healthy and gastrectomy groups using the non-parametric Kruskal–Wallis sum-rank test (P≤0.05). We modified the default calculation by controlling the multiple testing using Benjamini–Hochberg (BH) false discovery rate (FDR) correction procedure. All of the modification can be found in `LEfSe_Modif`.

    * **MaAsLin**

    The associations coefficient (variance explained) between the demographic data (groups, age, gender, smoking status, alcohol consumption status, and BMI) as explanatory variables and the detected microbial features as response were tested by the [MaAsLin R package](https://bitbucket.org/biobakery/maaslin/src/default/).  
To run the analysis simustaneusly we use....

>>>
packages:

input:

outputs

>>>


* **Statistical differences in demographic data**

The statistical differences were tested for the demographic data between healthy control and gastrectomy group to test the nature of the confounding factors. The two-sided MWU test (scipy.stats.mannwhitneyu version 0.18.1) was performed on numerical data (BMI, age, and dietary component information) and Fisher’s exact test (FisherExact 1.4.2) was performed on categorical data (gender, smoking status, and alcohol consumption status).
The script (writen in python 2.7.13) can be found in `Statistical_Analysis/Metadata_stats.py`.
>>>
packages:
* matplotlib
* pandas
* numpy
* seaborn
* statsmodel
* scipy
* FisherExact
* math
* collections

input:
* metadata table in .tsv format

outputs:
* statistical test (p value)
* boxplot
>>>

usage:
```python
python ./Statistical_Analysis/Metadata_stats.py metadata_table reference columns_to_test(e.g BMI, age,etc)
```

 




### Correlations 


## Command options for other tools

## Figure guide