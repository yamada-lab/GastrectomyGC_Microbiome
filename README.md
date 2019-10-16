# Gastrectomy Project Data Analysis

The **GastrectomyProject_Custom_Script** contains code to perform analysis for the gastrectomy paper.
Annotations table for microbial features (taxon, modules, and metabolites), reads profile table, subjects ID, and metadata table were availabe here `please insert the directory`. 
The raw sequencing data reported in this paper has been deposited at the DNA Data Bank of Japan (DDBJ) Sequence Read Archive (DRA), Tokyo, Japan under accession numbers **DRA007281, DRA008243, DRA006684, and DRA008156**. 
Bellow the example usage of sripts will be given. For several shell scripts (.sh), the code can be run in the terminal by ignoring the `qsub -cwd <options for qsub>`. Make sure to make script executable by using `chmod +x` command. 

## High quality reads process
The raw reads were subjected for the quality filtering before going to the annotations process. Please find the explanations of the process in **(Supplementary Material)**. 
The code were available in `High_quality_filtering/hq_reads.sh`

## Taxonomic annotations 
Taxonomic annotations were carried out using two different tools named [mOTU](https://motu-tool.org/) and [MetaPhlan2](https://bitbucket.org/biobakery/metaphlan2/src/default/)

### mOTU
The high-quality reads were subjected to mapping against the mOTU.v1.padded database using sequence identity and alignment cutoffs of 97% and 45 bp, respectively.
The annotations were run following the instruction in [mOTU](https://motu-tool.org/) tutorial.

### MetaPhlan2
High-quality reads were mapped to unique clade-specific marker genes identified from ~17,000 reference genomes (~13,500 bacterial and archaeal, ~3,500 viral, and ~110 eukaryotic). 
For annotations using MetaPhlan2, please see the HUMAnN2 part in the functional annotations.

## Functional annotations
Functional annotations were carried out in our in house pipeline (KEGG-based) and [HUMAnN2](https://bitbucket.org/biobakery/humann2/wiki/Home)

### KEGG-based

We annotated the functional modules using our *in house* pipeline. Please find the manuscript **(Supplementary material)** for detail processes. The code were available upon request. 

### HUMAnN2

Custom scripts for runing in SGE system:  `HUMAnN2_pipeline.sh`.
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
    *  **LEfSe**

    The differences in microbiome features (taxonomy, KO modules) relative abundance were calculated by linear discriminant analysis (LDA) effect size [LEfSe](https://bitbucket.org/biobakery/biobakery/wiki/lefse). LEfSe first identifies features that are statistically different between healthy and gastrectomy groups using the non-parametric Kruskal–Wallis sum-rank test (P≤0.05). We modified the default calculation by controlling the multiple testing using Benjamini–Hochberg (BH) false discovery rate (FDR) correction procedure. All of the modification can be found in `LEfSe_Modif`.
    The script for runing LEfSe presented in the manuscript together with the parameter were written in `LEfSe_Modif/LEfSeRun.sh`.
    usage example in SGE:
    ```
    qsub -cwd LEfSeRun.sh <columns for annotations samples> <input table contains the relative abundance> <Category analysis> <prefix for title>
    ```

    *  **MaAsLin**

    The associations coefficient (variance explained) between the possible confounding effects of clinical parameters (BMI, serum glucose and total cholesterol), demographic data (e.g age, gender) and medical history (e.g history of drug consumption and diseases) as explanatory variables and the detected microbial features as response were tested by the [MaAsLin R package](https://bitbucket.org/biobakery/maaslin/src/default/).

To run the analysis simultaneously we use `Data_Analysis/MaAsLinRun.sh`


>>>

input: 
* relative abundance tables (species, KEGG modules, metabolites)
* metadata files contains the samples ID and clinical informations (BMI, Age, Gender, etc.)
* selected samples ID
* statistic results from `LEfSeRun.sh`

outputs:
* microbial features associations with the metadata
* graph for associations coefficient (crude vs adjusted)  
>>>
  

The script for runing MaAsLin presented in the manuscript together with the parameter were written in `LEfSe_Modif/LEfSeRun.sh`.
usage example in SGE:
    ```
    qsub -cwd MaAsLin.sh <input table contains the relative abundance> <metadata tables> <prefix for title> <selected ID> <statistic tables from LEfSe>
    ```


* **Statistical differences in demographic data**

The statistical differences were tested for the demographic data between healthy control and gastrectomy group to test the nature of the confounding factors. The two-sided MWU test (scipy.stats.mannwhitneyu version 0.18.1) was performed on numerical data (BMI, age, and dietary component information) and Fisher’s exact test (FisherExact 1.4.2) was performed on categorical data (gender, smoking status, and alcohol consumption status).
The script (writen in python 2.7.13) can be found in `Data_Analysis/Metadata_stats.py`.
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
python ./Data_Analysis/Metadata_stats.py metadata_table reference columns_to_test(e.g BMI, age,etc)
```


### Correlations 
* **Species - species or Genus - Genus Correlation**

We estimated the microbial association in each of group using [SparCC](https://bitbucket.org/yonatanf/sparcc/src/default/)
Custom script for runing in SGE system were available at `Data_Analysis/SparCC_corr.sh`
usage:
```bash
qsub -cwd ./Data_Analysis/SparCC_corr.sh <Species/genus_relative abundance table> <List of selected species/genus>
```
After the calculations, the results were visualized using [igraph R packages](https://igraph.org/r/) and the script were available here `insert the script` 
* 

## Command options for other tools

## Figure guide