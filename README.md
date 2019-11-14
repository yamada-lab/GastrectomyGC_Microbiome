# Gastrectomy Project Data Analysis
**Author: Pande Putu Erawijantari**

**School of Life Science and Technology, Tokyo Institute of Technology**

**contact:**

**e-mail:erawijantari.p.aa@m.titech.ac.jp**

**twitter:@erawijantaript**

The **GastrectomyProject_Custom_Script** contains collection of scripts to perform analysis for the gastrectomy paper.
Annotations table for microbial features (taxon, modules, and metabolites), reads profile table, subjects ID, and metadata table were availabe in the `./Raw_Data`. 
The raw sequencing data reported in this paper has been deposited at the DNA Data Bank of Japan (DDBJ) Sequence Read Archive (DRA), Tokyo, Japan under accession numbers **DRA007281, DRA008243, DRA006684, and DRA008156**. 
We will give the example of runing the analysis and generating the figures in the explanation below.

## High quality reads process
The raw reads were subjected for the quality filtering before going to the annotations process. Please find the explanations of the process in the **(Supplementary Methods)**. 
The code were available in `./High_quality_filtering/hq_reads.sh`
The process need instalations of [Bowtie2 version 2.2.2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) and [cutadapt version 1.9.1](https://cutadapt.readthedocs.io/en/v1.9.1/guide.html).
For Phix control reads removal, we downloaded the reference [NC_001422.1](https://www.ncbi.nlm.nih.gov/genome/?term=NC_001422.1) (*Enterobacteria phage phiX174* sensu lato, complete genome) from the NCBI website and indexed it by Bowtie2.
For the removal of human genomes, we downloaded the reference [GRCh38](ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/GCA000001405.15_GRCh38_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/) from the NCBI website and indexed it by Bowtie2.
The reads profile before and after the process can be found in the **Supplementary Table S3**

## Taxonomic annotations 
Taxonomic annotations were carried out using two different tools named [mOTU](https://motu-tool.org/) and [MetaPhlan2](https://bitbucket.org/biobakery/metaphlan2/src/default/)
The annotations for each samples can be combine by utilizing the script in `./Data_Analysis/CombineMatrix.py`. For using the script the path that directed to each sample should be similar.
For example we put the all sample under the same directory `project/results/` then the usage for the script is
```
python ./Data_Analysis/CombineMatrix.py/ project/results/ <reference file to convert the header from sample name to sample ID> <suffix pattern of the results to be combine>
```

### mOTU
The high-quality reads were subjected to map against the mOTU.v1.padded database using sequence identity and alignment cutoffs of 97% and 45 bp, respectively.
The annotations were run following the instruction in [mOTU](https://motu-tool.org/) tutorial.
The table of species relative abundance annotated by mOTU can be found in `./Raw_Data/Raw_Species_mOTU.tsv`. 
Species with average relative abundance exceed 0.001% and appeared in at least 5% of samples number (five samples) were retrieved for downstream analysis can be found in `./Raw_Data/Species_mOTUAnn_5.1e-05.up.tsv` 
The downstream analysis for the taxon generated by this annotations includes: **PERMANOVA**,**LEfSe statistical test**,**MaAsLin associations test with metadata**, **phenotypic and oral microbes characterization**,**species and genus network by SparCC**, **genus and metabolite correlations**, and **Procrustest**.

### MetaPhlAn2
High-quality reads were mapped to unique clade-specific marker genes identified from ~17,000 reference genomes (~13,500 bacterial and archaeal, ~3,500 viral, and ~110 eukaryotic). 
For annotations using MetaPhlan2, please see the HUMAnN2 part in the functional annotations.
The table of species relative abundance annotated by MetaPhlAn2 can be found in `./Raw_Data/Raw_Species_MetaPhlan2.tsv`. 
Species with average relative abundance exceed 0.1% and appeared in at least 5% of samples number (five samples) were retrieved for downstream analysis can be found in `./Raw_Data/Species_MetaPhlan2.5.0.001.up.tsv` 
The downstream analysis for the taxon generated by this annotations includes: **PERMANOVA**,**LEfSe statistical test**,**MaAsLin associations test with metadata**, **phenotypic and oral microbes characterization**, and **Procrustest**.

## Functional annotations
Functional annotations were carried out in our in house pipeline (KEGG-based) and [HUMAnN2](https://bitbucket.org/biobakery/humann2/wiki/Home)
Both annotations generates the genes and KEGG Orthology table. We converted the KEGG Orthology tables into KEGG modules level by implementing [omixer-rpm](https://github.com/raeslab/omixer-rpm).
For module annotations we generated the reference database using the KEGG database 2017 version (available upon request).
The annotations for each samples can be combine by utilizing the script in `./Data_Analysis/CombineMatrix.py`. For using the script the path that directed to each sample should be similar.
For example we put the all sample under the same directory `project/results/` then the usage for the script is
```
python ./Data_Analysis/CombineMatrix.py/ project/results/ <reference file to convert the header from sample name to sample ID> <suffix pattern of the results to be combine>
```

### KEGG-based
We annotated the functional modules using our *in house* pipeline. Please find the manuscript **(Supplementary Methods)** for detail processes. The code were available upon request. 
The table of KEGG module relative abundance annotated by KEGG-based annotations can be found in `./Raw_Data/Raw_Module_KEGGbased.tsv`. 
The KEGG modules with average relative abundance of 0.0001% and appeared in at least 5% of samples number (five samples) were selected for the downstream analysis can be found in `./Raw_Data/Modules_InHouseAnn.5.1e-06.up.tsv` 
The downstream analysis for the KEGG modules generated by this annotations includes: **LEfSe statistical test** and **MaAsLin associations test with metadata**

### HUMAnN2

Custom scripts for runing in SGE system:  `./Data_Analysis/HUMAnN2_pipeline.sh`.

* `$1` : Sample ID
* `$2` : directory for the output
* `$3` : directory for the high quality reads

Example usage

```bash
qsub -l mem=4G -pe smp 8 ./Data_Analysis/HUMAnN2_pipeline.sh <sequence_id> <directory_out> <directory_in>
```
The table of KEGG module relative abundance annotated by HUMAnN2 annotations can be found in `./Raw_Data/Raw_Module_KEGGbased.tsv/Raw_Module_HUMAnN2.tsv`. 
The KEGG modules with average relative abundance of 0.0001% and appeared in at least 5% of samples number (five samples) were selected for the downstream analysis can be found in `./Raw_Data/Modules_Humann2Ann.5.1e-06.up.tsv`
HUMAnN2 also generated the KEGG Orthology profiles stratified by species that we used for MIMOSA (see below in the **Microbe-Metabolite (MIMOSA)** section). The file available upon request due to size limitations. 
The downstream analysis for the KEGG modules generated by this annotations includes: **LEfSe statistical test**, **MaAsLin associations test with metadata**, and **MIMOSA**

## Other annotations
### Oral microbes
We categorized the mOTU-annotated species and MetaPhlAn2-annotated species into oral microbes or others or others based on the [expanded Human Oral Microbiome Database (eHOMD)](http://www.homd.org/).
The tabularized file were downloaded from `http://www.homd.org/?name=HOMD`. If the species were matched with the species name in the table, we categorized the species as `oral`, while unmatched species were categorized as `other`.
We later calculate the total relative abundance of species labeled as `oral` in each group (Gastrectomy, n=50; Control, n=56) and performed statistical test (Mann Whitney-U) to compare the two groups.
Custom script for runing the categorizations and statistical test were available at `./Data_Analysis/Stat_OralOther.py`
Example usage:
```
python Stat_OralOther.py <species_table> <oral_spec_table>
```
The output of the script includes:
- ***OralOtherRatio.csv** : the calculation for each samples 
- ***calcOralOther.tsv** : statistical results
- ***Oral.pdf** : Boxplot of total relative abundance of species categorized as oral microbes (*use in Figure 3C,D*)
- ***OralOtherRatio.pdf**: Boxplot of ratio between total relative abundance of species categorized as oral microbes *versus* other

### Phenotype annotations (BugBase)
We predict the phenotypic properties of the species based by implementing [BugBase](https://bugbase.cs.umn.edu/).
BugBase were originally designed for OTU table from OTU picking for 16S sequencing reads. 
For OTU-picking from shotgun metagenomic, like our data we follow the workflow as define [here](https://github.com/knights-lab/BugBase) by using SHOGUN for BugBase.
After OTU-picking, we run the BugBase with default parameters as defined [here](https://github.com/knights-lab/BugBase). 
The results were shown in the **Fig 3E,F**. 

### Microbe-Metabolite (MIMOSA)
We implemented [MIMOSA](https://github.com/borenstein-lab/MIMOSA) to generate the community-specific metabolic network model and identify potential taxonomic and gene contributor.
The example of script and parameter that we used in the publications were available at `./Data_Analysis/Spec-Metab.sh`.
The input that we used for the analysis were generated by [HUMAnN2](https://bitbucket.org/biobakery/humann2/wiki/Home), for which we need to convert the format by executing `./Data_Analysis/convert_HUMANN2-PICRUST.R` script.
The analysis also need the link of KO to reactions file from KEGG which can be downloaded from `genes/ko/ko_reaction.list in the KEGG FTP download` and full KEGG reaction annotations from `ligand/reaction/reaction in the KEGG FTP download`.
For the manucript we used the 2017 version of database downloaded from KEGG.
The results were shown in the **Figure 5C**, **Supplementary Figure S7**, and **Supplementary Table 15**.

## Data analysis custom scripts

### Statistics
* **Univariate statistic for differential relative abundance analysis**
    *  **LEfSe**

    The differences in microbiome features (taxonomy, KO modules) relative abundance were calculated by linear discriminant analysis (LDA) effect size [LEfSe](https://bitbucket.org/biobakery/biobakery/wiki/lefse). LEfSe first identifies features that are statistically different between control and gastrectomy groups using the non-parametric Kruskal–Wallis sum-rank test (P≤0.05). We modified the default calculation by controlling the multiple testing using Benjamini–Hochberg (BH) false discovery rate (FDR) correction procedure. All of the modification can be found in `LEfSe_Modif`.
    The script for runing LEfSe presented in the manuscript together with the parameter were written in `LEfSe_Modif/LEfSeRun.sh`.
    
>>>
input: 
* relative abundance tables (species, KEGG modules, metabolites) e. g in ./Raw_Data/Raw_Metabolome.tsv
* reference file contains the sequence ID, sample ID and clinical informations needed for test. The format can be found in ./Raw_Data/Reference_name.txt
* selected samples ID, files contains the selected sample ID

outputs:
* ***StatsCombine*** : the statistical output contains the P value, q value, log10(LDA score), and enrichment shown in **Supplementary Table S5-14, S17, S18**
* ***.svg** : plot such as shown in **Figure 3A,B**, **Figure 4A,C,**, and **Supplementary Figure 4**
>>>    

usage example in SGE:
    ```
qsub -cwd LEfSeRun.sh <columns for annotations samples, e.g columns 1 for Control vs Gastrectomy, then type 1> <relative abundance table> <reference file> <midle pattern for title> <prefix for title> <list of ID to retrieve> <preference cut off> <relative abundance cut off>
```


* 
    * **MaAsLin**
    
    The associations coefficient (variance explained) between the possible confounding effects of clinical parameters (BMI, total cholesterol, diabetes medications, and gastric acids medications, age,gender), demographic data (e.g age, gender) and medical history (e.g history of drug consumption and diseases) as explanatory variables and the detected microbial features as response were tested by the [MaAsLin R package](https://bitbucket.org/biobakery/maaslin/src/default/).
    To run the analysis simultaneously we use `Data_Analysis/MaAsLinRun.sh`

>>>
input: 
* relative abundance tables (species, KEGG modules, metabolites) e. g in ./Raw_Data/Raw_Metabolome.tsv
* metadata files contains the samples ID and clinical informations (BMI, Age, Gender, etc.) e. g in `Metadata_Tables.txt`
* selected samples ID, files contains the selected sample ID
* statistic results from `LEfSeRun.sh` with pattern *StatsCombine in the title

outputs:
* microbial features associations with the metadata such as written in **Supplementary Table 10,11,12,14**
* graph for associations coefficient (crude vs adjusted) shown in **Supplementary Figure 2**
>>>
  

The script for runing MaAsLin presented in the manuscript together with the parameter were written in `LEfSe_Modif/LEfSeRun.sh`.
usage example in SGE:
    
    qsub -cwd MaAsLin.sh <input table contains the relative abundance> <metadata tables> <prefix for title> <selected ID> <statistic tables from LEfSe>


* **Statistical differences in demographic data**

The statistical differences were tested for the demographic data between control and gastrectomy group to test the nature of the confounding factors. The two-sided MWU test (scipy.stats.mannwhitneyu version 0.18.1) was performed on numerical data (BMI, age, and dietary component information) and Fisher’s exact test (FisherExact 1.4.2) was performed on categorical data (gender, smoking status, and alcohol consumption status).
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
* statistical test (p value) shown in **Table 1**,and **Supplementary Table S1**
* boxplot
>>>

usage:
```python
python ./Data_Analysis/Metadata_stats.py metadata_table reference columns_to_test(e.g BMI, age,etc)
```

* **PERMANOVA analysis on demographic data**

The overall species, KEGG modules and metabolites distribution were tested for the demographic data between different categorical or numerical grouping based on metadata. 

>>>
packages:
* phyloseq
* ggplot2
* pylr
* vegan
* tidyverse

input:
* relative abundance table
* metadata table

outputs:
* PERMANOVA table test shown in **Supplementary Table S4**

>>>
The script (writen in R and can be run at Rstudio) can be found in `Data_Analysis/PERMANOVA.R`. 



### Correlations 
* **Species - species or Genus - Genus Correlation**

We estimated the microbial association in each of group using [SparCC](https://bitbucket.org/yonatanf/sparcc/src/default/)
Custom script for runing in SGE system were available at `Data_Analysis/SparCC_corr.sh`
usage:
```bash
qsub -cwd ./Data_Analysis/SparCC_corr.sh <Species/genus_relative abundance table> <List of selected species/genus>
```
After the calculations, the results were visualized using [igraph R packages](https://igraph.org/r/) and the script were available in `./Data_Analysis/SparCC-Network.R` 


## Figure script's guide 

Below are the list of script and its usage to generate the raw version of figure's and supplementary figure's element.
- **Fig. 1A,B,D; Supplementary Figure S8 A**: *script*: ./Data_Analysis/PCoA_viz.R ; *usage*: `Rscript ./Data_Analysis/PCoA_viz.R`; the path were specified inside the script 
- **Fig. 1C**: ./Data_Analysis/BetaDiver_boxplot.py; *usage*: `python ./Data_Analysis/BetaDiver_boxplot.py <species relative abundance table>`
- **Fig. 2A,B,C,D**: ./Data_Analysis/alpha_diversity_boxplot.py; *usage*: `python ./Data_Analysis/alpha_diversity_boxplot.py <species relative abundance table> <metric e.g shannon,chao1>`
- **Fig. 2E,F**: ./Data_Analysis/alpha_phylum.py; *usage*: `python ./Data_Analysis/alpha_phylum.py <species relative abundance table with taxon level>`
- **Fig. 3A,B**: ./LEfSe_Modif/plot_cladogram.modif.py; *usage*: `python ./LEfSe_Modif/plot_cladogram.modif.py <results from LEfSe *.res run toward species relative abundance table with taxon level> <output_title> --dpi 300 --format svg`
- **Fig. 3C,D**: ./Data_Analysis/Stat_OralOther.py; *usage*: `python ./Data_Analysis/Stat_OralOther.py <species relative abudance table> <table of oral species downloaded from eHOMDB>`
- **Fig. 3E,F**: see the explanations in BugBase, the figures were generated by the tools
- **Fig 4A,C; Supplementary Figure S3; Supplementary Figure S8B,C,D,E,F** : this figures were combinations of boxplot generated by runing the script `./Data_Analysis/Boxplot_generator.py` and `./LefSe_Modif/plot_res.modif.py` and labeled based on the KEGG module definitions of 2017 version of KEGG database
- **Fig 4B; Supplementary Figure S4**: these figure were generated by HUMAnN2 utility's script called `humann2_barplot`; usage: `humann2_barplot --sort sum metadata --input <module stratified by species relative abudance table with metadata> --focal-feature <module of interest> --focal-metadatum <Metadata columns> --last-metadatum <Metadata colums> --output <output_title>.pdf -c <custom color for each species> --meta-colormap RdBu -x`
- **Fig 5A,B; Supplementary Figure S5**: ./Data_Analysis/SparCC-Network.R;  the path and customization were specified inside the script
- **Fig 5C**: ./Data_Analysis/Spearman_CorrGenus-Metab.py; please see the script that enable various customization
- **Supplementary Figure S2**:./Data_Analysis/ConfoundVarExp.py; usage: `python ConfoundVarExp.py <MaAsLin result adjusted> <MaAsLin result non-adjusted> -r <relative abundance table> -s <statistic table output from LefSe>`
- **Supplementary Figure S6**:./Data_Analysis/Procrustest.R
- **Supplementary Figure S9,S10**:./Data_Analysis/boxplot_gastrectomyType.py; usage: `python ./Data_Analysis/boxplot_gastrectomyType.py <relative abundance table> <metadata> <output_title> <feature_of_interest> -t boxplot`

