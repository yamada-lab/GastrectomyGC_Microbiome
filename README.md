# Gastrectomy Project Data Analysis

## Taxonomic annotations 
Taxonomic annotations were carried out using two different tools named [mOTU](https://motu-tool.org/) and [MetaPhlan2](https://bitbucket.org/biobakery/metaphlan2/src/default/)

### mOTU


### MetaPhlan2
For annotations using MetaPhlan2, please see the HUMAnN2 part in the Functional annotations

## Functional annotations
Functional annotations were carried out in our in house pipeline (KEGG-based) and [HUMAnN2](https://bitbucket.org/biobakery/humann2/wiki/Home)

### KEGG-based


### HUMAnN2

Custom scripts for runing in SGE system:  `HUMAnN2_pipeline.sh`.
For runing the scripts you need to specify
The format of the high quality reads should be `Sample_ID.samenames.fastq` or `Sample_ID.samenames.fastq.gz`
*`$1` : Sample ID
*`$2` : directory for the output
*`$3` : directory for the high quality reads

Example usage
```qsub -l mem=4G -pe smp 8 HUMAnN2_pipeline.sh <sequence_id> <directory_out> <directory_in>```


## Other annotations
### Oral microbes

### Phenotype annotations (BugBase)


### Microbe-Metabolite (MIMOSA)



## Data analysis custom scripts

### Statistics


### Correlations 


## Command options for other tools

## Figure guide