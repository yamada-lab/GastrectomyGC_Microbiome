#functions for repeated measure of PERMANOVA 
#Erawijantari-Titech-2019

library("phyloseq")
library("ggplot2")
library("plyr")
library("vegan")
library("ggplot2")
library("tidyverse")


#covariates for all
perindividual_covariates <- c(
  "Status","BMI","Age","Gender","SmokingStatus","DrinkingStatus","AlcoholConsumption",
  "complication","Glucose","Total_Cholesterol","HypertensionMed","Cholesterol",
  "DiabetesMed","GoutMed","GastricAcidMedication","Analgesic","Anticoagulant","Others",
  "LungCancer","LiverCancer","BreastCancer","UterineCancer","OtherCancers","Stroke",
  "CardiacInfarction","Angina","Hypertension","Diabetes","Dyslipidemia","Cataract","StomachUlcer",
  "StomachPolyps","DuodenalUlcers","ColorectalPolyps","ChronicHepatitis_LiverCirrhosis",
  "Gallstone","Ureteral_KidneyStones","Gout","HipFracture","Arm_wristFracture",
  "OtherDiseases"
)

#covariate for gastrectomy
gastrectomy_covariates <-c(
  "Surgery_Type","Reconstruction","Hypertension","Cholesterol","DiabetesMed",
  "GastricAcidMedication","Analgesic","Anticoagulant","Others","LungCancer",
  "BreastCancer","OtherCancers","Stroke","CardiacInfarction","Angina",
  "Hypertension","Diabetes","Dyslipidemia","Cataract","StomachUlcer","StomachPolyps",
  "DuodenalUlcers","ColorectalPolyps","Gallstone","Ureteral_KidneyStones","Gout",
  "Arm_wristFracture",
  "OtherDiseases","Dumping_Syndrome"
)


completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}

#permanova and get the value (R2 and p-val) --> create tables (for all subjects)
#covar contain the list of variables for permanova tes
#observe again like previous analysis is better
perm_all<-function(relab,metadata,covar){
  ###find alternative for subseting through loop
  #better do before converting to phyloseq
  dffinal=list()
  for(meta in covar){
    print(meta)
    set.seed(123)
    #metadata modif loop and relative abundance modif for subseting subject with NA
    meta_new<- completeFun(metadata, meta)
    #subset
    meta_new <- meta_new %>% select(meta)
    #select match samples
    relab_new <- relab[,colnames(relab) %in% rownames(meta_new)]
    sampledata = sample_data(meta_new)
    relabP = otu_table(relab_new, taxa_are_rows = TRUE)
    physeq1 = phyloseq(relabP,sampledata)
    #for metabolome
    physeq1  = transform_sample_counts(physeq1, function(x) x / sum(x) )
    #s_bray <- distance(physeq1, method = "bray")  ##same
    Dist <- phyloseq::distance(physeq1, method = "bray")
    sampledf <- data.frame(sample_data(physeq1))
    ad_result<-adonis(Dist~.,data=sampledf,permutations = 1000)
    outdf<-as.matrix(ad_result$aov.tab[1,])
    #row.names(outdf)<-meta
    dffinal<-rbind(dffinal,outdf) #not significant?
  }
  return(dffinal)
}


#example runing in Rstudio
#FROM LOCAL
#setwd("/Users/erawijantari/Documents/ProgTemp/PERMANOVA/")
#relabdf <-read.table("/Users/erawijantari/Documents/ProgTemp/PERMANOVA/2018-11-03_NCBI.species.abundances.5.1e-05.up.tsv", header = T, row.names = 1, sep = "\t")
#relabdf <-read.table("/Users/erawijantari/Documents/ProgTemp/PERMANOVA/2018-12-07_Metaphlan2_Species.5.0.001.up.tsv", header = T, row.names = 1, sep = "\t")
relabdf <-read.table("/Users/erawijantari/Documents/ProgTemp/PERMANOVA/190402_Metabolome.4.1e-06.up.tsv",header = T, row.names = 1, sep = "\t")
#clinical
sample_id <- read.delim("/Users/erawijantari/Documents/ProgTemp/PERMANOVA/190910_MetadataPhyloSeq",sep="\t",header=TRUE,row.names = 1,check.names = FALSE)
#for all
df_tot=perm_all(relabdf,sample_id,perindividual_covariates)
write.table(df_tot,"190930_PERMANOVA_All_Metabolome.tsv",sep="\t")

#for gastrectomy
#metagast <- sample_id %>% filter(Status=="Gastrectomy")
metagast <- subset(sample_id, Status=="Gastrectomy")
df_gas=perm_all(relabdf,metagast,gastrectomy_covariates)
write.table(df_gas,"190910_PERMANOVA_Gastrectomy_Metabolome.tsv",sep="\t")

