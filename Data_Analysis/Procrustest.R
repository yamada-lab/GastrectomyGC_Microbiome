library("phyloseq")
library("ggplot2")
library("plyr")
library("vegan")
library("ggplot2")
library("dplyr")

#mOTU
relab_mg1 <-read.table("/Volumes/Virtual/MountPoint/181103_Taxonomy/mOTU/2018-11-03_NCBI.species.abundances.5.1e-05.up.tsv",header = T, row.names = 1, sep = "\t")
#relab_mg1 <-read.table("/Volumes/Virtual/MountPoint/181103_Taxonomy/MetaPhlan2/2018-12-07_Metaphlan2_Species.5.0.001.up.tsv", header = T, row.names = 1, sep = "\t")
relab_mb <-read.table("/Volumes/Virtual/MountPoint/190402_Metabolome/190402_Metabolome.4.1e-06.up.tsv",header = T, row.names = 1, sep = "\t")

#head(relab_mb)
#hoge <- relab_mg1

#Sigonly

setwd("/Volumes/Virtual/MountPoint/190402_Metabolome/Procrustes")
relab_mg1 <-read.table("/Volumes/Virtual/MountPoint/190402_Metabolome/Procrustes/190919_mOTUSig.tsv",header = T, row.names = 1, sep = "\t")
#relab_mg1 <-read.table("/Volumes/Virtual/MountPoint/190402_Metabolome/Procrustes/190919_MetaPhlan2Sig", header = T, row.names = 1, sep = "\t")
relab_mb <-read.table("/Volumes/Virtual/MountPoint/190402_Metabolome/Procrustes/190919_MetaboliteSig.tsv",header = T, row.names = 1, sep = "\t")

#prop.table(as.matrix(relab_mb),margin = 2)
#a<-prcomp(t(hoge),scale. = F)
#plot(a$rotation[,1],a$rotation[,2])
##with_healthy
sample_id <- read.delim("/Volumes/Virtual/MountPoint/references/190127_MetadataPhyloseq.3",sep="\t",header=TRUE,row.names = 1,check.names = FALSE)

#converting to phyloseq format
#relab_mg <- subset(relab_mg, rownames(relab_mg) %in% rownames(relab_mb))

#select match samples
relab_mg <- relab_mg1[,colnames(relab_mg1) %in% colnames(relab_mb)]

#z-scoretransformations by features
relab_mg <- as.data.frame(t(scale(t(relab_mg))))
relab_mb <- as.data.frame(t(scale(t(relab_mb))))
#converting to phyloseq format
relabP_mg = otu_table(relab_mg, taxa_are_rows = TRUE)
metab= otu_table(relab_mb, taxa_are_rows = TRUE)
sampledata = sample_data(sample_id)
physeq_mg1 = phyloseq(relabP_mg,sampledata)
physeq_mb1 = phyloseq(metab,sampledata)

physeq_mb<-physeq_mb1
physeq_mg <- physeq_mg1

#PCA for metabolome and metagenome
ordu_mg1 = (ordinate(physeq_mg, method="RDA"))
ordu_mb1 = (ordinate(physeq_mb, method="RDA"))

cols <- c("Gastrectomy" = "#DC9238", "Healthy" = "#0DA4DA")

ord=plot_ordination(physeq_mg, ordu_mg1,type="samples",color="Status",axes = c(1,2,3)) 
ord + 
  stat_ellipse(type = "t") +
  theme_bw()+scale_colour_manual(values = cols)

ord=plot_ordination(physeq_mb, ordu_mb1,type="samples",color="Status",axes = c(1,2,3)) 
ord + 
  stat_ellipse(type = "t") +
  theme_bw()+scale_colour_manual(values = cols)

pro2<-protest(X=ordu_mg1,Y=ordu_mb1,scores="sites",symetric=FALSE)
pro2
#plot(pro2,kind=1)
#summary(pro2)

#ploting the procrustes results
plot(pro2$X,asp=1,pch=1,col=cols,scale=TRUE)
points(pro2$Yrot,col=cols,asp=1,pch=2)
segments(pro2$Yrot[,1], pro2$Yrot[,2], pro2$X[,1],pro2$X[,2],col=cols)