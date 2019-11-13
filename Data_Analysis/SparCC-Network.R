#network new
# load R graph library igraph
library(igraph)
library(corrplot)

#ordered by LDA score
#relabAll <-
#Healthy
NPS.corrH <- as.matrix(read.table("/Volumes/Virtual/MountPoint/181103_Taxonomy/mOTU/SparCC/190625_SparCCAllSpec/Healthy/2019-06-28_Species_CountSigOverlapmOTU.Healthy.sparcc.tsv", header = T, row.names = 1, sep = "\t"))
pathH="/Volumes/Virtual/MountPoint/181103_Taxonomy/mOTU/SparCC/190625_SparCCAllSpec/Healthy/Species_CountSigOverlapmOTU_pvals/two_sided5000.txt"
relabH <- as.matrix(read.table("/Volumes/Virtual/MountPoint/181103_Taxonomy/mOTU/SparCC/190625_SparCCAllSpec/Healthy/2019-06-28_Species_CountSigOverlapmOTU.Healthy.tsv", header = T, row.names = 1, sep = "\t"))

#Gastrectomy
NPS.corrG <- as.matrix(read.table("/Volumes/Virtual/MountPoint/181103_Taxonomy/mOTU/SparCC/190625_SparCCAllSpec/Gastrectomy/2019-06-28_Species_CountSigOverlapmOTU.Gastrectomy.sparcc.tsv", header = T, row.names = 1, sep = "\t"))
pathG="/Volumes/Virtual/MountPoint/181103_Taxonomy/mOTU/SparCC/190625_SparCCAllSpec/Gastrectomy/Species_CountSigOverlapmOTU_pvals/two_sided5000.txt"
relabG <- as.matrix(read.table("/Volumes/Virtual/MountPoint/181103_Taxonomy/mOTU/SparCC/190625_SparCCAllSpec/Gastrectomy/2019-06-28_Species_CountSigOverlapmOTU.Gastrectomy.tsv", header = T, row.names = 1, sep = "\t"))

 
##functions on the pval cuting
parsPCorr <- function(sparcc,pathpval,cut){
  pvals=read.table(pathpval,header=TRUE,sep="\t")
  pvals.mat=pvals[,2:ncol(pvals)]
  # set p-values of 0 to a non-zero, small p-value so we can take the logarithm
  pvals.mat[pvals.mat==0]=0.000000001
  # convert into significance
  sig.mat=-1*log10(pvals.mat)
  #Df parsing
  sparcc[sig.mat<=1.30103]=0
  sparcc<-as.matrix(sparcc)
  #select the correlation value
  sparcc[abs(sparcc)<cut]=0
  return(sparcc)
}

#graph viz
GraphViz<-function(corrtab,relab){
  #without modification
  g <- graph_from_adjacency_matrix(
    corrtab,
    mode="undirected",
    weighted=TRUE,
    diag=FALSE)
    means <- as.matrix(rowMeans(relab))
    #put the enrichment information genus
    #V(g)$enrichment <- c("Healthy", "Gastrectomy", "Gastrectomy", "Healthy", "Gastrectomy", "Gastrectomy",
    #                  "Gastrectomy", "Gastrectomy", "Healthy", "Gastrectomy","Gastrectomy","Gastrectomy",
    #                  "Healthy", "Gastrectomy", "Gastrectomy", "Healthy")
    
    #put the enrichment information species
    #V(g)$enrichment <- c("Gastrectomy", "Gastrectomy", "Gastrectomy", "Gastrectomy", "Gastrectomy",
    #                     "Gastrectomy","Gastrectomy", "Gastrectomy", "Gastrectomy", "Gastrectomy",
    #                     "Gastrectomy","Gastrectomy","Gastrectomy", "Gastrectomy", "Gastrectomy", 
    #                     "Gastrectomy","Gastrectomy","Gastrectomy", "Gastrectomy", "Gastrectomy", 
    #                     "Gastrectomy","Gastrectomy","Gastrectomy", "Gastrectomy", "Gastrectomy", 
    #                     "Gastrectomy","Gastrectomy", 
    #                     "Healthy", "Healthy","Healthy","Healthy","Healthy",
    #                    "Healthy","Healthy","Healthy","Healthy","Healthy",
    #                     "Healthy","Healthy")
    V(g)$relab <- c(means[,1])
    #Colour positive correlation edges as green
    E(g)[which(E(g)$weight<0)]$color <- "#3f7f93"
    #Colour positive correlation edges as red
    E(g)[which(E(g)$weight>0)]$color <- "#da3b46"
    x<-as.matrix(E(g)[which(E(g)$weight>0)])
    #Convert edge weights to absolute values
    E(g)$weight <- abs(E(g)$weight)
    #Change line width
    E(g)[which(E(g)$weight>0.2)]$width <- 0.5
    E(g)[which(E(g)$weight>0.4)]$width <- 0.5*10
    E(g)[which(E(g)$weight>0.6)]$width <- 0.5*20
    #E(g)[which(E(g)$weight>0.4)]
    #Assign names to the graph vertices (optional)
    #V(g)$name <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,
    #               17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39)
  return(g)
}

#union two graph
union<-function(graph1,graph2){
  #if else for the layout
  gUnion<-print_all(graph1 %u% graph2)
  return(gUnion)
}

#Layouting for smae positions
Layouteach<-function(graph1,graphS,layoutc){
  #if else for the layout
  if(layoutc=='fr'){
    #coords <- layout_with_fr(graphS)
    coords <- layout_with_kk(graphS)
  }else {
    coords <- layout_in_circle(graphS, order = V(graphS))
  }
  Colorvertex=c( '#0072BC','#F7941D')[1+(V(graph1)$enrichment=="Gastrectomy")]
  ## size of nodes based on degree
  #V(g)$vertex_degree <-  degree(g)*2.5
  ###size by relab
  rescale = function(x,a,b){b + x*a}
  V(graph1)$relab_val <-  rescale(V(graph1)$relab,0.01,5)
  plot(graph1, layout=coords, vertex.color=Colorvertex
       ,vertex.label.dist=1.0,vertex.size =V(graph1)$relab_val)
  ##implement the code below to save the figures
  #dev.copy(pdf,'myplot.pdf')
  #dev.off()
}

#Test_functions
A<-parsPCorr(NPS.corrH,pathH,0.2) #return a table
H<-GraphViz(A,relabH) #return a graph
B<-parsPCorr(NPS.corrG,pathH,0.2) #return a table
G<-GraphViz(B,relabG) #return a graph
V(G)$name
union1<-union(G,H)
#plot for gastrectomy
Layouteach(G,union1,'fr')
#plot for healthy
Layouteach(H,union1,'fr')

#Calculation network
pap.centz <- centr_degree(g)$centralization
pap.centz[2] <- centr_clo(g)$centralization
pap.centz[3] <- centr_betw(g)$centralization
pap.centz[4] <- centr_eigen(g)$centralization
# Table
tab.centz <- rbind(Healthy=pap.centz)
dimnames(tab.centz)[[2]] <- c("Degree", "Closeness", "Betweenness", "Eigenvector")
tab.centz


degree1 = as.matrix(degree(g))
strength1 = as.matrix(strength(g))
betweeness1 = as.matrix(betweenness(g))
matrixtotal =cbind(degree1, strength1, betweeness1) 
matrixtotal