#### RNA-seq_plot_2_conditions.R
#
# Usage is:
# Rscript RNA-seq_plot_2_conditions.R [folder name] [decimal symbol] [DEG Table 1] [Condition 1 name] [DEG Table 2] [Condition 2 name] [log2 FC threshold] [pval adj threshold] [organism] [output name]
#
# Table 1 and 2 must display the following column names and must be tab-delimited:
# Row.names external_gene_name  log2FoldChange  padj  chromosome_name description
#
# "organism" can be one of these, other options are available in biomaRt: celegans_gene_ensembl, dmelanogaster_gene_ensembl, hsapiens_gene_ensembl, mmusculus_gene_ensembl, scerevisiae_gene_ensembl
#
# For dependencies do in a terminal (remove # of course) :
#
# R
# install.packages('plotly')
# source("https://bioconductor.org/biocLite.R")
# biocLite("biomaRt")
# quit()

cat("Loading dependencies, if it fails do : head -15 RNA-seq_plot_2_conditions.R (", as.character(Sys.time()),")","\n")

library('biomaRt')
library('plotly')

args <- commandArgs(trailingOnly = TRUE)
setwd(args[1])
#setwd("~/Sukanya/")


cat("Done! Importing and formating data... (", as.character(Sys.time()),")","\n")

table1<-read.csv(file=args[3], sep="\t",header=T,dec=args[2])
table1<-data.frame(table1$Row.names,table1$external_gene_name,table1$log2FoldChange,table1$padj,table1$chromosome_name,table1$description)
colnames(table1)<-c("rows","genes","FC","pval","chr","desc")

table2<-read.csv(file=args[5], sep="\t",header=T,dec=args[2])
table2<-data.frame(table2$Row.names,table2$external_gene_name,table2$log2FoldChange,table2$padj,table2$chromosome_name,table2$description)
colnames(table2)<-c("rows","genes","FC","pval","chr","desc")

#table1<-read.csv(file="Mousebrain_CRE-MOF_KD.txt", sep="\t",header=T,dec=",")
#table2<-read.csv(file="Mousebrain_CRE-K3_KD.txt", sep="\t",header=T,dec=",")

#mof_DEG<-subset(mof,mof$padj<0.05 & abs(mof$log2FoldChange)>0.584962501)
#mof_DEG["DEG"]<-"TRUE"
#mof_DEG<-mof_DEG[,c(1,10)]
#mof<-merge(mof,mof_DEG,by="Row.names",all.x=T)

#k2_DEG<-subset(k2,k2$padj<0.05 & abs(k2$log2FoldChange)>0.584962501)
#k2_DEG["DEG"]<-"TRUE"
#k2_DEG<-k2_DEG[,c(1,10)]
#k2<-merge(k2,k2_DEG,by="Row.names",all.x=T)

#k3_DEG<-subset(k3,k3$padj<0.05 & abs(k3$log2FoldChange)>0.584962501)
#k3_DEG["DEG"]<-"TRUE"
#k3_DEG<-k3_DEG[,c(1,10)]
#k3<-merge(k3,k3_DEG,by="Row.names",all.x=T)

merged<-merge(table1,table2,by="rows")

name1=args[4]
name2=args[6]

cat("Done! Creating classes of DEG... (", as.character(Sys.time()),")","\n")

merged["group"]<-"Not DEG"

pval_thresh=args[7]

FC_thresh=as.numeric(args[8])
FC_thresh_minus= - as.numeric(paste(args[8]))

merged[which(merged['pval.x'] < pval_thresh & abs(merged['FC.x']) > FC_thresh |
merged['pval.y'] < pval_thresh & abs(merged['FC.y']) > FC_thresh),
"group"]<- "DEG"

merged[which(merged['pval.x'] < pval_thresh & (merged['FC.x']) > FC_thresh &
merged['pval.y'] < pval_thresh & (merged['FC.y']) > FC_thresh),
"group"]<- paste("Up", name1, "/", "Up", name2, sep=" ")

merged[which(merged['pval.x'] < pval_thresh & (merged['FC.x']) > FC_thresh &
merged['pval.y'] < pval_thresh & (merged['FC.y']) < FC_thresh_minus),
"group"]<- paste("Up", name1, "/", "Down", name2, sep=" ")

merged[which(merged['pval.x'] < pval_thresh & (merged['FC.x']) < FC_thresh_minus &
merged['pval.y'] < pval_thresh & (merged['FC.y']) < FC_thresh_minus),
"group"]<- paste("Down", name1, "/", "Down", name2, sep=" ")

merged[which(merged['pval.x'] < pval_thresh & (merged['FC.x']) < FC_thresh_minus &
merged['pval.y'] < pval_thresh & (merged['FC.y']) > FC_thresh),
"group"]<- paste("Down", name1, "/", "Up", name2, sep=" ")


cat("Done! (", as.character(Sys.time()),")","\n", "Number of genes in each category:","\n")
summary(as.factor(merged$group))



cat("Loading GO annotations, this will take some time... (", as.character(Sys.time()),")","\n")
##### GOID importation

ensembl=useMart("ensembl") #load ensembl db
#check for available datasets and their name
listDatasets(ensembl)
ensembl = useDataset(args[9],mart=ensembl) #import dmel datasets
attributes = listAttributes(ensembl) #list of attributes you can use for your query

# here we do the query to generate the raw GO annotation table
goids <- getBM(attributes = c('ensembl_gene_id','name_1006', 'go_linkage_type'),
mart = ensembl)

#  Merging of GO names and evidence code in the same cell
GO<-cbind(goids$ensembl_gene_id,paste(goids$name_1006,"(",goids$go_linkage_type,")",sep = ""))

# Aggregation of all GO information per gene, ",</br>" is used afterward by plot_ly for text display
GO<-aggregate(GO[,2]~GO[,1], GO, paste, collapse=",</br>")
GO[,1]<-as.factor(GO[,1])
colnames(GO)<-c('rows','GO')

merged<-merge(merged,GO,by="rows",all.x = T)

cat("Done! Plotting the data now... (", as.character(Sys.time()),")","\n")

##### Plot

# Font style for axis titles
f1 <- list(
family = "Arial, sans-serif",
size = 15,
color = "black"
)

# X axis
b <- list(
title = paste("Log2FC (", name1,")", sep=""),
titlefont = f1,
showticklabels = TRUE,
exponentformat = "E"
)

# Y axis
c <- list(
title = paste("Log2FC (", name2,")", sep=""),
titlefont = f1,
showticklabels = TRUE,
exponentformat = "E"
)


p1<-plot_ly(data = merged,
x = merged$FC.x,
y = merged$FC.y,
mode = "markers",
opacity = 0.6,
color = merged$group,
text =  ~paste('<b>Name: </b>', merged$genes.y,
'</br><b>Description: </b>', merged$desc.y,
'</br><b>p-value MOF: </b>', merged$pval.x,
'</br><b>p-value K2: </b>', merged$pval.y,
'</br><b>GOs: </b>', merged$GO)) %>%
layout(title = paste(name1, " vs ", name2, sep="")) %>%
layout(xaxis=b) %>%
layout(yaxis=c)

output=args[10]

htmlwidgets::saveWidget(as_widget(p1), paste(output, ".html", sep=""))

cat("Done! The graphic is stored in ", output,".html (", as.character(Sys.time()),")","\n", "Enjoy your data exploration now :)", "\n", sep="")
quit()
