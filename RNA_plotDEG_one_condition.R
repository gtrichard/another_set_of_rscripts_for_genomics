#### RNA-seq_plot_1_condition.R
## This script allows you to plot a nice interactive Volcano Plot.
## Gene names, GO term and gene description are aumotically retrieved 
## from the Ensembl Genes ID given as input, alongside their FC and 
## p-value after DEG Analysis.

# install.packages('plotly')
# source("https://bioconductor.org/biocLite.R")
# biocLite("biomaRt")
# quit()

library('biomaRt')
library('plotly')

setwd("")

table1<-read.csv(file="", sep="\t",header=T)
table1<-data.frame(table1$ID,table1$name,table1$log2FoldChange,table1$padj)
colnames(table1)<-c("ID", "genes","FC","pval")

pval_thresh=0.05
FC_thresh=as.numeric(1.5)
FC_thresh_minus= - as.numeric(paste(1.5))

merged<-table1
merged["group"]<-"Not DEG"
merged[which(merged['pval'] < pval_thresh),"group"]<- "p-value OK"
merged[which(abs(merged['FC']) > FC_thresh),"group"]<- "FC OK"
merged[which(merged['pval'] < pval_thresh & 
             abs(merged['FC']) > FC_thresh),"group"]<- "DEG"
merged[,"group"]<-as.factor(merged[,"group"])


cat("Done! (", as.character(Sys.time()),")","\n", "Number of genes in each category:","\n")
summary(as.factor(merged$group))

cat("Loading GO annotations, this will take some time... (", as.character(Sys.time()),")","\n")
##### GOID importation

ensembl=useMart("ensembl") #load ensembl db
#check for available datasets and their name
listDatasets(ensembl)
ensembl = useDataset("ecaballus_gene_ensembl",mart=ensembl) #import dmel datasets
attributes = listAttributes(ensembl) #list of attributes you can use for your query

# here we do the query to generate the raw GO annotation table
goids <- getBM(attributes = c('ensembl_gene_id','name_1006','namespace_1003'),mart = ensembl)
GO<-goids[which(goids$namespace_1003=="biological_process"),c(1,2)]

#  Merging of GO names and evidence code in the same cell
# GO<-cbind(goids$ensembl_gene_id,paste(goids$name_1006,"(",goids$go_linkage_type,")",sep = "" ))

# Aggregation of all GO information per gene, ",</br>" is used afterward by plot_ly for text display
GO<-aggregate(GO[,2]~GO[,1], GO, paste, collapse=", ")
GO[,1]<-as.factor(GO[,1])
colnames(GO)<-c("ID",'GO')

merged<-merge(merged,GO,by="ID",all.x = T)

# Add a line break every 100 characters to avoid display issues in plotly
merged["GO_new"]<-gsub("(.{100})", "\\1</br>", merged$GO) 

cat("Done! Plotting the data now... (", as.character(Sys.time()),")","\n")

##### Plot

# Font style for axis titles
#f1 <- list(
#family = "Arial, sans-serif",
#size = 15,
#color = "black"
#)

# X axis
#b <- list(
#title = paste("Log2FC (", name1,")", sep=""),
#titlefont = f1,
#showticklabels = TRUE,
#exponentformat = "E"
#)

# Y axis
#c <- list(
#title = paste("Log2FC (", name2,")", sep=""),
#titlefont = f1,
#showticklabels = TRUE,
#exponentformat = "E"
#)

library('plotly')

p1<-plot_ly(data = merged,
            domain=c(-2,2),
            x = merged$FC,
            y = -log(merged$pval,10),
            mode = "markers",
            trace = "scatter",
            opacity = 0.8,
            color = merged$complex,
#            hoverinfo= "x+y+text",
            textposition= "outside",
            colors = c("#999999", "#FF4C85", "#B7002A", "#FF4747", "#FFAC59", "#FF6B1C"),
#            colors = c("#999999", "#7C75FF", "#FF59D5", "#FFB638", "#51B6FF", "#FF4747"),
            text =  ~paste('<b>Name: </b>', merged$genes,
#                           '</br><b>Description: </b>', merged$desc,
#                           '</br><b>p-value: </b>', merged$pval,
#                           '</br><b>DEG-type: </b>', merged$group,
                           '</br><b>GOs: </b>', merged$GO_new))


line <- list()
lines <- list()
for (i in c(FC_thresh, FC_thresh_minus)) {
  line[["x0"]] <- i
  line[["x1"]] <- i
  line[["y0"]] <- 0
  line[["y1"]] <- max(-log(merged$pval,10))
  lines <- c(lines,list(line))
}

lines[[1]]$line$color="rgba(0,0,0,0.3)"
lines[[2]]$line$color="rgba(0,0,0,0.3)"
lines[[1]]$line$dash="dash"
lines[[2]]$line$dash="dash"

pf<-layout(p1,shapes=lines,
           xaxis = list(title = "Log2(KO/Control)"), 
           yaxis = list(title = "-Log10(adj. P-value)"))

pf #Display the plot within RStudio

htmlwidgets::saveWidget(as_widget(pf), paste(output, ".html", sep="")) #Save the plot in current working directory as .html file

cat("Done! The graphic is stored in ", output,".html (", as.character(Sys.time()),")","\n", "Enjoy your data exploration now :)", "\n", sep="")
quit()

# You can export that plot to your online plotly account. Create an account there and retrieve
# your API key to use that function

Sys.setenv("plotly_username"="blabla")
Sys.setenv("plotly_api_key"="lookituponline")
api_create(pf, filename = "KO_control_RNA-seq")
