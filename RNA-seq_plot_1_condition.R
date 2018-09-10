#### RNA-seq_plot_1_condition.R

# install.packages('plotly')
# source("https://bioconductor.org/biocLite.R")
# biocLite("biomaRt")
# quit()

library('biomaRt')
library('plotly')

setwd("")

table1<-read.csv(file="", sep="\t",header=T)
table1<-data.frame(table1$ID,table1$name,table1$log2FoldChange,table1$padj,table1$Complex)
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
goids <- getBM(attributes = c('ensembl_gene_id','name_1006', 'go_linkage_type'),
mart = ensembl)

#  Merging of GO names and evidence code in the same cell
GO<-cbind(goids$ensembl_gene_id,paste(goids$name_1006,"(",goids$go_linkage_type,")",sep = "" ))

# Aggregation of all GO information per gene, ",</br>" is used afterward by plot_ly for text display
GO<-aggregate(GO[,2]~GO[,1], GO, paste, collapse=",</br>")
GO[,1]<-as.factor(GO[,1])
colnames(GO)<-c("ID",'GO')

merged<-merge(merged,GO,by="ID",all.x = T)

cat("Done! Plotting the data now... (", as.character(Sys.time()),")","\n")

##### Plot

# Font style for axis titles
f1 <- list(
family = "Arial, sans-serif",
size = 15,
color = "black"
)

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
            x = merged$FC,
            y = -log(merged$pval,10),
            mode = "markers",
            opacity = 0.8,
            color = merged$group,
            text =  ~paste('<b>Name: </b>', merged$genes,
#  Gene description if you have it '</br><b>Description: </b>', merged$desc,
                           '</br><b>p-value: </b>', merged$pval,
                           '</br><b>DE group: </b>', merged$group,
                           '</br><b>GOs: </b>', merged$GO))


line <- list()
lines <- list()
for (i in c(FC_thresh, FC_thresh_minus)) {
  line[["x0"]] <- i
  line[["x1"]] <- i
  line[["y0"]] <- 0
  line[["y1"]] <- 250
  lines <- c(lines,list(line))
}

lines[[1]]$line$color="rgba(0,0,0,0.3)"
lines[[2]]$line$color="rgba(0,0,0,0.3)"
lines[[1]]$line$dash="dash"
lines[[2]]$line$dash="dash"

pf<-layout(p,shapes=lines,xaxis = list(title = "Log2(KO/Control)"), yaxis = list(title = "-Log10(adj. P-value)"))

chart_link = api_create(p, filename = "public-graph")
chart_link

Sys.setenv("plotly_username"="blabla")
Sys.setenv("plotly_api_key"="lookituponline")
api_create(pf, filename = "KO_control_RNA-seq")

htmlwidgets::saveWidget(as_widget(p1), paste(output, ".html", sep=""))

cat("Done! The graphic is stored in ", output,".html (", as.character(Sys.time()),")","\n", "Enjoy your data exploration now :)", "\n", sep="")
quit()
