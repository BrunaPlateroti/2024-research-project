## LOAD NECESSARY LIBRARY ##
library(data.table)
library(ggplot2)
library(ggrepel)

# Task 1 --> Identify GPCRs and TFs essential for the survival of cancer cells 
# Dataset for gene dependency and dataset for TFs and RTKs

## Loading the DepMapdata and gene name data ##
data_gene<- read.csv("/Users/bruna/Documents/University/Master/Research project/data/Selected_genesets.csv")
structure(data_gene)

DepMap_RNAi<- fread("/Users/bruna/Documents/University/Master/Research Project/data/DepMap_RNAi_2019_datasets.csv")
structure(DepMap_RNAi)

colnames(DepMap_RNAi) # columns = gene names
rownames(DepMap_RNAi) # rows = cell lines

# First - Get GPCR genes from data_gene table
gpcrs<- data_gene[data_gene$gene.set=="GPCRs",]$genes # Select only the GPCRs genes from the table
gpcrs<- unlist(strsplit(gpcrs, ";")) # Use unlist() to undo the list you exported 
#from the table and create a vector

# Select genes that contain data from the DepMap_RNAi table
gpcrs_with_data <- intersect(colnames(DepMap_RNAi), gpcrs) # with intersect() select 
# columns from the DepMap data that match data in gpcrs

# Create a table with genes from the gcprs_with_data
data_rnai_gpcrs<- DepMap_RNAi[,..gpcrs_with_data]

# Calculate mean, sd, and establish cutoff
DepMap_RNAi_copy <- DepMap_RNAi[,5:ncol(DepMap_RNAi)] # remove the first 4 columns

DepMap_RNAi_copy[is.na(DepMap_RNAi_copy)] <-0 # use is.na() to replace the NA values
# with the value zero

DepMap_RNAi_impact <- data.frame(below.minus.0.5=colSums(DepMap_RNAi_copy<(-0.5)),
                                 average.score= apply(DepMap_RNAi_copy,2,median))
# average.score - apply() avoid the loop and will apply the median to the matrix and 
# return the median of each columns
# colSums is the sum of the values in each column

DepMap_RNAi_impact[order(-DepMap_RNAi_impact$below.minus.0.5),] # ascending order

av_impact <- median(DepMap_RNAi_impact$average.score)
sd_impact <- sd(DepMap_RNAi_impact$average.score)

cutoff <- av_impact+(sd_impact*3)

#GPCR genes important for the survival of some cells
data_rnai_gpcrs[is.na(data_rnai_gpcrs)] <- 0

# cells that are kill/grow slow when the gene is deleted
data_rnai_impact_gpcrs <- data.frame(below.minus.0.5= colSums(data_rnai_gpcrs<(-cutoff)),
                                     average.score=apply(data_rnai_gpcrs,2,mean))
# colSums(data_rnai_gpcrs<(-cutoff) - calculates the count of elements in each columns 
# that meet the conditions
data_rnai_impact_gpcrs[order(-data_rnai_impact_gpcrs$below.minus.0.5),]

#add a column with all the values above the cutoff - Cells that can grow when a gene is deleted
data_rnai_impact_gpcrs$above.plus.0.5 <- colSums(data_rnai_gpcrs>cutoff)

data_rnai_impact_gpcrs <- data_rnai_impact_gpcrs[order(-data_rnai_impact_gpcrs$below.minus.0.5),]

# genes with significant impact
gpcrs_gene_with_significant_impact <- 
  rownames(data_rnai_impact_gpcrs[data_rnai_impact_gpcrs$below.minus.0.5>10 | 
                                    data_rnai_impact_gpcrs$above.plus.0.5>10,])

data_gpcrs_selected <- DepMap_RNAi[,..gpcrs_gene_with_significant_impact]

data_gpcrs_m <- reshape2::melt(data_gpcrs_selected)# melt() convert df with several measuraments
# columns in this canonical format - one row for every observed value
# Reshape2 package works with 2 functions: melt (take wide format data and melt into a long format data)
# and cast (the other way around)

# plot a box plot
par(mfrow=c(2,2)) # to divide the screen

plot_gpcrs <- ggplot(data_gpcrs_m, mapping=aes(x=variable,y=value))+ geom_boxplot()+
  geom_hline(yintercept = 0, linetype=1)+ geom_hline(yintercept = -0.5, linetype=2)+
  geom_hline(yintercept = 0.5, linetype=2)+ 
  theme(axis.text.x = element_text(angle = 90, vjust =0.5, hjust = 1 ))
plot_gpcrs

# Second - Get TFs genes from data_gene table
tfs<- data_gene[data_gene$gene.set=="TFs",]$genes # Select only the TFs genes from the table
tfs<- unlist(strsplit(tfs, ";")) # Use unlist() to undo the list you exported 
#from the table and create a vector

# Select genes that contain data from the DepMap_RNAi table
tfs_with_data <- intersect(colnames(DepMap_RNAi), tfs)

# Create a table with genes from the tfs_with_data
df_rnai_tfs<- DepMap_RNAi[,..tfs_with_data]

#TFs genes important for the survival of some cells
df_rnai_tfs[is.na(df_rnai_tfs)] <- 0

# cells that are kill/grow slow when the gene is deleted
df_rnai_impact_tfs <- data.frame(below.minus.0.5= colSums(df_rnai_tfs<(-cutoff)),
                                     average.score=apply(df_rnai_tfs,2,mean))

df_rnai_impact_tfs[order(-df_rnai_impact_tfs$below.minus.0.5),]

#add a column with all the values above the cutoff - Cells that can grow when a gene is deleted
df_rnai_impact_tfs$above.plus.0.5 <- colSums(df_rnai_tfs>cutoff)

df_rnai_impact_tfs <- df_rnai_impact_tfs[order(-df_rnai_impact_tfs$below.minus.0.5),]

# genes with significant impact
tfs_gene_with_significant_impact <- 
  rownames(df_rnai_impact_tfs[df_rnai_impact_tfs$below.minus.0.5>10 | 
                                    df_rnai_impact_tfs$above.plus.0.5>10,])

df_tfs_selected <- DepMap_RNAi[,..tfs_gene_with_significant_impact]

df_tfs_m <- reshape2::melt(df_tfs_selected)

# plot a box plot
plot_tfs <- ggplot(df_tfs_m, mapping=aes(x=variable,y=value))+ geom_boxplot()+
  geom_hline(yintercept = 0, linetype=1)+ geom_hline(yintercept = -0.5, linetype=2)+
  geom_hline(yintercept = 0.5, linetype=2)+ 
  theme(axis.text.x = element_text(angle = 90, vjust =0.5, hjust = 1 ))
plot_tfs

# Third - Get G protein genes from data_gene table
Gprotein<- data_gene[data_gene$gene.set=="G_proteins",]$genes # Select only the G_protein genes from the table
Gprotein<- unlist(strsplit(Gprotein, ";")) # Use unlist() to undo the list you exported 
#from the table and create a vector

# Select genes that contain data from the DepMap_RNAi table
G_protein_with_data <- intersect(colnames(DepMap_RNAi), Gprotein)

# Create a table with genes from the G_protein_with_data
G_protein_rnai_df<- DepMap_RNAi[,..G_protein_with_data]

# G_protein genes important for the survival of some cells
G_protein_rnai_df[is.na(G_protein_rnai_df)] <- 0

# cells that are kill/grow slow when the gene is deleted
G_protein_rnai_impact_df <- data.frame(below.minus.0.5= colSums(G_protein_rnai_df<(-cutoff)),
                                 average.score=apply(G_protein_rnai_df,2,mean))

G_protein_rnai_impact_df[order(-G_protein_rnai_impact_df$below.minus.0.5),]

#add a column with all the values above the cutoff - Cells that can grow when a gene is deleted
G_protein_rnai_impact_df$above.plus.0.5 <- colSums(G_protein_rnai_df>cutoff)

df_rnai_impact_tfs <- df_rnai_impact_tfs[order(-df_rnai_impact_tfs$below.minus.0.5),]

# genes with significant impact
G_protein_with_significant_impact <- 
  rownames(G_protein_rnai_impact_df[G_protein_rnai_impact_df$below.minus.0.5>10 | 
                                      G_protein_rnai_impact_df$above.plus.0.5>10,])

G_protein_selected <- DepMap_RNAi[,..G_protein_with_significant_impact]

G_protein_m <- reshape2::melt(G_protein_selected)

# plot a box plot
plot_G_protein <- ggplot(G_protein_m, mapping=aes(x=variable,y=value))+ geom_boxplot()+
  geom_hline(yintercept = 0, linetype=1)+ geom_hline(yintercept = -0.5, linetype=2)+
  geom_hline(yintercept = 0.5, linetype=2)+ 
  theme(axis.text.x = element_text(angle = 90, vjust =0.5, hjust = 1 ))
plot_G_protein

# example for HTR genes

