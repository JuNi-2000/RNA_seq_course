library(BiocManager)
library(tidyverse)
library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)



################# DE SEQ ###############################################

counts_data <- read.csv("feature_counts_formatted.tsv", sep = "\t", row.names = 1)
col_data <- read.csv("column_data.csv", sep = ";",  row.names = 1)



#construct DESeq dataset:
#coldata should be constructed with a single group with sample names rather than 
#two different groups with tissue and condition
dds <- DESeqDataSetFromMatrix(countData = counts_data, 
                              colData = col_data,
                              design = ~ Group)
#inspect dds object
dds


#prefiltering rows with less than 10 counts
keep <- rowSums(counts(dds)) >= 10
dds <-  dds[keep,]


#Quality control: Plot pca. Remove dependence on the variance (is automatically done by DEseq2 afterwards,
# we do it before to do exploratory data analysis because DEseq is not done yet.)
vst <- vst(dds, blind = TRUE)
plot_pca <- plotPCA(vst, intgroup = c("Group"), returnData = TRUE)

#make a nice plot
percentVar <- round(100 * attr(plot_pca, "percentVar"))
ggplot(plot_pca, aes(PC1, PC2, color=Group)) +
  geom_point(size=2) +
  ggtitle("Principal Component Analysis")+
  theme(plot.title = element_text(hjust = 0.5))+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()+
  theme(
    panel.background = element_rect(fill = "#BFD5E3", colour = "black",
                                    size = 1, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "white"), 
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                    colour = "white")
  )



#run DESeq
dds <- DESeq(dds)


######################## plot some normalized genes of interest (can be skipped for DE analysis) ##############################

GbP5 <- plotCounts(dds, 'ENSMUSG00000105504', intgroup = c('Group'), normalized = TRUE, returnData = TRUE)
boxplot(count ~ Group, data = GbP5, ylab = 'gene count', main = 'Normalized gene counts for GbP5')

#some descriptive stats about the gene 
GbP5.mean.lung.case<- mean(GbP5$count[GbP5$Group == 'Lung_Case'])
GbP5.mean.lung.control <- mean(GbP5$count[GbP5$Group == 'Lung_Control'])
GbP5.log2fold.Blood <- log2(GbP5.mean.lung.case/GbP5.mean.lung.control)
GbP5.sd.lung.case <- sd(GbP5$count[GbP5$Group == 'Lung_Case'])
GbP5.sd.lung.control <-sd(GbP5$count[GbP5$Group == 'Lung_Control'])

GbP5.mean.Blood.case<- mean(GbP5$count[GbP5$Group == 'Blood_Case'])
GbP5.mean.Blood.control <- mean(GbP5$count[GbP5$Group == 'Blood_Control'])
GbP5.log2fold.Blood <- log2(GbP5.mean.Blood.case/GbP5.mean.Blood.control)
GbP5.sd.Blood.case <- sd(GbP5$count[GbP5$Group == 'Blood_Case'])
GbP5.sd.Blood.control <-sd(GbP5$count[GbP5$Group == 'Blood_Control'])

#nicer plot for the gene counts:
ggplot(GbP5, aes(x = Group, y = count)) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "red") +
  ylab('gene count') +
  ggtitle('Normalized gene counts for GbP5')+
  scale_y_continuous(limits = c(0, 25000, 5000))+
  theme(
    panel.background = element_rect(fill = "#BFD5E3", colour = "black",
                                    size = 1, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "white"), 
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                    colour = "white")
  )



#Oas3 expression
Oas3 <- plotCounts(dds, 'ENSMUSG00000032661', intgroup = c('Group'), normalized = TRUE, returnData = TRUE)
boxplot(count ~ Group, data = Cxcl9, ylab = 'gene count', main = 'Normalized gene counts for Cxcl9')


#some descriptive stats about the gene 
Oas3.mean.lung.case<- mean(Oas3$count[Oas3$Group == 'Lung_Case'])
Oas3.mean.lung.control <- mean(Oas3$count[Oas3$Group == 'Lung_Control'])
Oas3.log2fold.lung <- log2(Oas3.mean.lung.case/Oas3.mean.lung.control)
Oas3.sd.lung.case <- sd(Oas3$count[Oas3$Group == 'Lung_Case'])
Oas3.sd.lung.control <-sd(Oas3$count[Oas3$Group == 'Lung_Control'])

Oas3.mean.Blood.case<- mean(Oas3$count[Oas3$Group == 'Blood_Case'])
Oas3.mean.Blood.control <- mean(Oas3$count[Oas3$Group == 'Blood_Control'])
Oas3.log2fold.blood <- log2(Oas3.mean.Blood.case/Oas3.mean.Blood.control)
Oas3.sd.Blood.case <- sd(Oas3$count[Oas3$Group == 'Blood_Case'])
Oas3.sd.Blood.control <-sd(Oas3$count[Oas3$Group == 'Blood_Control'])

#nicer plot for the gene counts:
ggplot(Oas3, aes(x = Group, y = count)) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "red") +
  ylab('gene count') +
  ggtitle('Normalized gene counts for Oas3')+
  scale_y_continuous(limits = c(0, 10000, 2000))+
  theme(
    panel.background = element_rect(fill = "#BFD5E3", colour = "black",
                                    size = 1, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "white"), 
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                    colour = "white")
  )




############################ RESULTS LUNG #############################

#extract the results from DESeq2
res_Lung <- results(dds, contrast = c("Group", "Lung_Case", "Lung_Control"))
res_Blood <- results(dds, contrast = c("Group", "Blood_Case", "Blood_Control"))


#convert to data frame for volcano plot
res_Lung <- as.data.frame(res_Lung)
res_Lung <- cbind(res_Lung, rownames(res_Lung))


#export the genes in a file so I can convert them to standard notation
write.csv(export, file = "genes_lung.txt", quote = FALSE, row.names = FALSE)

#import the new names
genes_lung <- read.csv("Lung_gene_names_results.csv", sep = ";", header = TRUE)
res_Lung <- cbind(res_Lung, gene_name = genes_lung$new)


# set default to Not diff expressed
res_Lung$diffexpressed <- FALSE
res_Lung$change <- "NO"

#mark upregulated genes
res_Lung$diffexpressed[res_Lung$log2FoldChange>1 & res_Lung$padj <= 0.05] <- TRUE
res_Lung$change[res_Lung$log2FoldChange>1 & res_Lung$padj <= 0.05] <- "UP"


#mark downregulated genes
res_Lung$diffexpressed[res_Lung$log2FoldChange<(-1) & res_Lung$padj <= 0.05] <- TRUE
res_Lung$change[res_Lung$log2FoldChange<(-1) & res_Lung$padj <= 0.05] <- "DOWN"


#which genes are significantly differnt?
padj_res_Lung <- res_Lung[which(res_Lung$padj <= 0.05),]

#How many genes are siginificantly differentially expressed?
padj_res_Lung[which(padj_res_Lung$diffexpressed == TRUE),]

#which genes are upregulated?
nrow(padj_res_Lung[which(padj_res_Lung$change == "UP"),])

#which genes are downregulated?
nrow(padj_res_Lung[which(padj_res_Lung$change == "DOWN"),])



#create the volcano plot for all the genes
EnhancedVolcano(padj_res_Lung,
                lab = padj_res_Lung$gene_name,
                x = 'log2FoldChange',
                y = 'pvalue',
                drawConnectors = TRUE,
                boxedLabels = TRUE)



###################### ENRICHMENT ANALYSIS LUNG #######################

library(org.Mm.eg.db)
library(clusterProfiler)
all_genes <- rownames(res_Lung) #lung and blood have the same gene list
significant_genes_lung <- res_Lung[res_Lung$padj <=0.05 & (res_Lung$diffexpressed == TRUE),]
significant_genes_lung <- rownames(na.omit(significant_genes_lung))


GO_results <- enrichGO(gene = significant_genes_lung,
                       universe = all_genes,
                       OrgDb = org.Mm.eg.db, 
                       keyType = "ENSEMBL",
                       ont = "BP",
                       readable = TRUE)


#convert to dataframe
as.data.frame(GO_results)

barplot(GO_results,
        showCategory = 15)+
  labs(title = "GO enrichment biological processes in Lung Tox vs WT")

#dotplot shows number in ascending order
dotplot(GO_results, showCategory=15) + ggtitle("Dotplot BP in Lung Tox vs WT")


############################ RESULTS BLOOD ###################

#convert to dataframe
res_Blood <- as.data.frame(res_Blood)

# set default to Not diff expressed
res_Blood$diffexpressed <- FALSE
res_Blood$change <- "NO"

#mark upregulated genes
res_Blood$diffexpressed[res_Blood$log2FoldChange>1 & res_Blood$padj <= 0.05] <- TRUE
res_Blood$change[res_Blood$log2FoldChange>1 & res_Blood$padj <= 0.05] <- "UP"

#mark downregulated genes
res_Blood$diffexpressed[res_Blood$log2FoldChange<(-1) & res_Blood$padj <= 0.05] <- TRUE
res_Blood$change[res_Blood$log2FoldChange<(-1) & res_Blood$padj <= 0.05] <- "DOWN"


#bind gene names to dataframe
res_Blood <- cbind(res_Blood, gene_name = genes_lung$new)


#significant genes
padj_res_Blood <- res_Blood[which(res_Blood$padj <= 0.05),]

#How many genes are siginificantly differentially expressed?
nrow(padj_res_Blood[which(padj_res_Blood$diffexpressed == TRUE),])

#how many genes are upregulated?
nrow(padj_res_Blood[which(padj_res_Blood$change == "UP"),])

#how many are downregulated?
nrow(padj_res_Blood[which(padj_res_Blood$change == "DOWN"),])



#create the volcano plot for all the genes
EnhancedVolcano(padj_res_Blood,
                lab = padj_res_Blood$gene_name,
                title = 'Results Blood Toxoplasma vs WT',
                x = 'log2FoldChange',
                y = 'pvalue',
                drawConnectors = TRUE,
                boxedLabels = TRUE)


#investigate single genes # exploratory, not necessary
EnhancedVolcano(res_Blood,
                lab = res_Blood$gene_name,
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'Results Blood Toxoplasma vs WT',
                selectLab = c('Cyp11a1'),
                xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = 10e-14,
                FCcutoff = 2.0,
                pointSize = 4.0,
                labSize = 6.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                colAlpha = 4/5,
                legendPosition = 'right',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black')




##################### ENRICHMENT ANALYISIS BLOOD ###################

#now for blood:
significant_genes_blood <- res_Blood[res_Blood$padj <=0.05 & (res_Blood$diffexpressed == TRUE),]
significant_genes_blood <- rownames(na.omit(significant_genes_blood))

GO_results_blood <- enrichGO(gene = significant_genes_blood,
                       universe = all_genes,
                       OrgDb = org.Mm.eg.db, 
                       keyType = "ENSEMBL",
                       ont = "BP",
                       readable = TRUE)

#look at the results
as.data.frame(GO_results_blood)

#plot:

barplot(GO_results_blood, showCategory = 15)+
  labs(title = "GO enrichment biological processes in blood Tox vs WT")

dotplot(GO_results_blood, showCategory=15) + ggtitle("Dotplot BP in Blood Tox vs WT")







