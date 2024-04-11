# run the differential gene expression analysis, 
# and the associated functional enrichment tests, 
# and then spend the rest of the lab interpreting the outputs 
# and determining how to frame your poster and what information to include


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

library("BiocManager")
BiocManager::install("clusterProfiler", force = TRUE)
BiocManager::install("org.Hs.eg.db")
nBiocManager::install("DOSE")
BiocManager::install("pathview")
BiocManager::install("enrichplot")
install.packages("ggupset")

library("clusterProfiler")
library("org.Hs.eg.db")
library("DOSE")
library("pathview")
library("dplyr")
library("enrichplot")
library("ggupset")



#------#------#------#------#------#------#------#------#------#------

# Load data
sampleFiles <- grep("CD",list.files(),value=TRUE)
all_files2 <- read.csv("lab9_all_counts.csv")
all_files2$X <- NULL
annotations_ahb <- read.csv("annotations_ahb.csv")

#------#------#------#------#------#------#------#------#------#------#------#------
# Get significance genes via DESeq2
# Assuming all_files2 is already loaded into your R session
# Extracting the gene names into a vector
library(DESeq2)
genes <- all_files2$CD_C1_H.Gene

# Removing the gene names from the count data to create a pure count matrix
counts <- all_files2[, -1]

# Setting the row names of the count matrix to be the gene names
rownames(counts) <- genes

# Creating a sample information data frame
# This should reflect the conditions each sample belongs to
sampleInfo <- data.frame(
  sampleName = colnames(counts),
  condition = rep(c("CD", "HS"), each = 3)
)



# Convert sample information to a DataFrame (DESeq2 class) for compatibility
sampleInfo <- DataFrame(sampleInfo)
# Creating the DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = sampleInfo,
                              design = ~ condition)
# Running the DESeq pipeline
dds <- DESeq(dds)
# Extracting results
res <- results(dds)

# Filtering for significant results, e.g., adjusted p-value < 0.05
sigGenes <- res[which(res$padj < 0.05), ]



# Create background dataset for hypergeometric testing using all genes tested for significance in the results
all_genes <- as.character(all_files2$CD_C1_H.Gene)

# Rename the genes so that they match the ENSEMBL database:
all_genes2<-as.character(annotations_ahb$gene_id[match(all_genes, annotations_ahb$gene_name)])




# Extract significant results (adjust the p-value to suit your data!)
signif_res <- res[res$padj < 0.05 & !is.na(res$padj), ]
signif_genes <- as.character(rownames(signif_res))

# Rename the significant genes to match ENSEMBL
signif_genes2 <- annotations_ahb$gene_id[match(signif_genes, annotations_ahb$gene_name)]


#------#------#------#------#------#------#------#------#------#------#------#------
#------#------#------#------#------#------#------#------#------#------#------#------



# Run GO enrichment analysis
ego <- enrichGO(gene = signif_genes2, 
                universe = all_genes2,
                keyType = "ENSEMBL",
                OrgDb = org.Hs.eg.db, 
                ont = "BP", 
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.05, 
                readable = TRUE)

# Output results from GO analysis to a table
## Output results from GO analysis to a table
cluster_summary <- data.frame(ego)
cluster_summary



#------#------#------#------#------#------#------#------#------#------#------#------
#------#------#------#------#------#------#------#------#------#------#------#------


mutate(ego, qscore = -log(p.adjust, base=10)) %>% 
  barplot(x="qscore")



#------#------#------#------#------#------#------#------#------#------#------#------
#------#------#------#------#------#------#------#------#------#------#------#------




dotplot(ego, showCategory=20)

dotplot(ego, showCategory=10)

#------#------#------#------#------#------#------#------#------#------#------#------
#------#------#------#------#------#------#------#------#------#------#------#------




## Enrichmap clusters the 20 most significant (by padj) GO terms to visualize relationships between terms
emapplot(pairwise_termsim(ego), showCategory = 10)

#------#------#------#------#------#------#------#------#------#------#------#------
#------#------#------#------#------#------#------#------#------#------#------#------




## To color genes by log2 fold changes, we need to extract the log2 fold changes from our results table creating a named vector
OE_foldchanges <- signif_res$log2FoldChange

names(OE_foldchanges) <- annotations_ahb$gene_id[match(rownames(signif_res), annotations_ahb$gene_name)]

## Cnetplot details the genes associated with one or more terms - by default gives the top 5 significant terms (by padj)
cnetplot(ego, 
         categorySize="pvalue", 
         showCategory = 5, 
         foldChange=OE_foldchanges, 
         vertex.label.font=6)

#------#------#------#------#------#------#------#------#------#------#------#------
#------#------#------#------#------#------#------#------#------#------#------#------


heatplot(ego, foldChange=OE_foldchanges, showCategory=20)

upsetplot(ego, showCategory=20)


