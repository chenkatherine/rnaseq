#' Performs normalization to housekeeping genes, differential gene expression 
#' (DGE), and visualizes DGE results with PCA plot and heatmap. Incorporates 
#' biomaRt to retrieve gene symbol and description of protein coding genes only 
#' for biological interpretation.

#' Load libraries
library(tidyverse)
library(DESeq2)
library(pheatmap)
library(openxlsx)
library(ggrepel)
library(ggplot2)
library(ggfortify)
library(limma)
library(RUVSeq)

#' Creates vector with genes of interest for project specific analyses
genes_of_interest <- c("ENSG00000248098.12", # branched chain keto acid dehydrogenase E1 
                       "ENSG00000112592.15", # TATA-box binding protein
                       "ENSG00000134852.15", #CLOCK
                       "ENSG00000179094.16", #PER1
                       "ENSG00000132326.13", #PER2
                       "ENSG00000049246.15", #PER3
                       "ENSG00000133794.20", #BMAL1
                       "ENSG00000008405.12", #CRY1
                       "ENSG00000121671.12", #CRY2
                       "ENSG00000008710.20", #PKD1
                       "ENSG00000118762.8" #PKD2
                      )

#' Creates vector with housekeeping genes for RUVSeq
housekeeping_genes <- c("ENSG00000248098.12", "ENSG00000112592.15") 
# branched chain keto acid dehydrogenase E1, TATA-box binding protein

#' Read in count matrix and annotations file.
#' Annotations file contains columns "sample", "condition, "type", "batch"
counts <- read.csv("procode_counts.csv")
annotations <- read.csv("annotations.csv")

#' Removes headers to equalize column count and row count
gene_ids <- counts[, 1]
counts_matrix <- as.matrix(counts[, -1])
rownames(counts_matrix) <- gene_ids

rownames(annotations) <- annotations$sample
annotations_sub <- annotations[, c("sample", "condition")]

#' Performs normalization with RUVSeq using the specified housekeeping genes
#' Extracts the normalized counts for DESeq
seq <- newSeqExpressionSet(counts = counts_matrix, 
                           phenoData = AnnotatedDataFrame(data = annotations_sub))

set <- RUVg(x = seq, 
            cIdx = housekeeping_genes, 
            k=1)

normalized_counts <- assayData(set)$normalizedCounts

#' Creates DESeqDataSet and runs DESeq2
#' Saves output to xlsx
dds <- DESeqDataSetFromMatrix(countData = normalized_counts,
                              colData = annotations,
                              design = ~ batch + condition)
dds <- DESeq(dds, minReplicatesForReplace = Inf)
res <- results(dds, cooksCutoff = FALSE, independentFiltering = TRUE)
write.xlsx(as.data.frame(res), "deseq_t500_procode_ruv.xlsx")

# Graphs PCA with normalized and batch-removed data
# Note: DESeq2 vignette suggests using removeBatchEffect only for visualization
# since ~ batch + condition param in dds already fixes that issue
vsd <- vst(dds, blind = FALSE)
vsd_corrected <- removeBatchEffect(assay(vsd), batch = annotations$batch)
pca <- prcomp(t(vsd_corrected))
autoplot(pca, data = annotations, color = "condition", label = TRUE)

#' Employs biomaRt to retrieve gene symbols and descriptions from Ensembl
#' Saves scraping result to xlsx for local matching
library(biomaRt)
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
attributes <- c("ensembl_gene_id", "hgnc_symbol", "description")
gene_ids <- gene_ids <- counts[, 1] |> as.character()
gene_ids_clean <- sub("\\..*", "", gene_ids)
filters <- "ensembl_gene_id"
gene_info <- getBM(attributes = attributes,
                   filters = filters,
                   values = gene_ids_clean,
                   mart = ensembl)
write.xlsx(as.data.frame(gene_info), "gene_descriptors.xlsx")

#' Creates heatmap to visualize differentially expressed genes
#' Loads in the top 500 differentially expressed genes by padj for analyses
#' File previously created by hand (extracting the top 500 after sorting by padj)
t500 <- as.data.frame(read.xlsx("deseq_t500_procode.xlsx"))
t500_gene_ids <- t500$GeneID

genes_to_graph <- t500_gene_ids[t500_gene_ids %in% rownames(vsd)]
mat <- vsd_corrected[genes_to_graph, ]
mat_scaled <- t(scale(t(mat)))

sample_annots <- data.frame(
  condition = annotations$condition,
  batch = annotations$batch
)


rownames(sample_annots) <- colnames(mat)
pheatmap(
  mat_scaled,
  annotation_col = sample_annots,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 8,
  scale = "none",
  main = "Heatmap of Genes of Interest (Batch-corrected VST)"
)
