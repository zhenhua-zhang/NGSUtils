#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE, ggrepel.max.overlaps = Inf)

#
## TODO list
#

#
## Output check list
#
# [x] Expression TMP values
# [x] Quality control
#     [x] Sample PCA
#     [x] Sample to sample distance (by rlog)
#     [x] SD plots, vst (variance stabilizing transformation)
#     [x] Count matrix (by rlog)
#     [x] Top various genes (TODO)
# [x] Differential expression gene list: 
#     [x] MAplot, shrunken by apeglm (please cite their paper)
#     [x] Timepoint d1, low vs high
#     [x] Timepoint d5, low vs high 
#     [x] Timepoint d13, low vs high
#     [x] Speed high, d5 vs d1
#     [x] Speed high, d13 vs d5
#     [x] Speed low, d5 vs d1
#     [x] Speed low, d13 vs d5
# [x] Annotation
#     [x] GO annotation and enrichement

#
## Packages (install it if missing)
#
# 1. Packages host on CRAN
if (!require(pheatmap)) { install.packages("pheatmap"); require(pheatmap) }
if (!require(ggplot2)) { install.packages("ggplot2"); require(ggplot2) }
if (!require(ggrepel)) { install.packages("ggrepel"); require(ggrepel) }
if (!require(ggupset)) { install.packages("ggupset"); require(ggupset) }
if (!require(stringi)) { install.packages("stringi"); require(stringi) }
# if (!require(WGCNA)) { install.packages("WGCNA"); require(WGCNA) }
if (!require(ape)) { install.packages("ape"); require(ape) }
if (!require(DBI)) { install.packages("DBI"); require(DBI)}

# 2. Packages hosted by Bioconductor
if (!requireNamespace("BiocManager")) { install.packages("BiocManager"); requireNamespace(BiocManager) }
if (!require(AnnotationForge)) { BiocManager::install("AnnotationForge"); require(AnnotationForge) }
if (!require(clusterProfiler)) { BiocManager::install("clusterProfiler"); require(clusterProfiler) }
if (!require(AnnotationHub)) { BiocManager::install("AnnotationHub"); require(AnnotationHub) }
if (!require(enrichplot)) { BiocManager::install("enrichplot"); require(enrichplot) }
if (!require(GOSemSim)) { BiocManager::install("GOSemSim"); require(GOSemSim) }
if (!require(DESeq2)) { BiocManager::install("DESeq2"); require(DESeq2) }
if (!require(vsn)) { BiocManager::install("vsn"); require(vsn) }

# 3. User's package. The GO annotation database built using AnnotationForge package.
if (!require(org.Csp2T21.eg.db)) { create_go_db <- TRUE }

#
## You should reset this work_dir variable to the path where your inputs file located.
#
work_dir <- "~/Documents/projects/wp_wyf_rnaseq/outputs/diff_exp_gene"
setwd(work_dir)

#
## Prepare annotation database
#
# Wether create the GO db?
create_go_db <- FALSE
if (create_go_db) {
  # Prepare data
  gff_file <- read.gff("~/Documents/projects/wp_wyf_rnaseq/inputs/ncbi/GCA_009194965.1_Conioc1_genomic.gff.gz")
  go_ann <- gff_file[apply(gff_file, 1, function(e) { ifelse(stri_detect(e["attributes"], fixed = 'GO'), TRUE, FALSE) }), ]
  go_ann <- apply(go_ann, 1, function(e) {
    parent_go_pro <- stri_match(e["attributes"], regex = c("Parent=rna-gnl\\|WGS:VSMA\\|(.*?);", "Ontology_term=(.*?;)", "product=(.*?);"))[1:3, 2]
    gene_id <- stri_replace(parent_go_pro[1], "", fixed="mRNA")
    go_terms <- stri_match(parent_go_pro[2], mode="all", regex="(GO:[0-9].*?)[,;%]")
    go_product <- parent_go_pro[3]
    if (length(go_terms) >= 1) go_terms <- paste0(go_terms[[1]][, 2], collapse=";")
    list(SYMBOL=gene_id, GENENAME=go_product, CHROMOSOME=e["seqid"], GO=go_terms, EVIDENCE="NCBI")})
  
  go_ann <- do.call(rbind.data.frame, unique(go_ann))  # Could be problem prone
  go_ann <- go_ann[apply(go_ann, 1, function(e) !any(e == "NA")), ]
  go_ann["GID"] <- as.character(100000:(99999+nrow(go_ann)))  # The GID has to be numeric
  
  # Generate database
  org_db_dir <- "~/Documents/projects/wp_wyf_rnaseq/inputs/godb"
  fSym <- go_ann[c("GID", "SYMBOL", "GENENAME")]
  if (nrow(fSym[duplicated(fSym), ]) != 0) { stop("Duplicated observations are not allowed!") }
  rownames(fSym) <- 1:nrow(fSym)

  fChr <- go_ann[c("GID", "CHROMOSOME")]
  rownames(fChr) <- 1:nrow(fChr)

  fGO <- base::do.call(rbind.data.frame, apply(go_ann[c("GID", "GO", "EVIDENCE")], 1, function(e) {
    go_terms <- sapply(stri_split(e["GO"], fixed = ";")[[1]], function(e) {
      paste0("GO:", stri_pad_left(sub("GO:", "", e), width=7, pad="0"))}) # GO id should be padded to width of 7 by zeros. E.g., GO:5886 -> GO:0005886
    n_go_terms <- length(go_terms)

    list(GID=as.character(rep(e["GID"], n_go_terms)), GO=go_terms, EVIDENCE=rep(e["EVIDENCE"], n_go_terms))
  }))
  rownames(fGO) <- 1:nrow(fGO)

  unlink(paste0(org_db_dir, "/org.Csp2T21.eg.db"), recursive=TRUE)
  makeOrgPackage(gene_info=fSym, chromosome=fChr, go=fGO,
                 version="0.1.0",
                 maintainer="Zhenhua Zhang <zhenhua.zhang@rug.nl>",
                 author="Zhenhua Zhang <zhenhua.zhang@rug.nl",
                 outputDir=org_db_dir,
                 tax_id="1571157",
                 genus="Coniochaeta",
                 species="sp2T21",
                 goTable="go")

  install.packages(paste0(org_db_dir, "/org.Csp2T21.eg.db"), repos=NULL)
}


#
## Read and prepare input data. Read count matrix, design matrix
#
# Load read counts from featureCounts
fcrc_file <- "fungi_readcounts.tsv"  # Ought to be specified accordingly
fcrc_dtfm <- read.csv(fcrc_file, sep = "\t", check.names = FALSE)

# Transform the input data.frame into matrix
non_count_cols <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")
count_cols <- colnames(fcrc_dtfm)[!colnames(fcrc_dtfm) %in% non_count_cols]
fcrc_mtrx <- as.matrix(fcrc_dtfm[, count_cols])

# Row names and colum names
rownames(fcrc_mtrx) <- fcrc_dtfm[, "Geneid"]  # Using Geneid columns as row.names
colnames(fcrc_mtrx) <- count_cols             # Ensure the col.names are given

# Load experiment design matrix
exds_file <- "experiment_design_matrix.tsv" # Ought to be specified accordingly
exds_dtfm <- read.csv(exds_file, sep="\t", row.names = 1)

# Convert the treatment into factors
my_sep <- "."
exds_dtfm$speed <- factor(exds_dtfm$speed, levels=c("High", "Low"))
exds_dtfm$timepoint <- factor(exds_dtfm$timepoint, levels=c("D1", "D5", "D13"))
exds_dtfm$condition <- factor(paste(exds_dtfm$speed, exds_dtfm$timepoint, sep=my_sep))

#
## Calculate TPM
#
# A = total-reads-mapped-to-transcript / transcript-length
# TPM = A / sum(A) * 10^6
calc_tpm <- TRUE
if (calc_tpm) {
  x <- fcrc_mtrx / fcrc_dtfm$Length
  tpm_dtfm <- as.data.frame(t(t(x) / colSums(x) * 1e6))
  tpm_dtfm[non_count_cols] <- fcrc_dtfm[non_count_cols]
  tpm_dtfm <- tpm_dtfm[c(non_count_cols, count_cols)]
  tpm_ofpath <- stri_replace(fcrc_file, "_tpm.csv", fixed = ".tsv")
  write.csv(tpm_dtfm, file=tpm_ofpath, quote=FALSE, row.names=FALSE)
}


#
## Checking read counts matrix, NOT designed for differential expressing gene (DGE) analysis
#
# Check the data quality after do DE analysis using design ~ condition
dds <- DESeqDataSetFromMatrix(countData=fcrc_mtrx, colData=exds_dtfm,
                              design=formula("~ speed + timepoint"))

# Remove gene of pretty low expression (e.g. 10 reads for all samples)
dds <- dds[rowSums(counts(dds)) >= 10, ]

# DESeq wrapper for all steps, including estimation of size factors and dispersions, negative
# bionomial Vald test
dds <- DESeq(dds, quiet=TRUE)

# Transformation by vst (variance stablizing transformation). The following visualization analysis
# will use data transformed by this method.
vsd <- vst(dds, blind=FALSE)

# Check if SD is positively correlated with mean
ggsave("qc_meanSdPlot.pdf", meanSdPlot(assay(vsd), rank=FALSE)$gg, width=10, height=10)

# Heatmap to check read count matrix (first 50, sorted by normalized row mean counts)
# The color of each cell doesn't represent the raw read counts, but transformed values by vst()
select_n <- 500
select <- order(rowMeans(counts(dds, normalized=TRUE)), decreasing=TRUE)[1:select_n]
annot <- as.data.frame(colData(dds)[, c("speed", "timepoint")])
pheatmap(assay(vsd)[select, ],
         color=colorRampPalette(c("white", "firebrick3"))(50),
         cluster_rows=FALSE,
         cluster_cols=FALSE,
         show_rownames=FALSE,
         annotation_col=annot,
         filename='qc_readCountMatrix.pdf',
         width=10,
         height=10)

# PCA to check the effects of speed and timepoint
pcadata <- plotPCA(vsd, intgroup=c("speed", "timepoint"), returnData=TRUE)
percent_var <- round(100 * attr(pcadata, "percentVar"))
ggsave('qc_samplePCA.pdf',
       ggplot(pcadata, aes_string("PC1", "PC2", color=exds_dtfm$speed, shape=exds_dtfm$timepoint)) +
         theme_bw() +
         geom_point(size=3) +
         xlab(paste0("PC1: " , percent_var[1], "% variance")) +
         ylab(paste0("PC2: " , percent_var[2], "% variance")),
       width=10,
       height=10,)

# Heatmap to check the sample-to-sample distances
samples_dists <- dist(t(assay(vsd)))  # dist() estimates distances row-wisely.
pheatmap(as.matrix(samples_dists),
         color=colorRampPalette(c("firebrick3", "white"))(50),
         clustering_distance_rows=samples_dists,
         clustering_distance_cols=samples_dists,
         filename="allCompPairs_sampleToSampleDistance.pdf",
         width=10,
         height=10)

# Gene clustering (FIXME: is this z-score of expression?)
# Check the gene expression divation from the mean expression across all samples (n=18)
top_n <- 1000
top_n_var_genes <- head(order(rowVars(assay(vsd)), decreasing=TRUE), top_n)
var_mat <- assay(vsd)[top_n_var_genes, ]
var_mat <- var_mat - rowMeans(var_mat)
anno_col <- as.data.frame(colData(vsd)[, c("speed", "timepoint")])
pheatmap(var_mat,
         annotation_col=anno_col,
         show_rownames=FALSE,
         filename=paste0("allCompPairs_geneClusteringForTop", top_n, "VariedGenes.pdf"),
         width=10,
         height=10)

#
## DEG analysis
#

# Check the data quality after do DE analysis using design ~ condition
dds <- DESeqDataSetFromMatrix(countData=fcrc_mtrx, colData=exds_dtfm,
                              design=formula("~ condition"))

# Remove gene of pretty low expression (e.g. 10 reads for all samples)
dds <- dds[rowSums(counts(dds)) >= 10, ]

# DESeq wrapper for all steps, including estimation of size factors and dispersions, negative
# bionomial Vald test
dds <- DESeq(dds, quiet=TRUE)

# All comparison pairs. NOTE: some of them are meaningless, but remove in the for-loop.
comp_pairs <- t(combn(rev(unique(as.character(exds_dtfm[, "condition"]))), 2, c))

# idx <- 2 # meant to debug and test
for(idx in 1:nrow(comp_pairs)) {
  comp <- comp_pairs[idx, 1:2]
  
  # The first part and the second part should be halfly equal.
  cp1 <- stri_split(comp[1], fixed=my_sep)[[1]]
  cp2 <- stri_split(comp[2], fixed=my_sep)[[1]]
  if (cp1[1] != cp2[1] && cp1[2] != cp2[2]) { next }
  
  # Decide the exact comparison pairs. E.g., comparing D1 and D13 under low speed
  base_comp <- cp1[2]
  exact_comp <- c(cp1[1], cp2[1])
  if (cp1[1] == cp2[1]){
    base_comp <- cp1[1]
    exact_comp <- c(cp1[2], cp2[2])}
  baseflnm <- paste(base_comp, paste(exact_comp, collapse="-vs-"), sep="_")
  basetitle <- paste(paste(exact_comp, collapse=" vs "), "at", base_comp)
  
  # The contrast to fecth results
  contrast <- c("condition", comp[1], comp[2])
  
  # Fetch DE analysis results using specific 
  sub_res <- results(dds, contrast=contrast, alpha=0.01, pAdjustMethod="bonferroni")
  sub_res <- sub_res[order(sub_res$padj), ]
  
  # Shrink the log2 fold change value for better visualization and ranking
  shrink <- "ashr"
  sub_res_lfcs <- lfcShrink(dds, contrast=contrast, type=shrink, quiet=TRUE)        # Shrink the log2FoldChange
  sub_res_lfcs <- sub_res_lfcs[order(sub_res_lfcs$padj), ]                          # Order the DE by p-values (padj)
  sub_res_lfcs[is.na(sub_res_lfcs$padj), c('pvalue', 'padj')] <- 1                  # Update NA p-values to 1
  sub_res_lfcs[sub_res_lfcs$padj == 0, c("pvalue", "padj")] <- .Machine$double.xmin # 0 p-values to double min
  
  # MA plot
  ylim_max <- max(abs(sub_res_lfcs$log2FoldChange))
  pdf(paste0(baseflnm, "_MAplot.pdf"))
  DESeq2::plotMA(sub_res_lfcs, ylim=c(-ylim_max, ylim_max),
                 xlab=paste0("Mean of normalized counts (", shrink, ")"),
                 main=basetitle)
  dev.off()
  
  # Draw the volcano plot.
  
  ## Get the data.frame for the volcano plot
  sub_res_lfcs_dtfm <- as.data.frame(sub_res_lfcs)
  sub_res_lfcs_dtfm["geneId"] <- gsub("gene-", "", rownames(sub_res_lfcs_dtfm))
  
  ## p-val and log2(FC) threshold
  max_pval <- 1e-6
  min_logfc <- 1
  
  ## Color per group
  sub_res_lfcs_dtfm["Groups"] <- "NS"
  
  select <- sub_res_lfcs_dtfm[, "padj"] < max_pval & abs(sub_res_lfcs_dtfm[, "log2FoldChange"]) > 1
  sub_res_lfcs_dtfm[select, "Groups"] <- "Log2(FC) and p-adj"
  
  select <- sub_res_lfcs_dtfm[, "padj"] < max_pval & abs(sub_res_lfcs_dtfm[, "log2FoldChange"]) <= 1
  sub_res_lfcs_dtfm[select, "Groups"] <- "P-adj"
  
  select <- sub_res_lfcs_dtfm[, "padj"] > max_pval & abs(sub_res_lfcs_dtfm[, "log2FoldChange"]) > 1
  sub_res_lfcs_dtfm[select, "Groups"] <- "Log2(FC)"
  
  ## Using factors to specific the order of colors
  sub_res_lfcs_dtfm[, "Groups"] <- factor(sub_res_lfcs_dtfm[, "Groups"], levels=c("Log2(FC) and p-adj", "P-adj", "Log2(FC)", "NS"))
  
  ## First top_n genes will be labele
  top_n <- 25
  
  ## Prepare a data.frame to draw labels
  sub_res_for_lab_dtfm <- head(sub_res_lfcs_dtfm, top_n)
  
  ## Draw and save the plot
  ggsave(paste0(baseflnm, '_volcanoPlot.pdf'),
         ggplot() +
           geom_point(aes_string(x="log2FoldChange", y="-log10(padj)", color="Groups"),
                      data=sub_res_lfcs_dtfm, alpha=0.5) +
           geom_text_repel(aes_string(x="log2FoldChange", y="-log10(padj)", label="geneId"),
                           data=sub_res_for_lab_dtfm, min.segment.length=0) +
           geom_vline(xintercept=c(-min_logfc, min_logfc), linetype="longdash", color="black") +
           geom_hline(yintercept=c(-log10(max_pval)), linetype="longdash", color="black") +
           xlab(bquote(~Log[2]~ "fold change")) +
           ylab(bquote(~-Log[10]~"(P-adj)")) +
           scale_color_manual(values=c("red", "blue", "darkgreen", "black")),
         width=13,
         height=12)
  
  # Functional analysis for DE genes
  cand_degenes <- sub_res_lfcs_dtfm[sub_res_lfcs_dtfm[, "Groups"] == "Log2(FC) and p-adj", ]
  go_level <- 2
  go_type <- "BP"
  for (go_type in c("CC", "BP", "MF")) {
    # Group GO terms
    gpgo <- groupGO(gene=cand_degenes[, "geneId"],
                   OrgDb=org.Csp2T21.eg.db,
                   ont=go_type,
                   level=go_level,
                   keyType="SYMBOL")
    write.csv(gpgo@result, paste0(baseflnm, "_DEGenes_GO", go_type, ".csv"), quote=FALSE, row.names=FALSE)
    
    # Enrich GO terms
    ergo <- enrichGO(gene=cand_degenes[, "geneId"],
                    universe=sub_res_lfcs_dtfm[, "geneId"],
                    OrgDb=org.Csp2T21.eg.db,
                    ont=go_type,
                    keyType="SYMBOL",
                    pAdjustMethod="BH",
                    pvalueCutoff=0.05,
                    qvalueCutoff=0.05)
    
    # Save the bubble plot and the network plot, if enriched to terms found
    if (nrow(ergo) > 0) {
      tryCatch({
        gg <- dotplot(ergo)
        ggsave(paste0(baseflnm, "_DEGenes_GO", go_type, "Enrich_bubble.pdf"), width=10, height=10)
      }, error=function(e) (print(e)))
      
      tryCatch({
        gg <- emapplot(ergo)
        ggsave(paste0(baseflnm, "_DEGenes_GO", go_type, "Enrich_network.pdf"), gg, width=10, height=10)
      }, error=function(e) {print(e)})
    }
  }

  # Write the results to the disk
  res_lfcs_colnames <- colnames(sub_res_lfcs_dtfm)
  sub_res_lfcs_dtfm <- sub_res_lfcs_dtfm[c("geneId", res_lfcs_colnames[!res_lfcs_colnames %in% c("geneId")])]
  sub_res_lfcs_dtfm <- sub_res_lfcs_dtfm[order(sub_res_lfcs_dtfm$geneId), ]
  write.csv(as.data.frame(sub_res_lfcs_dtfm), file=paste0(baseflnm, "_DEGenes_LfcShrink.csv"), quote=FALSE, row.names=FALSE)
}


# Top n genes (ordered by padj)
# top_n <- 50
# top_n_genes <- head(sub_res_lfcs, top_n)
# top_n_genes_betas <- coef(dds)[rownames(top_n_genes), -1]
# pheatmap(top_n_genes_betas,
#          cluster_cols=FALSE,
#          filename=paste0(baseflnm, '_heatmapTopGenes.pdf'),
#          width=10,
#          height=10)



# vim: set nowrap ai ft=r tw=100:
