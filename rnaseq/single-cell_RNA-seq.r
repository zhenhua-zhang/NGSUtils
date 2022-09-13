#!/usr/bin/env Rscript

# Here we used Seurat version 3.1.4 specifically.

if (!require(scales)) {
  install.packages("scales")
  require(scales)}

if (!require(Seurat)) {
  if (!require(devtools)) {
    install.packages("devtools")}
  
  devtools::install_version("Seurat", version="3.1.4", repos="https://cran.r-project.org")
  require(Seurat)}

if (!require(dplyr)) {
  install.packages("dplyr")
  require(dplyr)}

if (!require(stringi)) {
  install.packages("stringi")
  require(stringi)}

if (!require(ggplot2)) {
  install.packages("ggplot2")
  require(ggplot2)
}

if (!require(patchwork)) {
  install.packages("patchwork")
  require(patchwork) }

if (!require(ggrepel)) {
  install.packages("ggrepel")
  require(ggrepel)}

if (!requireNamespace("BiocManager")) {
  install.packages("BiocManager")
  requireNamespace("BiocManager")
}

# FIXME: to use SingleR and avoid the bug here (https://support.bioconductor.org/p/p132709/), downgrade dbplyr firstly.
if (!require(dbplyr)) {
  devtools::install_url("http://cran.r-project.org/src/contrib/Archive/dbplyr/dbplyr_1.3.0.tar.gz")
  require(dbplyr)
}

if (!require(SingleR)) {
  BiocManager::install("SingleR")
  require(SingleR)}

if (!require(scater)) {
  BiocManager::install("scater")
  require(scater)}


# Project dir as working dir
wkdir <- "~/Documents/projects/wp_sc_rnaseq"
setwd(wkdir)

use_forced <- TRUE
# Regular pipepline reads per cell. BD Rhapsody
reads_pc <- "inputs/workingdata/regular_pipeline/Combined_SchulteCov_DBEC_ReadsPerCell.csv"
tag_call <- "inputs/workingdata/regular_pipeline/SchulteCov_Sample_Tag_Calls.csv"
ana_tag <- "regular"

# Forched pipeline reads per cell. BD Rhapsody
if (use_forced) {
  reads_pc <- "inputs/workingdata/forced_count_pipeline/forced_Combined_SchulteCov_DBEC_ReadsPerCell.csv"
  tag_call <- "inputs/workingdata/forced_count_pipeline/forced_SchulteCov_Sample_Tag_Calls.csv"
  ana_tag <- "forced"
}

# Load readcounts and call results
covid_rpc <- read.csv(reads_pc, skip=7, row.names=1)
covid_meta <- read.csv(tag_call, skip=7)

# Remove calls of c("Multiplet", "Undetermined")
good_call <- covid_meta[!covid_meta[, "Sample_Tag"] %in% c("Multiplet", "Undetermined"), "Cell_Index"]
covid_rpc <- covid_rpc[as.character(good_call), ]

# Make a meta data frame
covid_meta <- covid_meta[covid_meta[, "Cell_Index"] %in% good_call, ]

# Gender
covid_meta[, "Sex"] <- "Male"
covid_meta[covid_meta[, "Sample_Name"] %in% c("S14", "S16", "S7"), "Sex"] <- "Female"

# Disease
covid_meta[, "Disease"] <- "Control"
covid_meta[covid_meta[, "Sample_Name"] %in% c("S2", "S7"), "Disease"] <- "Covid19"

# Age
covid_meta[, "Age"] <- 0
covid_meta[covid_meta[, "Sample_Name"] == "S2", "Age"] <- 88
covid_meta[covid_meta[, "Sample_Name"] == "S7", "Age"] <- 53
covid_meta[covid_meta[, "Sample_Name"] == "S14", "Age"] <- 46
covid_meta[covid_meta[, "Sample_Name"] == "S16", "Age"] <- 47

rownames(covid_meta) <- covid_meta[, "Cell_Index"]
covid_meta <- covid_meta[, c("Sample_Name", "Sex", "Disease", "Age")]

# Anti-body seq
abseq_cols <- colnames(covid_rpc)[stri_endswith(colnames(covid_rpc), fixed=".pAbO")]
covid_abseq_rpc <- t(covid_rpc[, abseq_cols, drop=FALSE])

# RNA-seq
rnaseq_cols <- colnames(covid_rpc)[!colnames(covid_rpc) %in% abseq_cols]
covid_rnaseq_rpc <- t(covid_rpc[, rnaseq_cols, drop=FALSE])

# Seurat starts
covid_obj <- CreateSeuratObject(counts=covid_rnaseq_rpc,
                                project="LncRNA-scRNAseq",
                                meta.data=covid_meta,
                                min.cells=3,
                                min.features=3)
covid_obj[["PRO"]] <- CreateAssayObject(counts=covid_abseq_rpc)


# QC Check list
# [ ] Number of unique genes detected in each cell
#    [ ] Low-quality cells of empty droplets
#    [ ] Cell doublets or multiplets may exhibit an aberrantly high gene counts
# [ ] Number of molecules detected within a cell
# [-] Percentage of reads tha mapp to the mitochodiral genome. No MT genes in the target panel

# VlnPlot(covid_obj, features=c("nFeature_RNA", "nCount_RNA"), pt.size=0.5, ncol=2) # Check unique feature per cell
# FeatureScatter(covid_obj, feature1="nCount_RNA", feature2="nFeature_RNA") # Check sensitivity of features per cell.

# Keep cells that have unique feature counts
# n_features <- 45
# covid_obj <- subset(covid_obj, subset=nFeature_RNA > n_features)

# Handle donor effects by integrating results per sample
covid_obj_list <- SplitObject(covid_obj, split.by="Sample_Name")
covid_obj_list <- lapply(X=covid_obj_list, FUN = function(e) {
  e <- NormalizeData(e)
  e <- FindVariableFeatures(e, selection.method = "vst")  # The default nfeatures = 2000
})

features <- SelectIntegrationFeatures(object.list=covid_obj_list)
covid_acr <- FindIntegrationAnchors(object.list=covid_obj_list, anchor.features = features)
covid_obj <- IntegrateData(anchorset=covid_acr, verbose=FALSE)

# Scaling the data, default assay is integrated
covid_obj <- ScaleData(covid_obj, features=rownames(covid_obj))

# Dimentional reduction: linear version by PCA
covid_obj <- RunPCA(covid_obj, features=VariableFeatures(covid_obj))
# VizDimLoadings(covid_obj, dims=1:2, reduction="pca") # Dimention loadings
# DimPlot(covid_obj, reduction="pca") # PCA plot
# DimHeatmap(covid_obj, dims=1, cells=500) # Heatmap plot

# Determine dimensionality
# covid_obj <- JackStraw(covid_obj, num.replicate=100, prop.freq=0.1)
# covid_obj <- ScoreJackStraw(covid_obj, dims=1:20)
# JackStrawPlot(covid_obj, dims=1:20)
# ElbowPlot(covid_obj)

# As the previous plots, the first top_dim dimentions were used for the downstream analysis
top_dim <- 1:15
covid_obj <- FindNeighbors(covid_obj, dims=top_dim)  # Neighbors
covid_obj <- FindClusters(covid_obj, resolution=0.4) # Find clusters
covid_obj <- RunUMAP(covid_obj, dims=top_dim)        # UMAP

# UMAP plot, grouped by different indicators
for (group in c("Sex", "Disease", "Age", "Sample_Name", "seurat_clusters")) {
  ggsave(paste0("./outputs/first_run/", ana_tag, "-", DefaultAssay(covid_obj), "-", group, "-umap.pdf"), 
         DimPlot(covid_obj, reduction="umap", group.by=group, label=TRUE, label.size=6, pt.size=0.1),
         width=10,
         height=10)
}

# Proteins
covid_obj <- NormalizeData(covid_obj, normalization.method="CLR", margin=2, assay="PRO")

# RNA-seq
covid_obj <- NormalizeData(covid_obj, normalization.method="CLR", margin=2, assay="RNA")

# Determine cluster biomarkers (differentially expressed features)
check_assay <- DefaultAssay(covid_obj)
markers <- FindAllMarkers(covid_obj, min.pct=0.2, logfc.threshold = 0.2, min.diff.pct = 0.1, assay=check_assay)
top_n <- markers %>% group_by(cluster) %>% top_n(n=15, wt=avg_logFC)
ggsave(paste0("./outputs/first_run/", ana_tag, "-", check_assay, "-heatmap.pdf"),
       DoHeatmap(covid_obj, features=top_n$gene, combine=TRUE),
       width=20,
       height=12)

#
## SingleR to assign cell types using reference dataset
#
# ref_marker_db <- DatabaseImmuneCellExpressionData()
db <- "MID"
if (db == "MID") {
  ref_marker_db <- BlueprintEncodeData()
} else if (db == "BED") {
  ref_marker_db <- MonacoImmuneData()
} else if (db == "DICED") {
  ref_marker_db <- DatabaseImmuneCellExpressionData()
} else if (db == "HPCA") {
  ref_marker_db <- HumanPrimaryCellAtlasData()
} else {
  db <- "MID"
  ref_marker_db <- BlueprintEncodeData()
}

covid_seurat_rc <- GetAssayData(covid_obj, slot="scale.data", assay="RNA") # GetAssayData to ensure using the same cells
pred <- SingleR(test=covid_seurat_rc, ref=ref_marker_db, labels=ref_marker_db$label.main, method="cluster",
                clusters=covid_obj@meta.data$seurat_clusters, assay.type.test=1)
cluster_names <- paste(as.character(0:(length(pred$labels)-1)), pred$labels, sep="-")
names(cluster_names) <- as.character(0:(length(pred$labels)-1))
save_name <- paste0("./outputs/first_run/",
                    paste(ana_tag, DefaultAssay(covid_obj), db, "seurat_clusters-umap.pdf", sep="-"))
ggsave(save_name,
       DimPlot(covid_obj, reduction="umap", label=TRUE, pt.size=0.5) +
         scale_color_discrete(name="Cell cluster", labels=cluster_names),
       width=10,
       height=10)


#
## Differential expression analysis
#
# TODO: rewrite as a function

#' Defferential expression analysis
#'
#'@author Zhenhua Zhang
#'
#'@param ob Seurat object
#'@param gb Group by
#'@param cc Cell cluster
#'@param pv P-value
#'@param lf Log2 fold-change
#'@param tn Top n genes
#'@param sn Save name
#'
#'@example 
#'
deanalysis <- function(ob, gb, cc, pv, lf, tn, sn) {
  observation <- "Disease"
  cell_cluster <- "4"
  de_markers <- FindMarkers(covid_obj, ident.1="Control", ident.2="Covid19",
                            group.by=observation, subset.ident=cell_cluster)
  
  ## Transform avg_logFC to avg_log2FC
  if ((!"avg_log2FC" %in% colnames(de_markers)) && ("avg_logFC" %in% colnames(de_markers))) {
    de_markers["avg_log2FC"] <- de_markers["avg_logFC"] * log2(exp(1))
  }
  
  ## Add gene labels
  de_markers["geneId"] <- rownames(de_markers)
  
  ## p-val and log2(FC) threshold
  max_pval <- 5e-2
  min_logfc <- 1
  
  ## Color per group
  de_markers["Groups"] <- "NS"
  
  select <- de_markers[, "p_val_adj"] < max_pval & abs(de_markers[, "avg_log2FC"]) > 1
  de_markers[select, "Groups"] <- "Log2(FC) and p-adj"
  
  select <- de_markers[, "p_val_adj"] < max_pval & abs(de_markers[, "avg_log2FC"]) <= 1
  de_markers[select, "Groups"] <- "P-adj"
  
  select <- de_markers[, "p_val_adj"] > max_pval & abs(de_markers[, "avg_log2FC"]) > 1
  de_markers[select, "Groups"] <- "Log2(FC)"
  
  ## Using factors to specific the order of colors
  de_markers[, "Groups"] <- factor(de_markers[, "Groups"], levels=c("Log2(FC) and p-adj", "P-adj", "Log2(FC)", "NS"))
  
  ## First top_n genes will be labele
  top_n <- 10
  
  ## Prepare a data.frame to draw labels
  select <- de_markers[, "Groups"] == "Log2(FC) and p-adj"
  sub_de_markers <- de_markers[head(which(select), top_n), ]
  
  save_name <- paste0("./outputs/first_run/",
                      paste(ana_tag, DefaultAssay(covid_obj), observation, cell_cluster,
                            "volcanoPlot.pdf", sep="-"))
  ggsave(save_name, ggplot() +
           geom_point(aes_string(x="avg_log2FC", y="-log10(p_val_adj)", color="Groups"),
                      data=de_markers, alpha=0.5, size=3) +
           geom_text_repel(aes_string(x="avg_log2FC", y="-log10(p_val_adj)", label="geneId"),
                           data=sub_de_markers, min.segment.length=0) +
           geom_vline(xintercept=c(-min_logfc, min_logfc), linetype="dashed", color="black") +
           geom_hline(yintercept=c(-log10(max_pval)), linetype="dashed", color="black") +
           xlab(bquote(~Log[2]~ "fold change")) +
           ylab(bquote(~-Log[10]~"(P-adj)")) +
           scale_color_manual(values=c("red", "blue", "darkgreen", "black")),
         width=8,
         height=7)
}

DefaultAssay(covid_obj) <- "RNA"
# Gene expression stacked violin plot.
features <- rownames(de_markers[de_markers$p_val_adj < 0.01, ])
cowplot::plot_grid(plotlist=base::lapply(VlnPlot(covid_obj, features=features, ncol=1, same.y.lims=TRUE,
                                                 pt.size=0, combine=FALSE),
                                         function(e, last) {
                                           gene_id <- e$labels$title
                                           e + labs(x=NULL, y=gene_id, title=NULL) +
                                             theme(legend.position="none",
                                                   plot.margin=unit(c(5, 5, 0, 5), units="pt"),
                                                   axis.ticks.x=element_blank(),
                                                   axis.text.x=element_blank())},
                                         last=features[length(features)]),
                   ncol=1, labels=NULL)


#
## LINC00211-controled genes
#
koup_oedw <- c("STC2", "MMP25", "PRAC2", "S100A9", "DBNDD1", "C10orf35", "TMEM176B", "LYZ", "S100A8", "TRIM8")
kodw_oeup <- c("PCYOX1L", "HSPA1A", "IRF5", "HIST1H2BH", "ITGAX", "C1orf186", "SLC43A3", "CHI3L1", "TNFSF14", "HIST1H3H")


# Save the object to avoid redoing the computationally intensive steps
saveRDS(covid_obj, "./outputs/first_run/lncRNA_scRNAseq-regular_pipeline.rds")

# Last but not least, show infomation of current session.
sessionInfo()
