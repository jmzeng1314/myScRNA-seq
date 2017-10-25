########################################
#                                      #
#    Step 1:install packages          #
#                                      #
########################################

install.packages("devtools")
source('https://bioconductor.org/biocLite.R')

install.packages("pheatmap")
biocLite("limma")
biocLite("scater",suppressUpdates=T) # scater用于自动化异常值检测
biocLite("RUVSeq",suppressUpdates=T)
biocLite("pcaMethods",suppressUpdates=T)
biocLite("SC3",suppressUpdates=T)
biocLite("M3Drop",suppressUpdates=T)
biocLite("TSCAN",suppressUpdates=T)
biocLite("monocle",suppressUpdates=T)
biocLite("destiny",suppressUpdates=T)

biocLite("edgeR",suppressUpdates=T)
biocLite("DESeq2",suppressUpdates=T)
biocLite("MAST",suppressUpdates=T)
biocLite("SummarizedExperiment",suppressUpdates=T)
biocLite("MultiAssayExperiment",suppressUpdates=T)
biocLite("scran",suppressUpdates=T)



install.packages("mvoutlier") # mvoutlier的依赖包
install.packages("statmod")
install.packages("ROCR")
install.packages('mnormt') 

library(devtools)
install_github("hemberg-lab/scRNA.seq.funcs")
devtools::install_github('satijalab/seurat')
devtools::install_github("JustinaZ/pcaReduce")
devtools::install_github('jw156605/SLICER')
devtools::install_github("hms-dbmi/scde", build_vignettes = FALSE)

library(limma)
library(pheatmap)
library(devtools)

library(edgeR)
library(DESeq2)
library(MAST)

library(MultiAssayExperiment)
library(SummarizedExperiment)

library(mvoutlier)
library(statmod)
library(ROCR)

library("scater") # scater用于自动化异常值检测
library("RUVSeq")
library("pcaMethods")
library("SC3")
library("M3Drop")
library("TSCAN")
library("monocle")
library("destiny")

library(mnormt)
library(seurat)
library(pcaReduce)
library(SLICER)

if (!require(animation)) install.packages("animation")
library(animation)
ani.options(interval = 0.1, nmax = 100)
par(mar = c(4, 4, 1, 0.5))
clt.ani()

 
########################################
#                                      #
#    Step 2: QC by scater package      #
#                                      #
########################################

library(scater, quietly = TRUE)
library(knitr)
options(stringsAsFactors = FALSE)

# 这个文件是表达矩阵，包括线粒体基因和 ERCC spike-ins 的表达量，可以用来做质控
molecules <- read.table("tung/molecules.txt", sep = "\t")

## 这个文件是表达矩阵涉及到的所有样本的描述信息，包括样本来源于哪个细胞，以及哪个批次。
anno <- read.table("tung/annotation.txt", sep = "\t", header = TRUE)

knitr::kable(
  head(molecules[ , 1:3]), booktabs = TRUE,
  caption = 'A table of the first 6 rows and 3 columns of the molecules table.'
)

knitr::kable(
  head(anno), booktabs = TRUE,
  caption = 'A table of the first 6 rows of the anno table.'
)


#First, create the scater SCESet classes:

pheno_data <- new("AnnotatedDataFrame", anno)
rownames(pheno_data) <- pheno_data$sample_id
dat <- scater::newSCESet(
  countData = molecules,
  phenoData = pheno_data
)

#Remove genes that are not expressed in any cell:

keep_feature <- rowSums(counts(dat) > 0) > 0
dat <- dat[keep_feature, ]

## 这里得到dat是一个SCESet对象

#Define control features (genes) - ERCC spike-ins and mitochondrial genes (provided by the authors):
## 线粒体基因和 ERCC spike-ins 都需要拿出来作为质控的指引。
## 这里可以找到89个ERCC spike-ins 序列
ercc <- featureNames(dat)[grepl("ERCC-", featureNames(dat))] 
## 线粒体基因序列，可以从GTF等注释文件提取
mt <- c("ENSG00000198899", "ENSG00000198727", "ENSG00000198888",
        "ENSG00000198886", "ENSG00000212907", "ENSG00000198786",
        "ENSG00000198695", "ENSG00000198712", "ENSG00000198804",
        "ENSG00000198763", "ENSG00000228253", "ENSG00000198938",
        "ENSG00000198840")


#Calculate the quality metrics:
## 这里非常重要，后续的分析都是基于次
dat_qc <- scater::calculateQCMetrics(
  dat,
  feature_controls = list(ERCC = ercc, MT = mt)
)


hist(
  dat_qc$total_counts,
  breaks = 100
)
abline(v = 25000, col = "red")
hist(
  dat_qc$total_features,
  breaks = 100
)
abline(v = 7000, col = "red")


scater::plotPhenoData(
  dat_qc,
  aes_string(x = "total_features",
             y = "pct_counts_feature_controls_MT",
             colour = "batch")
)

scater::plotPhenoData(
  dat_qc,
  aes_string(x = "total_features",
             y = "pct_counts_feature_controls_ERCC",
             colour = "batch")
)
## 可以看到整个 NA19098.r2都是有问题的，


filter_by_expr_features <- dat_qc$total_counts >25000 
table(filter_by_expr_features)  
filter_by_total_counts  <- dat_qc$total_features >7000
table(filter_by_total_counts) 
filter_by_ERCC   <- dat_qc$pct_counts_feature_controls_ERCC < 5 
table(filter_by_ERCC)
filter_by_MT<- dat_qc$pct_counts_feature_controls_MT <10
table(filter_by_MT)  


dat_qc$use <- (
  # sufficient features (genes)
  filter_by_expr_features &
    # sufficient molecules counted
    filter_by_total_counts &
    # sufficient endogenous RNA
    filter_by_ERCC &
    # remove cells with unusual number of reads in MT genes
    filter_by_MT
)
knitr::kable(
  as.data.frame(table(dat_qc$use)),
  booktabs = TRUE,
  row.names = FALSE,
  caption = 'The number of cells removed by manual filter (FALSE)'
)

dat_pca <- scater::plotPCA(dat_qc,
                           size_by = "total_features", 
                           shape_by = "use",
                           pca_data_input = "pdata",
                           detect_outliers = TRUE,
                           return_SCESet = TRUE)

knitr::kable(
  as.data.frame(table(dat_pca$outlier)),
  booktabs = TRUE,
  row.names = FALSE,
  caption = 'The number of cells removed by automatic filter (FALSE)'
)


scater::plotQC(dat_qc, type = "highest-expression")


filter_genes <- apply(counts(dat_qc[ , pData(dat_qc)$use]), 1, 
                      function(x) length(x[x > 1]) >= 2)
fData(dat_qc)$use <- filter_genes
knitr::kable(
  as.data.frame(table(filter_genes)),
  booktabs = TRUE,
  row.names = FALSE,
  caption = 'The number of genes removed by gene filter (FALSE)'
)

dim(dat_qc[fData(dat_qc)$use, pData(dat_qc)$use])
########################################
#                                      #
#    Step 3:  cluster by PCA or tSNE   #
#                                      #
########################################
set_exprs(dat_qc, "log2_counts") <- log2(counts(dat_qc) + 1)
dat_filter <-dat_qc[fData(dat_qc)$use, pData(dat_qc)$use]
endog_genes <- !fData(dat_filter)$is_feature_control

scater::plotPCA(dat_qc[endog_genes, ],
                ntop = 500,
                colour_by = "batch",
                size_by = "total_features",
                shape_by = "individual",
                exprs_values = "counts")
scater::plotPCA(dat_qc[endog_genes, ],
                ntop = 500,
                colour_by = "batch",
                size_by = "total_features",
                shape_by = "individual",
                exprs_values = "log2_counts")

scater::plotPCA(dat_filter[endog_genes, ],
                ntop = 500,
                colour_by = "batch",
                size_by = "total_features",
                shape_by = "individual",
                exprs_values = "log2_counts")

scater::plotTSNE(dat_qc[endog_genes, ],
                 ntop = 500,
                 perplexity = 130,
                 colour_by = "batch",
                 size_by = "total_features",
                 shape_by = "individual",
                 exprs_values = "log2_counts",
                 rand_seed = 123456)
scater::plotTSNE(dat_filter[endog_genes, ],
                 ntop = 500,
                 perplexity = 130,
                 colour_by = "batch",
                 size_by = "total_features",
                 shape_by = "individual",
                 exprs_values = "log2_counts",
                 rand_seed = 123456)


scater::plotQC(dat_filter[endog_genes, ],
               type = "find-pcs",
               variable = "total_features",
               exprs_values = "log2_counts")
scater::plotQC(dat_filter[endog_genes, ],
               type = "expl",
               exprs_values = "log2_counts",
               variables = c("total_features",
                             "total_counts",
                             "batch",
                             "individual",
                             "pct_counts_feature_controls_ERCC",
                             "pct_counts_feature_controls_MT"))
## normalization 
########################################
#                                      #
#    Step 4:  normalization            #
#                                      #
########################################

plotPCA(
  dat_filter[endog_genes, ],
  colour_by = "batch",
  size_by = "total_features",
  shape_by = "individual",
  exprs_values = "log2_counts"
)

## first for CPM 
plotPCA(
  dat_filter[endog_genes, ],
  colour_by = "batch",
  size_by = "total_features",
  shape_by = "individual",
  exprs_values = "exprs"
) 
## then for TMM
dat_filter <- normaliseExprs(
  dat_filter,
  method = "TMM",
  feature_set = endog_genes
)
plotPCA(
  dat_filter[endog_genes, ],
  colour_by = "batch",
  size_by = "total_features",
  shape_by = "individual",
  exprs_values = "norm_exprs"
)
## next for Size-factor (RLE)
dat_filter <- normaliseExprs(
  dat_filter,
  method = "RLE", 
  feature_set = endog_genes
)
plotPCA(
  dat_filter[endog_genes, ],
  colour_by = "batch",
  size_by = "total_features",
  shape_by = "individual",
  exprs_values = "norm_exprs"
)
## next for Upperquantile 
dat_filter <- normaliseExprs(
  dat_filter,
  method = "upperquartile", 
  feature_set = endog_genes,
  p = 0.99
)
plotPCA(
  dat_filter[endog_genes, ],
  colour_by = "batch",
  size_by = "total_features",
  shape_by = "individual",
  exprs_values = "norm_exprs"
)

## cluster 
########################################
#                                      #
#    Step 4:  clustering               #
#                                      #
########################################

library(pcaMethods)
library(pcaReduce)
library(SC3)
library(scater)
library(pheatmap)
set.seed(1234567)
pollen <- readRDS("pollen/pollen.rds")
pollen
table(pData(pollen)$cell_type1)
plotPCA(pollen, colour_by = "cell_type1")


pollen <- sc3_prepare(pollen, ks = 2:5) 
pollen <- sc3_estimate_k(pollen)
pollen@sc3$k_estimation
## SC3方法聚类的结果恰好是实验选取的细胞种类。
pollen <- sc3(pollen, ks = 11, biology = TRUE)
sc3_plot_consensus(pollen, k = 11, show_pdata = "cell_type1")
sc3_plot_silhouette(pollen, k = 11)
sc3_plot_expression(pollen, k = 11, show_pdata = "cell_type1")
sc3_plot_markers(pollen, k = 11, show_pdata = "cell_type1")
plotPCA(pollen, colour_by = "sc3_11_clusters")
#sc3_interactive(pollen)

