########################################
#                                      #
#    Step 1:install packages          #
#                                      #
########################################

install.packages("devtools")
source('https://bioconductor.org/biocLite.R')

install.packages("pheatmap")
biocLite("limma")
biocLite("scater") # scater用于自动化异常值检测
biocLite("RUVSeq",suppressUpdates=T)
biocLite("pcaMethods",suppressUpdates=T)
biocLite("SC3",suppressUpdates=T)
biocLite("M3Drop",suppressUpdates=T)
biocLite("TSCAN",suppressUpdates=T)
biocLite("monocle",suppressUpdates=T)
biocLite("destiny",suppressUpdates=T)

biocLite("edgeR")
biocLite("DESeq2")
biocLite("MAST")
biocLite("SummarizedExperiment")
biocLite("MultiAssayExperiment")



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
#    Step 1:install packages          #
#                                      #
########################################


