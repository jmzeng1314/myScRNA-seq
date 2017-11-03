# 写在前面

> 虽然会了那么多NGS组学分析，也一一写过全套教程了。但总躺在舒适区就不好了，还是得学点新东西哈。

就从这个scRNA-seq开始吧，这里面有一个 **old** 文件夹，是因为这个**scater**包更新的步子迈得太大了，以至于我2017/11/1 之前的学习笔记全部作废！


# scRNA-seq课程介绍

> 我会完整的学习完这个课程，并且记录自己的学习笔记，课程是：Analysis of single cell RNA-seq data course, Cambridge University, UK

只需要看下面两个文档即可：

* [gitbub](https://github.com/hemberg-lab/scRNA.seq.course)
* [文档](http://hemberg-lab.github.io/scRNA.seq.course)

# 环境配置

> 主要是需要安装R包啦。

如果是在自己的电脑里面就直接打开R，然后一个个包的安装吧，我把代码简单整理了一下,直接点击[installed_required_packages.R](installed_required_packages.R)

其实最方便的是用docker技术。
首先在自己的亚马逊云上面把他们的课程docker镜像下载下来，并且使用该镜像来创建一个容器。
命令如下：
```
docker run -it quay.io/hemberg-group/scrna-seq-course:latest R
```
等下载完成后，我们可以直接使用这个镜像来启动运行容器！

因为亚马逊是国外的，所以下载非常快，如果是腾讯云阿里云可能需要好几个小时。毕竟是3个多G呀!

这个docker镜像里面已经安装好了所有的软件环境和R包.

如果想了解[更多docker命令，可以点击学习](http://www.runoob.com/docker/docker-image-usage.html)

## 主要基于scater包

最新的文档如下：

| [HTML](http://bioconductor.org/packages/release/bioc/vignettes/scater/inst/doc/vignette-intro.html) | [R Script](http://bioconductor.org/packages/release/bioc/vignettes/scater/inst/doc/vignette-intro.R) | An introduction to the scater package    |
| ---------------------------------------- | ---------------------------------------- | ---------------------------------------- |
| [HTML](http://bioconductor.org/packages/release/bioc/vignettes/scater/inst/doc/vignette-dataviz.html) | [R Script](http://bioconductor.org/packages/release/bioc/vignettes/scater/inst/doc/vignette-dataviz.R) | Data visualisation methods in scater     |
| [HTML](http://bioconductor.org/packages/release/bioc/vignettes/scater/inst/doc/vignette-quantimport.html) | [R Script](http://bioconductor.org/packages/release/bioc/vignettes/scater/inst/doc/vignette-quantimport.R) | Expression quantification and import     |
| [HTML](http://bioconductor.org/packages/release/bioc/vignettes/scater/inst/doc/vignette-qc.html) | [R Script](http://bioconductor.org/packages/release/bioc/vignettes/scater/inst/doc/vignette-qc.R) | Quality control with scater              |
| [HTML](http://bioconductor.org/packages/release/bioc/vignettes/scater/inst/doc/vignette-transition.html) | [R Script](http://bioconductor.org/packages/release/bioc/vignettes/scater/inst/doc/vignette-transition.R) | Transition from SCESet to SingleCellExperiment |
| [PDF](http://bioconductor.org/packages/release/bioc/manuals/scater/man/scater.pdf) |                                          |                                          |

而且其GitHub的教程也更新了：http://hemberg-lab.github.io/scRNA.seq.course/

## 必须了解这个 `SingleCellExperiment` 对象

创建该对象代码如下：

```R
suppressPackageStartupMessages(library(scater))
data("sc_example_counts")
data("sc_example_cell_info")

## ----quickstart-make-sce, results='hide'-----------------------------------
gene_df <- DataFrame(Gene = rownames(sc_example_counts))
rownames(gene_df) <- gene_df$Gene
example_sce <- SingleCellExperiment(assays = list(counts = sc_example_counts), 
                                    colData = sc_example_cell_info, 
                                    rowData = gene_df)

example_sce <- normalise(example_sce)

## ----quickstart-add-exprs, results='hide'----------------------------------
exprs(example_sce) <- log2(
    calculateCPM(example_sce, use.size.factors = FALSE) + 1)

## ----filter-no-exprs-------------------------------------------------------
keep_feature <- rowSums(exprs(example_sce) > 0) > 0
example_sce <- example_sce[keep_feature,]

example_sceset <- calculateQCMetrics(example_sce, feature_controls = list(eg = 1:40)) 
 

colnames(colData(example_sceset))
colnames(rowData(example_sceset))
```

首先是基于样本的过滤，用 `colData(object)` 可以查看各个样本统计情况 

- `total_counts`: total number of counts for the cell (aka ‘library size’)

- `log10_total_counts`: total_counts on the log10-scale

- `total_features`: the number of features for the cell that have expression above the detection limit (default detection limit is zero)

- `filter_on_total_counts`: would this cell be filtered out based on its log10-total_counts being (by default) more than 5 median absolute deviations from the median log10-total_counts for the dataset?

- `filter_on_total_features`: would this cell be filtered out based on its total_features being (by default) more than 5 median absolute deviations from the median total_features for the dataset?

- `counts_feature_controls`: total number of counts for the cell that come from (a set of user-defined) control features. Defaults to zero if no control features are indicated.

- `counts_endogenous_features`: total number of counts for the cell that come from endogenous features (i.e. not control features). Defaults to `total_counts` if no control features are indicated.

- `log10_counts_feature_controls`: total number of counts from control features on the log10-scale. Defaults to zero (i.e. log10(0 + 1), offset to avoid infinite values) if no control features are indicated.

- `log10_counts_endogenous_features`: total number of counts from endogenous features on the log10-scale. Defaults to zero (i.e. log10(0 + 1), offset to avoid infinite values) if no control features are indicated.

- `n_detected_feature_controls`: number of defined feature controls that have expression greater than the threshold defined in the object. *`pct_counts_feature_controls`: percentage of all counts that come from the defined control features. Defaults to zero if no control features are defined.


然后是基于基因的过滤，用 `rowData(object)` 可以查看各个基因统计情况

- `mean_exprs`: the mean expression level of the gene/feature.
- `exprs_rank`: the rank of the feature’s expression level in the cell.
- `total_feature_counts`: the total number of counts mapped to that feature across all cells.
- `log10_total_feature_counts`: total feature counts on the log10-scale.
- `pct_total_counts`: the percentage of all counts that are accounted for by the counts mapping to the feature.
- `is_feature_control`: is the feature a control feature? Default is `FALSE` unless control features are defined by the user.
- `n_cells_exprs`: the number of cells for which the expression level of the feature is above the detection limit (default detection limit is zero).

所有的可视化都集中在了 `scater_gui` 这个函数产生的`shiny`网页里面：

- `plotScater`: a plot method exists for `SingleCellExperiment` objects, which gives an overview of expression across cells.
- `plotQC`: various methods are available for producing QC diagnostic plots.
- `plotPCA`: produce a principal components plot for the cells.
- `plotTSNE`: produce a t-distributed stochastic neighbour embedding (reduced dimension) plot for the cells.
- `plotDiffusionMap`: produce a diffusion map (reduced dimension) plot for the cells.
- `plotMDS`: produce a multi-dimensional scaling plot for the cells.
- `plotReducedDim`: plot a reduced-dimension representation of the cells.
- `plotExpression`: plot expression levels for a defined set of features.
- `plotPlatePosition`: plot cells in their position on a plate, coloured by cell metadata and QC metrics or feature expression level.
- `plotColData`: plot cell metadata and QC metrics.
- `plotRowData`: plot feature metadata and QC metrics.

可以充分的探索自己的数据，或者上面的每一个函数都可以对进行单独可视化，细致的探索自己的单细胞测序数据。

# 数据集的介绍

## Tung dataset

这个是芝加哥大学[ Yoav Gilad’s lab ](http://giladlab.uchicago.edu/)实验的[Tung et al. 2017)](http://hemberg-lab.github.io/scRNA.seq.course/exprs-qc.html#ref-Tung2017-ba)文章里面的数据。测了来源于hapmap计划的3个人的单细胞Single-cell RNA sequencing (scRNA-seq)。用的是single-cell Fluidigm C1 platform来做单细胞分选，材料是three human induced pluripotent stem cell (iPSC) lines of three Yoruba individuals (abbreviation: YRI) 

### 测序数据

数据可以直接去他们的[GitHub里面下载](http://jdblischak.github.io/singleCellSeq/analysis/),但是不建议从fastq开始，比较太大了一点,如果一定要下载，可以分成4个部分下载，log日志如下：
```
Downloaded: 1375 files, 129G in 4h 59m 1s (7.37 MB/s)
Downloaded: 1377 files, 101G in 5h 0m 56s (5.71 MB/s)
Downloaded: 1376 files, 112G in 6h 24m 35s (4.96 MB/s)
Downloaded: 1184 files, 96G in 4h 54m 6s (5.54 MB/s)
```
GEO数据库[(GSE77288)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE77288)里面也有作者上传的整理好的表达矩阵：
```
GSE77288_molecules-raw-single-per-lane.txt.gz	20.0 Mb	(ftp)(http)	TXT
GSE77288_molecules-raw-single-per-sample.txt.gz	7.9 Mb	(ftp)(http)	TXT
GSE77288_reads-raw-bulk-per-lane.txt.gz	1.6 Mb	(ftp)(http)	TXT
GSE77288_reads-raw-bulk-per-sample.txt.gz	323.1 Kb	(ftp)(http)	TXT
GSE77288_reads-raw-single-per-lane.txt.gz	30.6 Mb	(ftp)(http)	TXT
GSE77288_reads-raw-single-per-sample.txt.gz	12.4 Mb	(ftp)(http)	TXT
```

这个数据集里面既添加了[unique molecular identifiers (UMIs)](http://www.nature.com/nmeth/journal/v9/n1/full/nmeth.1778.html) 又添加了 ERCC spike-ins 来作为表达量的质量控制。

### 数据质控

C1平台是96孔板，所以96个样本是一个batch，细胞的平均测序量是6.3M，在0.4~11.2M之间。首先做 visual inspection 挑选真正的单孔单细胞，而且无污染。再各种PCA，tSNE聚类看它们的表现。

每一种细胞系测了3个96孔板，所以共有9个batch，近900个样本，经过一系列过滤之后的高质量的样本只有564个。

#### 再过滤一些细胞

- Only one cell observed per well.
- At least 1,556,255 mapped reads.
- Less than 36.4% unmapped reads.
- Less than 3.2% ERCC reads.
- More than 6,788 genes with at least one read.

After filtering, we maintained 564 high quality single cells (NA19098: 142, NA19101: 201, NA19239: 221).

#### 再过滤一些基因

- The quality control analyses were performed using all protein-coding genes (Ensembl GRCh37 release 82) with at least one observed read. 
- Using the high quality single cells, we further removed genes with low expression levels for downstream analyses. We removed all genes with a mean log2 cpm less than 2
- We also removed genes with molecule counts larger than 1,024 for the correction of collision probability. 

In the end we kept 13,058 endogenous genes and 48 ERCC spike-in genes.

接下来我们就重现一下这个分析。



#### step1 : 读取表达矩阵

todo

#### step2: 样本和基因分别过滤

todo

#### step3：表达矩阵的可视化

todo

