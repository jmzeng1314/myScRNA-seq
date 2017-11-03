# 写在前面

> 虽然会了那么多NGS组学分析，也一一写过全套教程了。但总躺在舒适区就不好了，还是得学点新东西哈。

就从这个scRNA-seq开始吧。


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

因为亚马逊是国外的，所以下载非常快，如果是腾讯云阿里云可能需要好几个小时。毕竟是6个多G呀!

这个docker镜像里面已经安装好了所有的软件环境和R包.

如果想了解[更多docker命令，可以点击学习](http://www.runoob.com/docker/docker-image-usage.html)

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

# 构建表达矩阵

> 可以选择STAR或者HISAT2，fastq格式的测序数据文件处理成bam格式的比对结果。

- To assess read quality, we ran FastQC  and observed a decrease in base quality at the 3′ end of the reads. 
- Thus we removed low quality bases from the 3′ end using sickle with default settings.
- To handle the UMI sequences at the 5′ end of each read, we used umitools to find all reads with a UMI of the pattern NNNNNGGG (reads without UMIs were discarded). 
- We then mapped reads to human genome hg19 (only including chromosomes 1–22, X, and Y, plus the ERCC sequences) with Subjunc , discarding non-uniquely mapped reads (option -u). 
- To obtain gene-level counts, we assigned reads to protein-coding genes (Ensembl GRCh37 release 82) and the ERCC spike-in genes using featureCounts . 

这里就得到了所有细胞的所有检测到的基因的表达量。上述流程里面包含了去除低质量reads，挑选标记着UMI的序列，subjunc工具来比对，featureCounts来做定量。

### 流程的shell代码示例

```
$<path_to_STAR>/STAR --runThreadN 1 --runMode alignReads
--readFilesIn reads1.fq.gz reads2.fq.gz --readFilesCommand zcat --genomeDir <path>
--parametersFiles FileOfMoreParameters.txt --outFileNamePrefix <outpath>/output
$<path_to_Salmon>/salmon quant -i salmon_transcript_index -1 reads1.fq.gz -2 reads2.fq.gz -p #threads -l A -g genome.gtf --seqBias --gcBias --posBias
python <RSeQCpath>/geneBody_coverage.py -i input.bam -r genome.bed -o output.txt
python <RSeQCpath>/bam_stat.py -i input.bam -r genome.bed -o output.txt
python <RSeQCpath>/split_bam.py -i input.bam -r rRNAmask.bed -o output.txt
# include multimapping
<featureCounts_path>/featureCounts -O -M -Q 30 -p -a genome.gtf -o outputfile input.bam
# exclude multimapping
<featureCounts_path>/featureCounts -Q 30 -p -a genome.gtf -o outputfile input.bam

```

这个流程代码可以更加多元化，自动化，上面只是一个简单的例子而已。

### 比对情况的简单QC

![比对的百分比](https://hemberg-lab.github.io/scRNA.seq.course/figures/Bergiers_exp1_mapping_by_cell.png)

还可以用[RSeQC](http://rseqc.sourceforge.net/)来查看转录组的各种比对情况：

![基因主体的覆盖偏差](https://hemberg-lab.github.io/scRNA.seq.course/figures/Exp1_RSEQC_geneBodyCoverage_plot_Combined.png) 

# 用scater包做QC

## 首先了解测序数据得到的表达矩阵


```r
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
```



Table: A table of the first 6 rows and 3 columns of the molecules table.

                   NA19098.r1.A01   NA19098.r1.A02   NA19098.r1.A03
----------------  ---------------  ---------------  ---------------
ENSG00000237683                 0                0                0
ENSG00000187634                 0                0                0
ENSG00000188976                 3                6                1
ENSG00000187961                 0                0                0
ENSG00000187583                 0                0                0
ENSG00000187642                 0                0                0

```r
knitr::kable(
    head(anno), booktabs = TRUE,
    caption = 'A table of the first 6 rows of the anno table.'
)
```



Table: A table of the first 6 rows of the anno table.

individual   replicate   well   batch        sample_id      
-----------  ----------  -----  -----------  ---------------
NA19098      r1          A01    NA19098.r1   NA19098.r1.A01 
NA19098      r1          A02    NA19098.r1   NA19098.r1.A02 
NA19098      r1          A03    NA19098.r1   NA19098.r1.A03 
NA19098      r1          A04    NA19098.r1   NA19098.r1.A04 
NA19098      r1          A05    NA19098.r1   NA19098.r1.A05 
NA19098      r1          A06    NA19098.r1   NA19098.r1.A06 

## 表达矩阵导入scater包做质控


```r
#First, create the scater SCESet classes:

pheno_data <- new("AnnotatedDataFrame", anno)
rownames(pheno_data) <- pheno_data$sample_id
dat <- scater::newSCESet(
    countData = molecules,
    phenoData = pheno_data
)

# 这个对象非常重要， pData(dat) 和  fData(dat) 分别记录了 表达矩阵的样本信息和基因信息。

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

dat_qc <- scater::calculateQCMetrics(
    dat,
    feature_controls = list(ERCC = ercc, MT = mt)
)
## 这里dat_qc非常重要，后续的分析都是基于此，还是看pData(dat_qc) 和  fData(dat_qc) 
dat_qc
```

```
## SCESet (storageMode: lockedEnvironment)
## assayData: 18726 features, 864 samples 
##   element names: counts, exprs 
## protocolData: none
## phenoData
##   rowNames: NA19098.r1.A01 NA19098.r1.A02 ... NA19239.r3.H12 (864
##     total)
##   varLabels: individual replicate ... is_cell_control (52 total)
##   varMetadata: labelDescription
## featureData
##   featureNames: ENSG00000237683 ENSG00000187634 ... ERCC-00171
##     (18726 total)
##   fvarLabels: mean_exprs exprs_rank ... is_feature_control (12
##     total)
##   fvarMetadata: labelDescription
## experimentData: use 'experimentData(object)'
## Annotation:
```

## 基于样本的过滤

### 基于文库大小和测到的基因数量来质控


```r
## 首先看每个样本的总reads数量
hist(
    dat_qc$total_counts,
    breaks = 100
)
abline(v = 25000, col = "red")
```

![](readme_files/figure-html/unnamed-chunk-2-1.png)<!-- -->

```r
## 然后看每个样本测到的基因数量
hist(
    dat_qc$total_features,
    breaks = 100
)
abline(v = 7000, col = "red")
```

![](readme_files/figure-html/unnamed-chunk-2-2.png)<!-- -->
画红线的左边的样本可能不合格，上面的threshold可能过于武断。

### 线粒体基因和 ERCC spike-ins 序列占比问题

单细胞测序想探究的是细胞里面基因的表达量，如果外源的ERCC spike-ins 序列占比过多，这样的数据也得抛弃，可能是细胞生存状态不佳，或者实验过程中的RNA降解等原因。

下面的两个QC图也是scater包的特色


```r
scater::plotPhenoData(
    dat_qc,
    aes_string(x = "total_features",
               y = "pct_counts_feature_controls_MT",
               colour = "batch")
)
```

![](readme_files/figure-html/unnamed-chunk-3-1.png)<!-- -->

```r
## 可以看到少量样本的线粒体序列占比超过了10%,可以考虑去除掉这些样本。

scater::plotPhenoData(
    dat_qc,
    aes_string(x = "total_features",
               y = "pct_counts_feature_controls_ERCC",
               colour = "batch")
)
```

![](readme_files/figure-html/unnamed-chunk-3-2.png)<!-- -->

```r
## 可以看到整个 NA19098.r2都是有问题的，外源的ERCC spike-ins含量太高了。可以考虑去除整个batch
```

### 选定过滤标准

```r
filter_by_expr_features <- dat_qc$total_counts >25000 
table(filter_by_expr_features)  
```

```
## filter_by_expr_features
## FALSE  TRUE 
##    46   818
```

```r
filter_by_total_counts  <- dat_qc$total_features >7000
table(filter_by_total_counts) 
```

```
## filter_by_total_counts
## FALSE  TRUE 
##   120   744
```

```r
filter_by_ERCC   <- dat_qc$pct_counts_feature_controls_ERCC < 5
## 这里我没有根据作者那样选择去除整个NA19098.r2批次的样本
# filter_by_ERCC <- reads$batch != "NA19098.r2" & reads$pct_counts_feature_controls_ERCC < 25
table(filter_by_ERCC)
```

```
## filter_by_ERCC
## FALSE  TRUE 
##   102   762
```

```r
filter_by_MT<- dat_qc$pct_counts_feature_controls_MT < 10
table(filter_by_MT)  
```

```
## filter_by_MT
## FALSE  TRUE 
##    31   833
```

### 得到最后的高质量单细胞测序数据


```r
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
```



Table: The number of cells removed by manual filter (FALSE)

Var1     Freq
------  -----
FALSE     207
TRUE      657

这4个标准合起来过滤掉了207个样本，最后剩下657个高质量测序数据。

### scater一站式过滤低质量样本

> scater包自己提供了一个基于PCA的QC标准，不需要自己根据文库大小，覆盖的基因数量，外源的ERCC spike-ins 含量以及线粒体DNA含量来进行人工过滤。

默认的标准如下：

- pct_counts_top100features
- total_features
- pct_counts_feature_controls
- n_detected_feature_controls
- log10_counts_endogenous_features
- log10_counts_feature_controls

一站式QC函数如下：


```r
dat_pca <- scater::plotPCA(dat_qc,
                  size_by = "total_features", 
                  shape_by = "use",
                  pca_data_input = "pdata",
                  detect_outliers = TRUE,
                  return_SCESet = TRUE)
```

```
## sROC 0.1-2 loaded
```

```
## The following cells/samples are detected as outliers:
## NA19098.r2.A01
## NA19098.r2.A02
## NA19098.r2.A06
## NA19098.r2.A09
## NA19098.r2.A10
## NA19098.r2.A12
## NA19098.r2.B01
## NA19098.r2.B03
## NA19098.r2.B04
## NA19098.r2.B05
## NA19098.r2.B07
## NA19098.r2.B11
## NA19098.r2.B12
## NA19098.r2.C01
## NA19098.r2.C02
## NA19098.r2.C03
## NA19098.r2.C04
## NA19098.r2.C05
## NA19098.r2.C06
## NA19098.r2.C07
## NA19098.r2.C08
## NA19098.r2.C09
## NA19098.r2.C10
## NA19098.r2.C11
## NA19098.r2.C12
## NA19098.r2.D01
## NA19098.r2.D02
## NA19098.r2.D03
## NA19098.r2.D04
## NA19098.r2.D07
## NA19098.r2.D08
## NA19098.r2.D09
## NA19098.r2.D10
## NA19098.r2.D12
## NA19098.r2.E01
## NA19098.r2.E02
## NA19098.r2.E03
## NA19098.r2.E04
## NA19098.r2.E05
## NA19098.r2.E06
## NA19098.r2.E07
## NA19098.r2.E12
## NA19098.r2.F01
## NA19098.r2.F02
## NA19098.r2.F07
## NA19098.r2.F08
## NA19098.r2.F09
## NA19098.r2.F10
## NA19098.r2.F11
## NA19098.r2.F12
## NA19098.r2.G01
## NA19098.r2.G02
## NA19098.r2.G03
## NA19098.r2.G05
## NA19098.r2.G06
## NA19098.r2.G08
## NA19098.r2.G09
## NA19098.r2.G10
## NA19098.r2.G11
## NA19098.r2.H01
## NA19098.r2.H02
## NA19098.r2.H03
## NA19098.r2.H04
## NA19098.r2.H05
## NA19098.r2.H06
## NA19098.r2.H07
## NA19098.r2.H08
## NA19098.r2.H10
## NA19098.r2.H12
## NA19101.r3.A02
## NA19101.r3.C12
## NA19101.r3.D01
## NA19101.r3.E08
## Variables with highest loadings for PC1 and PC2:
## 
##                                            PC1         PC2
## ---------------------------------  -----------  ----------
## pct_counts_top_100_features          0.4771343   0.3009332
## pct_counts_feature_controls          0.4735839   0.3309562
## n_detected_feature_controls          0.1332811   0.5367629
## log10_counts_feature_controls       -0.1427373   0.5911762
## total_features                      -0.5016681   0.2936705
## log10_counts_endogenous_features    -0.5081855   0.2757918
```

![](readme_files/figure-html/unnamed-chunk-6-1.png)<!-- -->

```r
knitr::kable(
  as.data.frame(table(dat_pca$outlier)),
  booktabs = TRUE,
  row.names = FALSE,
  caption = 'The number of cells removed by automatic filter (FALSE)'
)
```



Table: The number of cells removed by automatic filter (FALSE)

Var1     Freq
------  -----
FALSE     791
TRUE       73

## 基于基因的过滤

主要是看看表达量很高的那些基因是什么情况


```r
scater::plotQC(dat_qc, type = "highest-expression")
```

![](readme_files/figure-html/unnamed-chunk-7-1.png)<!-- -->

如图所示，表达排名前50个基因中有15个是外源的ERCC spike-ins，如果下次再做实验，可以考虑把外源的ERCC spike-ins稀释一下。

### 根据基因表达量来过滤

这个过滤的threshold取决于测序深度，需要仔细权衡。基于UMI技术的，标准可以是至少有两个细胞检测到了该基因的至少一个转录本。如果是基于reads的，那么至少有两个以上细胞检测到了该基因的至少5个reads。值得注意的是，要先过滤细胞，再过滤基因。


```r
filter_genes <- apply(counts(dat_qc[ , pData(dat_qc)$use]), 1, 
                      function(x) length(x[x > 1]) >= 2)
fData(dat_qc)$use <- filter_genes
knitr::kable(
  as.data.frame(table(filter_genes)),
  booktabs = TRUE,
  row.names = FALSE,
  caption = 'The number of genes removed by gene filter (FALSE)'
)
```



Table: The number of genes removed by gene filter (FALSE)

filter_genes     Freq
-------------  ------
FALSE            4660
TRUE            14066

有4660个基因不合格，被剔除掉了。

## 最终留下来的数据是


```r
dim(dat_qc[fData(dat_qc)$use, pData(dat_qc)$use])
```

```
## Features  Samples 
##    14066      657
```

657个样本的14066个基因的表达量。

PS:注意一下，上面的过滤基于的是UMI方法算出的表达矩阵，如果是纯粹的reads表达矩阵，结果会略微有一点不同。
the ERCC and MT filters are more strict for the reads-based analysis.

# 数据可视化-聚类

首先把reads  counts矩阵进行对数转换，如下：


```r
set_exprs(dat_qc, "log2_counts") <- log2(counts(dat_qc) + 1)
library(scater, quietly = TRUE)
options(stringsAsFactors = FALSE) 
dat_filter <-dat_qc[fData(dat_qc)$use, pData(dat_qc)$use]
endog_genes <- !fData(dat_filter)$is_feature_control
```

这里dat_qc是质控前的数据集，dat_filter是质控后的数据集，一般来说，我们做可视化只需要考虑细胞自身的基因表达量即可


### PCA
基于原始reads counts矩阵和对数化矩阵分开可视化

```r
scater::plotPCA(dat_qc[endog_genes, ],
                ntop = 500,
                colour_by = "batch",
                size_by = "total_features",
                shape_by = "individual",
                exprs_values = "counts")
```

![](readme_files/figure-html/unnamed-chunk-11-1.png)<!-- -->

```r
scater::plotPCA(dat_qc[endog_genes, ],
                ntop = 500,
                colour_by = "batch",
                size_by = "total_features",
                shape_by = "individual",
                exprs_values = "log2_counts")
```

![](readme_files/figure-html/unnamed-chunk-11-2.png)<!-- -->

```r
## 去除低质量细胞，和信号微弱的基因之后。
scater::plotPCA(dat_filter[endog_genes, ],
                ntop = 500,
                colour_by = "batch",
                size_by = "total_features",
                shape_by = "individual",
                exprs_values = "log2_counts")
```

![](readme_files/figure-html/unnamed-chunk-11-3.png)<!-- -->
很明显对数转换后的数据更适合把不同批次不同来源的细胞分开。它降低了第一主成分的可解释度，而且使得表达值的分布更趋近于正态性。但是仅仅是对数转换不足以去除细胞之间的测序技术误差，比如测序深度等等。真正的下游分析可以选择CPM的归一化方法。

### t-SNE

从算法的角度来说， tSNE (t-Distributed Stochastic Neighbor Embedding) combines dimensionality reduction (e.g. PCA) with random walks on the nearest-neighbour network  是优于PCA的。


```r
scater::plotTSNE(dat_qc[endog_genes, ],
                 ntop = 500,
                 perplexity = 130,
                 colour_by = "batch",
                 size_by = "total_features",
                 shape_by = "individual",
                 exprs_values = "log2_counts",
                 rand_seed = 123456)
```

![](readme_files/figure-html/unnamed-chunk-12-1.png)<!-- -->

```r
scater::plotTSNE(dat_filter[endog_genes, ],
                 ntop = 500,
                 perplexity = 130,
                 colour_by = "batch",
                 size_by = "total_features",
                 shape_by = "individual",
                 exprs_values = "log2_counts",
                 rand_seed = 123456)
```

![](readme_files/figure-html/unnamed-chunk-12-2.png)<!-- -->

# 剔除可能的 confounders, artifacts and biases




