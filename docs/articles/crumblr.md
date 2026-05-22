# Using crumblr in practice

### Introduction

Changes in cell type composition play an important role in health and
disease. Recent advances in single cell technology have enabled
measurement of cell type composition at increasing cell lineage
resolution across large cohorts of individuals. Yet this raises new
challenges for statistical analysis of these compositional data to
identify changes associated with a phenotype. We introduce `crumblr`, a
scalable statistical method for analyzing count ratio data using
precision-weighted linear models incorporating random effects for
complex study designs. Uniquely, `crumblr` performs tests of association
at multiple levels of the cell lineage hierarchy using multivariate
regression to increase power over tests of a single cell component. In
simulations, `crumblr` increases power compared to existing methods,
while controlling the false positive rate.

The `crumblr` package integrates with the
[`variancePartition`](https://www.bioconductor.org/packages/variancePartition/)
and [`dreamlet`](https://www.bioconductor.org/packages/dreamlet/)
packages in the Bioconductor ecosystem.

Here we consider counts for 8 cell types from quantified using single
cell RNA-seq data from unstimulated and interferon β stimulated PBMCs
from 8 subjects [(Kang, et al.,
2018)](https://www.nature.com/articles/nbt.4042).

The functions here incorporate the precision weights:

- [`variancePartition::fitExtractVarPartModel()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/fitExtractVarPartModel-method.md)
- [`variancePartition::dream()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/dream-method.md)
- [`limma::lmFit()`](https://rdrr.io/pkg/limma/man/lmFit.html)

## Installation

To install this package, start R and enter:

``` r

# 1) Make sure Bioconductor is installed
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# 2) Install crumblr and dependencies:
# From Bioconductor
BiocManager::install("crumblr")
```

## Analysis workflow

### Process data

Here we evaluate whether the observed cell proportions change in
response to interferon β. Given the results here, we cannot reject the
null hypothesis that interferon β does not affect the cell type
proportions.

``` r

library(crumblr)

# Load cell counts, clustering and metadata
# from Kang, et al. (2018) https://doi.org/10.1038/nbt.4042
data(IFNCellCounts)

# Apply crumblr transformation
# cobj is an EList object compatable with limma workflow
# cobj$E stores transformed values
# cobj$weights stores precision weights
#    corresponding to the regularized inverse variance
cobj <- crumblr(df_cellCounts)
```

### Variance partitioning

Decomposing the variance illustrates that more variation is explained by
subject than stimulation status.

``` r

library(variancePartition)

# Partition variance into components for Subject (i.e. ind)
#   and stimulation status, and residual variation
form <- ~ (1 | ind) + (1 | StimStatus)
vp <- fitExtractVarPartModel(cobj, form, info)

# Plot variance fractions
fig.vp <- plotPercentBars(vp)
fig.vp
```

![](crumblr_files/figure-html/vp-1.png)

### PCA

Performing PCA on the transformed cell counts indicates that the samples
cluster based on subject rather than stimulation status. Here, we use
[`standardize()`](http://DiseaseNeurogenomics.github.io/crumblr/reference/standardize-method.md)
so that each observation has approximately equal variance
(i.e. homoscedasticity) by dividing the CLR transformed frequencies by
their estimated sampling standard deviation. Transforming the data to be
approximately homoscedastic has been
[shown](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8)
to improve performance of PCA.

``` r

library(ggplot2)

# Perform PCA
# use crumblr::standardize() to get values with
# approximately equal sampling variance,
# which is a key property for downstream PCA and clustering analysis.
pca <- prcomp(t(standardize(cobj)))

# merge with metadata
df_pca <- merge(pca$x, info, by = "row.names")

# Plot PCA
#   color by Subject
#   shape by Stimulated vs unstimulated
ggplot(df_pca, aes(PC1, PC2, color = as.character(ind), shape = StimStatus)) +
  geom_point(size = 3) +
  theme_classic() +
  theme(aspect.ratio = 1) +
  scale_color_discrete(name = "Subject") +
  xlab("PC1") +
  ylab("PC2")
```

![](crumblr_files/figure-html/pca-1.png)

### Hierarchical clustering

The samples from the same subject also cluster together.

``` r

heatmap(cobj$E)
```

![](crumblr_files/figure-html/hclust-1.png)

### Differential testing

``` r

# Use variancePartition workflow to analyze each cell type
# Perform regression on each cell type separately
#  then use eBayes to shrink residual variance
# Also compatible with limma::lmFit()
fit <- dream(cobj, ~ StimStatus + ind, info)
fit <- eBayes(fit)

# Extract results for each cell type
topTable(fit, coef = "StimStatusstim", number = Inf)
```

    ##                         logFC    AveExpr          t     P.Value  adj.P.Val         B
    ## CD8 T cells       -0.25085170  0.0857175 -4.0787416 0.002436375 0.01949100 -1.279815
    ## Dendritic cells    0.37386979 -2.1849234  3.1619195 0.010692544 0.02738587 -2.638507
    ## CD14+ Monocytes   -0.10525402  1.2698117 -3.1226341 0.011413912 0.02738587 -2.709377
    ## B cells           -0.10478652  0.5516882 -3.0134349 0.013692935 0.02738587 -2.940542
    ## CD4 T cells       -0.07840101  2.0201947 -2.2318104 0.050869691 0.08139151 -4.128069
    ## FCGR3A+ Monocytes  0.07425165 -0.2567492  1.6647681 0.128337022 0.17111603 -4.935304
    ## NK cells           0.10270672  0.3797777  1.5181860 0.161321761 0.18436773 -5.247806
    ## Megakaryocytes     0.01377768 -1.8655172  0.1555131 0.879651456 0.87965146 -6.198336

#### Multivariate testing along a tree

We can gain power by jointly testing multiple cell types using a
multivariate statistical model, instead of testing one cell type at a
time. Here we construct a hierarchical clustering between cell types
based on gene expression from pseudobulk, and perform a multivariate
test for each internal node of the tree based on its leaf nodes. The
results for the leaves are the same as from `topTable()` above. At each
internal node
[`treeTest()`](http://DiseaseNeurogenomics.github.io/crumblr/reference/treeTest.md)
performs a fixed effects meta-analysis of the coefficients of the leaves
while modeling the covariance between coefficient estimates. In the
backend, this is implemented using
[`variancePartition::mvTest()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/mvTest-method.md)
and [remaCor](https://cran.r-project.org/package=remaCor) package.

Here the hierarchical clustering, `hcl`, is precomputed from pseudobulk
gene expression using `buildClusterTreeFromPB()`.

``` r

# Perform multivariate test across the hierarchy
res <- treeTest(fit, cobj, hcl, coef = "StimStatusstim")

# Plot hierarchy and testing results
plotTreeTest(res)
```

![](crumblr_files/figure-html/treeTest-1.png)

``` r

# Plot hierarchy and regression coefficients
plotTreeTestBeta(res)
```

![](crumblr_files/figure-html/treeTest-2.png)

##### Combined plotting

``` r

plotTreeTestBeta(res) +
  theme(legend.position = "bottom", legend.box = "vertical") |
  plotForest(res, hide = FALSE) |
  fig.vp
```

![](crumblr_files/figure-html/combined-1.png)

\

### Hierarchical clustering

The hierarchical clustering used by
[`treeTest()`](http://DiseaseNeurogenomics.github.io/crumblr/reference/treeTest.md)
can be computed a number of ways, depending on the available data and
biological question. For example, see
[Article](http://DiseaseNeurogenomics.github.io/crumblr/articles/integration.md)
for details about how hierarchical clustering was run in this dataset.

In general, hierarchical clustering can be computed from

- pseudobulked single cell gene expression

``` r

hcl <- buildClusterTreeFromPB(pb)
```

- cell type frequencies

``` r

# correlation matrix between all cell types
C <- cor(t(standardize(cobj)))

# convert to distance
dm <- as.dist(1 - abs(C))

# eval hierarchical clustering
hcl <- hclust(dm)
```

- [Newick formated](https://en.wikipedia.org/wiki/Newick_format) tree
  computed from external data

``` r

# Make sure packages are installed
# BiocManager::install(c("ctc", "ape", "phylogram"))
library(ape)
library(ctc)
library(phylogram)
library(tidyverse)

# Write tree to file, 
# edit manually
# then read back into R
# 

# Specify tree as text in Newick format
txt = "((CD14+ Monocytes,(B cells,(Dendritic cells,Megakaryocytes))),(CD8 T cells,(NK cells,(CD4 T cells,FCGR3A+ Monocytes))));"

# read from text
hcl_from_txt <- read.tree(text = txt) %>%
                  as.dendrogram.phylo %>%
                  as.hclust

# Alternatively, 
# write existing tree to file
# and edit manaully
write(hc2Newick(hcl),file='hcl.nwk')

hcl_from_txt2 <- read.tree(file = 'hcl.nwk') %>%
                  as.dendrogram.phylo %>%
                  as.hclust
```

## Considerations

### Computational scaling

The
[`crumblr()`](http://DiseaseNeurogenomics.github.io/crumblr/reference/crumblr.md)
workflow is very fast, especially compared to more demanding
differential expression analysies. Applying
[`crumblr()`](http://DiseaseNeurogenomics.github.io/crumblr/reference/crumblr.md)
requires \<1 sec, even for very large datsets. Differential testing with
[`dream()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/dream-method.md)
takes \< 10 seconds for fixed effects models and \< 30 seconds for mixed
models with typical sample sizes and number of cell types. Running
[`treeTest()`](http://DiseaseNeurogenomics.github.io/crumblr/reference/treeTest.md)
can be a little more demanding, but should finish in \< 30 seconds with
20 cell types.

### Complex study designs

The
[`crumblr()`](http://DiseaseNeurogenomics.github.io/crumblr/reference/crumblr.md)
workflow can handle complex study designs with repeated measures or
multiple random effects.
[`dream()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/dream-method.md)
uses [`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html) in the
backend to fit weighted linear mixed models. Considerations about
defining the regression model are described in documentation to
[`variancePartition`](https://diseaseneurogenomics.github.io/variancePartition)
or this
[book](https://people.math.ethz.ch/~maechler/MEMo-pages/lMMwR.pdf) by
the author of `lme4`.

### Tuning parameters

[`crumblr()`](http://DiseaseNeurogenomics.github.io/crumblr/reference/crumblr.md)
uses two tuning parameters accessable to the user. These are fixed to
default values in all simulations and data analysis. *The user is
strongly recommended to keep these dfault values*.

- In order to deal with zero counts,
  [`crumblr()`](http://DiseaseNeurogenomics.github.io/crumblr/reference/crumblr.md)
  uses a pseudocount (default: 0.5) added to all observed counts.

- For real data, the asymptotic variance formula can give weights that
  vary substantially across samples and give very high weights for a
  subset of samples. In order to address this, we regularize the weights
  to reduce the variation in the weights to have a maximum ratio
  (default of 5) between the maximum and specified quantile value
  (default of 5%).

## Session Info

    ## R version 4.5.1 (2025-06-13)
    ## Platform: aarch64-apple-darwin23.6.0
    ## Running under: macOS Sonoma 14.7.1
    ## 
    ## Matrix products: default
    ## BLAS/LAPACK: /opt/homebrew/Cellar/openblas/0.3.33/lib/libopenblasp-r0.3.33.dylib;  LAPACK version 3.12.0
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## time zone: America/New_York
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] muscat_1.24.0            crumblr_0.99.22          variancePartition_1.43.1
    ## [4] BiocParallel_1.44.0      limma_3.66.0             ggplot2_4.0.3           
    ## [7] BiocStyle_2.38.0        
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] splines_4.5.1               bitops_1.0-9                ggplotify_0.1.3            
    ##   [4] tibble_3.3.1                lpsymphony_1.38.0           lifecycle_1.0.5            
    ##   [7] Rdpack_2.6.6                edgeR_4.8.2                 doParallel_1.0.17          
    ##  [10] globals_0.19.1              lattice_0.22-9              MASS_7.3-65                
    ##  [13] backports_1.5.1             magrittr_2.0.5              sass_0.4.10                
    ##  [16] rmarkdown_2.31              jquerylib_0.1.4             yaml_2.3.12                
    ##  [19] otel_0.2.0                  sctransform_0.4.3           minqa_1.2.8                
    ##  [22] RColorBrewer_1.1-3          multcomp_1.4-30             abind_1.4-8                
    ##  [25] EnvStats_3.1.0              glmmTMB_1.1.14              GenomicRanges_1.62.1       
    ##  [28] purrr_1.2.2                 BiocGenerics_0.56.0         TH.data_1.1-5              
    ##  [31] yulab.utils_0.2.4           rappdirs_0.3.4              sandwich_3.1-1             
    ##  [34] gdtools_0.5.0               circlize_0.4.18             IRanges_2.44.0             
    ##  [37] S4Vectors_0.48.1            ggrepel_0.9.8               pbkrtest_0.5.5             
    ##  [40] irlba_2.3.7                 listenv_0.10.1              tidytree_0.4.7             
    ##  [43] parallelly_1.47.0           pkgdown_2.2.0               codetools_0.2-20           
    ##  [46] DelayedArray_0.36.1         scuttle_1.20.0              tidyselect_1.2.1           
    ##  [49] shape_1.4.6.1               aplot_0.2.9                 farver_2.1.2               
    ##  [52] ScaledMatrix_1.18.0         lme4_2.0-1                  viridis_0.6.5              
    ##  [55] matrixStats_1.5.0           stats4_4.5.1                Seqinfo_1.0.0              
    ##  [58] jsonlite_2.0.0              BiocNeighbors_2.4.0         GetoptLong_1.1.1           
    ##  [61] scater_1.38.1               survival_3.8-6              iterators_1.0.14           
    ##  [64] emmeans_2.0.3               systemfonts_1.3.2           foreach_1.5.2              
    ##  [67] progress_1.2.3              tools_4.5.1                 treeio_1.34.0              
    ##  [70] ragg_1.5.2                  Rcpp_1.1.1-1.1              blme_1.0-7                 
    ##  [73] glue_1.8.1                  gridExtra_2.3               SparseArray_1.10.10        
    ##  [76] xfun_0.57                   mgcv_1.9-4                  DESeq2_1.50.2              
    ##  [79] MatrixGenerics_1.22.0       dplyr_1.2.1                 withr_3.0.2                
    ##  [82] numDeriv_2016.8-1.1         BiocManager_1.30.27         fastmap_1.2.0              
    ##  [85] boot_1.3-32                 rsvd_1.0.5                  caTools_1.18.3             
    ##  [88] digest_0.6.39               estimability_1.5.1          R6_2.6.1                   
    ##  [91] gridGraphics_0.5-1          textshaping_1.0.5           colorspace_2.1-2           
    ##  [94] gtools_3.9.5                dichromat_2.0-0.1           RhpcBLASctl_0.23-42        
    ##  [97] tidyr_1.3.2                 generics_0.1.4              fontLiberation_0.1.0       
    ## [100] data.table_1.18.4           corpcor_1.6.10              prettyunits_1.2.0          
    ## [103] htmlwidgets_1.6.4           S4Arrays_1.10.1             pkgconfig_2.0.3            
    ## [106] gtable_0.3.6                ComplexHeatmap_2.26.1       S7_0.2.2                   
    ## [109] SingleCellExperiment_1.32.0 XVector_0.50.0              remaCor_0.0.20             
    ## [112] htmltools_0.5.9             fontBitstreamVera_0.1.1     bookdown_0.46              
    ## [115] TMB_1.9.21                  zigg_0.0.2                  clue_0.3-68                
    ## [118] scales_1.4.0                Biobase_2.70.0              png_0.1-9                  
    ## [121] fANCOVA_0.6-1               reformulas_0.4.4            ggfun_0.2.0                
    ## [124] knitr_1.51                  reshape2_1.4.5              rjson_0.2.23               
    ## [127] coda_0.19-4.1               nlme_3.1-169                nloptr_2.2.1               
    ## [130] cachem_1.1.0                zoo_1.8-15                  GlobalOptions_0.1.4        
    ## [133] stringr_1.6.0               KernSmooth_2.23-26          vipor_0.4.7                
    ## [136] parallel_4.5.1              desc_1.4.3                  pillar_1.11.1              
    ## [139] grid_4.5.1                  vctrs_0.7.3                 gplots_3.3.0               
    ## [142] slam_0.1-55                 BiocSingular_1.26.1         beachmat_2.26.0            
    ## [145] IHW_1.38.0                  xtable_1.8-8                cluster_2.1.8.2            
    ## [148] beeswarm_0.4.0              evaluate_1.0.5              mvtnorm_1.3-7              
    ## [151] cli_3.6.6                   locfit_1.5-9.12             compiler_4.5.1             
    ## [154] rlang_1.2.0                 crayon_1.5.3                future.apply_1.20.2        
    ## [157] labeling_0.4.3              fdrtool_1.2.18              ggbeeswarm_0.7.3           
    ## [160] plyr_1.8.9                  fs_2.1.0                    ggiraph_0.9.6              
    ## [163] stringi_1.8.7               viridisLite_0.4.3           lmerTest_3.2-1             
    ## [166] lazyeval_0.2.3              aod_1.3.3                   fontquiver_0.2.1           
    ## [169] Matrix_1.7-5                hms_1.1.4                   patchwork_1.3.2            
    ## [172] future_1.70.0               statmod_1.5.2               SummarizedExperiment_1.40.0
    ## [175] rbibutils_2.4.1             Rfast_2.1.5.2               broom_1.0.13               
    ## [178] RcppParallel_5.1.11-2       bslib_0.11.0                ggtree_4.0.5               
    ## [181] ape_5.8-1
