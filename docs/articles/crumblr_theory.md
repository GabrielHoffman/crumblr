# Normal approximation vs. empirical simulation

### Asymptotic normal approximation

Let the vector {\bf p} be the true fractions across D categories.
Consider C total counts sampled from a Dirichlet-multinomial (DMN)
distribution with overdispersion \tau, where \tau=1 reduces to the
multinomial distribution. The [centered log
ratio](https://rdrr.io/cran/compositions/man/clr.html) (CLR) of the
i^{th} estimated fraction, \hat p_i is

\text{clr}(\hat p_i) = \log(\hat p_i) - \frac{1}{D}\sum\_{j=1}^D
\log(\hat p_j) , and we show that the sampling variance is well
approximated by

\text{var}\[\text{clr}(\hat p_i)\] = \frac{\tau}{C} \left\[
\frac{1}{\hat p_i} - \frac{2}{ D \hat p_i} + \frac{1}{D^2}\sum\_{j=1}^D
\frac{1}{\hat p_j} \right\] .

### Simulations

The sampling variance is derived from asymototic theory, so we examine
its behavior for finite total counts. Here we evaluate the empirical
variance from 1,000 draws from a Dirichlet-multinomial distribution
while varying D, \tau, and C. A pseudocount of 0.5 is added to the
observed counts since the asymptotic theory is not defined for counts of
zero.

Here we plot the standard deviation after CLR transform from the
empirical DMN and the asymptotic normal approximation under a range of
conditions. Results are shown for instances with at least 2 counts.

#### D=2 categories

![](crumblr_theory_files/figure-html/plot1-1.png)

#### D=15 categories

![](crumblr_theory_files/figure-html/plot2-1.png)

#### Interpretation

The asymptotic standard deviation shows good agreement with the
empirical results even for small values of C, *when at least 2 counts
are observed*. In practice, it is often reasonable to assume a
sufficient number of counts before a variable is included in an
analysis. Importantly, with less than 2 counts the asymptotic theory
gives a *larger* standard deviation than the emprical results (results
not shown). Therefore, this approach is conservative and should not
underestimate the true amount of variation. The asymptotic normal
approximation is most accurate for large total counts C, large
proportions p, and small overdispersion \tau.

##### Consideration of overdispersion

Based on Equation (2), the variance of the CLR-transformed proportions
is a *linear* function of \tau. Importantly, downstream analysis of the
CLR-transformed proportions with a precision-weighted linear (mixed)
model or a variance stabilizing transform depends only on the *relative*
variances. Since relative variances are invariant to the scale of \tau,
for these applications the value of \tau can be set to 1 instead of
being estimated from the data.

For other applications, `crumblr` can estimate \tau from the data by
using `crumblr(counts, tau=NULL)`. This calls `dmn.mle()` to estimate
the parameters of the DMN distribution and is substantially faster than
alternatives.

But note that due the theoretical properties of the variance estimate of
CLR-transformed proportions, the precision weights are invariant to the
scale of \tau. We can see this empirically:

``` r

library(crumblr)

data(IFNCellCounts)

counts <- df_cellCounts

# run crumblr with different tau values
# show part of the weights matrix
crumblr(counts, tau=1)$weights[1:3, 1:3]
```

    ##                  ctrl101 ctrl1015 ctrl1016
    ## B cells         2.704599 5.000000 3.094438
    ## CD14+ Monocytes 1.822824 4.143929 2.716263
    ## CD4 T cells     2.044461 3.615231 2.433429

``` r

crumblr(counts, tau=5)$weights[1:3, 1:3]
```

    ##                  ctrl101 ctrl1015 ctrl1016
    ## B cells         2.704599 5.000000 3.094438
    ## CD14+ Monocytes 1.822824 4.143929 2.716263
    ## CD4 T cells     2.044461 3.615231 2.433429

``` r

crumblr(counts, tau=NULL)$weights[1:3, 1:3]
```

    ##                  ctrl101 ctrl1015 ctrl1016
    ## B cells         2.704599 5.000000 3.094438
    ## CD14+ Monocytes 1.822824 4.143929 2.716263
    ## CD4 T cells     2.044461 3.615231 2.433429

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
    ## [1] parallel  stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] muscat_1.24.0            variancePartition_1.43.1 BiocParallel_1.44.0     
    ##  [4] limma_3.66.0             lubridate_1.9.5          forcats_1.0.1           
    ##  [7] stringr_1.6.0            dplyr_1.2.1              purrr_1.2.2             
    ## [10] readr_2.2.0              tidyr_1.3.2              tibble_3.3.1            
    ## [13] tidyverse_2.0.0          glue_1.8.1               HMP_2.0.1               
    ## [16] dirmult_0.1.3-5          crumblr_1.2.0            ggplot2_4.0.3           
    ## [19] BiocStyle_2.38.0        
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] fs_2.1.0                    matrixStats_1.5.0           bitops_1.0-9               
    ##   [4] RColorBrewer_1.1-3          doParallel_1.0.17           numDeriv_2016.8-1.1        
    ##   [7] sctransform_0.4.3           tools_4.5.1                 backports_1.5.1            
    ##  [10] R6_2.6.1                    vegan_2.7-3                 lazyeval_0.2.3             
    ##  [13] mgcv_1.9-4                  GetoptLong_1.1.1            permute_0.9-10             
    ##  [16] withr_3.0.2                 prettyunits_1.2.0           gridExtra_2.3              
    ##  [19] fdrtool_1.2.18              cli_3.6.6                   Biobase_2.70.0             
    ##  [22] textshaping_1.0.5           sandwich_3.1-1              labeling_0.4.3             
    ##  [25] slam_0.1-55                 sass_0.4.10                 mvtnorm_1.3-7              
    ##  [28] S7_0.2.2                    blme_1.0-7                  pkgdown_2.2.0              
    ##  [31] systemfonts_1.3.2           yulab.utils_0.2.4           dichromat_2.0-0.1          
    ##  [34] scater_1.38.1               parallelly_1.47.0           generics_0.1.4             
    ##  [37] gridGraphics_0.5-1          shape_1.4.6.1               gtools_3.9.5               
    ##  [40] Matrix_1.7-5                ggbeeswarm_0.7.3            S4Vectors_0.48.1           
    ##  [43] abind_1.4-8                 lifecycle_1.0.5             multcomp_1.4-30            
    ##  [46] yaml_2.3.12                 edgeR_4.8.2                 SummarizedExperiment_1.40.0
    ##  [49] gplots_3.3.0                SparseArray_1.10.10         grid_4.5.1                 
    ##  [52] crayon_1.5.3                lattice_0.22-9              beachmat_2.26.0            
    ##  [55] pillar_1.11.1               knitr_1.51                  ComplexHeatmap_2.26.1      
    ##  [58] GenomicRanges_1.62.1        rjson_0.2.23                boot_1.3-32                
    ##  [61] estimability_1.5.1          corpcor_1.6.10              future.apply_1.20.2        
    ##  [64] codetools_0.2-20            ggiraph_0.9.6               ggfun_0.2.0                
    ##  [67] fontLiberation_0.1.0        data.table_1.18.4           vctrs_0.7.3                
    ##  [70] png_0.1-9                   treeio_1.34.0               Rdpack_2.6.6               
    ##  [73] gtable_0.3.6                cachem_1.1.0                zigg_0.0.2                 
    ##  [76] xfun_0.57                   rbibutils_2.4.1             S4Arrays_1.10.1            
    ##  [79] Rfast_2.1.5.2               Seqinfo_1.0.0               coda_0.19-4.1              
    ##  [82] reformulas_0.4.4            survival_3.8-6              SingleCellExperiment_1.32.0
    ##  [85] iterators_1.0.14            statmod_1.5.2               TH.data_1.1-5              
    ##  [88] nlme_3.1-169                pbkrtest_0.5.5              ggtree_4.0.5               
    ##  [91] fontquiver_0.2.1            progress_1.2.3              EnvStats_3.1.0             
    ##  [94] bslib_0.11.0                TMB_1.9.21                  irlba_2.3.7                
    ##  [97] vipor_0.4.7                 KernSmooth_2.23-26          otel_0.2.0                 
    ## [100] rpart_4.1.27                colorspace_2.1-2            BiocGenerics_0.56.0        
    ## [103] DESeq2_1.50.2               tidyselect_1.2.1            emmeans_2.0.3              
    ## [106] compiler_4.5.1              BiocNeighbors_2.4.0         desc_1.4.3                 
    ## [109] fontBitstreamVera_0.1.1     DelayedArray_0.36.1         bookdown_0.46              
    ## [112] scales_1.4.0                caTools_1.18.3              remaCor_0.0.20             
    ## [115] rappdirs_0.3.4              digest_0.6.39               minqa_1.2.8                
    ## [118] rmarkdown_2.31              aod_1.3.3                   XVector_0.50.0             
    ## [121] RhpcBLASctl_0.23-42         htmltools_0.5.9             pkgconfig_2.0.3            
    ## [124] lme4_2.0-1                  MatrixGenerics_1.22.0       lpsymphony_1.38.0          
    ## [127] fastmap_1.2.0               rlang_1.2.0                 GlobalOptions_0.1.4        
    ## [130] htmlwidgets_1.6.4           farver_2.1.2                jquerylib_0.1.4            
    ## [133] IHW_1.38.0                  zoo_1.8-15                  jsonlite_2.0.0             
    ## [136] BiocSingular_1.26.1         magrittr_2.0.5              scuttle_1.20.0             
    ## [139] ggplotify_0.1.3             patchwork_1.3.2             Rcpp_1.1.1-1.1             
    ## [142] ape_5.8-1                   viridis_0.6.5               gdtools_0.5.0              
    ## [145] stringi_1.8.7               MASS_7.3-65                 plyr_1.8.9                 
    ## [148] listenv_0.10.1              ggrepel_0.9.8               splines_4.5.1              
    ## [151] hms_1.1.4                   circlize_0.4.18             locfit_1.5-9.12            
    ## [154] reshape2_1.4.5              stats4_4.5.1                ScaledMatrix_1.18.0        
    ## [157] evaluate_1.0.5              RcppParallel_5.1.11-2       rpart.plot_3.1.4           
    ## [160] BiocManager_1.30.27         nloptr_2.2.1                tzdb_0.5.0                 
    ## [163] foreach_1.5.2               future_1.70.0               clue_0.3-68                
    ## [166] rsvd_1.0.5                  broom_1.0.13                xtable_1.8-8               
    ## [169] fANCOVA_0.6-1               tidytree_0.4.7              viridisLite_0.4.3          
    ## [172] ragg_1.5.2                  lmerTest_3.2-1              aplot_0.2.9                
    ## [175] glmmTMB_1.1.14              beeswarm_0.4.0              IRanges_2.44.0             
    ## [178] cluster_2.1.8.2             globals_0.19.1              timechange_0.4.0
