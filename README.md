### <u>C</u>ount <u>r</u>atio <u>u</u>ncertainty <u>m</u>odeling <u>b</u>ased <u>l</u>inear <u>r</u>egression <img src="man/figures/logo.png" align="right" alt="" width="120"/>

The `crumblr` package enables analysis of count ratio data using weighted linear models.  `crumblr`'s fast approximation to generalized linear models (GLM's) and generalized linear mixed models (GLMM's) allows used standard workflows to analyize count ratio data.  This can be easier, faster, and more flexible than other models of count data.

## Example usage:

```r
library(crumblr)

# Load cell counts from Kang, et al. (2018)
#  https://doi.org/10.1038/nbt.4042
data(IFNCellCounts)

# Apply crumblr transformation 
# cobj is an EList object compatable with limma workflow
# cobj$E stores transformed values 
# cobj$weights stores precision weights
cobj = crumblr(cellCounts)

# Use limma work flow to analyze each cell
# Perform regression on each cell type separately
#  then use eBayes to shrink residual variance
library(limma)

design = model.matrix(~StimStatus, info)
fit = lmFit(cobj, design)
fit = eBayes(fit)

# extract results for each cell type
topTable(fit, coef="StimStatusstim", number=Inf, sort.by="none")
#>                         logFC    AveExpr          t   P.Value adj.P.Val
#> B cells           -0.09882673  0.5516882 -0.6989167 0.4943409 0.7706781
#> CD14+ Monocytes   -0.12323960  1.2698117 -0.9641347 0.3489000 0.7706781
#> CD4 T cells       -0.08092309  2.0201947 -0.3882013 0.7028311 0.7706781
#> CD8 T cells       -0.14896362  0.0857175 -0.3589998 0.7241489 0.7706781
#> Dendritic cells    0.32314289 -2.1849234  1.1280014 0.2754583 0.7706781
#> FCGR3A+ Monocytes  0.07153536 -0.2567492  0.3637101 0.7206938 0.7706781
#> Megakaryocytes     0.06001827 -1.8655172  0.2963275 0.7706781 0.7706781
#> NK cells           0.08379008  0.3797777  0.5846539 0.5666914 0.7706781
```



## Install
```r
# repo is currently private, so need to include your userid and password
devtools::install_github("GabrielHoffman/crumblr", auth_token=XXXXX)
```
