# Forest plot

Forest plot

Forest plot of effect size estimates at the leaves of the tree

## Usage

``` r
plotForest(x, ...)

# S4 method for class 'treedata'
plotForest(x, ..., hide = FALSE)
```

## Arguments

- x:

  result from
  [`treeTest()`](http://DiseaseNeurogenomics.github.io/crumblr/reference/treeTest.md)

- ...:

  other arguments

- hide:

  hide rownames and legend

## Value

ggplot2 object

## Examples

``` r
library(variancePartition)

# Load cell counts, clustering and metadata
# from Kang, et al. (2018) https://doi.org/10.1038/nbt.4042
data(IFNCellCounts)

# Apply crumblr transformation
cobj <- crumblr(df_cellCounts)

# Use dream workflow to analyze each cell separately
fit <- dream(cobj, ~ StimStatus + ind, info)
fit <- eBayes(fit)

# Perform multivariate test across the hierarchy
res <- treeTest(fit, cobj, hcl, coef = "StimStatusstim")

# Plot log fold changes from coef
plotForest(res)

```
