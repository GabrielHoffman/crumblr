# Plot tree coefficients from multivariate testing

Plot tree coefficients from multivariate testing at each node. Only
applicable top fixed effect tests

## Usage

``` r
plotTreeTestBeta(
  tree,
  low = "blue",
  mid = "white",
  high = "red",
  xmax.scale = 1.5,
  fdr.cutoff = 0.05
)
```

## Arguments

- tree:

  phylo object storing tree

- low:

  low color on gradient

- mid:

  mid color on gradient

- high:

  high color on gradient

- xmax.scale:

  expand the x-axis by this factor so leaf labels fit in the plot

- fdr.cutoff:

  value used as the FDR cutoff for significance annotation. Defaults to
  0.05.

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

# Plot hierarchy, no tests are significant at FDR < 0.05
plotTreeTestBeta(res)

```
