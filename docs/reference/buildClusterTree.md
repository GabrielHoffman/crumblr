# Perform hierarchical clustering on reducedDim

Perform hierarchical clustering dimension reduction from single cell
expression data

## Usage

``` r
buildClusterTree(
  sce,
  reduction,
  labelCol,
  method.dist = c("cosine", "euclidean", "maximum", "manhattan", "canberra", "binary",
    "minkowski"),
  method.hclust = c("complete", "ward.D", "ward.D2")
)
```

## Arguments

- sce:

  `SingleCellExperiment` object

- reduction:

  field of reducedDims(sce) to use

- labelCol:

  column in `SingleCellExperiment` storing cell type annotations

- method.dist:

  method for `dist(..,method=method.dist)`

- method.hclust:

  method for `hclust(..,method=method.hclust)`

## Value

hierarchical clustering computed by
[`hclust()`](https://rdrr.io/r/stats/hclust.html)

## Examples

``` r
library(muscat)

data(example_sce)

hcl_test = buildClusterTree(example_sce, "TSNE", "cluster_id")
```
