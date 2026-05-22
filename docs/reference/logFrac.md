# Log fractions and precision weights

Compute log fractions and precision weights from matrix of c ounts,
where columns are variables and rows are samples

## Usage

``` r
logFrac(counts, pseudocount = 0.5, max.ratio = 5, quant = 0.05)
```

## Arguments

- counts:

  count data with samples as rows and variables are columns

- pseudocount:

  added to counts to avoid issues with zeros

- max.ratio:

  regularize estimates of the weights to have a maximum ratio of
  `max.ratio` between the maximum and `quant` quantile value

- quant:

  quantile value used for `max.ratio`

## Value

An `EList` object with the following components:

- E: :

  numeric matrix of log transformed counts

- weights: :

  numeric matrix of observation-level inverse-variance weights

## Details

For real data, the asymptotic variance formula can give weights that
vary substantially across samples and give very high weights for a
subset of samples. In order to address this, we regularize the weights
to reduce the variation in the weights to have a maximum ratio of
`max.ratio` between the maximum and `quant` quantile value.

## See also

[`limma::voom()`](https://rdrr.io/pkg/limma/man/voom.html),
[`variancePartition::dream()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/dream-method.md)

## Examples

``` r
# set probability of each category
prob <- c(0.1, 0.2, 0.3, 0.5)

# number of total counts
countsTotal <- 300

# number of samples
n_samples <- 100

# simulate info for each sample
info <- data.frame(Age = rgamma(n_samples, 50, 1))
rownames(info) <- paste0("sample_", 1:n_samples)

# simulate counts from multinomial
counts <- t(rmultinom(n_samples, size = countsTotal, prob = prob))
colnames(counts) <- paste0("cat_", 1:length(prob))
rownames(counts) <- paste0("sample_", 1:n_samples)

# run logFrac on counts
cobj <- logFrac(counts)

# run standard variancePartition analysis on crumblr results
library(variancePartition)

fit <- dream(cobj, ~ Age, info)
fit <- eBayes(fit)

topTable(fit, coef = "Age", sort.by = "none")
#>               logFC    AveExpr          t   P.Value adj.P.Val         B
#> cat_1  0.0004874472 -2.4063982  0.1499668 0.8810907 0.9119836 -8.856395
#> cat_2 -0.0013536071 -1.6940978 -0.7022378 0.4841504 0.9119836 -8.562192
#> cat_3  0.0008978300 -1.3177560  0.5717690 0.5687521 0.9119836 -8.642022
#> cat_4 -0.0001136616 -0.7907723 -0.1108149 0.9119836 0.9119836 -8.791745
```
