# Count ratio uncertainty modeling based linear regression

Count ratio uncertainty modeling based linear regression (crumblr)
returns CLR-transformed counts and observation-level inverse-variance
weights for use in weighted linear models.

## Usage

``` r
crumblr(
  counts,
  pseudocount = 0.5,
  method = c("clr", "clr_2class"),
  tau = 1,
  max.ratio = 5,
  quant = 0.05
)

# S4 method for class 'matrix'
crumblr(
  counts,
  pseudocount = 0.5,
  method = c("clr", "clr_2class"),
  tau = 1,
  max.ratio = 5,
  quant = 0.05
)

# S4 method for class 'data.frame'
crumblr(
  counts,
  pseudocount = 0.5,
  method = c("clr", "clr_2class"),
  tau = 1,
  max.ratio = 5,
  quant = 0.05
)
```

## Arguments

- counts:

  count data with samples as rows and variables are columns

- pseudocount:

  added to counts to avoid issues with zeros

- method:

  `"clr"` computes standard centered log ratio and precision weights
  based on the delta approximation. `"clr_2class"` computes the
  [`clr()`](http://DiseaseNeurogenomics.github.io/crumblr/reference/clr.md)
  transform for category `i` using 2 classes: 1) counts in category i,
  and 2) counts \_not\_ in category i.

- tau:

  overdispersion parameter for Dirichlet multinomial. If `NULL`,
  estimate from observed counts.

- max.ratio:

  regularize estimates of the weights to have a maximum ratio of
  `max.ratio` between the maximum and `quant` quantile value

- quant:

  quantile value used for `max.ratio`

## Value

An `EList` object with the following components:

- E: :

  numeric matrix of CLR transformed counts

- weights: :

  numeric matrix of observation-level inverse-variance weights

## Details

Evaluate the centered log ratio (CLR) transform of the count matrix, and
the asymptotic theoretical variances of each transformed observation.
The asymptotic normal approximation is increasingly accurate for small
overdispersion \\\tau\\, large total counts \\C\\, and large proportions
\\p\\, but shows good agreement with the empirical results in most
situations. In practice, it is often reasonable to assume a sufficient
number of counts before a variable is included in an analysis anyway.
But the feasibility of this assumption is up to the user to determine.

Given the array `p` storing proportions for one sample across all
categories, the delta approximation uses the term `1/p`. This can be
unstable for small values of `p`, and the estimated variances can be
sensitive to small changes in the proportions. To address this, the
`"clr_2class"` method computes the
[`clr()`](http://DiseaseNeurogenomics.github.io/crumblr/reference/clr.md)
transform for category `i` using 2 classes: 1) counts in category i, and
2) counts \_not\_ in category i. Since class (2) now sums counts across
all other categories, the small proportions are avoided and the variance
estimates are more stable.

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

# run crumblr on counts
cobj <- crumblr(counts)

# run standard variancePartition analysis on crumblr results
library(variancePartition)
#> Loading required package: limma
#> Loading required package: BiocParallel
#> 
#> Attaching package: ‘variancePartition’
#> The following objects are masked from ‘package:limma’:
#> 
#>     eBayes, topTable

fit <- dream(cobj, ~Age, info)
fit <- eBayes(fit)

topTable(fit, coef = "Age", sort.by = "none")
#>               logFC    AveExpr          t    P.Value  adj.P.Val         B
#> cat_1 -0.0007423261 -0.8290861 -0.3516740 0.72577236 0.72577236 -8.787818
#> cat_2 -0.0031282260 -0.1728203 -2.0052997 0.04744828 0.09489656 -6.773709
#> cat_3  0.0032081485  0.2454492  2.2337843 0.02757133 0.09489656 -6.282402
#> cat_4  0.0006850339  0.7564572  0.4775859 0.63391632 0.72577236 -8.611400
```
