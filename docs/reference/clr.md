# Centered log ratio transform

Compute the centered log ratio (CLR) transform of a count matrix.

## Usage

``` r
clr(counts, pseudocount = 0.5)
```

## Arguments

- counts:

  count data with samples as rows and variables are columns

- pseudocount:

  added to counts to avoid issues with zeros

## Value

matrix of CLR transformed counts

## Details

The CLR of a vector `x` of counts in `D` categories is defined as
`clr(x) = log(x) - mean(log(x))`. For details see van den Boogaart and
Tolosana-Delgado (2013).

## References

Van den Boogaart, K. Gerald, and Raimon Tolosana-Delgado. Analyzing
compositional data with R. Vol. 122. Berlin: Springer, 2013.

## See also

`compositions::clr()`

## Examples

``` r
# set probability of each category
prob <- c(0.1, 0.2, 0.3, 0.5)

# number of total counts
countsTotal <- 300

# number of samples
n_samples <- 5

# simulate info for each sample
info <- data.frame(Age = rgamma(n_samples, 50, 1))
rownames(info) <- paste0("sample_", 1:n_samples)

# simulate counts from multinomial
counts <- t(rmultinom(n_samples, size = countsTotal, prob = prob))
colnames(counts) <- paste0("cat_", 1:length(prob))
rownames(counts) <- paste0("sample_", 1:n_samples)

# centered log ratio
clr(counts)
#>               cat_1      cat_2     cat_3     cat_4
#> sample_1 -0.8764783 -0.1372400 0.2221340 0.7915844
#> sample_2 -0.5899417 -0.2534694 0.1863204 0.6570906
#> sample_3 -0.5751955 -0.1740688 0.1114370 0.6378272
#> sample_4 -0.9834305 -0.1240479 0.3538077 0.7536706
#> sample_5 -0.7635900 -0.2205036 0.2788277 0.7052659
```
