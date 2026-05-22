# Inverse of Centered log ratio transform

Compute the inverse centered log ratio (CLR) transform of a count
matrix.

## Usage

``` r
clrInv(x)
```

## Arguments

- x:

  CLR transform values

## Value

matrix of fractions

## Details

Given the CLR transformed values, compute the original fractions

## References

Van den Boogaart, K. Gerald, and Raimon Tolosana-Delgado. Analyzing
compositional data with R. Vol. 122. Berlin: Springer, 2013.

## See also

`compositions::clrInv()`

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

# Fractions
counts / rowSums(counts)
#>               cat_1     cat_2     cat_3     cat_4
#> sample_1 0.10333333 0.1700000 0.2633333 0.4633333
#> sample_2 0.10666667 0.1766667 0.2700000 0.4466667
#> sample_3 0.07666667 0.2266667 0.2233333 0.4733333
#> sample_4 0.08666667 0.2100000 0.2833333 0.4200000
#> sample_5 0.06000000 0.1966667 0.3133333 0.4300000

# centered log ratio, with zero pseudocount
clr(counts, 0)
#>               cat_1       cat_2      cat_3     cat_4
#> sample_1 -0.7334465 -0.23560802 0.20201420 0.7670403
#> sample_2 -0.7163433 -0.21178728 0.21236996 0.7157606
#> sample_3 -0.9933862  0.09062731 0.07581222 0.8269467
#> sample_4 -0.9119446 -0.02690638 0.27261015 0.6662408
#> sample_5 -1.2023823 -0.01521665 0.45054069 0.7670583

# recover fractions from CLR transformed values
clrInv(clr(counts, 0))
#>               cat_1     cat_2     cat_3     cat_4
#> sample_1 0.10333333 0.1700000 0.2633333 0.4633333
#> sample_2 0.10666667 0.1766667 0.2700000 0.4466667
#> sample_3 0.07666667 0.2266667 0.2233333 0.4733333
#> sample_4 0.08666667 0.2100000 0.2833333 0.4200000
#> sample_5 0.06000000 0.1966667 0.3133333 0.4300000
```
