# MLE for Dirichlet Multinomial

MLE for Dirichlet Multinomial

## Usage

``` r
dmn_mle(counts, ...)
```

## Arguments

- counts:

  matrix with rows as samples and columns as categories

- ...:

  additional arguments passed to
  [`optim()`](https://rdrr.io/r/stats/optim.html)

## Value

list storing `alpha` parameter estimates, logLik, and details about
convergence

- `alpha`:

  estimated \\alpha\\ parameters

- `overdispersion`:

  Overdispersion value \\1 + \rho^2(n-1)\\ compared to multinomial

- `logLik`:

  value of function

- `scale`:

  scaling of \\\alpha\\ parameters computed in a second optimization
  step

- `evals`:

  number of function evaluations in step 1

- `convergence`:

  convergence details from step 1

## Details

Maximize Dirichlet Multinomial (DMN) log-likelihood with
[`optim()`](https://rdrr.io/r/stats/optim.html) using log likelihood
function and its gradient. This method uses a second round of
optimization to estimate the scale of \\\alpha\\ parameters, which is
necessary for accurate estimation of overdispersion metric.

The covariance between counts in each category from DMN distributed data
is \\n(diag(p) - pp^T) (1 + \rho^2(n-1))\\ for \\n\\ total counts, and
vector of proportions \\p\\, where \\\rho^2 = 1 / (a_0 + 1)\\ and \\a_0
= \sum_i \alpha_i\\. The count data is overdispersed by a factor of
\\1 + \rho^2(n-1)\\ compared to a multinomial (MN) distribution. As
\\a_0\\ increases, the DMN converges to the MN.

See
<https://en.wikipedia.org/wiki/Dirichlet-multinomial_distribution#Matrix_notation>

## See also

Other functions also estimate DMN parameters. `MGLM::MGLMfit()` and
[`dirmult::dirmult()`](https://rdrr.io/pkg/dirmult/man/dirmult.html)
give good parameter estimates but are slower.
[`Rfast::dirimultinom.mle()`](https://rdrr.io/pkg/Rfast/man/multinom.mle.html)
often fails to converge

## Examples

``` r
library(HMP)
#> Loading required package: dirmult
#> 
#> Attaching package: ‘HMP’
#> The following object is masked from ‘package:dirmult’:
#> 
#>     weirMoM

set.seed(1)

n_samples <- 1000
n_counts <- 5000
alpha <- c(500, 1000, 2000)

# Dirichlet.multinomial
counts <- Dirichlet.multinomial(rep(n_counts, n_samples), alpha)

fit <- dmn_mle(counts)

fit
#> $alpha
#>    Taxa 1    Taxa 2    Taxa 3 
#>  506.3946 1015.2996 2027.6421 
#> 
#> $overdispersion
#> [1] 2.408036
#> 
#> $logLik
#> [1] -4777957
#> 
#> $scale
#> [1] 0.7098813
#> 
#> $evals
#> function gradient 
#>        2        2 
#> 
#> $convergence
#> [1] 0
#> 

# overdispersion: true value
a0 <- sum(alpha)
rhoSq <- 1 / (a0 + 1)
1 + rhoSq * (n_counts - 1)
#> [1] 2.427878

# multinomial, so overdispersion is 1
counts <- t(rmultinom(n_samples, n_counts, prob = alpha / sum(alpha)))

dmn_mle(counts)
#> $alpha
#> [1] 2165173253 4324433613 8651164924
#> 
#> $overdispersion
#> [1] 1
#> 
#> $logLik
#> [1] -4779160
#> 
#> $scale
#> [1] 70011.04
#> 
#> $evals
#> function gradient 
#>       33       33 
#> 
#> $convergence
#> [1] 0
#> 
#
#
```
