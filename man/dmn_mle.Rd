% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dmn_mle.R
\name{dmn_mle}
\alias{dmn_mle}
\title{MLE for Dirichlet Multinomial}
\usage{
dmn_mle(counts, ...)
}
\arguments{
\item{counts}{matrix with rows as samples and columns as categories}

\item{...}{additional arguments passed to \code{optim()}}
}
\value{
list storing \code{alpha} parameter estimates, logLik, and details about convergence
\describe{
 \item{\code{alpha}}{estimated \eqn{alpha} parameters}
 \item{\code{overdispersion}}{Overdispersion value \eqn{1 + \rho^2(n-1)} compared to multinomial}
 \item{\code{logLik}}{value of function}
 \item{\code{scale}}{scaling of \eqn{\alpha} parameters computed in a second optimization step}
 \item{\code{evals}}{number of function evaluations in step 1}
 \item{\code{convergence}}{convergence details from step 1}
}
}
\description{
MLE for Dirichlet Multinomial
}
\details{
Maximize Dirichlet Multinomial (DMN) log-likelihood with \code{optim()} using log likelihood function and its gradient.  This method uses a second round of optimization to estimate the scale of \eqn{\alpha} parameters, which is necessary for accurate estimation of overdispersion metric.

The covariance between counts in each category from DMN distributed data is \eqn{n(diag(p) - pp^T) (1 + \rho^2(n-1))} for \eqn{n} total counts, and vector of proportions  \eqn{p}, where \eqn{\rho^2 = 1 / (a_0 + 1)} and \eqn{a_0 = \sum_i \alpha_i}.  The count data is overdispersed by a factor of \eqn{1 + \rho^2(n-1)} compared to a multinomial (MN) distribution.  As \eqn{a_0} increases, the DMN converges to the MN.

See \url{https://en.wikipedia.org/wiki/Dirichlet-multinomial_distribution#Matrix_notation}
}
\examples{
library(HMP)

set.seed(1)

n_samples <- 1000
n_counts <- 5000
alpha <- c(500, 1000, 2000)

# Dirichlet.multinomial
counts <- Dirichlet.multinomial(rep(n_counts, n_samples), alpha)

fit <- dmn_mle(counts)

fit

# overdispersion: true value
a0 <- sum(alpha)
rhoSq <- 1 / (a0 + 1)
1 + rhoSq * (n_counts - 1)

# multinomial, so overdispersion is 1
counts <- t(rmultinom(n_samples, n_counts, prob = alpha / sum(alpha)))

dmn_mle(counts)
#
#
}
\seealso{
Other functions also estimate DMN parameters.  \code{MGLM::MGLMfit()} and \code{dirmult::dirmult()} give good parameter estimates but are slower.  \code{Rfast::dirimultinom.mle()} often fails to converge
}
