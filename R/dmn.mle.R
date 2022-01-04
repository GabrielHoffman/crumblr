# Gabriel Hoffman
# Dec 26, 2021

#' MLE for Dirichlet Multinomial
#'
#' MLE for Dirichlet Multinomial
#' 
#' @param counts matrix with rows as samples and columns as categories
#' @param ... additional arguments passed to \code{optim()} 
#' 
#' @return list storing \code{alpha} parameter estimates, logLik, and details about convergence
#' \itemize{
#'  \item{\code{alpha}}{estimated \eqn{alpha} parameters}
#'  \item{\code{overdispersion}}{Overdispersion value \eqn{1 + \rho^2(n-1)} compared to multinomial}
#'  \item{\code{logLik}}{value of function}
#'  \item{\code{scale}}{scaling of \eqn{\alpha} parameters computed in a second optimization step}
#'  \item{\code{evals}}{number of function evaluations in step 1}
#'  \item{\code{convergence}}{convergence details from step 1}
#' }
#' 
#' @examples
#' 
#' library(HMP)
#' 
#' set.seed(1)
#' 
#' n_samples = 1000
#' n_counts = 5000
#' alpha = c(500, 1000, 2000)
#' 
#' # Dirichlet.multinomial
#' counts = Dirichlet.multinomial(rep(n_counts, n_samples), alpha)
#' 
#' fit = dmn.mle(counts)
#' 
#' fit
#' 
#' # overdispersion: true value
#' a0 = sum(alpha)
#' rhoSq = 1 / (a0 + 1)
#' 1 + rhoSq*(n_counts-1) 
#' 
#' # multinomial
#' counts = t(rmultinom(n_samples, n_counts, prob=alpha / sum(alpha)))
#'
#' dmn.mle(counts)
#'
# # where n is number of counts
# n = rowsums(counts)[1]
# a0 = sum(alpha)
# p = alpha / a0
#
# # a0 = (1-rhoSq)/rhoSq
# rhoSq = 1 / (a0 + 1)
#
# colMeans(counts)
# n *(alpha/a0)
#
# cov(counts)
# n * (diag(p) - tcrossprod(p)) *                                    
#'
#' @details Maximize Dirichlet Multinomial (DMN) log-likelihood with \code{optim()} using log likelihood function and its gradient.  This method sses a second round of optimization to estimate the scale of \eqn{\alpha} parameters, which is necessary for accurate estimation of overdispersion metric.
#' 
#' The covariance between counts in each category from DMN distributed data is \eqn{n(diag(p) - pp^T) (1 + \rho^2(n-1))} for \eqn{n} total counts, and vector of proportions  \eqn{p}, where \eqn{\rho^2 = 1 / (a_0 + 1)} and \eqn{a_0 = \sum_i \alpha_i}.  The count data is overdispersed by a factor of \eqn{1 + \rho^2(n-1)} compared to a multinomial (MN) distribution.  As \eqn{a_0} increases, the DMN converges to the MN.
#' 
#' See \url{https://en.wikipedia.org/wiki/Dirichlet-multinomial_distribution#Matrix_notation}
#' 
#' @seealso Other functions also estimate DMN parameters.  \code{MGLM::MGLMfit()} and \code{dirmult::dirmult()} give good parameter estimates but are slower.  \code{Rfast::dirimultinom.mle()} often fails to converge
#' @importFrom Rfast Lgamma Digamma colmeans rowsums
#' @importFrom stats optim optimize 
#' @export
dmn.mle = function(counts,...){

  # First round
  ##############

  # negative log-likelihood
  # adapted from https://github.com/RfastOfficial/Rfast/blob/1f14fbb7a0729b6c4122fc3ce12701ad4801f1a3/R/discrete_mle.R#L157
  f = function(a1, x, rs){
    # a1 = exp(a1)
    y <- x + a1
    sa <- sum(a1)
    value = n * Lgamma( sa ) - sum( Lgamma( rs + sa ) ) - n * sum( Lgamma( a1 ) ) + sum( Lgamma( y ) )
    -1*value
  }

  # gradient
  # Note that including the gradient in optim() substantially increases speed
  gr = function(a1, x, rs){
    # a1 = exp(a1)
    y <- x + a1
    sa <- sum(a1)
    value = n * Digamma(sa) - sum( Digamma(rs + sa) ) - n * Digamma(a1) + rowsums( Digamma(y) )
    -1*value
  }

  # pre-compute values
  n <- nrow(counts)  ## sample size
  rs <- rowsums(counts)
  # init <- log(colmeans(counts) )
  init <- colmeans(counts) 
  
  # Maximize likelihood
  # # use parameters in log-space to avoid constraint at zero
  fit = optim(init, f, gr=gr, x=t(counts), rs=rs,..., method="L-BFGS-B", lower = 1e-04)

  # convergence check
  if( fit$convergence != 0 ){
    warning("Optimization did not converge: ", fit$message)
  }

  # Second round
  ##############

  # the first round does a goode job at estimating relative values of alpha
  # but has issues estimating the scale of alpha for large values
  # Here, keep relative values constant and do 1D optimization of scale

  f2 = function(s, a1, x, rs){ 
    a1 = a1*s
    y <- x + a1
    sa <- sum(a1)
    value = n * Lgamma( sa ) - sum( Lgamma( rs + sa ) ) - n * sum( Lgamma( a1 ) ) + sum( Lgamma( y ) )
    -1*value
  }
  
  x = seq(-1, 5, length.out=10)
  df = lapply( 2:length(x), function(i){
    interval = 10^c(x[i-1], x[i])
    fit2 = optimize(f2, interval=interval, a1 = fit$par, x=t(counts), rs=rs  )
    data.frame(fit2)
  })
  df = do.call(rbind, df)

  # get best scale
  s = df[which.min(df$objective),'minimum']

  # set parameter names and apply scale
  # alpha = exp(fit$par)
  alpha = fit$par * s
  names(alpha) = colnames(counts)

  # DMN covariance is 
  # n * (diag(p) - tcrossprod(p)) * (1 + rhoSq*(n-1))
  # so overdispersion is measured as 1 + rhoSq*(n-1)
  # based on https://en.wikipedia.org/wiki/Dirichlet-multinomial_distribution#Matrix_notation

  # a0 = (1-rhoSq)/rhoSq
  a0 = sum(alpha)
  rhoSq = 1 / (a0 + 1)
  overdispersion = 1 + rhoSq*(mean(rs)-1) 

  list( alpha       = alpha, 
        overdispersion = overdispersion,    
        logLik      = -1*fit$value, 
        scale       = s,
        evals       = fit$counts, 
        convergence = fit$convergence)
}












