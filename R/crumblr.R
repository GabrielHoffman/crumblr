# Gabriel Hoffman
# October 26, 2021

#' Centered log ratio transform
#'
#' Compute the centered log ratio (CLR) transform of a count matrix.
#'
#' @param counts count data with samples as rows and variables are columns
#' @param pseudocount added to counts to avoid issues with zeros
#'
#' @details The CLR of a vector \code{x} of counts in \code{D} categories is defined as
#' \code{clr(x) = log(x) - mean(log(x))}. For details see van den Boogaart and  Tolosana-Delgado (2013).
#'
#' @references{
#'   \insertRef{van2013analyzing}{crumblr}
#' }
#'
#' @return matrix of CLR transformed counts
#'
#' @examples
#' # set probability of each category
#' prob <- c(0.1, 0.2, 0.3, 0.5)
#'
#' # number of total counts
#' countsTotal <- 300
#'
#' # number of samples
#' n_samples <- 100
#'
#' # simulate info for each sample
#' info <- data.frame(Age = rgamma(n_samples, 50, 1))
#' rownames(info) <- paste0("sample_", 1:n_samples)
#'
#' # simulate counts from multinomial
#' counts <- t(rmultinom(n_samples, size = n_samples, prob = prob))
#' colnames(counts) <- paste0("cat_", 1:length(prob))
#' rownames(counts) <- paste0("sample_", 1:n_samples)
#'
#' # centered log ratio
#' clr(counts)[1:4, ]
#'
#' @import Rdpack
#' @seealso \code{compositions::clr}
#' @export
clr <- function(counts, pseudocount = 0.5) {
  if (!is.data.frame(counts) & !is.matrix(counts)) {
    counts <- matrix(counts, nrow = 1)
  }

  log(counts + pseudocount) - rowMeans(log(counts + pseudocount))
}


#' Count ratio uncertainty modeling based linear regression
#'
#' Count ratio uncertainty modeling based linear regression (crumblr) returns CLR-transformed counts and observation-level inverse-variance weights for use in weighted linear models.
#'
#' @param counts count data with samples as rows and variables are columns
#' @param pseudocount added to counts to avoid issues with zeros
#' @param method \code{"clr"} computes standard centered log ratio and precision weights based on the delta approximation. \code{"clr_2class"} computes the \code{clr()} transform for category \code{i} using 2 classes: 1) counts in category i, and 2) counts _not_ in category i.
#' @param tau overdispersion parameter for Dirichlet multinomial.  If \code{NULL}, estimate from observed counts.
#' @param max.ratio regularize estimates of the weights to have a maximum ratio of \code{max.ratio} between the maximum and \code{quant} quantile value
#' @param quant quantile value used for \code{max.ratio}
#'
#' @return  An \code{EList} object with the following components:
#' \describe{
#'  \item{E: }{numeric matrix of CLR transformed counts}
#'  \item{weights: }{numeric matrix of observation-level inverse-variance weights}
#' }
#'
#' @details
#' Evalute the centered log ratio (CLR) transform of the count matrix, and the asymptotic theoretical variances of each transformed observation.  The asymptotic normal approximation is increasingly accurate for small overdispersion \eqn{\tau}, large total counts \eqn{C}, and large proportions \eqn{p}, but shows good agreement with the empirical results in most situtations. In practice, it is often reasonable to assume a sufficient number of counts before a variable is included in an analysis anyway.  But the feasability of this assumption is up to the user to determine.
#'
#' Given the array \code{p} storing proportions for one sample across all categories, the delta approximation uses the term \code{1/p}.  This can be unstable for small values of \code{p}, and the estimated variances can be sensitive to small changes in the proprtions.  To address this, the \code{"clr_2class"} method computes the \code{clr()} transform for category \code{i} using 2 classes: 1) counts in category i, and 2) counts _not_ in category i. Since class (2) now sums counts across all other categories, the small proportions are avoided and the variance estimates are more stable.
#' 
#' For real data, the asymptotic variance formula can give weights that vary substantially across samples and give very high weights for a subset of samples.  In order to address this, we regularize the weights to reduce the variation in the weights to have a maximum ratio of \code{max.ratio} between the maximum and \code{quant} quantile value.  
#' 
#' @examples
#' # set probability of each category
#' prob <- c(0.1, 0.2, 0.3, 0.5)
#'
#' # number of total counts
#' countsTotal <- 300
#'
#' # number of samples
#' n_samples <- 100
#'
#' # simulate info for each sample
#' info <- data.frame(Age = rgamma(n_samples, 50, 1))
#' rownames(info) <- paste0("sample_", 1:n_samples)
#'
#' # simulate counts from multinomial
#' counts <- t(rmultinom(n_samples, size = countsTotal, prob = prob))
#' colnames(counts) <- paste0("cat_", 1:length(prob))
#' rownames(counts) <- paste0("sample_", 1:n_samples)
#'
#' # run crumblr on counts
#' cobj <- crumblr(counts)
#'
#' # run standard variancePartition analysis on crumblr results
#' library(variancePartition)
#'
#' fit <- dream(cobj, ~Age, info)
#' fit <- eBayes(fit)
#'
#' topTable(fit, coef = "Age", sort.by = "none")
#' @seealso \code{limma::voom}, \code{variancePartition::dream}
#' @export
setGeneric(
  "crumblr",
  function(counts, pseudocount = 0.5, method = c("clr", "clr_2class"), tau = 1, max.ratio=5, quant=0.05) {
    standardGeneric("crumblr")
  }
)



#' @rdname crumblr
#' @aliases crumblr,matrix-method
#' @importFrom methods new
setMethod(
  "crumblr", "matrix",
  function(counts, pseudocount = 0.5, method = c("clr", "clr_2class"), tau = 1, max.ratio=5, quant=0.05) {
    method <- match.arg(method)

    if (method == "clr") {
      res <- .crumblr(counts, pseudocount, tau, max.ratio, quant)
    } else {
      res <- .clr_2class(counts, pseudocount)
    }
    res
  }
)

#' @rdname crumblr
#' @aliases crumblr,data.frame-method
#' @importFrom methods new
setMethod(
  "crumblr", "data.frame",
  function(counts, pseudocount = 0.5, method = c("clr", "clr_2class"), tau = 1, max.ratio=5, quant=0.05) {
    method <- match.arg(method)

    if (method == "clr") {
      res <- .crumblr(as.matrix(counts), pseudocount, tau, max.ratio, quant)
    } else {
      res <- .clr_2class(as.matrix(counts), pseudocount, tau, max.ratio, quant)
    }
    res
  }
)

# @param max.ratio regularize estimates of the weights to have a maximum ratio of \code{max.ratio} between the maximum and \code{quant} qauntile value
# @param quant qauntile value used for \code{max.ratio}
#' @importFrom stats quantile
.cap_ratio = function(W, max.ratio=5, quant=0.05){
  t(apply(W, 1, function(x){
    x = x / quantile(x, quant)
    pmin(x,max.ratio)
    }))
}

# @param max.ratio regularize estimates of the weights to have a maximum ratio of \code{max.ratio} between the maximum and \code{quant} qauntile value
# @param quant qauntile value used for \code{max.ratio}
.crumblr <- function(counts, pseudocount = 0.5, tau = 1, max.ratio=5, quant=0.05) {
  D <- ncol(counts)

  # estimate overdispersion from observed counts
  if (is.null(tau)) {
    tau <- dmn.mle(counts)$overdispersion
  }

  # Compute asymptotic variance for each observation
  # var_asymp = tau * (1/p - 2/(p*D) + sum(1/p)/D^2) / C
  var_asymp <- apply(counts + pseudocount, 1, function(x) {
    C <- sum(x) # total counts
    p <- x / C # fractions
    tau * (1 / p - 2 / (p * D) + sum(1 / p) / D^2) / C
  })

  # regularise the variance estimates
  weights = 1 / var_asymp
  weights = .cap_ratio( weights, max.ratio, quant)

  # get CLR values
  Y_clr <- t(clr(counts, pseudocount))

  # return clr transformed counts and the corresponding observation level precision
  # (i.e. inverse variances)
  new("EList", list(E = Y_clr, 
                    weights = weights))
}


# for each category, perform clr on two class coded as
# 1) counts from category i
# 2) sum of all other counts
# This avoids the 1/p term for small values of p in the
# delta approximation for rare cell types
.clr_2class <- function(counts, pseudocount = 0.5, tau = 1, max.ratio=5, quant=0.05) {
  # sunm counts for each row
  rs <- rowSums(counts)

  # for each category
  res <- lapply(colnames(counts), function(x) {
    # create a 2 category count data.frame:
    # 1) counts of target
    # 2) sum of all other counts
    df <- data.frame(
      positive = counts[, x],
      other = rs - counts[, x]
    )
    colnames(df)[1] <- x

    # run clr with crumblr weights on 2 class counts
    obj <- .crumblr(df, pseudocount, tau, max.ratio, quant)

    # extract positive class
    obj[1, ]
  })

  # combine results
  cobj <- do.call(rbind, res)

  cobj
}
