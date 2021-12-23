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
#' 
#' # set probability of each category
#' prob = c(0.1, 0.2, 0.3, 0.5)
#' 
#' # number of total counts
#' countsTotal = 300
#' 
#' # number of samples
#' n_samples = 100
#' 
#' # simulate info for each sample
#' info = data.frame(Age = rgamma(n_samples, 50, 1))
#' rownames(info) = paste0("sample_", 1:n_samples)
#' 
#' # simulate counts from multinomial
#' counts = t(rmultinom(n_samples, size = n_samples, prob = prob))
#' colnames(counts) = paste0("cat_", 1:length(prob))
#' rownames(counts) = paste0("sample_", 1:n_samples)
#' 
#' # centered log ratio
#' clr(counts)
#'
#' @import Rdpack
#' @seealso \code{compositions::clr}
#' @export
clr = function(counts, pseudocount = 0.5){

	if( !is.data.frame(counts)  & !is.matrix(counts) ){
		counts = matrix(counts, nrow=1)
	}

	log(counts + pseudocount) - rowMeans(log(counts + pseudocount))
}

#' Count ratio uncertainty modeling based linear regression
#' 
#' Count ratio uncertainty modeling based linear regression (crumblr) returns CLR-transformed counts and observation-level inverse-variance weights for use in weighted linear models.
#' 
#' @param counts count data with samples as rows and variables are columns
#' @param pseudocount added to counts to avoid issues with zeros
#' 
#' @return  An \code{EList} object with the following components:
#' \itemize{
#'  \item{E: }{numeric matrix of CLR transformed counts}
#'  \item{weights: }{numeric matrix of observation-level inverse-variance weights}
#' }
#' 
#' @details 
#' Evalute the centered log ratio (CLR) transform of the count matrix, and the asymptotic theoretical variances of each transformed observation.  The asymptotic theoretical variances show remarkable agreement with the emprical results even for small total (i.e. row) counts, for observations of at least 2 counts.  In practice, it is often reasonable to assume a sufficient number of counts before a variable is included in an analysis anyway.  But the feasability of this assumption is up to the user to determine.     
#' 
#' Importantly, with less than 2 counts the asymptotic theory gives a *larger* variance  than the emprical results.  Therefore, this approach is conservative and will not underestimate the true amount of variation.  
#'
#' @examples
#' 
#' # set probability of each category
#' prob = c(0.1, 0.2, 0.3, 0.5)
#' 
#' # number of total counts
#' countsTotal = 300
#' 
#' # number of samples
#' n_samples = 100
#' 
#' # simulate info for each sample
#' info = data.frame(Age = rgamma(n_samples, 50, 1))
#' rownames(info) = paste0("sample_", 1:n_samples)
#' 
#' # simulate counts from multinomial
#' counts = t(rmultinom(n_samples, size = n_samples, prob = prob))
#' colnames(counts) = paste0("cat_", 1:length(prob))
#' rownames(counts) = paste0("sample_", 1:n_samples)
#' 
#' # run crumblr on counts
#' cobj = crumblr(counts)
#' 
#' # run standard limma analysis on crumblr results
#' library(limma)
#' 
#' design = model.matrix(~Age, info)
#' fit = lmFit(cobj, design)
#' fit = eBayes(fit)
#' 
#' topTable(fit, coef="Age", sort.by="none")
#' 
#' @seealso \code{limma::voom}, \code{limma::lmFit}
#' @importFrom methods new
#' @import limma
#' @export
crumblr = function(counts, pseudocount = 0.5){

	D = ncol(counts)
	
	# Compute asymptotic variance for each observation
	# var_asymp = (1/p - 2/(p*D) + sum(1/p)/D^2) / n
	var_asymp = apply(counts + pseudocount, 1, function(x){
		n = sum(x) # total counts
		p = x / n # fractions
		(1/p - 2/(p*D) + sum(1/p)/D^2) / n
		})
	
	Y_clr = t(clr(counts, pseudocount))

	# return clr transformed counts and the corresponding observation level precision
	# (i.e. inverse variances)
	new("EList", list(E = Y_clr, weights = 1/var_asymp))
}

























