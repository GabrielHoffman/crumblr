# Gabriel Hoffman
# Dec 6, 2021


#' Variance stabilizing transform from precision weights
#'
#' Compute variance stabilizing transform (VST) from precision weights by scaling each observation by their respective weights
#'
#' @param x object storing data to be transformed
#'
#' @details A variance stabilizing transform is usually described in terms of a parametric model of the observed data.  Instead, here inverse variance of each observation are stored in \code{x$weight} and the VST divides the observed data by the variances.
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
#' fit = lmFit(vst(cobj), design)
#' fit = eBayes(fit)
#' 
#' topTable(fit, coef="Age", sort.by="none")
#' 
#' @export
#' @docType methods
#' @rdname vst-method
setGeneric("vst", signature="x",
	function( x )
      standardGeneric("vst")
)

#' @export
#' @rdname vst-method
#' @aliases vst,EList-method
setMethod("vst", "EList",
	function( x ){
 		x$E * x$weights	
 	}
)

