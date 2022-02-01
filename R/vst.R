# Gabriel Hoffman
# Dec 6, 2021


#' Variance stabilizing transform from precision weights
#'
#' Compute variance stabilizing transform (VST) from precision weights by scaling each observation by their respective weights
#'
#' @param x object storing data to be transformed
#'
#' @details A variance stabilizing transform is usually described in terms of a parametric model of the observed data.  Instead, here inverse variance of each observation are stored in \code{x$weight} and the VST divides the observed data by the scaled standard deviations
#'
#' @examples
#' # set probability of each category
#' prob = c(0.1, 0.2, 0.3, 0.5)
#' 
#' # number of total counts
#' countsTotal = 300
#' 
#' # number of samples
#' n_samples = 100
#' 
#' # simulate counts from multinomial
#' counts = t(rmultinom(n_samples, size = countsTotal, prob = prob))
#' colnames(counts) = paste0("cat_", 1:length(prob))
#' rownames(counts) = paste0("sample_", 1:n_samples)
#' 
#' # run crumblr on counts
#' cobj = crumblr(counts)
#' 
#' # apply variance stabilizing transform (vst)
#' df_vst = vst(cobj)
#' 
#' # Perform PCA on VST transformed data
#' pca = prcomp(t(df_vst))
#' df_pca = as.data.frame(pca$x)
#' 
#' ggplot(df_pca, aes(PC1, PC2)) + geom_point() + theme_classic() + theme(aspect.ratio=1)
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
 		# x$E * x$weights	

 		# scale by sd instead of var
 		# set mean scale to 1
 		with(x, E * (sqrt(weights)/rowMeans(sqrt(weights))))	
 	}
)

