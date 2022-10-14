# Gabriel Hoffman
# Dec 6, 2021


#' Studentize observations using precision weights
#'
#' Compute studentized observations from precision weights by scaling each observation by their respective standard error
#'
#' @param x object storing data to be transformed
#'
#' @details A Student transform scales the observed value but its standard error.  Here \code{x$weight} stores the inverse variance of each observation.
#'
#' @return matrix of studentized values
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
#' # apply student transform 
#' df_student = studentize(cobj)
#' 
#' # Perform PCA on student transformed data
#' pca = prcomp(t(df_student))
#' df_pca = as.data.frame(pca$x)
#' 
#' ggplot(df_pca, aes(PC1, PC2)) + geom_point() + theme_classic() + theme(aspect.ratio=1)
#' @export
#' @docType methods
#' @rdname studentize-method
setGeneric("studentize", signature="x",
	function( x )
      standardGeneric("studentize")
)

#' @export
#' @rdname studentize-method
#' @aliases studentize,EList-method
setMethod("studentize", "EList",
	function( x ){
 		
 		# W = with(x, weights / rowMeans(weights))
 		# scale by sd instead of var
 		# set mean scale to 1
 		with(x, E * (sqrt(weights)/rowMeans(sqrt(weights))))	
 	}
)














