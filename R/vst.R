# Gabriel Hoffman
# Dec 6, 2021


#' Variance stabilizing transform from precision weights
#'
#' Compute variance stabilizing transform (VST) from precision weights by scaling each observation by their respective weights
#'
#' @param x object storing data to be transformed
#'
#' @details A variance stabilizing transform is usually described in terms of a parametric model of the observed data.  Instead, here inverse variance of each observation are stored in \code{x$weight} and the VST divides the observed data by the variances's.
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
 		x$E / x$weights	
 	}
)

