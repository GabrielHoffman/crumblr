#' Cell counts following interferon treatment
#'
#' Counts are from single cell RNA-seq data from treated and untreated samples from \insertCite{kang2018multiplexed}{crumblr}.
#'
#' @docType data
#'
#' @usage data(IFNCellCounts)
#'
#' @format
#' \itemize{
#' 	\item \code{info} is metadata for each sample
#' 	\item \code{df_cellCounts} data.frame of counts for each sample
#' 	\item \code{hcl} cluster of cell types based on pseudobulk expression
#' }
#'
#' @keywords datasets
#'
#' @references{
#'   \insertRef{kang2018multiplexed}{crumblr}
#' }
#' @name IFNCellCounts
NULL


#' @docType data
#' @name IFNCellCounts
#' @keywords datasets
"info"

#' @docType data
#' @name IFNCellCounts
#' @keywords datasets
"df_cellCounts"

#' @docType data
#' @name IFNCellCounts
#' @keywords datasets
"hcl"
