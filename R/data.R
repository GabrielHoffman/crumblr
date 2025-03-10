#' Cell counts following interferon treatment
#'
#' Counts are from single cell RNA-seq data from treated and untreated samples from Kang, et al (2018).
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
#' @references
#'   Kang, Hyun Min, et al. "Multiplexed droplet single-cell RNA-sequencing using natural genetic variation." Nature Biotechnology 36.1 (2018): 89-94.
#' 
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
