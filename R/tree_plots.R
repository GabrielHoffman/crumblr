#' Plot tree with results from multivariate testing
#'
#' Plot tree with results from multivariate testing
#'
#' @param tree phylo object storing tree
#' @param low low color on gradient
#' @param mid mid color on gradient
#' @param high high color on gradient
#' @param xmax.scale expand the x-axis by this factor so leaf labels fit in the plot
#'
#' @examples
#' library(variancePartition)
#'
#' # Load cell counts, clustering and metadata
#' # from Kang, et al. (2018) https://doi.org/10.1038/nbt.4042
#' data(IFNCellCounts)
#'
#' # Apply crumblr transformation
#' cobj <- crumblr(df_cellCounts)
#'
#' # Use dream workflow to analyze each cell separately
#' fit <- dream(cobj, ~ StimStatus + ind, info)
#' fit <- eBayes(fit)
#'
#' # Perform multivariate test across the hierarchy
#' res <- treeTest(fit, cobj, hcl, coef = "StimStatusstim")
#'
#' # Plot hierarchy and testing results
#' plotTreeTest(res)
#'
#' # Extract results for first 3 nodes
#' res[1:3, ]
#' @import ggtree ggplot2
#' @export
plotTreeTest <- function(tree, low = "grey90", mid = "red", high = "darkred", xmax.scale = 1.5) {
  # PASS R check
  isTip <- label <- node <- FDR <- NULL

  fig <- ggtree(tree, branch.length = "none") +
    geom_tiplab(color = "black", size = 3, hjust = 0, offset = .2) +
    geom_point2(aes(label = node, color = pmin(4, -log10(FDR)), size = pmin(4, -log10(FDR)))) +
    scale_color_gradient2(name = bquote(-log[10] ~ FDR), limits = c(0, 4), low = low, mid = mid, high = high, midpoint = -log10(0.01)) +
    scale_size_area(name = bquote(-log[10] ~ FDR), limits = c(0, 4)) +
    geom_text2(aes(label = "+", subset = FDR < 0.05), color = "white", size = 6, vjust = .3, hjust = .5) +
    theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))

  # get default max value of x-axis
  xmax <- layer_scales(fig)$x$range$range[2]

  # increase x-axis width
  fig + xlim(0, xmax * xmax.scale)
}



#' Plot tree coefficients from multivariate testing
#'
#' Plot tree coefficients from multivariate testing at each node.  Only applicable top fixed effect tests
#'
#' @param tree phylo object storing tree
#' @param low low color on gradient
#' @param mid mid color on gradient
#' @param high high color on gradient
#' @param xmax.scale expand the x-axis by this factor so leaf labels fit in the plot
#'
#' @examples
#' library(variancePartition)
#'
#' # Load cell counts, clustering and metadata
#' # from Kang, et al. (2018) https://doi.org/10.1038/nbt.4042
#' data(IFNCellCounts)
#'
#' # Apply crumblr transformation
#' cobj <- crumblr(df_cellCounts)
#'
#' # Use dream workflow to analyze each cell separately
#' fit <- dream(cobj, ~ StimStatus + ind, info)
#' fit <- eBayes(fit)
#'
#' # Perform multivariate test across the hierarchy
#' res <- treeTest(fit, cobj, hcl, coef = "StimStatusstim")
#'
#' # Plot hierarchy, no tests are significant
#' plotTreeTestBeta(res)
#' @import ggtree ggplot2
#' @export
plotTreeTestBeta <- function(tree, low = "blue", mid = "white", high = "red", xmax.scale = 1.5) {
  # PASS R check
  isTip <- label <- node <- FDR <- NULL

  # comparison only works for fixed effects models
  if (!all(tree@data$method %in% c("FE", "FE.empirical"))) {
    stop("tree1 must be evaluated with a fixed effect model")
  }

  beta_max <- tree %>%
    as_tibble() %>%
    pull(beta) %>%
    abs() %>%
    max()

  fig <- ggtree(tree, branch.length = "none") +
    geom_tiplab(color = "black", size = 3, hjust = 0, offset = .2) +
    geom_point2(aes(label = node, color = beta, size = pmin(4, -log10(FDR)))) +
    scale_color_gradient2(name = bquote(beta), low = low, mid = mid, high = high, midpoint = 0, limits = c(-beta_max, beta_max)) +
    scale_size_area(name = bquote(-log[10] ~ FDR), limits = c(0, 4)) +
    geom_text2(aes(label = "+", subset = FDR < 0.05), color = "white", size = 6, vjust = .3, hjust = .5) +
    theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))

  # get default max value of x-axis
  xmax <- layer_scales(fig)$x$range$range[2]

  # increase x-axis width
  fig + xlim(0, xmax * xmax.scale)
}

#' Forest plot
#'
#' Forest plot
#'
#' @param x object to be plotted
#' @param ... other arguments
#'
#' @rdname plotForest-methods
#' @export
setGeneric("plotForest", function(x, ...) {
  standardGeneric("plotForest")
})




.plotForest <- function(tree, low = "blue", mid = "grey70", high = "red", hide = FALSE) {
  # PASS R check
  n_features <- label <- se <- NA

  # get results at nodes
  tab <- tree %>%
    as_tibble() %>%
    as_tibble() %>%
    filter(n_features == 1)

  # get plot of tree for node order
  fig.tree <- plotTreeTestBeta(tree)

  # Same order for tree and forest plot
  lvls <- rev(get_taxa_name(fig.tree))
  tab$label <- factor(tab$label, lvls)

  beta_max <- max(abs(tab$beta))

  fig <- ggplot(tab, aes(label, beta)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", linewidth = 1) +
    geom_errorbar(aes(ymin = beta - 1.96 * se, ymax = beta + 1.96 * se), width = 0) +
    geom_point(aes(color = beta)) +
    theme_classic() +
    coord_flip() +
    xlab("") +
    ylab("Effect size") +
    scale_color_gradient2(name = bquote(beta), low = low, mid = mid, high = high, midpoint = 0, limits = c(-beta_max, beta_max))

  if (hide) {
    fig <- fig +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "bottom"
      )
  }
  fig
}



#' Forest plot of effect size estimates
#'
#' Forest plot of effect size estimates at the leaves of the tree
#'
#' @param x result from \code{treeTest()}
#' @param ... other arguments
#' @param hide hide rownames and legend
#'
#' @examples
#' library(variancePartition)
#'
#' # Load cell counts, clustering and metadata
#' # from Kang, et al. (2018) https://doi.org/10.1038/nbt.4042
#' data(IFNCellCounts)
#'
#' # Apply crumblr transformation
#' cobj <- crumblr(df_cellCounts)
#'
#' # Use dream workflow to analyze each cell separately
#' fit <- dream(cobj, ~ StimStatus + ind, info)
#' fit <- eBayes(fit)
#'
#' # Perform multivariate test across the hierarchy
#' res <- treeTest(fit, cobj, hcl, coef = "StimStatusstim")
#'
#' # Plot log fold changes from coef
#' plotForest(res)
#'
#' @rdname plotForest-methods
#' @aliases plotForest,#' -method
#' @export
setMethod(
  "plotForest", signature(x = "treedata"),
  function(x, ..., hide = FALSE) {
    .plotForest(x, ..., hide = hide)
  }
)
