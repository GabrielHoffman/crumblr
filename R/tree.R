# Gabriel Hoffman
# October 7, 2022
#
# Perform mvTest() on a hierarchical tree



# For each node and tip return the set of offspring tips

#' @import tidytree
getTips <- function(hc) {
  # convert to table
  tb <- as_tibble(as.phylo(hc))

  tipLabels <- hc$labels

  # for each node, return the tips below it
  res <- lapply(tb$node, function(id) {
    # given node, identify offspring
    nodes <- offspring(tb, id)$node

    # if there are no offspring
    if (length(nodes) == 0) {
      # return the tip label
      lab <- tb$label[match(id, tb$node)]
    } else {
      # get labels of offspring
      labels <- tb$label[match(nodes, tb$node)]

      # kepe only tip labels
      lab <- labels[labels %in% tipLabels]
    }
    lab
  })
  names(res) <- tb$node
  res
}

# Assign names to internal nodes by concatenating
# names from the tips
NameInternalNodes <- function(tb) {
  # for each internal node
  # assign name based on tips
  df <- lapply(tb$node[is.na(tb$label)], function(id) {
    labs <- offspring(tb, id)$label
    txt <- paste0(labs[!is.na(labs)], collapse = "/")

    data.frame(node = id, label = txt)
  })
  df <- do.call(rbind, df)

  # assign new labels
  i <- match(df$node, tb$node)
  tb$label[i] <- df$label

  tb
}





#' Perform multivariate testing along a hierarchy
#'
#' Perform multivariate testing using \code{mvTest()} along the nodes of tree
#'
#' @param fit \code{MArrayLM} object return by \code{lmFit()} or \code{dream()}
#' @param obj \code{EList} object returned by \code{voom()}
#' @param hc hierarchical clustering as an \code{hclust} object
#' @param coef name of coefficient to be extracted
#' @param method statistical method used to perform multivariate test.  See details.  \code{'FE'} is a fixed effect test that models the covariance between coefficients. \code{'FE.empirical'} use compute empirical p-values by sampling from the null distribution and fitting with a gamma. \code{'RE2C'} is a random effect test of heterogeneity of the estimated coefficients that models the covariance between coefficients, and also incorporates a fixed effects test too. \code{'tstat'} combines the t-statistics and models the covariance between coefficients. \code{'sidak'} returns the smallest p-value and accounting for the number of tests. \code{'fisher'} combines the p-value using Fisher's method assuming independent tests.
#' @param shrink.cov shrink the covariance matrix between coefficients using the Schafer-Strimmer method
#'
#' @details See package \code{remaCor} for details about the \code{remaCor::RE2C()} test, and see \code{remaCor::LS()} for details about the fixed effect test.  When only 1 feature is selected, the original t-statistic and p-value are returned.
#'
#' @seealso \code{variancePartition::mvTest}
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
#' @importFrom tidytree as_tibble left_join as.treedata
#' @importFrom variancePartition mvTest
#' @importFrom stats p.adjust
#' @importFrom dplyr bind_rows
#' @export
treeTest <- function(fit, obj, hc, coef, method = c("FE.empirical", "FE", "RE2C", "tstat", "sidak", "fisher"), shrink.cov = TRUE) {
  method <- match.arg(method)

  # get tips for each node
  testSets <- getTips(hc)

  uniqTips <- unique(unlist(testSets))

  if (!identical(sort(uniqTips), sort(rownames(fit)))) {
    stop("Tree and fit have non-shared labels")
  }

  # for each node
  res <- lapply(seq(length(testSets)), function(i) {
    labels <- testSets[[i]]

    # peform a multivariate test based on labels
    res <- mvTest(fit, obj, labels, coef, method, shrink.cov)

    tibble(node = as.integer(names(testSets)[i]), res)
  })
  res <- bind_rows(res)

  # convert to table
  tb <- as_tibble(as.phylo(hc))

  # join to combine with test results
  tb <- left_join(tb, res, by = "node")

  # assign names to internal nodes based on tips
  tb <- NameInternalNodes(tb)

  tb$FDR <- p.adjust(tb$pvalue, "BH")

  as.treedata(tb)
}

#' Perform hierarchical clustering on reducedDim
#'
#' Perform hierarchical clustering dimension reduction from single cell expression data
#'
#' @param sce \code{SingleCellExperiment} object
#' @param reduction field of reducedDims(sce) to use
#' @param labelCol column in \code{SingleCellExperiment} storing cell type annotations
#' @param method.dist method for \code{dist(..,method=method.dist)}
#' @param method.hclust method for \code{hclust(..,method=method.hclust)}
#'
#' @importFrom SingleCellExperiment reducedDim colData
#' @importFrom stats dist hclust
#' @export
buildClusterTree <- function(sce, reduction, labelCol, method.dist = c("cosine", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"), method.hclust = c("complete", "ward.D", "ward.D2")) {
  method.dist <- match.arg(method.dist)
  method.hclust <- match.arg(method.hclust)

  if (!labelCol %in% colnames(colData(sce))) {
    stop("labelCol must be a column in colData(sce)")
  }

  # extract low dimensional embeddings
  embeddings <- reducedDim(sce, reduction)

  # extract unique cell labels
  lvls <- colData(sce)[[labelCol]]
  if (is.factor(lvls)) {
    lvls <- levels(lvls)
  } else {
    lvls <- unique(lvls)
  }

  # for each cell type, return mean of each dimension
  df <- lapply(lvls, function(x) {
    idx <- colData(sce)[[labelCol]] == x

    colMeans(embeddings[idx, , drop = FALSE])
  })
  df <- do.call(rbind, df)
  rownames(df) <- lvls

  # compute pairwise distances
  if (method.dist == "cosine") {
    dst <- dist.cosine(df)
  } else {
    dst <- dist(df, method = method.dist)
  }

  # perform clustering
  hclust(dst, method = method.hclust)
}


#' Compare difference in estimates between two trees
#'
#' Compare difference in cofficient estimates between two trees.  For node \code{i}, the test evaluates \code{tree1[i] - tree2[i] = 0}.
#'
#' @param tree1 object of type \code{treedata} from \code{treeTest()}
#' @param tree2 object of type \code{treedata} from \code{treeTest()}
#'
#' @details When a fixed effect test is performed at each node using \code{treeTest()} with \code{method = "FE.empirical"} or \code{method = "FE"}, a coefficient estimate and standard error are estimated for each node based on the children.  This function performs a two-sample z-test to test if a given coefficient from \code{tree1} is significantly different from the corresponding coefficient in \code{tree2}.
#'
#' @examples
#' library(variancePartition)
#'
#' # Load cell counts, clustering and metadata
#' # from Kang, et al. (2018) https://doi.org/10.1038/nbt.4042
#' data(IFNCellCounts)
#'
#' # Simulate a factor with 2 levels called DiseaseRand
#' set.seed(123)
#' info$DiseaseRand <- sample(LETTERS[seq(2)], nrow(info), replace = TRUE)
#' info$DiseaseRand <- factor(info$DiseaseRand, LETTERS[seq(2)])
#'
#' # Apply crumblr transformation
#' cobj <- crumblr(df_cellCounts)
#'
#' # Use dream workflow to analyze each cell separately
#' fit <- dream(cobj, ~ StimStatus + ind, info)
#' fit <- eBayes(fit)
#'
#' # Perform multivariate test across the hierarchy
#' res1 <- treeTest(fit, cobj, hcl, coef = "StimStatusstim")
#'
#' # Perform same test, but on DiseaseRand
#' fit2 <- dream(cobj, ~DiseaseRand, info)
#' fit2 <- eBayes(fit2)
#' res2 <- treeTest(fit2, cobj, hcl, coef = "DiseaseRandB")
#'
#' # Compare the coefficient estimates at each node
#' # Test if res1 - res2 is significantly different from zero
#' resDiff <- diffTree(res1, res2)
#'
#' resDiff
#'
#' plotTreeTest(resDiff)
#'
#' plotTreeTestBeta(resDiff)
#' @importFrom dplyr inner_join
#' @export
diffTree <- function(tree1, tree2) {
  df1 <- as_tibble(tree1)
  df2 <- as_tibble(tree2)

  # comparison only works for fixed effects models
  if (!all(df1$method %in% c("FE", "FE.empirical"))) {
    stop("tree1 must be evaluated with a fixed effect model")
  }

  if (!all(df2$method %in% c("FE", "FE.empirical"))) {
    stop("tree2 must be evaluated with a fixed effect model")
  }

  # join data.frames from each tree
  df_merge <- inner_join(df1, df2, by = c("parent", "node", "label"))

  # initialize
  df_out <- df1

  # perform two-sample z-test
  df_out$beta <- with(df_merge, beta.x - beta.y)
  df_out$se <- with(df_merge, sqrt(se.x^2 + se.y^2))
  df_out$stat <- with(df_out, beta / se)
  df_out$pvalue <- with(df_out, 2 * pnorm(abs(stat), lower.tail = FALSE))
  df_out$FDR <- with(df_out, p.adjust(pvalue, "BH"))

  as.treedata(df_out)
}
