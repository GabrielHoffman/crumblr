# Gabriel Hoffman
# October 7, 2022
#
# Perform mvTest() on a hierarchical tree



# For each node and tip return the set of offspring tips

#' @import tidytree 
getTips = function(hc){

	# convert to table
	tb = as_tibble(as.phylo(hc))

	tipLabels = hc$labels

	# for each node, return the tips below it
	res = lapply(tb$node, function(id){
		# given node, identify offspring
		nodes = offspring(tb, id)$node

		# if there are no offspring
		if( length(nodes) == 0){
			# return the tip label
			lab = tb$label[match(id, tb$node)]
		}else{
			# get labels of offspring
			labels = tb$label[match(nodes, tb$node)]

			# kepe only tip labels
			lab = labels[labels %in% tipLabels]
		}
		lab
	})
	names(res) = tb$node
	res
}

# Assign names to internal nodes by concatenating 
# names from the tips
NameInternalNodes = function(tb){

	# for each internal node
	# assign name based on tips
	df = lapply(tb$node[is.na(tb$label)], function(id){

		labs = offspring(tb, id)$label
		txt = paste0(labs[!is.na(labs)], collapse='/')

		data.frame(node = id, label = txt)
	})
	df = do.call(rbind, df)

	# assign new labels
	i = match(df$node, tb$node)
	tb$label[i] = df$label

	tb
}


#' Plot tree with results from multivariate testing
#' 
#' Plot tree with results from multivariate testing
#' 
#' @param tree phylo object storing tree
#' @param low low color on gradient
#' @param mid mid color on gradient
#' @param high high color on gradient
#' 
#' @examples
#' library(variancePartition)
#' 
#' # Load cell counts from Kang, et al. (2018)
#' #  https://doi.org/10.1038/nbt.4042
#' data(IFNCellCounts)
#' 
#' # Apply crumblr transformation 
#' cobj = crumblr(cellCounts)
#' 
#' # Use dream workflow to analyze each cell separately
#' fit = dream(cobj, ~ StimStatus + ind, info)
#' fit = eBayes(fit)
#' 
#' # Create a hierarchical cluster of cell types
#' # NOTE: for example only
#' # Create clustering using prior knowledge 
#' # or single cell data
#' hc = hclust(dist(t(cellCounts)))
#' 
#' # Perform multivariate test across the hierarchy
#' res = treeTest( fit, cobj, hc, coef="StimStatusstim", method="RE2C")
#' 
#' # Plot hierarchy and testing results
#' # Adjust xlim() until text fits in window 
#' plotTreeTest(res) + xlim(0, 7)
#' 
#' @import ggtree ggplot2
#' @export
plotTreeTest = function(tree, low="grey90", mid = "red", high="darkred"){

	# PASS R check
	isTip = label = node = FDR = NULL

	ggtree(tree, branch.length = "none") +
	    # geom_text2(aes(label = paste0('     ', label), subset=isTip), color = "black", size=3, hjust=0) + 
	    geom_tiplab(color = "black", size=3, hjust=0, offset=.2) +
	    geom_point2(aes(label = node, color=pmin(4,-log10(FDR)), size=pmin(4,-log10(FDR)))) + 
	    scale_color_gradient2(name = bquote(-log[10]~FDR), limits=c(0,4), low=low, mid=mid, high=high, midpoint=-log10(0.01)) +
	    scale_size_area(name = bquote(-log[10]~FDR), limits=c(0,4)) +
	    geom_text2(aes(label = '+', subset=FDR < 0.05), color = "white", size=6, vjust=.3, hjust=.5) +
	    theme(legend.position="bottom", aspect.ratio=1)
}


#' Perform multivariate testing along a hierarchy
#' 
#' Perform multivariate testing using \code{mvTest()} along the nodes of tree
#' 
#' @param fit \code{MArrayLM} object return by \code{lmFit()} or \code{dream()}
#' @param obj \code{EList} object returned by \code{voom()}
#' @param hc hierarchical clustering as an \code{hclust} object
#' @param coef name of coefficient to be extracted
#' @param method statistical method used to perform multivariate test.  See details. \code{'RE2C'} is a random effect test of heterogeneity of the estimated coefficients that models the covariance between coefficients, and also incorporates a fixed effects test too. \code{'FE'} is a fixed effect test that models the covariance between coefficients.  \code{'tstat'} combines the t-statistics and models the covariance between coefficients. \code{'sidak'} returns the smallest p-value and accounting for the number of tests. \code{'fisher'} combines the p-value using Fisher's method assuming independent tests.
#'  
#' @details See package \code{remaCor} for details about the \code{remaCor::RE2C()} test, and see \code{remaCor::LS()} for details about the fixed effect test.  When only 1 feature is selected, the original t-statistic and p-value are returned.
#'
#' @seealso \code{variancePartition::mvTest}
#' @examples
#' library(variancePartition)
#' 
#' # Load cell counts from Kang, et al. (2018)
#' #  https://doi.org/10.1038/nbt.4042
#' data(IFNCellCounts)
#' 
#' # Apply crumblr transformation 
#' cobj = crumblr(cellCounts)
#' 
#' # Use dream workflow to analyze each cell separately
#' fit = dream(cobj, ~ StimStatus + ind, info)
#' fit = eBayes(fit)
#' 
#' # Create a hierarchical cluster of cell types
#' # NOTE: for example only
#' # Create clustering using prior knowledge 
#' # or single cell data
#' hc = hclust(dist(t(cellCounts)))
#' 
#' # Perform multivariate test across the hierarchy
#' res = treeTest( fit, cobj, hc, coef="StimStatusstim", method="RE2C")
#' 
#' # Plot hierarchy and testing results
#' # Adjust xlim() until text fits in window 
#' plotTreeTest(res) + xlim(0, 7)
#' 
#' @importFrom tidytree as_tibble left_join as.treedata
#' @importFrom variancePartition mvTest
#' @importFrom stats p.adjust
#' @export
treeTest = function(fit, obj, hc, coef, method = c("RE2C", "FE", "tstat", "sidak", "fisher")){

	# get tips for each node
	testSets = getTips(hc)

	uniqTips = unique(unlist(testSets))

	if( ! identical(sort(uniqTips), sort(rownames(fit))) ){
		stop("Tree and fit have non-shared labels")
	}

	# for each node
	res = lapply(seq(length(testSets)), function(i){

		labels = testSets[[i]]

		# peform a multivariate test based on labels
		res = mvTest(fit, obj, labels, coef, method)

		tibble(node = as.integer(names(testSets)[i]), res)
	})
	res = do.call(rbind, res)

	# convert to table
	tb = as_tibble(as.phylo(hc))

	# join to combine with test results 
	tb = left_join(tb, res, by="node")

	# assign names to internal nodes based on tips
	tb = NameInternalNodes(tb)

	tb$FDR = p.adjust(tb$pvalue, "BH")

	as.treedata(tb)
}

#' Perform hierarchical clustering on reducedDim
#'
#' Perform hierarchical clustering dimension reduction from single cell expression data 
#' 
#' @param sce \code{SingleCellExperiment} object
#' @param reduction field of reducedDims(sce) to use
#' @param labelCol column in \code{SingleCellExperiment} storing cell type annotations
#' 
#' @importFrom SingleCellExperiment reducedDim colData
#' @importFrom stats dist hclust
#' @export
buildClusterTree <- function(sce, reduction, labelCol){
   
	if( ! labelCol %in% colnames(colData(sce)) ){
		stop("labelCol must be a column in colData(sce)")
	}

   	# extract low dimensional embeddings
    embeddings <- reducedDim(sce, reduction)

    # extract unique cell labels
    lvls = colData(sce)[[labelCol]]
    if( is.factor(lvls) ){
    	lvls = levels(lvls)
    }else{
    	lvls = unique(lvls)
    }

    # for each cell type, return mean of each dimension
    df <- lapply( lvls, function(x) {
        	idx = colData(sce)[[labelCol]] == x

        	colMeans(embeddings[idx,,drop=FALSE])
          })
    df = do.call(rbind, df)
    rownames(df) = lvls

    # compute pairwise distances
    dst <- dist( df )
 
 	# perform clustering
 	hclust( dst )
}







