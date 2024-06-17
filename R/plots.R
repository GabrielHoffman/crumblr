# Compute two dimensional density
#' @importFrom MASS kde2d
get_density <- function(x, y, n = 250) {
  dens <- kde2d(x = x, y = y, n = n)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

#' Scatter plot with 2D density using viridis colors
#'
#' Scatter plot with 2D density using viridis colors
#'
#' @param x the x-coordinates of points in the plot
#' @param y the y-coordinates of points in the plot
#' @param size size of point
#'
#' @return plot from ggplot2
#'
#' @examples
#' # simulate data
#' M <- Rfast::rmvnorm(1000, mu = c(0, 0), sigma = diag(1, 2))
#'
#' # create 2D density plot
#' plotScatterDensity(M[, 1], M[, 2])
#'
#' @import ggplot2
#' @importFrom viridis scale_color_viridis
#' @export
plotScatterDensity <- function(x, y, size = 1) {
  # pass R CMD check
  density <- NULL

  # convert two vectors to a data frame
  df <- data.frame(cbind(x, y))

  # determine limits of the plot
  lim <- with(df, max(abs(c(x, y))))

  # Compute 2D density
  df$density <- get_density(df$x, df$y, n = 100)

  # Scatter plot colored by density
  ggplot(df, aes(x, y, color = density)) +
    geom_point(size = size) +
    theme_classic() +
    theme(aspect.ratio = 1, legend.position = "bottom", plot.title = element_text(hjust = 0.5)) +
    scale_color_viridis(n.breaks = 3) +
    guides(fill = guide_colourbar(barwidth = 0.5))
}



#' Plot row standard deviations versus rank of row means
#'
#' Diagnositic plot for homoscedasticity across variables
#'
#' @param x data matrix
#'
#' @return plot from ggplot2
#'
#' @details
#' Plot the sd versus rank mean of each row like \code{vsn::meanSdPlot}.  Also show the coefficient of variation of the variances.  A lower value indicates stronger variance stabilization
#'
#' @examples
#' # set probability of each category
#' prob <- runif(300)
#'
#' # number of samples
#' n_samples <- 1000
#'
#' # number of counts
#' nCounts <- 3000
#'
#' # simulate counts from multinomial
#' counts <- t(rmultinom(n_samples, size = nCounts, prob = prob))
#' colnames(counts) <- paste0("cat_", 1:length(prob))
#' rownames(counts) <- paste0("sample_", 1:n_samples)
#'
#' # keep categories with at least 5 counts in at least 10 samples
#' keep <- colSums(counts > 5) > 10
#'
#' # run crumblr on counts
#' cobj <- crumblr(counts[, keep], max.ratio=10)
#'
#' # Plot for CLR
#' # For each sample, plot rank of mean vs sd
#' fig1 = meanSdPlot(cobj$E) + ggtitle("CLR")
#'
#' # run crumblr::standardize()
#' df_std <- standardize(cobj)
#'
#' # Standardized crumblr 
#' fig2 = meanSdPlot(df_std) + ggtitle("Standardized crumblr")
#'
#' # Standardizing the crumblr results better stabilizes
#' # the variances across variables
#' fig1 | fig2
#' @seealso \code{vsn::meanSdPlot}
#' @import ggplot2
#' @importFrom stats var
#' @export
meanSdPlot <- function(x) {
  df <- data.frame(
    var = apply(x, 1, var),
    mean = apply(x, 1, mean)
  )

  label <- paste0("cv: ", round(cv(df$var), digits = 3))

  ymax <- sqrt(max(df$var))
  ymin <- sqrt(min(df$var))
  # ypos = (ymax-ymin) * .99 + ymin
  # xpos = nrow(df)*.01
  ypos <- (ymax - ymin) * .05 + ymin
  xpos <- nrow(df)

  plotScatterDensity(rank(df$mean), sqrt(df$var)) + xlab("rank(mean)") + ylab("Standard deviation") + annotate("text", label = label, x = xpos, y = ypos, color = "black", hjust = 1)
}



# coefficient of variation
#' @importFrom stats sd
cv <- function(x, ...) {
  sd(x, ...) / mean(x, ...)
}
