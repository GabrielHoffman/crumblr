% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/logFrac.R
\name{logFrac}
\alias{logFrac}
\title{Log fractions and precision weights}
\usage{
logFrac(counts, pseudocount = 0.5, max.ratio = 5, quant = 0.05)
}
\arguments{
\item{counts}{count data with samples as rows and variables are columns}

\item{pseudocount}{added to counts to avoid issues with zeros}

\item{max.ratio}{regularize estimates of the weights to have a maximum ratio of \code{max.ratio} between the maximum and \code{quant} quantile value}

\item{quant}{quantile value used for \code{max.ratio}}
}
\value{
An \code{EList} object with the following components:
\describe{
 \item{E: }{numeric matrix of log transformed counts}
 \item{weights: }{numeric matrix of observation-level inverse-variance weights}
}
}
\description{
Compute log fractions and precision weights from matrix of c ounts, where columns are variables and rows are samples
}
\details{
For real data, the asymptotic variance formula can give weights that vary substantially across samples and give very high weights for a subset of samples.  In order to address this, we regularize the weights to reduce the variation in the weights to have a maximum ratio of \code{max.ratio} between the maximum and \code{quant} quantile value.
}
\examples{
# set probability of each category
prob <- c(0.1, 0.2, 0.3, 0.5)

# number of total counts
countsTotal <- 300

# number of samples
n_samples <- 100

# simulate info for each sample
info <- data.frame(Age = rgamma(n_samples, 50, 1))
rownames(info) <- paste0("sample_", 1:n_samples)

# simulate counts from multinomial
counts <- t(rmultinom(n_samples, size = countsTotal, prob = prob))
colnames(counts) <- paste0("cat_", 1:length(prob))
rownames(counts) <- paste0("sample_", 1:n_samples)

# run logFrac on counts
cobj <- logFrac(counts)

# run standard variancePartition analysis on crumblr results
library(variancePartition)

fit <- dream(cobj, ~ Age, info)
fit <- eBayes(fit)

topTable(fit, coef = "Age", sort.by = "none")
}
\seealso{
\code{limma::voom()}, \code{variancePartition::dream()}
}
