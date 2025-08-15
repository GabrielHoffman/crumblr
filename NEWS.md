# crumblr 0.99.22
 - August 11, 2025
 - update docs
 - add `fdr.cutoff` argument to `plotTreeTest()` and `plotTreeTestBeta()

# crumblr 0.99.21
 - March 10, 2025
 - address Bioc comments

# crumblr 0.99.20
 - March 10, 2025
 - address Bioc comments

# crumblr 0.99.19
 - Feb 12, 2025
 - address Bioc comments: https://github.com/Bioconductor/Contributions/issues/3699#issuecomment-2640806359

# crumblr 0.99.17
 - Feb 6, 2025
 - address Bioc comments: https://github.com/Bioconductor/Contributions/issues/3699#issuecomment-2640806359



# crumblr 0.99.13
 - Jan 14, 2025
 - Bioc submission

# crumblr 0.99.11
 - June 14, 2024
 - improve documentation 
 - add `plotForest()`
 - add `standardize()`


# crumblr 0.99.10
 - June 12, 2024
 - add `diffTree()` and `plotTreeTestBeta()`
 - regularize weight estimates

# crumblr 0.99.9
 - Feb 1, 2024
 - add `crumblr(..,method=method)` for `method %in% c("clr", "clr_2class")`

# crumblr 0.99.8
 - May 12, 2023
 - in `treeTest()` add argument `shrink.cov = TRUE`.  This will give slightly different p-values, but protects against false positives when the number of features get large compared to the sample size

# crumblr 0.99.7
 - April 5, 2023
 - fix plotting bug in `plotTreeTest()`

# crumblr 0.99.6
 - Dec 7, 2022
 - update docs, add vignette, dependencies

# crumblr 0.99.5
 - Nov 29, 2022
 - bug fix and update docs

# crumblr 0.99.4
 - enforce `variancePartition` version dependency

# crumblr 0.99.3
 - `buildClusterTree()` uses cosine distance by default
 - `treeTest()` now uses `method="treeTest"` by default

# crumblr 0.99.2
 - update `plotTreeTest()`
 - Use variancePartition dependency on GitHub

# crumblr 0.99.1
 - add functions to perform `mvTest()` on trees

# crumblr 0.99.0
 - Jun 24, 2022
 - Fix version number for compatibility with Bioconductor

# crumblr 1.0.7
 - Feb 1, 2022
 - fix `vst()`
 - add plots `plotScatterDensity()` and `meanSdPlot()` for VST
 - add vignette vst.Rmd

# crumblr 1.0.6
 - Jan 31, 2022
 - add logo, example data

# crumblr 1.0.5
 - Dec 27, 2021
 - add `dmn.mle()`
 - use generic `crumblr()` function

# crumblr 1.0.4
 - Dec 23, 2021
 - Update vignette
 - export `clr()`

# crumblr 1.0.3
 - Dec 20, 2021
 - fix `vst` bug

# crumblr 1.0.2
 - Dec 20, 2021
 - add `vst`

# crumblr 1.0.1
 - Oct 29, 2021
 - handle `data.frame`

# crumblr 1.0.0
 - Oct 28, 2021
 - Initial version