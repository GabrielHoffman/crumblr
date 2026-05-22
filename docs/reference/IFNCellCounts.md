# Cell counts following interferon treatment

Counts are from single cell RNA-seq data from treated and untreated
samples from Kang, et al (2018).

## Usage

``` r
data(IFNCellCounts)

data(info)

data(df_cellCounts)

data(hcl)
```

## Format

- `info` is metadata for each sample

- `df_cellCounts` data.frame of counts for each sample

- `hcl` cluster of cell types based on pseudobulk expression

An object of class `data.frame` with 16 rows and 4 columns.

An object of class `matrix` (inherits from `array`) with 16 rows and 8
columns.

An object of class `hclust` of length 7.

## References

Kang, Hyun Min, et al. "Multiplexed droplet single-cell RNA-sequencing
using natural genetic variation." Nature Biotechnology 36.1 (2018):
89-94.
