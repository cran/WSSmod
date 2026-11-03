
<!-- README.md is generated from README.Rmd. Please edit that file -->

# WSSmod <a href="https://Longo-Lab.github.io/WSSmod/"><img src="man/figures/logo.png" align="right" height="138" alt="WSSmod website" /></a>

<!-- badges: start -->

[![Version](https://img.shields.io/badge/version-0.1.0-blue.svg)](https://github.com/Longo-Lab/WSSmod/blob/main/DESCRIPTION)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check](https://github.com/Longo-Lab/WSSmod/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Longo-Lab/WSSmod/actions/workflows/R-CMD-check.yaml)
[![codecov](https://codecov.io/github/Longo-Lab/WSSmod/graph/badge.svg?token=YU275M964D)](https://app.codecov.io/github/Longo-Lab/WSSmod)
<!-- badges: end -->

WSSmod provides tools to construct Weighted Sum Score (WSS) module-level
summaries from proteomic data. Proteins are grouped into STRING network
modules, sign-flipped to align with the direction of core biomarker
association, and combined into a single module score per sample using a
weighted average based on network connectivity (alpha centrality) and
biomarker association strength.

This package accompanies the Longo Lab preprint, Butler et al. (2026),
[Alzheimer’s disease blood biomarkers reveal proteomic modules of
disease
progression](https://www.medrxiv.org/content/10.64898/2026.08.14.26360394v1.full-text).

## Installation

You can install the development version of WSSmod from GitHub with:

``` r
# install.packages("pak")
pak::pak("Longo-Lab/WSSmod")
```

## Example

Using a custom module/weighting table:

``` r
library(WSSmod)

expr_matrix <- matrix(
  c(1, 2, 3, 4, 5, 6, 7, 8),
  nrow = 2,
  dimnames = list(c("sample1", "sample2"), c("geneA", "geneB", "geneC", "geneD"))
)

result <- data.table::data.table(
  module = c("mod1", "mod1", "mod1", "mod2"),
  symbol = c("geneA", "geneB", "geneC", "geneD"),
  mean_beta = c(1, -1, 1, 1),
  mean_alpha_scaled = c(1, 1, 2, 1)
)

WSSmod::calculate_WSS(expr_matrix, result, min_module_size = 1)
#> $scores
#>         mod1 mod2
#> sample1  2.0    7
#> sample2  2.5    8
#> 
#> $summary_stats
#>    module   Min    Q1 Median  Mean    Q3   Max
#>    <char> <num> <num>  <num> <num> <num> <num>
#> 1:   mod1     2 2.125   2.25  2.25 2.375   2.5
#> 2:   mod2     7 7.250   7.50  7.50 7.750   8.0
#> 
#> $module_sd
#>      mod1      mod2 
#> 0.3535534 0.7071068 
#> 
#> $module_proteins
#> $module_proteins$mod1
#> [1] "geneC" "geneA" "geneB"
#> 
#> $module_proteins$mod2
#> [1] "geneD"
```

WSSmod also ships prebuilt module/weighting reference sets so you don’t
have to build your own module table from scratch. `calculate_WSS()` uses
`prebuilt = "core_AD_plasma_biomarkers"` by default when `result =
NULL`:

``` r
list_prebuilt_wss()

# expr_matrix here has samples as rows and protein/gene symbols as columns,
# matching the symbols used by the prebuilt reference set
WSSmod::calculate_WSS(expr_matrix)

# Human-readable module term labels for the same reference set
load_prebuilt_terms("core_AD_plasma_biomarkers")
```

## Predicting module scores from core biomarkers

`predict_WSS()` predicts module scores for a prebuilt reference set from
a small set of covariates, using a fitted `joinet` model, without
needing the full proteomic panel. For `"core_AD_plasma_biomarkers"` this
predicts all 75 module scores from Age, Gender, and the 8 core plasma
biomarkers used in the original analysis. This requires the `glmnet` and
`joinet` packages.

The 8 biomarker columns must already be rank-based inverse-normal
transformed before calling `predict_WSS()`. `normalize_wss_biomarkers()`
builds these columns (keeping the same column names) from raw biomarker
values using one of two methods:

  - `method = "project"` (default) — project each raw value onto a
    bundled reference distribution via `project_rank_norm()`. Works for
    a single new patient or a small/differently-distributed cohort,
    since each value is scored independently against the fixed
    reference.
  - `method = "self"` — rank-normalize (`RNOmni::RankNorm()`) within
    your own batch, if you have a sizeable cohort of your own with a
    distribution you’re comfortable normalizing against directly instead
    of the bundled reference. Requires multiple samples, since a rank
    transform needs multiple values to rank against.

<!-- end list -->

``` r
# raw_biomarkers: one row per sample, the 8 raw (un-normalized) biomarker columns.
# Works for a single patient (method = "project" is the default):
raw_biomarkers <- data.frame(
  PlasmaPTau181 = 1.5, PlasmaAB142P = 25, PlasmaAB140P = 300, PlasmaABRatio = 0.09,
  PlasmapTau217 = 0.3, PlasmapTau217_AB42Ratio = 0.012, PlasmaGFAP = 60, PlasmaNfL = 22
)
normalized <- normalize_wss_biomarkers(raw_biomarkers)

newdata <- cbind(Age = 72, Gender = 1, normalized) # Gender: 1 = male, 0 = not male
pred <- predict_WSS(newdata)
pred$meta
```

## License

The WSSmod code is released under the MIT license. The bundled reference
files are derived from the data and results of the preprint above, which
is licensed under [CC BY-NC 4.0
International](https://creativecommons.org/licenses/by-nc/4.0/); please
observe its terms, including non-commercial use, when reusing them.
