## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

# Keep data.table within CRAN's thread limits
data.table::setDTthreads(2)

## ----load---------------------------------------------------------------------
library(WSSmod)

## ----list---------------------------------------------------------------------
list_prebuilt_wss()

## ----scores-------------------------------------------------------------------
weights <- load_prebuilt_wss()
dim(weights)
names(weights)

## ----scores-head--------------------------------------------------------------
head(weights[, c("module", "symbol", "term", "mean_beta", "mean_alpha_scaled")])

## ----duplicates---------------------------------------------------------------
c(
  rows = nrow(weights),
  proteins = length(unique(weights$symbol)),
  modules = length(unique(weights$module))
)

# Proteins that are in more than one module
multi <- weights[, .N, by = symbol][N > 1]
nrow(multi)
weights[symbol == multi$symbol[1], c("module", "symbol", "term", "mean_alpha_scaled")]

## ----dedup--------------------------------------------------------------------
res <- calculate_WSS(wss_example_expression)
c(
  modules_scored = ncol(res$scores),
  proteins_scored = length(unique(unlist(res$module_proteins)))
)

## ----terms--------------------------------------------------------------------
terms <- load_prebuilt_terms()
dim(terms)
head(terms)

## ----terms-match--------------------------------------------------------------
modules <- colnames(res$scores)[1:3]
terms$term[match(modules, terms$name)]

## ----terms-group--------------------------------------------------------------
table(terms$group)

## ----model--------------------------------------------------------------------
model_info <- load_prebuilt_wss_model()
names(model_info)
model_info$variant
model_info$x_cols
head(model_info$outcomes)
length(model_info$outcomes)

## ----reference----------------------------------------------------------------
reference <- load_prebuilt_biomarker_reference()
head(reference, 3)
reference[, .(n = .N, median = signif(median(value), 3)), by = biomarker]

## ----simulated----------------------------------------------------------------
dim(wss_example_expression)
dim(wss_example_biomarkers)

