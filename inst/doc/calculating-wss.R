## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

# Keep data.table within CRAN's thread limits
data.table::setDTthreads(2)

# RNOmni is a suggested package; chunks that need it are skipped without it
has_rnomni <- requireNamespace("RNOmni", quietly = TRUE)

## ----load---------------------------------------------------------------------
library(WSSmod)

## ----mock-raw-----------------------------------------------------------------
# Mock raw values: positive and right-skewed, like intensities
raw <- exp(wss_example_expression)

## ----z-log2-------------------------------------------------------------------
# z-scaled log2 values, one protein per column
z_log2 <- scale(log2(raw))
dim(z_log2)

## ----int, eval = has_rnomni---------------------------------------------------
# Rank-based inverse normal transformation, one protein per column
int_values <- apply(raw, 2, RNOmni::RankNorm)
dim(int_values)

## ----default------------------------------------------------------------------
res <- calculate_WSS(wss_example_expression)
names(res)

## ----scores-------------------------------------------------------------------
dim(res$scores)
round(head(res$scores[, 1:4], 3), 2)

## ----summary------------------------------------------------------------------
head(res$summary_stats, 3)
round(head(res$module_sd, 3), 2)

## ----module-proteins----------------------------------------------------------
# Number of proteins per module
summary(lengths(res$module_proteins))
head(res$module_proteins[["Merged.M17"]])

## ----boxplot, fig.width = 6, fig.height = 3.5, fig.alt = "Horizontal boxplots of WSS scores for the first eight modules, labeled by module term, across the 50 simulated samples."----
# Spread of the first eight module scores (simulated data, no biological
# meaning), labeled by module term
terms <- load_prebuilt_terms("core_AD_plasma_biomarkers")
labels <- terms$term[match(colnames(res$scores)[1:8], terms$name)]
par(mar = c(4, 14, 1, 1))
boxplot(res$scores[, 1:8], names = labels, horizontal = TRUE, las = 1,
        cex.axis = 0.7, xlab = "WSS score")

## ----custom-table-------------------------------------------------------------
proteins <- colnames(wss_example_expression)[1:12]

custom <- data.table::data.table(
  module = c(rep("modA", 5), rep("modB", 6), rep("modC", 2)),
  symbol = c(proteins[1:5], proteins[6:10], proteins[1], proteins[11:12]),
  mean_beta = c(rep(1, 5), rep(-1, 6), rep(1, 2)),
  mean_alpha_scaled = c(0.9, 0.8, 0.7, 0.6, 0.5, 0.9, 0.8, 0.7, 0.6, 0.5, 0.3,
                        0.8, 0.8),
  specificity = c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.9,
                  0.5, 0.5)
)
custom

## ----custom-default-----------------------------------------------------------
custom_res <- calculate_WSS(wss_example_expression, custom)
lengths(custom_res$module_proteins)

## ----custom-dedup-------------------------------------------------------------
by_specificity <- calculate_WSS(
  wss_example_expression,
  custom,
  dedup_by = "specificity"
)
lengths(by_specificity$module_proteins)

## ----custom-size--------------------------------------------------------------
with_small <- calculate_WSS(
  wss_example_expression,
  custom,
  min_module_size = 2
)
lengths(with_small$module_proteins)

## ----partial------------------------------------------------------------------
# Keep only the first 300 proteins; the warnings are suppressed here
partial <- suppressWarnings(calculate_WSS(wss_example_expression[, 1:300]))

# Modules that could be scored, and modules with no proteins in the matrix
length(partial$module_proteins)
sum(colSums(is.na(partial$scores)) == nrow(partial$scores))

## ----missing-values-----------------------------------------------------------
# One missing value in a protein that a module uses
with_na <- wss_example_expression
with_na[1, 2] <- NA

# calculate_WSS() stops and names the affected protein
tryCatch(calculate_WSS(with_na), error = function(e) conditionMessage(e))

## ----projection, eval = has_rnomni--------------------------------------------
reference_rows <- 1:40
new_rows <- 41:50

# Rank-normalize the reference cohort in-sample, one protein per column
reference_int <- apply(raw[reference_rows, ], 2, RNOmni::RankNorm)

# Project each new sample onto the frozen reference, protein by protein
new_int <- vapply(
  colnames(raw),
  function(protein) {
    project_rank_norm(raw[reference_rows, protein], raw[new_rows, protein])
  },
  numeric(length(new_rows))
)
rownames(new_int) <- rownames(raw)[new_rows]

## ----projection-scores, eval = has_rnomni-------------------------------------
reference_scores <- calculate_WSS(reference_int)
new_scores <- calculate_WSS(new_int)

dim(reference_scores$scores)
dim(new_scores$scores)

