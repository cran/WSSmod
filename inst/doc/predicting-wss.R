## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

# Keep data.table within CRAN's thread limits
data.table::setDTthreads(2)

# RNOmni, glmnet and joinet are suggested packages; chunks that need them are
# skipped without them
has_rnomni <- requireNamespace("RNOmni", quietly = TRUE)
has_model <- requireNamespace("glmnet", quietly = TRUE) &&
  requireNamespace("joinet", quietly = TRUE)

## ----load---------------------------------------------------------------------
library(WSSmod)

## ----inputs-------------------------------------------------------------------
dim(wss_example_biomarkers)
head(wss_example_biomarkers, 3)

## ----split--------------------------------------------------------------------
# The eight raw biomarker columns, without Age and Gender
biomarker_cols <- setdiff(names(wss_example_biomarkers), c("Age", "Gender"))
raw_biomarkers <- wss_example_biomarkers[, biomarker_cols]

## ----project-one--------------------------------------------------------------
# One patient
one_patient <- normalize_wss_biomarkers(raw_biomarkers[1, ])
round(one_patient, 2)

## ----project-cohort-----------------------------------------------------------
# A cohort of 20 patients
projected <- normalize_wss_biomarkers(raw_biomarkers)
all.equal(projected[1, ], one_patient)

## ----self, eval = has_rnomni--------------------------------------------------
self_normalized <- normalize_wss_biomarkers(raw_biomarkers, method = "self")
round(head(self_normalized, 3), 2)

## ----self-one, eval = has_rnomni----------------------------------------------
tryCatch(
  normalize_wss_biomarkers(raw_biomarkers[1, ], method = "self"),
  error = function(e) conditionMessage(e)
)

## ----predict, eval = has_model------------------------------------------------
newdata <- cbind(
  wss_example_biomarkers[, c("Age", "Gender")],
  projected
)
pred <- predict_WSS(newdata)
names(pred)

## ----predict-shape, eval = has_model------------------------------------------
dim(pred$meta)
round(pred$meta[1:3, 1:4], 2)

## ----predict-one, eval = has_model--------------------------------------------
single <- predict_WSS(cbind(wss_example_biomarkers[1, c("Age", "Gender")],
                            one_patient))
all.equal(single$meta[1, ], pred$meta[1, ])

## ----spread, eval = has_model-------------------------------------------------
# Standard deviation of the predicted scores, first five modules
round(apply(pred$meta, 2, sd)[1:5], 2)

## ----boxplot, eval = has_model, fig.width = 6, fig.height = 3.5, fig.alt = "Horizontal boxplots of predicted WSS scores for the first eight modules, labeled by module term, across 20 simulated patients."----
# Spread of the first eight predicted module scores (simulated data, no
# biological meaning), labeled by module term
terms <- load_prebuilt_terms("core_AD_plasma_biomarkers")
labels <- terms$term[match(colnames(pred$meta)[1:8], terms$name)]
par(mar = c(4, 14, 1, 1))
boxplot(pred$meta[, 1:8], names = labels, horizontal = TRUE, las = 1,
        cex.axis = 0.7, xlab = "Predicted WSS score")

