suppressPackageStartupMessages({
  library(data.table)
  library(SuperLearner)
  library(glmnet)
  library(pROC)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) {
  stop(
    "Usage: ensemble_superlearner.R predictors phenotype phenotype_column covariates trait_type output_prefix seed",
    call. = FALSE
  )
}

predictor_path <- args[[1L]]
phenotype_path <- args[[2L]]
phenotype_column <- args[[3L]]
covariates <- if (nzchar(args[[4L]])) strsplit(args[[4L]], ",", fixed = TRUE)[[1L]] else character()
trait_type <- args[[5L]]
output_prefix <- args[[6L]]
seed <- as.integer(args[[7L]])
set.seed(seed)

predictors <- fread(predictor_path)
phenotype <- fread(phenotype_path)
if (!("IID" %in% names(predictors))) stop("Predictor table must contain IID", call. = FALSE)
if (!(phenotype_column %in% names(phenotype))) {
  stop(sprintf("Phenotype column is absent: %s", phenotype_column), call. = FALSE)
}

phenotype_id <- names(phenotype)[[1L]]
phenotype_columns <- unique(c(phenotype_id, phenotype_column, covariates))
missing_covariates <- setdiff(covariates, names(phenotype))
if (length(missing_covariates)) {
  stop(sprintf("Missing phenotype covariates: %s", paste(missing_covariates, collapse = ",")), call. = FALSE)
}

data <- merge(
  phenotype[, ..phenotype_columns],
  predictors,
  by.x = phenotype_id,
  by.y = "IID"
)
data <- as.data.frame(na.omit(data))
if (nrow(data) < 4L) stop("At least four complete samples are required for the ensemble", call. = FALSE)

# Preserve the published-development workflow's two-step split, including the
# shared boundary row between the second train/test split.
data <- data[ceiling(nrow(data) / 2):nrow(data), , drop = FALSE]
if (length(covariates) == 0L || trait_type == "binary") {
  eta <- data[[phenotype_column]]
} else {
  formula <- as.formula(sprintf("%s ~ %s", phenotype_column, paste(covariates, collapse = " + ")))
  eta <- summary(lm(formula, data = data))$residuals
}

predictor_columns <- setdiff(names(predictors), "IID")
split_index <- ceiling(nrow(data) / 2)
train_rows <- 1:split_index
test_rows <- split_index:nrow(data)
data_train <- data[train_rows, predictor_columns, drop = FALSE]
data_test <- data[test_rows, predictor_columns, drop = FALSE]
eta_train <- eta[train_rows]
eta_test <- eta[test_rows]

individual_metrics <- vapply(seq_along(predictor_columns), function(index) {
  values <- as.numeric(data_test[[index]])
  if (trait_type == "quantitative") {
    cor(eta_test, values)^2
  } else {
    as.numeric(auc(suppressMessages(roc(data[[phenotype_column]][test_rows], values))))
  }
}, numeric(1L))

best_index <- which.max(individual_metrics)
best_metric <- individual_metrics[[best_index]]
best_prediction <- as.numeric(data_test[[best_index]])
prediction_source <- predictor_columns[[best_index]]
superlearner_metric <- NA_real_

if (trait_type == "quantitative") {
  model <- SuperLearner(
    Y = eta_train,
    X = data_train,
    family = gaussian(),
    SL.library = c("SL.glmnet", "SL.ridge")
  )
  superlearner_prediction <- as.numeric(predict(model, data_test, onlySL = TRUE)$pred)
  superlearner_metric <- cor(eta_test, superlearner_prediction)^2
  if (!is.na(superlearner_metric) && superlearner_metric >= best_metric) {
    best_metric <- superlearner_metric
    best_prediction <- superlearner_prediction
    prediction_source <- "SuperLearner(SL.glmnet,SL.ridge)"
  }
}

predictions <- data.frame(
  IID = data[[phenotype_id]][test_rows],
  observed = data[[phenotype_column]][test_rows],
  prediction = best_prediction,
  source = prediction_source,
  stringsAsFactors = FALSE
)
metrics <- data.frame(
  metric = if (trait_type == "quantitative") "R2" else "AUC",
  value = best_metric,
  selected_source = prediction_source,
  superlearner_metric = superlearner_metric,
  seed = seed,
  stringsAsFactors = FALSE
)

fwrite(predictions, paste0(output_prefix, ".predictions.tsv"), sep = "\t")
fwrite(metrics, paste0(output_prefix, ".metrics.tsv"), sep = "\t")
