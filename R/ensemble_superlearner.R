suppressPackageStartupMessages({
  library(data.table)
  library(SuperLearner)
  library(glmnet)
  library(nnet)
  library(pROC)
})

SL.lasso.HERCULES <- function(Y, X, newX, family, obsWeights, id, ...) {
  SuperLearner::SL.glmnet(
    Y = Y, X = X, newX = newX, family = family,
    obsWeights = obsWeights, id = id, alpha = 1, ...
  )
}

SL.ridge.HERCULES <- function(Y, X, newX, family, obsWeights, id, ...) {
  SuperLearner::SL.glmnet(
    Y = Y, X = X, newX = newX, family = family,
    obsWeights = obsWeights, id = id, alpha = 0, ...
  )
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 9L) {
  stop(
    paste(
      "Usage: ensemble_superlearner.R validation_predictors validation_phenotype",
      "test_predictors test_phenotype_or_empty phenotype_column covariates",
      "trait_type output_prefix seed"
    ),
    call. = FALSE
  )
}

validation_predictor_path <- args[[1L]]
validation_phenotype_path <- args[[2L]]
test_predictor_path <- args[[3L]]
test_phenotype_path <- args[[4L]]
phenotype_column <- args[[5L]]
covariates <- if (nzchar(args[[6L]])) strsplit(args[[6L]], ",", fixed = TRUE)[[1L]] else character()
trait_type <- args[[7L]]
output_prefix <- args[[8L]]
seed <- as.integer(args[[9L]])
set.seed(seed)

score_columns <- c("target_stage1_score", "calibrated_stage2_score")
read_scores <- function(path, label) {
  values <- fread(path)
  expected <- c("IID", score_columns)
  if (!identical(names(values), expected)) {
    stop(
      sprintf("%s score table must contain exactly: %s", label, paste(expected, collapse = ", ")),
      call. = FALSE
    )
  }
  if (anyDuplicated(values$IID)) {
    stop(sprintf("%s score table contains duplicate IID values", label), call. = FALSE)
  }
  if (any(!is.finite(as.matrix(values[, ..score_columns])))) {
    stop(sprintf("%s score table contains non-finite predictors", label), call. = FALSE)
  }
  values
}

validation_scores <- read_scores(validation_predictor_path, "Validation")
test_scores <- read_scores(test_predictor_path, "Test")
overlap <- intersect(validation_scores$IID, test_scores$IID)
if (length(overlap)) {
  stop(
    sprintf("Validation and test IIDs must be disjoint; overlap includes %s", overlap[[1L]]),
    call. = FALSE
  )
}

validation_phenotype <- fread(validation_phenotype_path)
required_validation <- c("IID", phenotype_column, covariates)
missing_validation <- setdiff(required_validation, names(validation_phenotype))
if (length(missing_validation)) {
  stop(
    sprintf("Validation phenotype is missing columns: %s", paste(missing_validation, collapse = ",")),
    call. = FALSE
  )
}
if (anyDuplicated(validation_phenotype$IID)) {
  stop("Validation phenotype contains duplicate IID values", call. = FALSE)
}

validation_rows <- match(validation_scores$IID, validation_phenotype$IID)
if (anyNA(validation_rows)) {
  stop("Validation phenotype does not cover every scored validation IID", call. = FALSE)
}
validation <- as.data.frame(cbind(
  validation_scores,
  as.data.frame(validation_phenotype[validation_rows, ..required_validation][, IID := NULL])
))
if (!all(complete.cases(validation[, c(phenotype_column, covariates), drop = FALSE]))) {
  stop("Validation phenotype or covariates contain missing values", call. = FALSE)
}
if (nrow(validation) < 10L) {
  stop("At least ten complete target validation samples are required", call. = FALSE)
}

validation_x <- validation[, score_columns, drop = FALSE]
validation_y <- validation[[phenotype_column]]
covariate_model <- NULL
if (trait_type == "quantitative" && length(covariates)) {
  formula <- as.formula(
    sprintf("%s ~ %s", phenotype_column, paste(covariates, collapse = " + "))
  )
  covariate_model <- lm(formula, data = validation)
  validation_y <- residuals(covariate_model)
}

if (trait_type == "quantitative") {
  learner_library <- c("SL.lasso.HERCULES", "SL.ridge.HERCULES", "SL.nnet")
  meta_method <- "method.NNLS"
  family_object <- gaussian()
} else if (trait_type == "binary") {
  unique_y <- sort(unique(validation_y))
  if (!all(unique_y %in% c(0, 1)) || length(unique_y) != 2L) {
    stop("Binary validation phenotype must contain both 0 and 1", call. = FALSE)
  }
  learner_library <- c("SL.lasso.HERCULES", "SL.nnet")
  meta_method <- "method.AUC"
  family_object <- binomial()
} else {
  stop("trait_type must be quantitative or binary", call. = FALSE)
}

model <- SuperLearner(
  Y = validation_y,
  X = validation_x,
  family = family_object,
  SL.library = learner_library,
  method = meta_method
)

test_x <- as.data.frame(test_scores[, ..score_columns])
test_prediction <- as.numeric(predict(model, newdata = test_x, onlySL = TRUE)$pred)
if (trait_type == "binary") {
  test_prediction <- pmin(1, pmax(0, test_prediction))
}

# SuperLearner records wall-clock timings in the fitted object. They are not
# part of the fitted prediction function and make otherwise identical model
# serializations differ across installations, so omit them from the portable
# reproducible representation.
model$times <- NULL

model_bundle <- list(
  superlearner = model,
  learner_library = learner_library,
  meta_method = meta_method,
  trait_type = trait_type,
  phenotype_column = phenotype_column,
  covariates = covariates,
  covariate_model = covariate_model,
  seed = seed,
  score_columns = score_columns
)
saveRDS(model_bundle, paste0(output_prefix, ".model.rds"))

predictions <- data.frame(
  IID = test_scores$IID,
  prediction = test_prediction,
  stringsAsFactors = FALSE
)

metadata <- data.frame(
  trait_type = trait_type,
  learner_library = paste(learner_library, collapse = ","),
  meta_method = meta_method,
  seed = seed,
  stringsAsFactors = FALSE
)
fwrite(metadata, paste0(output_prefix, ".model-metadata.tsv"), sep = "\t")

if (nzchar(test_phenotype_path)) {
  test_phenotype <- fread(test_phenotype_path)
  required_test <- c("IID", phenotype_column, covariates)
  missing_test <- setdiff(required_test, names(test_phenotype))
  if (length(missing_test)) {
    stop(
      sprintf("Test phenotype is missing columns: %s", paste(missing_test, collapse = ",")),
      call. = FALSE
    )
  }
  if (anyDuplicated(test_phenotype$IID)) {
    stop("Test phenotype contains duplicate IID values", call. = FALSE)
  }

  test_rows <- match(predictions$IID, test_phenotype$IID)
  if (anyNA(test_rows)) {
    stop("Test phenotype does not cover every scored test IID", call. = FALSE)
  }
  evaluated <- cbind(
    predictions,
    as.data.frame(test_phenotype[test_rows, ..required_test][, IID := NULL])
  )
  observed <- evaluated[[phenotype_column]]
  if (trait_type == "quantitative" && !is.null(covariate_model)) {
    observed <- observed - as.numeric(predict(covariate_model, newdata = evaluated))
  }
  predictions$observed <- observed

  if (trait_type == "quantitative") {
    metric_name <- "R2"
    metric_value <- cor(observed, predictions$prediction)^2
  } else {
    metric_name <- "AUC"
    metric_value <- as.numeric(
      auc(suppressMessages(roc(observed, predictions$prediction, quiet = TRUE)))
    )
  }
  metrics <- data.frame(
    metric = metric_name,
    value = metric_value,
    seed = seed,
    stringsAsFactors = FALSE
  )
  fwrite(metrics, paste0(output_prefix, ".metrics.tsv"), sep = "\t")
}

fwrite(predictions, paste0(output_prefix, ".predictions.tsv"), sep = "\t")
