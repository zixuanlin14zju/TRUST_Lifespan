# Shared execution, data validation and GAMLSS prediction utilities.
# This is analysis code, not an installable R package.

if (file.exists("config.local.R")) sys.source("config.local.R", envir = environment())

oef_input <- function(name = "") {
  mapping <- getOption("oef.input_map", list())
  if (nzchar(name) && !is.null(mapping[[name]])) return(mapping[[name]])
  file.path(getOption("oef.data_dir", "data"), name)
}

oef_output <- function(name = "") {
  file.path(getOption("oef.output_dir", "outputs"), name)
}

oef_start_step <- function(step) {
  options(oef.current_step = step)
  root <- oef_output("logs")
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  writeLines(c(paste("Step:", step), paste("Started:", Sys.time()),
               "Session at entry (runner also records the final session):",
               capture.output(sessionInfo())),
             file.path(root, paste0(step, "_", stamp, "_start_session.txt")))
  invisible(step)
}

oef_require_columns <- function(data, columns) {
  missing <- setdiff(columns, names(data))
  if (length(missing)) stop("Required columns missing: ", paste(missing, collapse = ", "))
  invisible(TRUE)
}

oef_binary <- function(x, label) {
  missing <- is.na(x) | trimws(as.character(x)) %in% c("", "NA", "N/A", "NaN")
  y <- suppressWarnings(as.numeric(as.character(x)))
  bad <- !missing & (!is.finite(y) | !(y %in% c(0, 1)))
  if (any(bad)) stop(label, " must contain 0, 1 or missing values; unexpected: ",
                     paste(unique(as.character(x[bad])), collapse = ", "))
  y[missing] <- NA_real_
  as.integer(y)
}

# Freeze scalar tuning parameters embedded in formulas. Without this, a saved
# model created inside a loop can depend on a later value of df_mu/df_sigma.
# Data variables and spline basis construction are not changed.
oef_freeze_formula <- function(f, data_names) {
  if (!inherits(f, "formula")) return(f)
  env <- environment(f)
  variables <- setdiff(all.vars(f), data_names)
  constants <- list()
  for (nm in variables) {
    if (exists(nm, envir = env, inherits = TRUE)) {
      value <- get(nm, envir = env, inherits = TRUE)
      if (is.atomic(value) && length(value) == 1L && !is.na(value)) constants[[nm]] <- value
    }
  }
  rewrite <- function(expr) {
    if (is.symbol(expr) && as.character(expr) %in% names(constants)) {
      return(constants[[as.character(expr)]])
    }
    if (is.call(expr)) return(as.call(lapply(as.list(expr), rewrite)))
    expr
  }
  result <- f
  for (i in seq.int(2L, length(result))) result[[i]] <- rewrite(result[[i]])
  environment(result) <- env
  result
}

oef_fit_gamlss <- function(formula, ...) {
  args <- list(formula = formula, ...)
  if (is.null(args$data)) stop("oef_fit_gamlss requires an explicit data argument.")
  dat <- as.data.frame(args$data)
  for (nm in intersect(c("formula", "sigma.fo", "nu.fo", "tau.fo"), names(args))) {
    args[[nm]] <- oef_freeze_formula(args[[nm]], names(dat))
  }
  args$data <- dat
  fit <- do.call(gamlss::gamlss, args)
  # predict.gamlss otherwise re-evaluates call$data, which can refer to a local
  # variable that no longer exists after a fit helper or a loop has returned.
  attr(fit, "oef_training_data") <- dat
  fit$call$data <- quote(.oef_training_data)
  fit
}

oef_predict_all <- function(object, newdata, data = NULL, random = "zero", ...) {
  if (is.null(data)) data <- attr(object, "oef_training_data", exact = TRUE)
  if (is.null(data)) {
    stop("Training data are not stored in this model. Supply data= explicitly, ",
         "or refit with oef_fit_gamlss().")
  }
  pred <- gamlss::predictAll(object, newdata = newdata, data = data,
                             random = random, ...)
  for (nm in intersect(c("mu", "sigma", "nu", "tau"), names(pred))) {
    if (length(pred[[nm]]) != nrow(newdata) || any(!is.finite(pred[[nm]]))) {
      stop("Invalid newdata prediction for parameter ", nm,
           ". No scalar replacement, truncation, or padding is permitted.")
    }
  }
  pred
}
