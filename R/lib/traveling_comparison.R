# Paired site comparisons using the traveling-subject table.
traveling_site_comparison <- function(dat) {
  required <- c("name", "OEF_SJTU", "OEF_ZJU")
  if (!all(required %in% names(dat))) {
    stop("Required columns: ", paste(required, collapse = ", "))
  }
  ids <- trimws(as.character(dat$name))
  if (anyNA(ids) || any(!nzchar(ids))) stop("Participant names must be non-missing.")
  if (anyDuplicated(ids)) stop("Duplicate participant names; resolve before analysis.")

  metric_columns <- list(OEF = c("OEF_SJTU", "OEF_ZJU"),
                         R2 = c("R2_SJTU", "R2_ZJU"))
  results <- list()
  paired_values <- list()
  availability <- list()
  for (metric in names(metric_columns)) {
    columns <- metric_columns[[metric]]
    present <- columns %in% names(dat)
    if (!any(present)) {
      availability[[metric]] <- data.frame(Metric = metric, Status = "COLUMNS_NOT_PROVIDED",
                                           N_pairs = 0L, stringsAsFactors = FALSE)
      next
    }
    if (!all(present)) stop("Supply both ", paste(columns, collapse = " and "), ".")
    if (!all(vapply(dat[columns], is.numeric, logical(1)))) {
      stop(metric, " columns must be numeric.")
    }
    values <- unlist(dat[columns], use.names = FALSE)
    if (any(!is.na(values) & !is.finite(values))) stop(metric, " contains non-finite values.")
    if (metric == "OEF" && any(values < 0 | values > 1, na.rm = TRUE)) {
      stop("OEF must be a fraction between 0 and 1.")
    }
    if (metric == "R2" && any(values <= 0, na.rm = TRUE)) stop("R2 must be positive s^-1.")

    ok <- is.finite(dat[[columns[1]]]) & is.finite(dat[[columns[2]]])
    pairs <- data.frame(name = ids[ok], SJTU = dat[[columns[1]]][ok],
                        ZJU = dat[[columns[2]]][ok], stringsAsFactors = FALSE)
    if (nrow(pairs) < 3L) stop("At least three paired participants are required for ", metric, ".")
    tt <- stats::t.test(pairs$SJTU, pairs$ZJU, paired = TRUE)
    ct <- stats::cor.test(pairs$SJTU, pairs$ZJU, method = "pearson")
    results[[metric]] <- data.frame(
      Metric = metric, N_pairs = nrow(pairs), SJTU_mean = mean(pairs$SJTU),
      SJTU_SD = stats::sd(pairs$SJTU), ZJU_mean = mean(pairs$ZJU), ZJU_SD = stats::sd(pairs$ZJU),
      Mean_difference_SJTU_minus_ZJU = mean(pairs$SJTU - pairs$ZJU),
      CI_low = unname(tt$conf.int[1]), CI_high = unname(tt$conf.int[2]),
      Paired_t = unname(tt$statistic), Paired_df = unname(tt$parameter), Paired_P = tt$p.value,
      Pearson_r = unname(ct$estimate), Pearson_P = ct$p.value,
      stringsAsFactors = FALSE)
    paired_values[[metric]] <- data.frame(Metric = metric, pairs, row.names = NULL)
    availability[[metric]] <- data.frame(Metric = metric, Status = "ANALYZED",
                                         N_pairs = nrow(pairs), stringsAsFactors = FALSE)
  }
  list(summary = do.call(rbind, results), paired_values = do.call(rbind, paired_values),
       availability = do.call(rbind, availability))
}
