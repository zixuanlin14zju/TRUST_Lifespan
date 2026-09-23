# Threshold-only delta-R2 quality-control input.
# Load R/lib/oef_utils.R before this file for oef_binary().
# 1 = strictly above the named threshold; 0 = at or below it.

oef_validate_dr2_flags <- function(data) {
  columns <- c("dR2_gt5", "dR2_gt10")
  missing_columns <- setdiff(columns, names(data))
  if (length(missing_columns)) {
    stop("Missing QC flag columns: ", paste(missing_columns, collapse = ", "),
         call. = FALSE)
  }
  flags <- data.frame(
    dR2_gt5 = oef_binary(data[["dR2_gt5"]], "dR2_gt5"),
    dR2_gt10 = oef_binary(data[["dR2_gt10"]], "dR2_gt10")
  )
  if (any(xor(is.na(flags$dR2_gt5), is.na(flags$dR2_gt10)))) {
    stop("dR2_gt5 and dR2_gt10 must be missing together: both flags derive from the same measurement.",
         call. = FALSE)
  }
  if (any(flags$dR2_gt10 == 1L & flags$dR2_gt5 == 0L, na.rm = TRUE)) {
    stop("Inconsistent QC flags: dR2_gt10 = 1 requires dR2_gt5 = 1.",
         call. = FALSE)
  }
  flags
}
