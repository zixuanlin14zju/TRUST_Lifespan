# Step15_Paired_CBF_Cognitive_Aging.R
# Public analysis entry point. Start from the repository root, or use:
#   Rscript run_step.R 15
source("R/lib/oef_utils.R")
oef_start_step("Step15")

library(readxl)
library(dplyr)
library(broom)
library(openxlsx)

InputFile <- oef_input("trust_cognitive_aging_cbf.xlsx")
if (!file.exists(InputFile)) stop("Required private input is missing: ", InputFile)
OutputDir <- file.path(oef_output(), "paired_cbf_cognitive")
dir.create(OutputDir, recursive = TRUE, showWarnings = FALSE)

D0 <- read_excel(
  InputFile,
  .name_repair = "unique"
)

cat("\n================ COLUMN NAMES ================\n")
print(names(D0))

# ============================================================
# 03) HELPER FUNCTIONS
# ============================================================

find_col <- function(dat, target) {

  nm <- names(dat)

  hit <- which(
    tolower(trimws(nm)) ==
      tolower(trimws(target))
  )

  if (length(hit) == 0) {
    stop(
      paste0(
        "Cannot find column: ",
        target
      )
    )
  }

  return(
    nm[hit[1]]
  )
}

find_col_optional <- function(dat, targets) {

  nm <- names(dat)

  for (target in targets) {

    hit <- which(
      tolower(trimws(nm)) ==
        tolower(trimws(target))
    )

    if (length(hit) > 0) {
      return(nm[hit[1]])
    }
  }

  return(NA_character_)
}

to_numeric <- function(x) {

  suppressWarnings(
    as.numeric(
      as.character(x)
    )
  )
}

to_binary <- function(x) {

  xc <- tolower(
    trimws(
      as.character(x)
    )
  )

  xn <- suppressWarnings(
    as.numeric(xc)
  )

  out <- rep(
    NA_real_,
    length(xc)
  )

  idx_numeric <-
    !is.na(xn) &
    xn %in% c(0, 1)

  out[idx_numeric] <-
    xn[idx_numeric]

  out[
    xc %in% c(
      "yes",
      "y",
      "true",
      "positive",
      "pos"
    )
  ] <- 1

  out[
    xc %in% c(
      "no",
      "n",
      "false",
      "negative",
      "neg"
    )
  ] <- 0

  return(out)
}

# ============================================================
# 04) IDENTIFY COLUMNS
# ============================================================

DiagnosisCol <- find_col(
  D0,
  "Diagnosis"
)

AgeCol <- find_col(
  D0,
  "Age"
)

SexCol <- find_col(
  D0,
  "Sex"
)

CBFCol <- find_col(
  D0,
  "CBF"
)

OEFzCol <- find_col(
  D0,
  "OEF zscore"
)

MoCACol <- find_col(
  D0,
  "MoCA"
)

BMIcodeCol <- find_col(
  D0,
  "BMI_Code"
)

HTcodeCol <- find_col(
  D0,
  "HT_Code"
)

HLcodeCol <- find_col(
  D0,
  "HL_Code"
)

DBcodeCol <- find_col(
  D0,
  "DB_Code"
)

# ------------------------------------------------------------
# Optional columns for CMRO2 analysis
# ------------------------------------------------------------
# If the source file already contains CMRO2, it will be used.
# Otherwise CMRO2 will be derived from CBF, raw OEF and Hct.

CMRO2Col <- find_col_optional(
  D0,
  c(
    "CMRO2",
    "CMRO_2",
    "CMRO₂"
  )
)

HctCol <- find_col_optional(
  D0,
  c(
    "Hct",
    "HCT",
    "Hematocrit",
    "Haematocrit"
  )
)

YaCol <- find_col_optional(
  D0,
  c(
    "Ya",
    "SaO2",
    "SaO2_percent",
    "Arterial oxygen saturation"
  )
)

cat("\n================ OPTIONAL CMRO2 INPUTS ================\n")
cat("CMRO2 column:", ifelse(is.na(CMRO2Col), "not found", CMRO2Col), "\n")
cat("Hct column:", ifelse(is.na(HctCol), "not found", HctCol), "\n")
cat("Ya column:", ifelse(is.na(YaCol), "not found", YaCol), "\n")

# ============================================================
# 05) IDENTIFY RAW OEF COLUMN
# ============================================================

OEF_candidates <- names(D0)[
  grepl(
    "^OEF($|\\.\\.\\.[0-9]+$)",
    names(D0),
    ignore.case = TRUE
  )
]

OEF_candidates <- OEF_candidates[
  !grepl(
    "zscore",
    OEF_candidates,
    ignore.case = TRUE
  )
]

if (length(OEF_candidates) == 0) {
  stop(
    "No raw OEF column detected."
  )
}

OEF_nonmissing <- sapply(
  OEF_candidates,
  function(x) {

    sum(
      !is.na(
        to_numeric(
          D0[[x]]
        )
      )
    )
  }
)

OEFcheck <- data.frame(
  Column = OEF_candidates,
  NonMissing = OEF_nonmissing
)

cat(
  "\n================ RAW OEF CANDIDATES ================\n"
)

print(OEFcheck)

RawOEFCol <- getOption("oef.pku_raw_oef_column", NULL)
if (is.null(RawOEFCol)) {
  if (length(OEF_candidates) != 1L) {
    stop("Multiple raw OEF columns: ", paste(OEF_candidates, collapse = ", "),
         ". Set options(oef.pku_raw_oef_column = 'EXACT_REPAIRED_COLUMN_NAME') in config.local.R. ",
         "The column with most observations is no longer selected automatically.")
  }
  RawOEFCol <- OEF_candidates[[1L]]
}
if (!(RawOEFCol %in% OEF_candidates)) stop("Configured raw OEF column is not among the raw OEF candidates.")

cat(
  "\nSelected raw OEF column:",
  RawOEFCol,
  "\n"
)

# ============================================================
# 06) CLEAN DATA
# ============================================================

D <- D0 %>%

  mutate(
    SourceRow = seq_len(nrow(D0)),

    # --------------------------------------------------------
    # Diagnosis
    # --------------------------------------------------------

    DiagnosisCode =
      to_numeric(
        .data[[DiagnosisCol]]
      ),

    DiagnosisF =
      factor(
        DiagnosisCode,
        levels = c(
          0,
          1,
          2
        ),
        labels = c(
          "Normal",
          "MCI",
          "Dementia"
        )
      ),

    # --------------------------------------------------------
    # Demographics
    # --------------------------------------------------------

    Age =
      to_numeric(
        .data[[AgeCol]]
      ),

    SexF =
      factor(
        trimws(
          as.character(
            .data[[SexCol]]
          )
        )
      ),

    # --------------------------------------------------------
    # Physiology
    # --------------------------------------------------------

    RawOEF =
      to_numeric(
        .data[[RawOEFCol]]
      ),

    OEF_z =
      to_numeric(
        .data[[OEFzCol]]
      ),

    CBF =
      to_numeric(
        .data[[CBFCol]]
      ),

    Existing_CMRO2 =
      if (!is.na(CMRO2Col)) {
        to_numeric(.data[[CMRO2Col]])
      } else {
        NA_real_
      },

    Measured_Hct =
      if (!is.na(HctCol)) {
        to_numeric(.data[[HctCol]])
      } else {
        NA_real_
      },

    Measured_Ya =
      if (!is.na(YaCol)) {
        to_numeric(.data[[YaCol]])
      } else {
        NA_real_
      },

    # --------------------------------------------------------
    # Cognition
    # --------------------------------------------------------

    
    MoCA =
      to_numeric(
        .data[[MoCACol]]
      ),

    # --------------------------------------------------------
    # Vascular risk
    # --------------------------------------------------------

    BMI_Code2 =
      to_binary(
        .data[[BMIcodeCol]]
      ),

    HT_Code2 =
      to_binary(
        .data[[HTcodeCol]]
      ),

    HL_Code2 =
      to_binary(
        .data[[HLcodeCol]]
      ),

    DB_Code2 =
      to_binary(
        .data[[DBcodeCol]]
      )
  )

# ============================================================
# 07) CALCULATE VRS
# ============================================================

D <- D %>%

  mutate(

    VRS_complete =
      complete.cases(
        BMI_Code2,
        HT_Code2,
        HL_Code2,
        DB_Code2
      ),

    VRS =
      ifelse(
        VRS_complete,

        BMI_Code2 +
          HT_Code2 +
          HL_Code2 +
          DB_Code2,

        NA_real_
      )
  )

# ============================================================
# 07B) CALCULATE / IMPORT CMRO2
# ============================================================
# CMRO2 is used only as a complementary reviewer analysis.
# Priority:
#   1) use an existing CMRO2 column if present;
#   2) otherwise derive CMRO2 from CBF + raw OEF + Hct.
#
# Adult reference Hct values match the primary TRUST analysis:
#   male = 0.42, female = 0.40.
# If individual Hct is available, it is used instead.
#
# IMPORTANT: if Sex is numerically coded differently in the
# source file, edit MaleSexValues / FemaleSexValues below.
# ============================================================

MaleSexValues <- c(
  "male",
  "m",
  "1",
  "男"
)

FemaleSexValues <- c(
  "female",
  "f",
  "2",
  "女"
)

Primary_Ya <- 0.98

# Oxygen carrying capacity coefficient used in the existing
# CBF-CMRO2 analyses: Ch = 20.3863636364 * Hct_fraction.
Ch_per_Hct <- 20.3863636364

D <- D %>%

  mutate(

    Sex_clean =
      tolower(
        trimws(
          as.character(SexF)
        )
      ),

    Hct_measured_fraction =
      case_when(
        is.na(Measured_Hct) ~ NA_real_,
        Measured_Hct > 1.5 ~ Measured_Hct / 100,
        TRUE ~ Measured_Hct
      ),

    Hct_reference =
      case_when(
        Sex_clean %in% MaleSexValues ~ 0.42,
        Sex_clean %in% FemaleSexValues ~ 0.40,
        TRUE ~ NA_real_
      ),

    Hct_for_CMRO2 =
      ifelse(
        !is.na(Hct_measured_fraction),
        Hct_measured_fraction,
        Hct_reference
      ),

    RawOEF_fraction =
      case_when(
        is.na(RawOEF) ~ NA_real_,
        abs(RawOEF) > 1.5 ~ RawOEF / 100,
        TRUE ~ RawOEF
      ),

    # Keep the same arterial saturation assumption as the
    # primary analysis for comparability.
    Ya_for_CMRO2 = Primary_Ya,

    Ch =
      Ch_per_Hct * Hct_for_CMRO2,

    Derived_CMRO2 =
      CBF *
      Ch *
      Ya_for_CMRO2 *
      RawOEF_fraction,

    CMRO2 =
      ifelse(
        !is.na(Existing_CMRO2),
        Existing_CMRO2,
        Derived_CMRO2
      )
  )

cat("\n================ CMRO2 CONSTRUCTION ================\n")
cat("Existing CMRO2 N =", sum(!is.na(D$Existing_CMRO2)), "\n")
cat("Derived CMRO2 N =", sum(!is.na(D$Derived_CMRO2)), "\n")
cat("Final CMRO2 N =", sum(!is.na(D$CMRO2)), "\n")
cat("Hct measured N =", sum(!is.na(D$Hct_measured_fraction)), "\n")
cat("Hct reference N =", sum(is.na(D$Hct_measured_fraction) & !is.na(D$Hct_reference)), "\n")

if (
  all(is.na(D$Existing_CMRO2)) &&
  any(!is.na(D$CBF) & !is.na(D$RawOEF_fraction) & is.na(D$Hct_for_CMRO2))
) {

  cat(
    "WARNING: Some subjects cannot receive reference Hct because Sex coding was not recognized.\n"
  )

  cat("Observed Sex values:\n")
  print(table(D$SexF, useNA = "ifany"))
}

cat(
  "\n================ DIAGNOSIS COUNTS ================\n"
)

print(
  table(
    D$DiagnosisF,
    useNA = "ifany"
  )
)

cat(
  "\n================ VRS DISTRIBUTION ================\n"
)

print(
  table(
    D$VRS,
    useNA = "ifany"
  )
)

# ============================================================
# 08) ORIGINAL OEF-ELIGIBLE COHORT
#
# IMPORTANT:
# CBF IS NOT REQUIRED HERE.
# ============================================================

D_OEF <- D %>%

  filter(
    !is.na(RawOEF),
    !is.na(OEF_z),
    !is.na(Age),
    !is.na(SexF),
    !is.na(DiagnosisF)
  ) %>%

  droplevels()

cat(
  "\n================ OEF COHORT ================\n"
)

cat(
  "OEF cohort N =",
  nrow(D_OEF),
  "\n"
)

print(
  table(
    D_OEF$DiagnosisF
  )
)

# ------------------------------------------------------------

# --- Model-specific cohorts; preserve the original missing-data rules ---
D_CBF <- D_OEF %>%
  filter(!is.na(CBF), DiagnosisF %in% c("Normal", "MCI")) %>%
  droplevels()
D_CBF_VRS <- D_CBF %>% filter(!is.na(VRS)) %>% droplevels()

D_MoCA_MCI_VRS <- D_OEF %>%
  filter(DiagnosisF == "MCI", !is.na(MoCA), !is.na(VRS)) %>%
  droplevels()
D_CBF_MoCA_MCI_VRS <- D_CBF %>%
  filter(DiagnosisF == "MCI", !is.na(MoCA), !is.na(VRS)) %>%
  droplevels()

# Preserve the original S3 cohort: do NOT additionally require RawOEF/OEF_z.
# CMRO2 was either supplied or derived above using the original rules.
D_CMRO2_MoCA_MCI_VRS <- D %>%
  filter(DiagnosisF == "MCI", !is.na(MoCA), !is.na(CMRO2),
         !is.na(Age), !is.na(SexF), !is.na(VRS)) %>%
  droplevels()

ModelResults <- list()
ModelInfo <- list()

fit_model <- function(
    ModelName,
    Formula,
    Data) {

  Data <- droplevels(
    Data
  )

  
  fit <- tryCatch(

    lm(
      Formula,
      data = Data
    ),

    error = function(e) {

      cat(
        "\nMODEL FAILED:",
        ModelName,
        "\n"
      )

      cat(
        e$message,
        "\n"
      )

      return(NULL)
    }
  )

  
  if (is.null(fit)) {

    ModelInfo[[ModelName]] <<-

      data.frame(

        Model =
          ModelName,

        Formula =
          paste(
            deparse(Formula),
            collapse = ""
          ),

        N =
          nrow(Data),

        R2 =
          NA,

        Adj_R2 =
          NA,

        AIC =
          NA,

        BIC =
          NA,

        Status =
          "FAILED"
      )

    return(NULL)
  }

  
  tmp <- broom::tidy(
    fit,
    conf.int = TRUE
  )

  tmp$Model <-
    ModelName

  
  ModelResults[[ModelName]] <<-

    tmp %>%

    select(
      Model,
      everything()
    )

  
  g <- broom::glance(
    fit
  )

  
  ModelInfo[[ModelName]] <<-

    data.frame(

      Model =
        ModelName,

      Formula =
        paste(
          deparse(Formula),
          collapse = ""
        ),

      N =
        nobs(fit),

      R2 =
        g$r.squared,

      Adj_R2 =
        g$adj.r.squared,

      AIC =
        AIC(fit),

      BIC =
        BIC(fit),

      Status =
        "OK"
    )

  
  return(
    fit
  )
}

# --- Six retained fits, with original model names for result matching ---

fit_model(

  "N1_CBF_Dx",

  CBF ~
    DiagnosisF +
    Age +
    SexF,

  D_CBF
)

fit_model(

  "O2_OEF_Dx_plusCBF",

  OEF_z ~
    DiagnosisF +
    CBF +
    Age +
    SexF,

  D_CBF
)

fit_model(

  "O4_OEF_Dx_CBF_plusVRS",

  OEF_z ~
    DiagnosisF +
    CBF +
    Age +
    SexF +
    VRS,

  D_CBF_VRS
)

fit_model(
  "G_MoCA_MCI_3_plusVRS",
  MoCA ~ OEF_z + Age + SexF + VRS,
  D_MoCA_MCI_VRS
)

fit_model(

  "R6_MoCA_MCI_CBFplusVRS",

  MoCA ~
    OEF_z +
    CBF +
    Age +
    SexF +
    VRS,

  D_CBF_MoCA_MCI_VRS
)

fit_model(

  "S3_MoCA_MCI_CMRO2_plusVRS",

  MoCA ~
    CMRO2 +
    Age +
    SexF +
    VRS,

  D_CMRO2_MoCA_MCI_VRS
)

# --- Full coefficients are retained only for these six fits ---
AllModelResults <- bind_rows(ModelResults)
AllModelInfo <- bind_rows(ModelInfo)
if (nrow(AllModelInfo) != 6L || any(AllModelInfo$Status != "OK")) {
  stop("At least one model failed. See model failure messages above; ",
       "no workbook has been written.")
}

# Seven reported effects from six fits. The reported CBF--MoCA coefficient is
# from R6 (the same jointly adjusted model), NOT a newly fitted CBF-only model.
ReportedTerms <- data.frame(
  Model = c("N1_CBF_Dx", "O2_OEF_Dx_plusCBF", "O4_OEF_Dx_CBF_plusVRS",
            "G_MoCA_MCI_3_plusVRS", "R6_MoCA_MCI_CBFplusVRS",
            "R6_MoCA_MCI_CBFplusVRS", "S3_MoCA_MCI_CMRO2_plusVRS"),
  term = c("DiagnosisFMCI", "DiagnosisFMCI", "DiagnosisFMCI", "OEF_z",
           "OEF_z", "CBF", "CMRO2"),
  Analysis = c("MCI vs Normal: CBF adjusted for age and sex",
               "MCI vs Normal: OEFz adjusted for CBF, age, and sex",
               "MCI vs Normal: OEFz adjusted for CBF, age, sex, and VRS",
               "MCI MoCA: OEFz adjusted for age, sex, and VRS; CBF not required",
               "MCI MoCA: OEFz adjusted for CBF, age, sex, and VRS",
               "MCI MoCA: CBF term from the same OEFz/CBF/age/sex/VRS model",
               "MCI MoCA: CMRO2 adjusted for age, sex, and VRS"),
  stringsAsFactors = FALSE
)
KeyEffects <- ReportedTerms %>%
  left_join(AllModelResults, by = c("Model", "term")) %>%
  left_join(AllModelInfo %>% select(Model, N), by = "Model") %>%
  select(Analysis, Model, N, term, estimate, std.error, statistic,
         conf.low, conf.high, p.value)
if (nrow(KeyEffects) != 7L || any(!is.finite(KeyEffects$estimate)) ||
    any(!is.finite(KeyEffects$p.value))) {
  stop("A reported coefficient is missing/non-finite; check group coding and model rank.")
}

# Descriptive sample/row auditing, not additional hypothesis tests.
cohort_list <- list(
  Paired_Normal_MCI = D_CBF,
  Paired_Normal_MCI_VRS = D_CBF_VRS,
  MCI_MoCA_VRS_CBF_not_required = D_MoCA_MCI_VRS,
  MCI_MoCA_VRS_CBF_available = D_CBF_MoCA_MCI_VRS,
  MCI_MoCA_VRS_CMRO2_available = D_CMRO2_MoCA_MCI_VRS
)
SampleFlow <- bind_rows(lapply(names(cohort_list), function(nm) {
  x <- cohort_list[[nm]]
  data.frame(Cohort = nm, N = nrow(x),
             Normal_N = sum(x$DiagnosisF == "Normal"),
             MCI_N = sum(x$DiagnosisF == "MCI"))
}))
AnalysisRows <- bind_rows(lapply(names(cohort_list), function(nm) {
  x <- cohort_list[[nm]]
  data.frame(Cohort = rep(nm, nrow(x)), SourceRow = x$SourceRow)
}))

TimeStamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
OutputWorkbook <- file.path(OutputDir,
  paste0("Cognitive_OEF_CBF_MoCA_", TimeStamp, ".xlsx"))
wb <- createWorkbook()
output_sheets <- list(
  KeyEffects = KeyEffects,
  AllModels = AllModelResults,
  ModelInfo = AllModelInfo,
  SampleFlow = SampleFlow,
  AnalysisRows = AnalysisRows,
  OEFcheck = OEFcheck
)
header_style <- createStyle(textDecoration = "bold", border = "Bottom")
for (sheet_name in names(output_sheets)) {
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, output_sheets[[sheet_name]])
  addStyle(wb, sheet_name, header_style, rows = 1,
           cols = seq_len(ncol(output_sheets[[sheet_name]])), gridExpand = TRUE)
  setColWidths(wb, sheet_name, cols = seq_len(ncol(output_sheets[[sheet_name]])),
               widths = 18)
  setColWidths(wb, sheet_name, cols = 1, widths = 58)
  freezePane(wb, sheet_name, firstRow = TRUE)
}
saveWorkbook(wb, OutputWorkbook, overwrite = FALSE)
cat("\nSAMPLE FLOW\n")
print(SampleFlow)
cat("\nKEY EFFECTS\n")
print(KeyEffects, row.names = FALSE)
cat("\nResults saved to:\n", OutputWorkbook, "\n")
