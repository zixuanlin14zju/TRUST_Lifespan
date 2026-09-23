# Copy to repository root as config.local.R. That file is ignored by git.
# Paths below are examples, not the author's private paths.
options(oef.data_dir = "data", oef.output_dir = "outputs")
# Override an individual canonical input without renaming a private workbook:
# options(oef.input_map = list(
#   "trust_hc_blood_sensitivity.xlsx" = "/your/private/input/hc_final.xlsx",
#   "trust_cognitive_aging_cbf.xlsx" = "/your/private/input/cognitive_cbf.xlsx"
# ))
# In the cognitive-aging workbook, readxl may repair duplicate OEF headers to
# OEF...6 and OEF...10. Choose the analytically correct column explicitly:
# options(oef.pku_raw_oef_column = "OEF...6")
