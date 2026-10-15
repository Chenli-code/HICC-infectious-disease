###############################################################################
## Supplementary materials | human-induced climate change and infectious
## disease burden
##
## This script assembles supplementary figures and tables from the final,
## quality-controlled outputs of Figures 1-5. It does not refit the health
## models or recompute the climate attribution.
##
## Reference strategy retained throughout:
##   * Respiratory infectious diseases: disease-specific MMT
##   * Gastrointestinal/enterovirus diseases: P50 curve reference
##   * Vector-borne diseases: P50 curve reference
## Disease-relevant attribution for the latter two groups uses each city's
## four hottest consecutive months.
##
## Outputs:
##   Supplementary/Figures  Supplementary Figures S1-S5
##   Supplementary/Tables   Supplementary Tables S1-S8
##   Supplementary/V01-Supplementary-figures.pptx
##   Supplementary/V01-Supplementary-tables.xlsx
##   Supplementary/V01-Supplementary-index.csv
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(scales)
})

ROOT <- Sys.getenv(
  "ATTR_PROJECT_ROOT",
  getwd()
)

source(file.path(ROOT, "03.0-user-route-config.R"))

source(USER_CORE)

# -------------------------------- paths --------------------------------------
SUPP_DIR <- file.path(USER_ROOT, "Supplementary")

SUPP_FIG_DIR <- file.path(SUPP_DIR, "Figures")

SUPP_TAB_DIR <- file.path(SUPP_DIR, "Tables")

dir.create(SUPP_FIG_DIR, recursive = TRUE, showWarnings = FALSE)

dir.create(SUPP_TAB_DIR, recursive = TRUE, showWarnings = FALSE)

FIG2_DIR <- file.path(USER_ROOT, "Figure2")

FIG4_DIR <- file.path(USER_ROOT, "Figure4")

FIG5_DIR <- file.path(USER_ROOT, "Figure5")

DIRECT_GCM_FILE <- file.path(
  USER_OUT,
  "direct_counterfactual_af",
  "V01-direct-ACC-AF-national-by-gcm.csv"
)

ANNUAL_FILE <- file.path(
  USER_OUT,
  "04.1-annual-major-class-components-summary.csv"
)

ANNUAL_SLOPE_FILE <- file.path(
  USER_OUT,
  "04.1-annual-major-class-trend-summary.csv"
)

CITY_FILE <- file.path(
  FIG4_DIR,
  "V01-Figure4-direct-city-map-data.csv"
)

REGIONAL_FILE <- file.path(
  FIG4_DIR,
  "V01-Figure4-direct-regional-forest-summary.csv"
)

DISEASE_GCM_FILE <- file.path(
  FIG5_DIR,
  "V01-Figure5-disease-specific-direct-attribution-by-gcm.csv"
)

DISEASE_STATUS_FILE <- file.path(
  FIG5_DIR,
  "V01-Figure5-disease-specific-analysis-status.csv"
)

REFERENCE_FILE <- file.path(
  FIG2_DIR,
  "V01-Figure2-model-reference-temperatures.csv"
)

SELECTED_RR_FILE <- file.path(
  FIG2_DIR,
  "V01-Figure2-selected-percentile-RR.csv"
)

COVERAGE_FILES <- setNames(
  file.path(
    USER_OUT,
    paste0(
      "erf_pooled_",
      c(
        "respiratory",
        "gastrointestinal",
        "vectorborne"
      ),
      "_coverage.csv"
    )
  ),
  c(
    "respiratory",
    "gastrointestinal",
    "vectorborne"
  )
)

required_files <- c(
  DIRECT_GCM_FILE,
  ANNUAL_FILE,
  ANNUAL_SLOPE_FILE,
  CITY_FILE,
  REGIONAL_FILE,
  DISEASE_GCM_FILE,
  DISEASE_STATUS_FILE,
  REFERENCE_FILE,
  SELECTED_RR_FILE,
  unname(COVERAGE_FILES)
)

if (any(!file.exists(required_files))) {
  stop(
    "Required final-analysis files are missing:\n",
    paste(
      required_files[!file.exists(required_files)],
      collapse = "\n"
    )
  )
}

GCM_BOOT_N <- 10000L

GCM_BOOT_SEED <- 20260724L

disease_order <- c(
  "respiratory",
  "gastrointestinal",
  "vectorborne"
)

disease_labels <- c(
  respiratory = "Respiratory infectious diseases",
  gastrointestinal =
    "Gastrointestinal and enterovirus diseases",
  vectorborne = "Vector-borne diseases"
)

supp_bootstrap_indices <- acc_bootstrap_indices(
  length(USER_GCMS),
  repetitions = GCM_BOOT_N,
  seed = GCM_BOOT_SEED
)

bootstrap_median <- function(x) {
  result <- acc_bootstrap_median(x, supp_bootstrap_indices)
  c(
    estimate = unname(result["estimate"]),
    low = unname(result["lower"]),
    high = unname(result["upper"]),
    gcm_min = unname(result["gcm_min"]),
    gcm_max = unname(result["gcm_max"])
  )
}

write_table <- function(x, name) {
  fwrite(
    as.data.table(x),
    file.path(
      SUPP_TAB_DIR,
      paste0(name, ".csv")
    )
  )
}

s1_target <- file.path(
  SUPP_FIG_DIR,
  "V01-SFigure1-regions-and-scenario-temperatures.tif"
)

# ------------- Supplementary Figure S2: GCM-specific direct AF ---------------
direct_gcm <- fread(DIRECT_GCM_FILE)

direct_gcm[, tag := factor(
  tag,
  levels = disease_order
)]

direct_gcm[, gcm := factor(
  gcm,
  levels = rev(USER_GCMS)
)]

direct_gcm[, disease := factor(
  unname(disease_labels[as.character(tag)]),
  levels = unname(disease_labels[disease_order])
)]

disease_status <- fread(DISEASE_STATUS_FILE)

# -------------------------- Supplementary tables -----------------------------
# Table S1: model definitions and data coverage
coverage <- rbindlist(
  lapply(
    names(COVERAGE_FILES),
    function(tag_i) {
      x <- fread(COVERAGE_FILES[[tag_i]])
      x[, tag := tag_i]
      x
    }
  ),
  fill = TRUE
)

references <- fread(REFERENCE_FILE)

table_s1 <- merge(
  coverage,
  references[, .(
    tag,
    disease,
    reference_temperature,
    reference_type,
    common_temperature_min,
    common_temperature_max
  )],
  by = "tag",
  all.x = TRUE,
  sort = FALSE
)

table_s1[, `:=`(
  study_start = as.character(USER_START_DATE),
  study_end = as.character(USER_END_DATE),
  gcm_count = length(USER_GCMS)
)]

setcolorder(
  table_s1,
  c(
    "tag",
    "disease",
    "lag",
    "reference_type",
    "reference_temperature",
    "common_temperature_min",
    "common_temperature_max",
    "health_cities",
    "total_cases",
    "observed_climate_cities",
    "climate_imputed_cities",
    "gcm_count",
    "study_start",
    "study_end"
  )
)

# Table S2: selected percentile RRs
table_s2 <- fread(
  SELECTED_RR_FILE,
  sep = ",",
  header = TRUE
)

# Table S3: GCM-specific national estimates
table_s3 <- copy(direct_gcm)

table_s3[, `:=`(
  tag = as.character(tag),
  gcm = as.character(gcm),
  disease = as.character(disease)
)]

table_s3[, `:=`(
  tag_order = match(tag, disease_order),
  gcm_order = match(gcm, USER_GCMS)
)]

setorder(
  table_s3,
  tag_order,
  gcm_order
)

table_s3[, c("tag_order", "gcm_order") := NULL]

# Table S4: coherent factual-versus-counterfactual direct attribution
required_direct_columns <- c(
  "tag", "gcm", "total_cases", "counterfactual_cases",
  "an_acc_direct", "af_acc_direct"
)

if (any(!required_direct_columns %in% names(direct_gcm))) {
  stop(
    "The direct-ACC file is outdated. Rerun ",
    "V01-07-calculate-direct-counterfactual-AF.R."
  )
}

table_s4 <- direct_gcm[, {
  factual <- bootstrap_median(total_cases)
  counterfactual <- bootstrap_median(counterfactual_cases)
  af_acc <- bootstrap_median(af_acc_direct)
  an_acc <- bootstrap_median(an_acc_direct)
  list(
    disease = unname(disease_labels[tag[1L]]),
    factual_cases = factual["estimate"],
    factual_cases_low = factual["low"],
    factual_cases_high = factual["high"],
    counterfactual_cases = counterfactual["estimate"],
    counterfactual_cases_low = counterfactual["low"],
    counterfactual_cases_high = counterfactual["high"],
    direct_acc_cases = an_acc["estimate"],
    direct_acc_cases_low = an_acc["low"],
    direct_acc_cases_high = an_acc["high"],
    direct_acc_af = af_acc["estimate"],
    direct_acc_af_low = af_acc["low"],
    direct_acc_af_high = af_acc["high"]
  )
}, by = tag]

table_s4[, tag_order := match(
  tag,
  disease_order
)]

setorder(table_s4, tag_order)

table_s4[, tag_order := NULL]

# Table S5: regional direct AF and cases
table_s5 <- fread(REGIONAL_FILE)

# Table S6: disease-specific included results
table_s6 <- disease_status[
  eligible == TRUE
]

setorder(table_s6, disease_order)

# Table S7: complete disease-level eligibility and exclusion audit
table_s7 <- copy(disease_status)

setorder(table_s7, disease_order)

# Table S8: annual trend slopes
table_s8 <- fread(ANNUAL_SLOPE_FILE)

# ----------------------------- Excel workbook --------------------------------
tables_for_workbook <- list(
  S1_model_coverage = table_s1,
  S2_percentile_RR = table_s2,
  S3_GCM_estimates = table_s3,
  S4_national_summary = table_s4,
  S5_regional = table_s5,
  S6_disease_results = table_s6,
  S7_eligibility_audit = table_s7,
  S8_annual_trends = table_s8
)

workbook_file <- file.path(
  SUPP_DIR,
  "V01-Supplementary-tables.xlsx"
)

wb <- openxlsx::createWorkbook()

header_style <- openxlsx::createStyle(
  textDecoration = "bold",
  fgFill = "#E8ECEE",
  border = "Bottom",
  halign = "center",
  valign = "center",
  wrapText = TRUE
)

for (sheet_i in names(tables_for_workbook)) {
  sheet_data <- as.data.frame(
    tables_for_workbook[[sheet_i]]
  )
  openxlsx::addWorksheet(
    wb,
    sheet_i,
    gridLines = FALSE
  )
  openxlsx::writeData(
    wb,
    sheet_i,
    sheet_data,
    headerStyle = header_style,
    withFilter = TRUE
  )
  openxlsx::freezePane(
    wb,
    sheet_i,
    firstRow = TRUE
  )
  openxlsx::setColWidths(
    wb,
    sheet_i,
    cols = seq_len(ncol(sheet_data)),
    widths = "auto"
  )
}

openxlsx::saveWorkbook(
  wb,
  workbook_file,
  overwrite = TRUE
)

# -------------------------- supplement index --------------------------------
supp_index <- data.table(
  type = c(
    rep("Figure", 5),
    rep("Table", 8)
  ),
  id = c(
    paste0("Figure S", 1:5),
    paste0("Table S", 1:8)
  ),
  title = c(
    "Regional classification and factual versus hist-nat temperatures",
    "GCM-specific direct anthropogenic climate-attributable fractions",
    "Annual variation in direct anthropogenic climate-attributable fractions",
    "City-level distributions of direct AF by region",
    "Disease-specific direct AF across the 10 GCMs",
    "Model specification and data coverage",
    "Relative risks at selected temperature percentiles",
    "Major-class direct attribution estimates by GCM",
    "Coherent factual-versus-counterfactual direct attribution",
    "Regional direct attribution estimates",
    "Disease-specific included results",
    "Disease-level eligibility and exclusion audit",
    "Annual attribution trend slopes"
  ),
  file = c(
    basename(s1_target),
    "V01-SFigure2-GCM-specific-direct-AF",
    "V01-SFigure3-annual-direct-AF-trends",
    "V01-SFigure4-city-AF-distributions-by-region",
    "V01-SFigure5-disease-by-GCM-heatmap",
    "V01-STable1-model-specification-and-coverage.csv",
    "V01-STable2-selected-temperature-percentile-RRs.csv",
    "V01-STable3-major-class-estimates-by-GCM.csv",
    "V01-STable4-national-ensemble-estimates.csv",
    "V01-STable5-regional-direct-attribution.csv",
    "V01-STable6-disease-specific-results.csv",
    "V01-STable7-disease-eligibility-audit.csv",
    "V01-STable8-annual-attribution-trend-slopes.csv"
  )
)

fwrite(
  supp_index,
  file.path(
    SUPP_DIR,
    "V01-Supplementary-index.csv"
  )
)

caption_lines <- c(
  "# Supplementary figures and tables",
  "",
  paste0(
    "All attribution estimates use 10 equally weighted GCMs. ",
    "Point estimates are GCM medians and 95% uncertainty intervals ",
    "are obtained from 10,000 GCM-bootstrap resamples. Respiratory ",
    "diseases use the fitted disease-specific minimum-risk temperature. ",
    "Gastrointestinal/enterovirus and vector-borne curves use P50; their ",
    "disease-relevant attribution uses each city's four hottest consecutive ",
    "months."
  ),
  "",
  "## Supplementary figures",
  "",
  paste0(
    "**Figure S1. Regional classification and temperatures under the ",
    "factual and natural-only counterfactual climates.** The panels ",
    "show the seven geographical regions and the spatial distributions ",
    "of factual temperature, hist-nat temperature and their difference. ",
    "City-level temperature summaries are restricted to unique study cities ",
    "with health outcome data."
  ),
  "",
  paste0(
    "**Figure S2. GCM-specific direct anthropogenic climate-attributable ",
    "fractions.** Points show estimates from each GCM and solid vertical ",
    "lines show the 10-GCM median. Positive values indicate cases ",
    "avoidable in the absence of anthropogenic warming; negative values ",
    "indicate a reduction in burden under the factual climate."
  ),
  "",
  paste0(
    "**Figure S3. Annual variation in direct anthropogenic ",
    "climate-attributable fractions.** Lines show annual point estimates ",
    "and shaded areas show 95% GCM-bootstrap uncertainty intervals."
  ),
  "",
  paste0(
    "**Figure S4. City-level distributions of direct anthropogenic ",
    "climate-attributable fractions by region.** Violin plots show the ",
    "city distribution and embedded boxplots show medians and ",
    "interquartile ranges. Values are winsorised at disease-specific ",
    "P1 and P99 for display only; analyses use the original estimates."
  ),
  "",
  paste0(
    "**Figure S5. Disease-specific attribution across the 10 GCMs.** ",
    "Cells show the direct anthropogenic climate-attributable fraction ",
    "for each disease-GCM combination. Colours are centred at zero and ",
    "clipped at the 97.5th percentile for visual comparability."
  ),
  "",
  "## Supplementary tables",
  "",
  "**Table S1.** Model specification, reference temperature and data coverage.",
  "**Table S2.** Relative risks at selected observed temperature percentiles.",
  "**Table S3.** Major-class direct attribution estimates from each GCM.",
  paste0(
    "**Table S4.** Observed factual cases, predicted natural-only ",
    "counterfactual cases and their direct paired difference."
  ),
  "**Table S5.** Regional direct attribution estimates.",
  "**Table S6.** Disease-specific results included in the main analysis.",
  "**Table S7.** Complete disease-level eligibility and exclusion audit.",
  "**Table S8.** Annual attribution trend slopes."
)

writeLines(
  caption_lines,
  file.path(
    SUPP_DIR,
    "V01-Supplementary-captions.md"
  ),
  useBytes = TRUE
)
