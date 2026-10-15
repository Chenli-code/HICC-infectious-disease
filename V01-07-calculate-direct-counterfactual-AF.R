###############################################################################
## Direct paired attribution for the three major infectious-disease classes
##
## For every GCM, city and day:
##   counterfactual risk ratio = RR(hist-nat) / RR(factual)
##   direct ACC cases          = y * (1 - counterfactual risk ratio)
##
## The final point estimate is the median of the 10 GCM-specific estimates.
## The 95% uncertainty interval is obtained by resampling only these 10 GCM
## point estimates with replacement 10,000 times. No health-model coefficient,
## city or day is resampled.
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(dlnm)
})

ROOT <- Sys.getenv(
  "ATTR_PROJECT_ROOT",
  getwd()
)
source(file.path(ROOT, "03.0-user-route-config.R"))
source(USER_CORE)

# -------------------------------- settings ----------------------------------

BASELINE <- Sys.getenv("ACC_CASE_BASELINE", "observed_daily")
if (!BASELINE %in% c("observed_daily", "day_of_year")) {
  stop("ACC_CASE_BASELINE must be 'observed_daily' or 'day_of_year'.")
}
REUSE_CACHE <- tolower(
  Sys.getenv("ACC_REUSE_CACHE", "true")
) == "true"

OUT_DIR <- file.path(USER_OUT, "direct_counterfactual_af")
CACHE_DIR <- file.path(OUT_DIR, "by_gcm_cache_v2")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(CACHE_DIR, recursive = TRUE, showWarnings = FALSE)

CITY_GCM_FILE <- file.path(
  OUT_DIR, "V01-direct-ACC-AF-by-city-and-gcm.csv"
)
NATIONAL_GCM_FILE <- file.path(
  OUT_DIR, "V01-direct-ACC-AF-national-by-gcm.csv"
)
NATIONAL_SUMMARY_FILE <- file.path(
  OUT_DIR, "V01-direct-ACC-AF-national-summary.csv"
)
AUDIT_FILE <- file.path(
  OUT_DIR, "V01-direct-ACC-AF-identity-audit.csv"
)

fit_files <- setNames(
  file.path(USER_OUT, paste0("erf_pooled_", USER_SPECS$tag, ".rds")),
  USER_SPECS$tag
)
stopifnot(
  all(file.exists(fit_files)),
  all(file.exists(vapply(USER_GCMS, user_pair_file, character(1))))
)

cache_file <- function(gcm) {
  file.path(
    CACHE_DIR,
    paste0("direct_ACC_", user_safe_gcm(gcm), ".rds")
  )
}

# --------------------------- one-GCM calculation ----------------------------

calculate_one_gcm <- function(gcm) {
  climate_file <- user_pair_file(gcm)
  cache_path <- cache_file(gcm)
  metadata <- list(
    formula_version = ACC_FORMULA_VERSION,
    baseline = BASELINE,
    climate_md5 = acc_file_md5(climate_file),
    fit_md5 = vapply(fit_files, acc_file_md5, character(1))
  )

  if (
    REUSE_CACHE &&
      acc_cache_is_current(cache_path, metadata)
  ) {
    message("Reusing verified direct-ACC cache: ", gcm)
    return(readRDS(cache_path)$result)
  }

  message("Calculating direct paired ACC: ", gcm)
  pair <- readRDS(climate_file)
  disease_results <- vector("list", length(fit_files))

  for (ii in seq_along(fit_files)) {
    tag <- names(fit_files)[ii]
    fit <- readRDS(fit_files[[ii]])
    cases <- as.data.table(fit$crossover)[, .(
      code = as.character(code),
      date = as.Date(date),
      n = as.numeric(n)
    )]
    case_weights <- acc_prepare_case_weights(
      cases = cases,
      dates = USER_DATES,
      baseline = BASELINE
    )
    codes <- sort(unique(case_weights$code))
    climate <- acc_complete_climate_pair(
      pair = pair,
      codes = codes,
      dates = USER_DATES
    )
    analysis_data <- merge(
      climate,
      case_weights,
      by = c("code", "date"),
      all.x = TRUE
    )
    analysis_data[is.na(n), n := 0]
    setorder(analysis_data, code, date)

    city_rows <- analysis_data[
      ,
      {
        result <- acc_direct_city(
          factual_temperature = tf,
          histnat_temperature = tn,
          cases = n,
          fit = fit
        )
        result
      },
      by = code
    ]
    city_rows[, `:=`(
      tag = tag,
      classification = unname(USER_DISEASE_LABELS[tag]),
      gcm = gcm,
      province_code = substr(code, 1L, 2L),
      lag = as.integer(fit$lag),
      baseline = BASELINE,
      formula_version = ACC_FORMULA_VERSION
    )]
    acc_check_identity(
      city_rows$total_cases,
      city_rows$counterfactual_cases,
      city_rows$an_acc_direct,
      city_rows$af_acc_direct
    )
    disease_results[[ii]] <- city_rows
  }

  result <- rbindlist(disease_results, use.names = TRUE, fill = TRUE)
  saveRDS(
    list(metadata = metadata, result = result),
    cache_path,
    compress = TRUE
  )
  result
}

city_gcm <- rbindlist(
  lapply(USER_GCMS, calculate_one_gcm),
  use.names = TRUE,
  fill = TRUE
)
city_gcm[, gcm_order__ := match(gcm, USER_GCMS)]
setorder(city_gcm, tag, gcm_order__, code)
city_gcm[, gcm_order__ := NULL]

# ------------------------- national aggregation -----------------------------

national_gcm <- city_gcm[, .(
  total_cases = sum(total_cases),
  included_forward_cases = sum(included_forward_cases),
  an_acc_direct = sum(an_acc_direct),
  counterfactual_cases = sum(counterfactual_cases)
), by = .(
  tag, classification, gcm, lag, baseline, formula_version
)]
national_gcm[, af_acc_direct := (
  100 * an_acc_direct / total_cases
)]
acc_check_identity(
  national_gcm$total_cases,
  national_gcm$counterfactual_cases,
  national_gcm$an_acc_direct,
  national_gcm$af_acc_direct
)

national_summary <- acc_summarise_gcm_metrics(
  data = national_gcm,
  grouping_columns = c(
    "tag", "classification", "lag",
    "baseline", "formula_version"
  ),
  metric_columns = c(
    "an_acc_direct",
    "af_acc_direct",
    "counterfactual_cases"
  ),
  gcm_levels = USER_GCMS
)
national_summary[, gcm_count := length(USER_GCMS)]

audit <- national_gcm[, .(
  case_identity_error = (
    total_cases - counterfactual_cases - an_acc_direct
  ),
  af_identity_error = (
    af_acc_direct - 100 * an_acc_direct / total_cases
  )
), by = .(tag, gcm, baseline, formula_version)]

# -------------------------------- outputs -----------------------------------

fwrite(city_gcm, CITY_GCM_FILE)
fwrite(national_gcm, NATIONAL_GCM_FILE)
fwrite(national_summary, NATIONAL_SUMMARY_FILE)
fwrite(audit, AUDIT_FILE)

message("Direct paired ACC results written to: ", OUT_DIR)
print(national_summary[, .(
  classification,
  af = af_acc_direct_estimate,
  lower = af_acc_direct_lower,
  upper = af_acc_direct_upper,
  cases = an_acc_direct_estimate
)])

