###############################################################################
## Disease-relevant attribution | City-specific four hottest consecutive
## months
##
## Purpose
##   Replace the former percentile-threshold attribution definition with a
##   calendar-based warm season. Respiratory disease retains its full
##   non-MMT temperature-attributable burden; gastrointestinal/enterovirus and
##   vector-borne diseases use positive temperature-attributable burden during
##   each city's four hottest consecutive months.
##
## Warm-season definition
##   For every city, the four consecutive calendar months with the highest
##   mean observed ERA5 temperature during 2008-2019 are selected. Consecutive
##   windows are evaluated cyclically across the end of the calendar year
##   (for example, November-February is an eligible window). The same fixed
##   city-specific months are then applied to the paired factual and hist-nat
##   climates for every GCM. The selection therefore cannot be affected by the
##   anthropogenic climate signal being estimated.
##
## Primary estimand
##   Respiratory:
##     full non-MMT AN in the factual and hist-nat worlds.
##   Gastrointestinal/enterovirus and vector-borne:
##                         AN_factual,warm4 - AN_histnat,warm4
##     warm4 ACC share = 100 --------------------------------------
##                                  AN_factual,warm4
##
##   AN_factual,warm4 and AN_histnat,warm4 are positive
##   temperature-attributable cases accumulated during the same city-specific
##   four hottest consecutive months in the two climate worlds. Positive parts
##   are used so that protective deficits (RR < 1) cannot shrink the
##   denominator and inflate the percentage.
##
## Secondary estimand
##   The reference-invariant direct paired ACC AF among all cases during the
##   same four months is retained as a validation analysis.
##
## Uncertainty
##   Point estimates are ratios of the mean paired burden across 10 equally
##   weighted GCMs. The 95% uncertainty intervals use 10,000 bootstrap
##   resamples of the 10 GCM point estimates only.
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
source(file.path(ROOT, "V01-12-reference-AC-core.R"))

# ------------------------------- settings ------------------------------------

ANALYSIS_VERSION <- paste0(
  "disease-relevant-burden-city-warmest-",
  "four-consecutive-months-v4"
)
BASELINE <- "day_of_year"
WARM_MONTH_COUNT <- USER_WARM_MONTH_COUNT
REUSE_CACHE <- tolower(
  Sys.getenv("WARM4_AF_REUSE_CACHE", "true")
) == "true"

OUT_DIR <- file.path(
  USER_OUT,
  "warm4_attribution"
)
CACHE_DIR <- file.path(OUT_DIR, "cache")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(CACHE_DIR, recursive = TRUE, showWarnings = FALSE)

CITY_MONTH_FILE <- file.path(
  OUT_DIR,
  "V01-warmest-four-months-by-city.csv"
)
CITY_GCM_FILE <- file.path(
  OUT_DIR,
  "V01-warm4-temperature-burden-city-by-gcm.csv"
)
NATIONAL_GCM_FILE <- file.path(
  OUT_DIR,
  "V01-warm4-temperature-burden-national-by-gcm.csv"
)
SUMMARY_FILE <- file.path(
  OUT_DIR,
  "V01-warm4-temperature-burden-national-summary.csv"
)

fit_files <- setNames(
  file.path(
    USER_OUT,
    paste0("erf_pooled_", USER_SPECS$tag, ".rds")
  ),
  USER_SPECS$tag
)
stopifnot(
  all(file.exists(fit_files)),
  all(file.exists(vapply(
    USER_GCMS,
    user_pair_file,
    character(1)
  ))),
  file.exists(USER_PATH_OBS)
)
fit_reference_current <- vapply(
  seq_along(fit_files),
  function(ii) {
    fit_i <- readRDS(fit_files[[ii]])
    expected_i <- USER_SPECS$reference[
      match(names(fit_files)[ii], USER_SPECS$tag)
    ]
    identical(
      tolower(as.character(fit_i$reference)),
      tolower(as.character(expected_i))
    )
  },
  logical(1)
)
if (!all(fit_reference_current)) {
  stop(
    "One or more pooled ERFs use an outdated reference. Rerun ",
    "03.1-fit-user-route-pooled-erf.R before this script."
  )
}

# --------------- city-specific four hottest consecutive months --------------

message(
  "Deriving city-specific four hottest consecutive months ",
  "from observed climate..."
)

health_codes <- sort(unique(unlist(lapply(
  fit_files,
  function(path) {
    fit <- readRDS(path)
    as.character(fit$crossover$code)
  }
))))

observed <- as.data.table(readRDS(USER_PATH_OBS))
observed[, `:=`(
  code = user_pad_code(id),
  date = as.Date(date),
  tmean = as.numeric(tmean)
)]
observed <- observed[
  date >= USER_START_DATE &
    date <= USER_END_DATE,
  .(code, date, tmean)
]
observed[, province_code := substr(code, 1L, 2L)]

province_daily <- observed[, .(
  tmean_province = mean(tmean, na.rm = TRUE)
), by = .(province_code, date)]
national_daily <- observed[, .(
  tmean_national = mean(tmean, na.rm = TRUE)
), by = date]

observed_panel <- CJ(
  code = health_codes,
  date = USER_DATES,
  unique = TRUE
)
observed_panel[, province_code := substr(code, 1L, 2L)]
observed_panel <- merge(
  observed_panel,
  observed[, .(code, date, tmean)],
  by = c("code", "date"),
  all.x = TRUE
)
observed_panel <- merge(
  observed_panel,
  province_daily,
  by = c("province_code", "date"),
  all.x = TRUE
)
observed_panel <- merge(
  observed_panel,
  national_daily,
  by = "date",
  all.x = TRUE
)
observed_panel[, tmean := fcoalesce(
  tmean,
  tmean_province,
  tmean_national
)]
if (any(!is.finite(observed_panel$tmean))) {
  stop("Observed-temperature completion failed for one or more city-days.")
}
observed_panel[, month := as.integer(format(date, "%m"))]

monthly_temperature <- observed_panel[, .(
  observed_monthly_mean_temperature = mean(tmean)
), by = .(code, month)]

if (
  monthly_temperature[, uniqueN(month), by = code][
    ,
    any(V1 != 12L)
  ] ||
    any(!is.finite(
      monthly_temperature$observed_monthly_mean_temperature
    ))
) {
  stop(
    "Twelve finite observed monthly temperatures are required for every city."
  )
}

# Evaluate all 12 cyclic four-month windows for every city. The window with
# the highest mean observed temperature is retained. If two windows have the
# same mean, the one with the earlier starting month is selected
# deterministically.
window_membership <- CJ(
  code = unique(monthly_temperature$code),
  warm_start_month = seq_len(12L),
  warm_rank = seq_len(WARM_MONTH_COUNT),
  unique = TRUE
)
window_membership[, month := (
  warm_start_month + warm_rank - 2L
) %% 12L + 1L]
window_membership <- merge(
  window_membership,
  monthly_temperature,
  by = c("code", "month"),
  all.x = TRUE
)

candidate_windows <- window_membership[, .(
  warm_window_mean_temperature = mean(
    observed_monthly_mean_temperature
  ),
  warm_end_month = month[which.max(warm_rank)]
), by = .(code, warm_start_month)]
setorder(
  candidate_windows,
  code,
  -warm_window_mean_temperature,
  warm_start_month
)
selected_windows <- candidate_windows[, .SD[1L], by = code]

warm_months <- merge(
  window_membership,
  selected_windows[, .(
    code,
    warm_start_month,
    warm_end_month,
    warm_window_mean_temperature
  )],
  by = c("code", "warm_start_month"),
  all = FALSE
)
setorder(warm_months, code, warm_rank)

warm_month_audit <- warm_months[, .(
  warm_start_month = first(warm_start_month),
  warm_start_month_name = month.abb[first(warm_start_month)],
  warm_end_month = first(warm_end_month),
  warm_end_month_name = month.abb[first(warm_end_month)],
  warm_months = paste(month, collapse = ","),
  warm_month_names = paste(
    month.abb[month],
    collapse = ", "
  ),
  warm_window_mean_temperature = first(
    warm_window_mean_temperature
  ),
  coldest_selected_month_mean = min(
    observed_monthly_mean_temperature
  ),
  hottest_selected_month_mean = max(
    observed_monthly_mean_temperature
  ),
  selected_month_count = .N
), by = code]

if (
  any(warm_month_audit$selected_month_count != WARM_MONTH_COUNT) ||
    anyDuplicated(warm_months[, .(code, month)])
) {
  stop("Each city must have exactly four unique selected months.")
}
fwrite(warm_month_audit, CITY_MONTH_FILE)

warm_month_list <- split(
  warm_months$month,
  warm_months$code
)

rm(
  observed,
  province_daily,
  national_daily,
  observed_panel,
  monthly_temperature,
  window_membership,
  candidate_windows,
  selected_windows
)
invisible(gc())

# -------------------------- city-level estimator -----------------------------

warm4_direct_city <- function(
    factual_temperature,
    histnat_temperature,
    cases,
    dates,
    warm_month_numbers,
    fit
) {
  total_cases <- sum(cases, na.rm = TRUE)
  forward_cases <- acc_forward_weights(cases, fit$lag)
  x_factual <- acc_one_basis(
    factual_temperature,
    fit$argvar
  )
  x_histnat <- acc_one_basis(
    histnat_temperature,
    fit$argvar
  )
  x_delta <- x_factual - x_histnat
  x_reference <- as.numeric(acc_one_basis(
    as.numeric(fit$cen),
    fit$argvar
  ))
  date_month <- as.integer(format(as.Date(dates), "%m"))

  valid <- is.finite(forward_cases) &
    is.finite(factual_temperature) &
    is.finite(histnat_temperature) &
    rowSums(!is.finite(x_delta)) == 0L

  empty_result <- function() {
    data.table(
      total_cases = total_cases,
      included_forward_cases = 0,
      warm4_case_weight = 0,
      warm4_city_days = 0L,
      an_factual_temperature_warm4 = 0,
      an_histnat_temperature_warm4 = 0,
      an_acc_temperature_warm4 = 0,
      ac_temperature_warm4 = NA_real_,
      an_factual_temperature_warm4_signed = 0,
      an_histnat_temperature_warm4_signed = 0,
      an_acc_direct_warm4 = 0,
      af_acc_direct_warm4 = NA_real_,
      af_acc_direct_warm4_all_cases = 0
    )
  }
  if (!any(valid) || total_cases <= 0) {
    return(empty_result())
  }

  forward_total <- sum(forward_cases[valid], na.rm = TRUE)
  if (!is.finite(forward_total) || forward_total <= 0) {
    return(empty_result())
  }

  weights <- forward_cases[valid] * total_cases / forward_total
  warm <- date_month[valid] %in% warm_month_numbers
  warm4_case_weight <- sum(weights[warm], na.rm = TRUE)
  if (!any(warm) || warm4_case_weight <= 0) {
    return(empty_result())
  }

  delta_eta <- drop(
    x_delta[valid, , drop = FALSE] %*%
      fit$reduced_coef
  )
  centred_factual <- sweep(
    x_factual[valid, , drop = FALSE],
    2L,
    x_reference,
    "-"
  )
  centred_histnat <- sweep(
    x_histnat[valid, , drop = FALSE],
    2L,
    x_reference,
    "-"
  )
  eta_factual <- drop(
    centred_factual %*% fit$reduced_coef
  )
  eta_histnat <- drop(
    centred_histnat %*% fit$reduced_coef
  )
  daily_an_factual <- weights * (
    1 - exp(-eta_factual)
  )
  daily_an_histnat <- weights * (
    1 - exp(-eta_histnat)
  )
  an_factual_temperature <- sum(
    pmax(daily_an_factual[warm], 0),
    na.rm = TRUE
  )
  an_histnat_temperature <- sum(
    pmax(daily_an_histnat[warm], 0),
    na.rm = TRUE
  )
  an_acc_temperature <-
    an_factual_temperature -
      an_histnat_temperature
  daily_direct_cases <- weights * (
    1 - exp(-delta_eta)
  )
  an_acc_direct_warm4 <- sum(
    daily_direct_cases[warm],
    na.rm = TRUE
  )

  data.table(
    total_cases = total_cases,
    included_forward_cases = sum(weights),
    warm4_case_weight = warm4_case_weight,
    warm4_city_days = sum(warm),
    an_factual_temperature_warm4 =
      an_factual_temperature,
    an_histnat_temperature_warm4 =
      an_histnat_temperature,
    an_acc_temperature_warm4 =
      an_acc_temperature,
    ac_temperature_warm4 = if (
      an_factual_temperature >
        sqrt(.Machine$double.eps)
    ) {
      100 * an_acc_temperature /
        an_factual_temperature
    } else {
      NA_real_
    },
    an_factual_temperature_warm4_signed = sum(
      daily_an_factual[warm],
      na.rm = TRUE
    ),
    an_histnat_temperature_warm4_signed = sum(
      daily_an_histnat[warm],
      na.rm = TRUE
    ),
    an_acc_direct_warm4 = an_acc_direct_warm4,
    af_acc_direct_warm4 =
      100 * an_acc_direct_warm4 / warm4_case_weight,
    af_acc_direct_warm4_all_cases =
      100 * an_acc_direct_warm4 / total_cases
  )
}

# ---------------------------- GCM calculation --------------------------------

cache_file <- function(gcm) {
  file.path(
    CACHE_DIR,
    paste0(
      "warm4_",
      "temperature_burden_",
      user_safe_gcm(gcm),
      ".rds"
    )
  )
}

calculate_one_gcm <- function(gcm) {
  climate_file <- user_pair_file(gcm)
  cache_path <- cache_file(gcm)
  metadata <- list(
    analysis_version = ANALYSIS_VERSION,
    baseline = BASELINE,
    warm_month_count = WARM_MONTH_COUNT,
    warm_month_file_md5 = acc_file_md5(CITY_MONTH_FILE),
    climate_md5 = acc_file_md5(climate_file),
    fit_md5 = vapply(
      fit_files,
      acc_file_md5,
      character(1)
    )
  )

  if (
    REUSE_CACHE &&
      acc_cache_is_current(cache_path, metadata)
  ) {
    message("Reusing warm-season cache: ", gcm)
    return(readRDS(cache_path)$result)
  }

  message("Calculating warm-season direct AF: ", gcm)
  pair <- readRDS(climate_file)
  disease_results <- vector(
    "list",
    length(fit_files)
  )

  for (ii in seq_along(fit_files)) {
    tag_i <- names(fit_files)[ii]
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

    city_rows <- analysis_data[, {
      months_i <- warm_month_list[[code[1L]]]
      if (
        is.null(months_i) ||
          length(months_i) != WARM_MONTH_COUNT
      ) {
        stop("Missing four-month definition for city ", code[1L])
      }
      warm_i <- warm4_direct_city(
        factual_temperature = tf,
        histnat_temperature = tn,
        cases = n,
        dates = date,
        warm_month_numbers = months_i,
        fit = fit
      )
      reference_i <- reference_ac_city(
        factual_temperature = tf,
        histnat_temperature = tn,
        cases = n,
        fit = fit
      )
      primary_i <- select_disease_relevant_ac(
        reference_result = reference_i,
        warm4_result = data.table(
          an_factual_warm4 =
            warm_i$an_factual_temperature_warm4,
          an_histnat_warm4 =
            warm_i$an_histnat_temperature_warm4
        ),
        class_tag = tag_i
      )
      cbind(
        warm_i,
        primary_i[, .(
          primary_burden_definition,
          an_factual_primary,
          an_histnat_primary,
          an_acc_primary,
          ac_primary
        )]
      )
    }, by = code]

    city_rows <- merge(
      city_rows,
      warm_month_audit,
      by = "code",
      all.x = TRUE
    )
    city_rows[, `:=`(
      tag = tag_i,
      classification =
        unname(USER_DISEASE_LABELS[tag_i]),
      gcm = gcm,
      lag = as.integer(fit$lag),
      baseline = BASELINE,
      warm_month_definition =
        paste0(
          "City-specific four consecutive calendar months with the highest ",
          "mean observed temperature"
        ),
      analysis_version = ANALYSIS_VERSION
    )]
    disease_results[[ii]] <- city_rows
  }

  result <- rbindlist(
    disease_results,
    use.names = TRUE,
    fill = TRUE
  )
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

national_gcm <- city_gcm[, .(
  total_cases = sum(total_cases),
  included_forward_cases = sum(included_forward_cases),
  warm4_case_weight = sum(warm4_case_weight),
  warm4_city_days = sum(warm4_city_days),
  an_factual_temperature_warm4 =
    sum(an_factual_temperature_warm4),
  an_histnat_temperature_warm4 =
    sum(an_histnat_temperature_warm4),
  an_factual_temperature_warm4_signed =
    sum(an_factual_temperature_warm4_signed),
  an_histnat_temperature_warm4_signed =
    sum(an_histnat_temperature_warm4_signed),
  an_acc_direct_warm4 = sum(an_acc_direct_warm4),
  an_factual_primary = sum(an_factual_primary),
  an_histnat_primary = sum(an_histnat_primary),
  primary_burden_definition =
    first(primary_burden_definition)
), by = .(
  tag,
  classification,
  gcm,
  lag,
  baseline,
  warm_month_definition,
  analysis_version
)]
national_gcm[, `:=`(
  an_acc_temperature_warm4 =
    an_factual_temperature_warm4 -
      an_histnat_temperature_warm4,
  an_acc_primary =
    an_factual_primary - an_histnat_primary,
  ac_temperature_warm4 = fifelse(
    an_factual_temperature_warm4 >
      sqrt(.Machine$double.eps),
    100 * (
      an_factual_temperature_warm4 -
        an_histnat_temperature_warm4
    ) / an_factual_temperature_warm4,
    NA_real_
  ),
  af_acc_direct_warm4 =
    100 * an_acc_direct_warm4 / warm4_case_weight,
  af_acc_direct_warm4_all_cases =
    100 * an_acc_direct_warm4 / total_cases
)]
national_gcm[, ac_primary := fifelse(
  abs(an_factual_primary) >
    sqrt(.Machine$double.eps),
  100 * an_acc_primary / an_factual_primary,
  NA_real_
)]

if (
  any(!is.finite(national_gcm$ac_temperature_warm4)) ||
  any(!is.finite(national_gcm$af_acc_direct_warm4)) ||
    any(national_gcm$warm4_case_weight <= 0)
) {
  stop("Warm-season national AF contains invalid values.")
}

# --------------------------- GCM-only uncertainty ----------------------------

bootstrap_indices <- acc_bootstrap_indices(
  gcm_count = length(USER_GCMS),
  repetitions = ACC_BOOTSTRAP_REPETITIONS,
  seed = ACC_BOOTSTRAP_SEED + 13L
)

summarise_warm4 <- function(
    acc_cases,
    warm_cases,
    total_cases
) {
  acc_cases <- as.numeric(acc_cases)
  warm_cases <- as.numeric(warm_cases)
  total_cases <- as.numeric(total_cases)
  if (
    length(acc_cases) != ncol(bootstrap_indices) ||
      any(!is.finite(acc_cases)) ||
      any(!is.finite(warm_cases)) ||
      any(!is.finite(total_cases))
  ) {
    stop("One finite warm-season estimate is required per GCM.")
  }

  derive <- function(index) {
    acc_i <- mean(acc_cases[index])
    warm_i <- mean(warm_cases[index])
    total_i <- mean(total_cases[index])
    c(
      acc_cases = acc_i,
      warm_case_weight = warm_i,
      af_warm4 = 100 * acc_i / warm_i,
      af_all_cases = 100 * acc_i / total_i
    )
  }

  point <- derive(seq_along(acc_cases))
  simulated <- t(apply(
    bootstrap_indices,
    1L,
    derive
  ))
  interval <- apply(
    simulated,
    2L,
    quantile,
    probs = c(0.025, 0.975),
    na.rm = TRUE,
    names = FALSE
  )
  list(
    point = point,
    lower = interval[1L, ],
    upper = interval[2L, ]
  )
}

national_summary <- national_gcm[, {
  temperature_burden <- reference_ac_bootstrap_summary(
    factual = an_factual_temperature_warm4,
    histnat = an_histnat_temperature_warm4,
    total_cases = total_cases,
    bootstrap_indices = bootstrap_indices
  )
  primary_burden <- reference_ac_bootstrap_summary(
    factual = an_factual_primary,
    histnat = an_histnat_primary,
    total_cases = total_cases,
    bootstrap_indices = bootstrap_indices
  )
  result <- summarise_warm4(
    acc_cases = an_acc_direct_warm4,
    warm_cases = warm4_case_weight,
    total_cases = total_cases
  )
  list(
    ensemble_method =
      "Ratio of mean paired burdens across 10 GCMs",
    an_factual_temperature_warm4 =
      unname(temperature_burden$point["factual"]),
    an_factual_temperature_warm4_low =
      unname(temperature_burden$lower["factual"]),
    an_factual_temperature_warm4_high =
      unname(temperature_burden$upper["factual"]),
    an_histnat_temperature_warm4 =
      unname(temperature_burden$point["histnat"]),
    an_histnat_temperature_warm4_low =
      unname(temperature_burden$lower["histnat"]),
    an_histnat_temperature_warm4_high =
      unname(temperature_burden$upper["histnat"]),
    an_acc_temperature_warm4 =
      unname(temperature_burden$point["acc"]),
    an_acc_temperature_warm4_low =
      unname(temperature_burden$lower["acc"]),
    an_acc_temperature_warm4_high =
      unname(temperature_burden$upper["acc"]),
    ac_temperature_warm4 =
      unname(temperature_burden$point["ac"]),
    ac_temperature_warm4_low =
      unname(temperature_burden$lower["ac"]),
    ac_temperature_warm4_high =
      unname(temperature_burden$upper["ac"]),
    ac_temperature_warm4_gcm_median =
      median(ac_temperature_warm4),
    ac_temperature_warm4_gcm_min =
      min(ac_temperature_warm4),
    ac_temperature_warm4_gcm_max =
      max(ac_temperature_warm4),
    primary_burden_definition =
      first(primary_burden_definition),
    an_factual_primary =
      unname(primary_burden$point["factual"]),
    an_factual_primary_low =
      unname(primary_burden$lower["factual"]),
    an_factual_primary_high =
      unname(primary_burden$upper["factual"]),
    an_histnat_primary =
      unname(primary_burden$point["histnat"]),
    an_histnat_primary_low =
      unname(primary_burden$lower["histnat"]),
    an_histnat_primary_high =
      unname(primary_burden$upper["histnat"]),
    an_acc_primary =
      unname(primary_burden$point["acc"]),
    an_acc_primary_low =
      unname(primary_burden$lower["acc"]),
    an_acc_primary_high =
      unname(primary_burden$upper["acc"]),
    ac_primary =
      unname(primary_burden$point["ac"]),
    ac_primary_low =
      unname(primary_burden$lower["ac"]),
    ac_primary_high =
      unname(primary_burden$upper["ac"]),
    an_acc_direct_warm4 =
      unname(result$point["acc_cases"]),
    an_acc_direct_warm4_low =
      unname(result$lower["acc_cases"]),
    an_acc_direct_warm4_high =
      unname(result$upper["acc_cases"]),
    warm4_case_weight =
      unname(result$point["warm_case_weight"]),
    af_acc_direct_warm4 =
      unname(result$point["af_warm4"]),
    af_acc_direct_warm4_low =
      unname(result$lower["af_warm4"]),
    af_acc_direct_warm4_high =
      unname(result$upper["af_warm4"]),
    af_acc_direct_warm4_gcm_median =
      median(af_acc_direct_warm4),
    af_acc_direct_warm4_gcm_min =
      min(af_acc_direct_warm4),
    af_acc_direct_warm4_gcm_max =
      max(af_acc_direct_warm4),
    af_acc_direct_warm4_all_cases =
      unname(result$point["af_all_cases"]),
    af_acc_direct_warm4_all_cases_low =
      unname(result$lower["af_all_cases"]),
    af_acc_direct_warm4_all_cases_high =
      unname(result$upper["af_all_cases"]),
    gcm_count = .N
  )
}, by = .(
  tag,
  classification,
  lag,
  baseline,
  warm_month_definition,
  analysis_version
)]

# Add the existing all-year direct AF only as a comparator. It is not used in
# the warm-season calculation.
all_year_file <- file.path(
  USER_OUT,
  "direct_counterfactual_af",
  "V01-direct-ACC-AF-national-summary.csv"
)
if (file.exists(all_year_file)) {
  all_year <- fread(all_year_file)[, .(
    tag,
    all_year_direct_acc_af = af_acc_direct_estimate,
    all_year_direct_acc_af_low = af_acc_direct_lower,
    all_year_direct_acc_af_high = af_acc_direct_upper
  )]
  national_summary <- merge(
    national_summary,
    all_year,
    by = "tag",
    all.x = TRUE
  )
}

national_summary[, tag_order__ := match(
  tag,
  USER_SPECS$tag
)]
setorder(national_summary, tag_order__)
national_summary[, tag_order__ := NULL]

# -------------------------------- outputs ------------------------------------

fwrite(city_gcm, CITY_GCM_FILE)
fwrite(national_gcm, NATIONAL_GCM_FILE)
fwrite(national_summary, SUMMARY_FILE)

message("\nWarm-season temperature-burden attribution results:")
print(
  national_summary[, .(
    classification,
    warm4_month_definition = warm_month_definition,
    factual_temperature_cases =
      round(an_factual_temperature_warm4),
    anthropogenic_temperature_cases =
      round(an_acc_temperature_warm4),
    anthropogenic_share_of_temperature_burden = sprintf(
      "%+.1f%% (%+.1f to %+.1f)",
      ac_primary,
      ac_primary_low,
      ac_primary_high
    ),
    validation_direct_acc_af = sprintf(
      "%+.1f%% (%+.1f to %+.1f)",
      af_acc_direct_warm4,
      af_acc_direct_warm4_low,
      af_acc_direct_warm4_high
    )
  )]
)
message("\nSaved:")
message("  ", CITY_MONTH_FILE)
message("  ", CITY_GCM_FILE)
message("  ", NATIONAL_GCM_FILE)
message("  ", SUMMARY_FILE)
