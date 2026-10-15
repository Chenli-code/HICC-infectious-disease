###############################################################################
## Figure 5 | Disease-specific burdens attributable to human-induced climate
## change
##
## Each disease is analysed with its own previously fitted pooled ERF and
## disease-specific lag. All anthropogenic effects use the same direct paired
## factual-versus-hist-nat contrast as Figures 3 and 4:
##
##   AF_ACC(t) = 1 - RR_hist-nat(t) / RR_factual(t)
##
## Point estimates are medians across 10 equally weighted GCMs. Percentile
## 95% GCM-bootstrap intervals use 10,000 resamples of the 10 GCM estimates.
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(dlnm)
  library(scales)
})

ROOT <- Sys.getenv(
  "ATTR_PROJECT_ROOT",
  getwd()
)

source(file.path(ROOT, "10.0-Figure3-disease-specification.R"))

source(USER_CORE)

source(file.path(ROOT, "V01-12-reference-AC-core.R"))

# -------------------------------- paths --------------------------------------
FIGURE5_DIR <- file.path(USER_ROOT, "Figure5")

FIGURE5_CACHE <- file.path(FIGURE5_DIR, "direct_by_gcm_cache_v4")

dir.create(FIGURE5_DIR, recursive = TRUE, showWarnings = FALSE)

dir.create(FIGURE5_CACHE, recursive = TRUE, showWarnings = FALSE)

BY_GCM_FILE <- file.path(
  FIGURE5_DIR,
  "V01-Figure5-disease-specific-direct-attribution-by-gcm.csv"
)

SUMMARY_FILE <- file.path(
  FIGURE5_DIR,
  "V01-Figure5-disease-specific-direct-attribution-summary.csv"
)

STATUS_FILE <- file.path(
  FIGURE5_DIR,
  "V01-Figure5-disease-specific-analysis-status.csv"
)

WARM4_MONTH_FILE <- file.path(
  USER_OUT,
  "warm4_attribution",
  "V01-warmest-four-months-by-city.csv"
)

if (!file.exists(WARM4_MONTH_FILE)) {
  stop(
    "Missing city-specific warm-season definition. Rerun ",
    "V01-13-hottest-four-months-AF.R."
  )
}

warm_month_audit <- fread(
  WARM4_MONTH_FILE,
  colClasses = list(character = "code")
)

warm_month_audit[, code := user_pad_code(code)]

warm_month_list <- setNames(
  lapply(
    strsplit(warm_month_audit$warm_months, ",", fixed = TRUE),
    as.integer
  ),
  warm_month_audit$code
)

REUSE_DIRECT <- tolower(
  Sys.getenv("FIGURE5_REUSE_DIRECT", "true")
) == "true"

GCM_BOOT_N <- 10000L

GCM_BOOT_SEED <- 20260723L

MIN_TOTAL_CASES <- 100L

# Disease-level estimates are omitted if the paired direct AF exceeds +/-100%
# in any GCM. Raw results remain available for auditing.
MAX_ABS_AF <- 100

gcm_bootstrap_indices <- acc_bootstrap_indices(
  length(USER_GCMS),
  repetitions = GCM_BOOT_N,
  seed = GCM_BOOT_SEED
)

bootstrap_median <- function(x) {
  result <- acc_bootstrap_median(x, gcm_bootstrap_indices)
  c(
    point = unname(result["estimate"]),
    low = unname(result["lower"]),
    high = unname(result["upper"])
  )
}

cache_file <- function(gcm) {
  file.path(
    FIGURE5_CACHE,
    paste0(
      "V01-direct-disease-ACC-",
      "direct-paired-warm4-ac-v7-",
      user_safe_gcm(gcm),
      ".rds"
    )
  )
}

evaluate_disease_direct <- function(
  fit,
  climate,
  disease,
  disease_tag,
  classification,
  class_tag,
  gcm
) {
  counts <- as.data.table(fit$sparse_counts)[, .(
    code = as.character(code),
    date = as.Date(date),
    n = as.numeric(n)
  )]
  dat <- climate[code %in% as.character(fit$health_codes)]
  dat <- merge(
    dat,
    counts,
    by = c("code", "date"),
    all.x = TRUE
  )
  dat[is.na(n), n := 0]
  setorder(dat, code, date)
  city_results <- dat[, {
    acc_direct_city(
      factual_temperature = tf,
      histnat_temperature = tn,
      cases = n,
      fit = fit
    )
  }, by = code]
  total_cases <- sum(city_results$total_cases)
  included_forward_cases <- sum(
    city_results$included_forward_cases
  )
  an_acc_direct <- sum(city_results$an_acc_direct)
  counterfactual_total_cases <- sum(
    city_results$counterfactual_cases
  )
  af_acc_direct <- 100 * an_acc_direct / total_cases
  acc_check_identity(
    total_cases,
    counterfactual_total_cases,
    an_acc_direct,
    af_acc_direct
  )

  # The disease-relevant AC uses a day-of-year climatological baseline.
  # Respiratory diseases retain AN over the full curve relative to the fitted
  # MMT. Gastrointestinal and vector-borne diseases use positive AN during each
  # city's four hottest consecutive months.
  reference_weights <- acc_prepare_case_weights(
    cases = counts,
    dates = USER_DATES,
    baseline = "day_of_year"
  )
  reference_data <- merge(
    climate[code %in% as.character(fit$health_codes)],
    reference_weights,
    by = c("code", "date"),
    all.x = TRUE
  )
  reference_data[is.na(n), n := 0]
  setorder(reference_data, code, date)
  reference_city <- reference_data[, {
    months_i <- warm_month_list[[code[1L]]]
    if (
      is.null(months_i) ||
        length(months_i) != USER_WARM_MONTH_COUNT
    ) {
      stop(
        "Missing four-month definition for city ",
        code[1L]
      )
    }
    reference_i <- reference_ac_city(
      factual_temperature = tf,
      histnat_temperature = tn,
      cases = n,
      fit = fit
    )
    warm4_i <- warm4_ac_city(
      factual_temperature = tf,
      histnat_temperature = tn,
      cases = n,
      dates = date,
      warm_month_numbers = months_i,
      fit = fit
    )
    select_disease_relevant_ac(
      reference_result = reference_i,
      warm4_result = warm4_i,
      class_tag = class_tag
    )
  }, by = code]
  reference_national <- reference_city[, .(
    total_cases_reference = sum(total_cases),
    an_factual_reference = sum(an_factual_primary),
    an_histnat_reference = sum(an_histnat_primary),
    an_acc_reference = sum(an_acc_primary),
    reference_burden_definition =
      first(primary_burden_definition)
  )]
  reference_national[, ac_reference := fifelse(
    abs(an_factual_reference) >
      sqrt(.Machine$double.eps),
    100 * an_acc_reference / an_factual_reference,
    NA_real_
  )]

  data.table(
    disease = disease,
    disease_tag = disease_tag,
    classification = classification,
    class_tag = class_tag,
    lag = as.integer(fit$lag),
    reference_temperature =
      as.numeric(fit$cen),
    reference_type =
      fifelse(
        tolower(as.character(fit$reference)) ==
          "mrt",
        "MMT",
        toupper(as.character(fit$reference))
      ),
    gcm = gcm,
    status = "estimated",
    total_cases = total_cases,
    included_forward_cases = included_forward_cases,
    an_acc_direct = an_acc_direct,
    af_acc_direct = af_acc_direct,
    total_cases_reference =
      reference_national$total_cases_reference,
    an_factual_reference =
      reference_national$an_factual_reference,
    an_histnat_reference =
      reference_national$an_histnat_reference,
    an_acc_reference =
      reference_national$an_acc_reference,
    ac_reference =
      reference_national$ac_reference,
    reference_burden_definition =
      reference_national$reference_burden_definition,
    reference_baseline = "day_of_year",
    counterfactual_total_cases =
      counterfactual_total_cases,
    formula_version = ACC_FORMULA_VERSION
  )
}

# ----------------------------- direct calculation ----------------------------
fit_index <- lapply(
  seq_len(nrow(FIG3_DISEASE_SPECS)),
  function(ii) {
    spec <- FIG3_DISEASE_SPECS[ii]
    path <- figure3_fit_file(
      spec$disease_tag
    )
    fit <- if (file.exists(path)) {
      tryCatch(
        readRDS(path),
        error = function(e) NULL
      )
    } else {
      NULL
    }
    list(
      spec = spec,
      fit = fit,
      path = path
    )
  }
)

fit_reference_current <- vapply(
  fit_index,
  function(z) {
    if (is.null(z$fit)) {
      return(TRUE)
    }
    identical(
      tolower(as.character(z$fit$reference)),
      tolower(as.character(z$spec$reference))
    )
  },
  logical(1)
)

if (!all(fit_reference_current)) {
  stop(
    "One or more disease-specific ERFs use an outdated reference. Rerun ",
    "10.1-fit-Figure3-disease-erf.R before Figure 5."
  )
}

all_health_codes <- sort(unique(unlist(
  lapply(fit_index, function(z) {
    if (is.null(z$fit)) {
      character()
    } else {
      as.character(z$fit$health_codes)
    }
  })
)))

fit_md5 <- vapply(
  fit_index,
  function(z) acc_file_md5(z$path),
  character(1)
)

for (gg in seq_along(USER_GCMS)) {
  gcm_i <- USER_GCMS[gg]
  cache_i <- cache_file(gcm_i)
  cache_metadata <- list(
    formula_version = ACC_FORMULA_VERSION,
    climate_md5 = acc_file_md5(user_pair_file(gcm_i)),
    fit_md5 = fit_md5,
    direct_baseline = "observed_daily",
    reference_baseline = "day_of_year",
    warm4_month_file_md5 =
      acc_file_md5(WARM4_MONTH_FILE),
    reference_ac_version =
      "disease-relevant-warm4-v3"
  )
  if (
    REUSE_DIRECT &&
      acc_cache_is_current(cache_i, cache_metadata)
  ) {
    message(sprintf(
      "[%02d/%02d] Reusing direct disease attribution: %s",
      gg, length(USER_GCMS), gcm_i
    ))
    next
  }

  message(sprintf(
    "[%02d/%02d] Loading paired climate: %s",
    gg, length(USER_GCMS), gcm_i
  ))
  pair <- readRDS(user_pair_file(gcm_i))
  climate <- acc_complete_climate_pair(
    pair = pair,
    codes = all_health_codes,
    dates = USER_DATES
  )
  rm(pair)
  invisible(gc())

  result <- vector(
    "list",
    length(fit_index)
  )
  for (ii in seq_along(fit_index)) {
    z <- fit_index[[ii]]
    spec <- z$spec
    fit <- z$fit
    if (is.null(fit)) {
      result[[ii]] <- data.table(
        disease = spec$disease,
        disease_tag = spec$disease_tag,
        classification = spec$classification,
        class_tag = spec$class_tag,
        lag = spec$lag,
        gcm = gcm_i,
        status = "not estimable",
        error = "Disease-specific ERF unavailable"
      )
      next
    }

    message(sprintf(
      "  %02d/%02d %s",
      ii,
      length(fit_index),
      spec$disease
    ))
    result[[ii]] <- tryCatch(
      evaluate_disease_direct(
        fit = fit,
        climate = climate,
        disease = spec$disease,
        disease_tag = spec$disease_tag,
        classification = spec$classification,
        class_tag = spec$class_tag,
        gcm = gcm_i
      ),
      error = function(e) {
        data.table(
          disease = spec$disease,
          disease_tag = spec$disease_tag,
          classification = spec$classification,
          class_tag = spec$class_tag,
          lag = spec$lag,
          gcm = gcm_i,
          status = "failed",
          error = conditionMessage(e)
        )
      }
    )
  }

  result_gcm <- rbindlist(
    result,
    fill = TRUE
  )
  saveRDS(
    list(
      metadata = cache_metadata,
      result = result_gcm
    ),
    cache_i,
    compress = TRUE
  )
  rm(result, result_gcm, climate)
  invisible(gc())
}

available_cache <- vapply(
  USER_GCMS,
  function(gcm_i) {
    cache_i <- cache_file(gcm_i)
    metadata_i <- list(
      formula_version = ACC_FORMULA_VERSION,
      climate_md5 = acc_file_md5(user_pair_file(gcm_i)),
      fit_md5 = fit_md5,
      direct_baseline = "observed_daily",
      reference_baseline = "day_of_year",
      warm4_month_file_md5 =
        acc_file_md5(WARM4_MONTH_FILE),
      reference_ac_version =
        "disease-relevant-warm4-v3"
    )
    acc_cache_is_current(cache_i, metadata_i)
  },
  logical(1)
)

if (!all(available_cache)) {
  stop(
    "Direct disease-attribution cache is incomplete for: ",
    paste(
      USER_GCMS[!available_cache],
      collapse = ", "
    )
  )
}

by_gcm <- rbindlist(
  lapply(
    USER_GCMS,
    function(gcm_i) {
      readRDS(cache_file(gcm_i))$result
    }
  ),
  fill = TRUE
)

by_gcm[
  reference_type == "MRT",
  reference_type := "MMT"
]

by_gcm[, disease_order := match(
  disease,
  FIG3_DISEASE_SPECS$disease
)]

by_gcm[, gcm_order := match(
  gcm,
  USER_GCMS
)]

setorder(
  by_gcm,
  disease_order,
  gcm_order
)

by_gcm[, gcm_order := NULL]

fwrite(by_gcm, BY_GCM_FILE)

# -------------------------- summaries and eligibility ------------------------
estimated <- by_gcm[
  status == "estimated" &
    is.finite(af_acc_direct) &
    is.finite(an_acc_direct) &
    is.finite(an_factual_reference) &
    is.finite(an_histnat_reference) &
    is.finite(ac_reference) &
    is.finite(counterfactual_total_cases)
]

set.seed(GCM_BOOT_SEED)

summary_direct <- estimated[, {
  boot_af <- bootstrap_median(
    af_acc_direct
  )
  boot_cases <- bootstrap_median(
    an_acc_direct
  )
  boot_counterfactual <- bootstrap_median(
    counterfactual_total_cases
  )
  reference_boot <- reference_ac_bootstrap_summary(
    factual = an_factual_reference,
    histnat = an_histnat_reference,
    total_cases = total_cases_reference,
    bootstrap_indices = gcm_bootstrap_indices
  )
  list(
    class_tag = first(class_tag),
    lag = first(lag),
    reference_temperature =
      first(reference_temperature),
    reference_type =
      first(reference_type),
    total_cases = first(total_cases),
    gcm_count = uniqueN(gcm),
    min_counterfactual_cases = min(
      counterfactual_total_cases,
      na.rm = TRUE
    ),
    max_abs_direct_af = max(
      abs(af_acc_direct),
      na.rm = TRUE
    ),
    factual_cases = median(total_cases),
    counterfactual_cases =
      boot_counterfactual["point"],
    counterfactual_cases_low =
      boot_counterfactual["low"],
    counterfactual_cases_high =
      boot_counterfactual["high"],
    af_acc_direct = boot_af["point"],
    af_acc_direct_low = boot_af["low"],
    af_acc_direct_high = boot_af["high"],
    an_acc_direct = boot_cases["point"],
    an_acc_direct_low = boot_cases["low"],
    an_acc_direct_high = boot_cases["high"],
    positive_gcm = sum(
      af_acc_direct > 0,
      na.rm = TRUE
    ),
    negative_gcm = sum(
      af_acc_direct < 0,
      na.rm = TRUE
    ),
    reference_burden_definition =
      first(reference_burden_definition),
    reference_baseline =
      first(reference_baseline),
    an_factual_reference =
      unname(reference_boot$point["factual"]),
    an_histnat_reference =
      unname(reference_boot$point["histnat"]),
    an_acc_reference =
      unname(reference_boot$point["acc"]),
    an_acc_reference_low =
      unname(reference_boot$lower["acc"]),
    an_acc_reference_high =
      unname(reference_boot$upper["acc"]),
    ac_reference =
      unname(reference_boot$point["ac"]),
    ac_reference_low =
      unname(reference_boot$lower["ac"]),
    ac_reference_high =
      unname(reference_boot$upper["ac"])
  )
}, by = .(
  disease,
  disease_tag,
  classification,
  disease_order
)]

status <- merge(
  FIG3_DISEASE_SPECS[, .(
    disease,
    disease_tag,
    classification,
    disease_order
  )],
  summary_direct,
  by = c(
    "disease",
    "disease_tag",
    "classification",
    "disease_order"
  ),
  all.x = TRUE,
  sort = FALSE
)

status[, exclusion_reason := fifelse(
  is.na(gcm_count),
  "Disease-specific model or attribution unavailable",
  fifelse(
    gcm_count < length(USER_GCMS),
    paste0(
      "Only ", gcm_count,
      " of ", length(USER_GCMS),
      " GCM estimates available"
    ),
    fifelse(
      total_cases < MIN_TOTAL_CASES,
      paste0(
        "Fewer than ",
        MIN_TOTAL_CASES,
        " observed cases"
      ),
      fifelse(
        !is.finite(min_counterfactual_cases) |
          min_counterfactual_cases <= 0,
        "Non-positive counterfactual total cases in at least one GCM",
        fifelse(
          !is.finite(max_abs_direct_af) |
            max_abs_direct_af > MAX_ABS_AF,
          paste0(
            "Numerically unstable disease-level estimate: ",
            "|AF| > ", MAX_ABS_AF,
            "% in at least one GCM"
          ),
          ""
        )
      )
    )
  )
)]

status[, eligible := exclusion_reason == ""]

fwrite(status, STATUS_FILE)

plot_data <- status[
  eligible == TRUE
]

if (!nrow(plot_data)) {
  stop(
    "No disease passed the prespecified Figure 5 eligibility rules."
  )
}

setorder(plot_data, disease_order)

# Disease-relevant contributions are reported only when the estimate and both
# uncertainty bounds lie within +/-100%; otherwise "not estimable".
plot_data[, reference_ac_eligible :=
  is.finite(ac_reference) &
    is.finite(ac_reference_low) &
    is.finite(ac_reference_high) &
    abs(ac_reference) <= 100 &
    abs(ac_reference_low) <= 100 &
    abs(ac_reference_high) <= 100]

fwrite(plot_data, SUMMARY_FILE)
