###############################################################################
## Figure 4-2 | Geography of disease-relevant attributable contribution (AC)
##
## Respiratory:
##   AC = (AN_factual - AN_histnat) / AN_factual
##   where AN includes the full temperature curve relative to the fitted MMT.
##
## Gastrointestinal/enterovirus and vector-borne:
##   AC uses positive temperature-attributable AN accumulated during each
##   city's four hottest consecutive months.
##
## City and regional percentages are always derived after aggregating factual
## and hist-nat AN. Percentages are never averaged across cities. The ensemble
## point estimate is the ratio of the two 10-GCM mean burdens; 95% uncertainty
## intervals use 10,000 bootstrap resamples of the 10 paired GCM estimates.
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(sf)
  library(scales)
})

ROOT <- Sys.getenv(
  "ATTR_PROJECT_ROOT",
  getwd()
)

source(file.path(ROOT, "03.0-user-route-config.R"))

source(USER_CORE)

source(file.path(ROOT, "V01-12-reference-AC-core.R"))

# -------------------------------- paths --------------------------------------
FIGURE4_DIR <- file.path(USER_ROOT, "Figure4")

dir.create(FIGURE4_DIR, recursive = TRUE, showWarnings = FALSE)

REFERENCE_AC_DIR <- file.path(
  USER_OUT,
  "warm4_attribution"
)

CITY_GCM_FILE <- file.path(
  REFERENCE_AC_DIR,
  "V01-warm4-temperature-burden-city-by-gcm.csv"
)

GCM_BOOT_N <- 10000L

GCM_BOOT_SEED <- 20260725L

disease_order <- c(
  "respiratory", "gastrointestinal", "vectorborne"
)

disease_short <- c(
  respiratory = "Respiratory",
  gastrointestinal = "Gastrointestinal /\nenterovirus",
  vectorborne = "Vector-borne"
)

region_order <- c(
  "Northeast", "North China", "East China", "Central China",
  "South China", "Southwest", "Northwest"
)

region_short <- c(
  "Northeast" = "Northeast (NE)",
  "North China" = "North China (NC)",
  "East China" = "East China (EC)",
  "Central China" = "Central China (CC)",
  "South China" = "South China (SC)",
  "Southwest" = "Southwest (SW)",
  "Northwest" = "Northwest (NW)",
  "National overall" = "National overall"
)

province_to_region <- function(province_code) {
  fcase(
    province_code %chin% c("21", "22", "23"), "Northeast",
    province_code %chin% c("11", "12", "13", "14", "15"), "North China",
    province_code %chin% c(
      "31", "32", "33", "34", "35", "36", "37"
    ), "East China",
    province_code %chin% c("41", "42", "43"), "Central China",
    province_code %chin% c("44", "45", "46"), "South China",
    province_code %chin% c("50", "51", "52", "53", "54"), "Southwest",
    province_code %chin% c("61", "62", "63", "64", "65"), "Northwest",
    default = NA_character_
  )
}

signed_label <- function(x, accuracy = 0.1) {
  paste0(
    ifelse(x > 0, "+", ""),
    number(x, accuracy = accuracy)
  )
}

# ------------------------------- analysis data -------------------------------
city_gcm <- fread(CITY_GCM_FILE)

required_columns <- c(
  "tag", "code", "gcm", "total_cases",
  "an_factual_primary", "an_histnat_primary",
  "an_acc_primary", "primary_burden_definition"
)

if (any(!required_columns %in% names(city_gcm))) {
  stop(
    "Warm-season city file is outdated. Rerun ",
    "V01-13-hottest-four-months-AF.R."
  )
}

city_gcm[, `:=`(
  code = user_pad_code(code),
  province_code = substr(user_pad_code(code), 1L, 2L),
  region = province_to_region(
    substr(user_pad_code(code), 1L, 2L)
  )
)]

city_gcm <- city_gcm[tag %chin% disease_order]

city_gcm[, gcm_order__ := match(gcm, USER_GCMS)]

setorder(city_gcm, tag, code, gcm_order__)

city_gcm[, gcm_order__ := NULL]

bootstrap_indices <- acc_bootstrap_indices(
  gcm_count = length(USER_GCMS),
  repetitions = GCM_BOOT_N,
  seed = GCM_BOOT_SEED
)

# City point estimates are ratios of 10-GCM mean burdens.
city <- city_gcm[, .(
  classification = first(classification),
  province_code = first(province_code),
  region = first(region),
  gcm_count = uniqueN(gcm),
  total_cases = mean(total_cases),
  an_factual_reference = mean(an_factual_primary),
  an_histnat_reference = mean(an_histnat_primary)
), by = .(tag, code)]

city[, `:=`(
  an_acc_reference =
    an_factual_reference - an_histnat_reference,
  ac_reference = fifelse(
    abs(an_factual_reference) >
      sqrt(.Machine$double.eps),
    100 * (
      an_factual_reference - an_histnat_reference
    ) / an_factual_reference,
    NA_real_
  )
)]

# Region-level and national percentages are calculated after summing city AN.
region_gcm <- city_gcm[!is.na(region), .(
  total_cases = sum(total_cases),
  an_factual_reference = sum(an_factual_primary),
  an_histnat_reference = sum(an_histnat_primary)
), by = .(tag, gcm, region)]

region_gcm[, `:=`(
  an_acc_reference =
    an_factual_reference - an_histnat_reference,
  ac_reference = fifelse(
    abs(an_factual_reference) >
      sqrt(.Machine$double.eps),
    100 * (
      an_factual_reference - an_histnat_reference
    ) / an_factual_reference,
    NA_real_
  )
)]

national_gcm <- city_gcm[, .(
  total_cases = sum(total_cases),
  an_factual_reference = sum(an_factual_primary),
  an_histnat_reference = sum(an_histnat_primary)
), by = .(tag, gcm)]

national_gcm[, `:=`(
  region = "National overall",
  an_acc_reference =
    an_factual_reference - an_histnat_reference
)]

national_gcm[, ac_reference := fifelse(
  abs(an_factual_reference) >
    sqrt(.Machine$double.eps),
  100 * an_acc_reference / an_factual_reference,
  NA_real_
)]

forest_gcm <- rbindlist(
  list(region_gcm, national_gcm),
  use.names = TRUE,
  fill = TRUE
)

forest_gcm[, gcm_order__ := match(gcm, USER_GCMS)]

setorder(forest_gcm, tag, region, gcm_order__)

forest_gcm[, gcm_order__ := NULL]

forest_summary <- forest_gcm[, {
  result <- reference_ac_bootstrap_summary(
    factual = an_factual_reference,
    histnat = an_histnat_reference,
    total_cases = total_cases,
    bootstrap_indices = bootstrap_indices
  )
  list(
    estimate = unname(result$point["ac"]),
    ci_low = unname(result$lower["ac"]),
    ci_high = unname(result$upper["ac"]),
    an_acc = unname(result$point["acc"]),
    factual_an = unname(result$point["factual"]),
    histnat_an = unname(result$point["histnat"])
  )
}, by = .(tag, region)]

forest_summary[, `:=`(
  disease = factor(
    unname(disease_short[tag]),
    levels = unname(disease_short[disease_order])
  ),
  region_label = factor(
    unname(region_short[region]),
    levels = c(
      "National overall",
      "Northwest (NW)", "Southwest (SW)",
      "South China (SC)", "Central China (CC)",
      "East China (EC)", "North China (NC)",
      "Northeast (NE)"
    )
  ),
  is_national = region == "National overall"
)]

forest_summary[, value_label := paste0(
  signed_label(estimate, 0.1), "% (",
  signed_label(ci_low, 0.1), " to ",
  signed_label(ci_high, 0.1), ")"
)]

STUDY_YEARS <- as.integer(format(USER_END_DATE, "%Y")) -
  as.integer(format(USER_START_DATE, "%Y")) + 1L

region_cases <- region_gcm[, .(
  annual_acc_cases = mean(an_acc_reference) / STUDY_YEARS,
  annual_gcm_low = min(an_acc_reference) / STUDY_YEARS,
  annual_gcm_high = max(an_acc_reference) / STUDY_YEARS
), by = .(tag, region)]

region_cases[, abs_total := sum(
  abs(annual_acc_cases),
  na.rm = TRUE
), by = tag]

region_cases[, composition_pct := fifelse(
  abs_total > 0,
  100 * abs(annual_acc_cases) / abs_total,
  NA_real_
)]

region_cases[, contribution_label := fifelse(
  composition_pct >= 5,
  paste0(
    fifelse(annual_acc_cases >= 0, "+", "\u2212"),
    number(composition_pct, accuracy = 1),
    "%"
  ),
  ""
)]

region_cases[, `:=`(
  disease = factor(
    unname(disease_short[tag]),
    levels = unname(disease_short[disease_order])
  ),
  region = factor(region, levels = region_order)
)]

fwrite(
  city,
  file.path(
    FIGURE4_DIR,
    "V01-Figure4-2-city-AC-map-data.csv"
  )
)

fwrite(
  forest_summary,
  file.path(
    FIGURE4_DIR,
    "V01-Figure4-2-regional-AC-summary.csv"
  )
)

fwrite(
  region_cases,
  file.path(
    FIGURE4_DIR,
    "V01-Figure4-2-regional-reference-ACC-cases.csv"
  )
)
