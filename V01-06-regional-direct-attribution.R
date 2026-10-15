###############################################################################
## Figure 4 | Geographical heterogeneity in anthropogenic climate-attributable
## infectious disease burden
##
## Direct factual-versus-hist-nat estimand:
##   AF_ACC = sum[(1 - RR_histnat / RR_factual) * forward cases] /
##            total observed cases * 100
##
## A-C City-level direct ACC AF among all observed cases
## D   Regional direct ACC AF (10-GCM median and 95% GCM-bootstrap interval)
## E   Regional composition of direct ACC-attributable cases, normalised to
##     100% within each disease group
##
## No panel divides by a reference-temperature-based net burden.
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

# -------------------------------- paths --------------------------------------
FIGURE4_DIR <- file.path(USER_ROOT, "Figure4")

dir.create(FIGURE4_DIR, recursive = TRUE, showWarnings = FALSE)

DIRECT_DIR <- file.path(USER_OUT, "direct_counterfactual_af")

CITY_GCM_FILE <- file.path(
  DIRECT_DIR, "V01-direct-ACC-AF-by-city-and-gcm.csv"
)

NATIONAL_GCM_FILE <- file.path(
  DIRECT_DIR, "V01-direct-ACC-AF-national-by-gcm.csv"
)

disease_order <- c(
  "respiratory", "gastrointestinal", "vectorborne"
)

disease_labels <- c(
  respiratory = "Respiratory infectious diseases",
  gastrointestinal = "Gastrointestinal and enterovirus diseases",
  vectorborne = "Vector-borne diseases"
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

national_gcm <- fread(NATIONAL_GCM_FILE)

STUDY_YEARS <- as.integer(
  format(USER_END_DATE, "%Y")
) - as.integer(
  format(USER_START_DATE, "%Y")
) + 1L

GCM_BOOT_N <- 10000L

GCM_BOOT_SEED <- 20260723L

# Percentile bootstrap based only on the 10 GCM point estimates.
gcm_bootstrap_indices <- acc_bootstrap_indices(
  length(USER_GCMS),
  repetitions = GCM_BOOT_N,
  seed = GCM_BOOT_SEED
)

gcm_median_bootstrap <- function(
  x
) {
  result <- acc_bootstrap_median(x, gcm_bootstrap_indices)
  list(
    estimate = unname(result["estimate"]),
    ci_low = unname(result["lower"]),
    ci_high = unname(result["upper"]),
    gcm_min = unname(result["gcm_min"]),
    gcm_max = unname(result["gcm_max"]),
    gcm_count = length(x)
  )
}

city_gcm[, `:=`(
  code = user_pad_code(code),
  province_code = sprintf(
    "%02d", as.integer(province_code)
  ),
  region = province_to_region(
    sprintf("%02d", as.integer(province_code))
  )
)]

city_gcm <- city_gcm[
  tag %chin% disease_order
]

# City point estimates are medians of paired GCM-specific direct AFs.
city <- city_gcm[, .(
  classification = first(classification),
  province_code = first(province_code),
  region = first(region),
  total_cases = first(total_cases),
  gcm_count = uniqueN(gcm),
  an_acc_direct = median(
    an_acc_direct, na.rm = TRUE
  ),
  an_acc_direct_min = min(
    an_acc_direct, na.rm = TRUE
  ),
  an_acc_direct_max = max(
    an_acc_direct, na.rm = TRUE
  ),
  af_acc_direct = median(
    af_acc_direct, na.rm = TRUE
  ),
  af_acc_direct_min = min(
    af_acc_direct, na.rm = TRUE
  ),
  af_acc_direct_max = max(
    af_acc_direct, na.rm = TRUE
  ),
  positive_gcm = sum(an_acc_direct > 0, na.rm = TRUE),
  negative_gcm = sum(an_acc_direct < 0, na.rm = TRUE)
), by = .(tag, code)]

city[, direction := fifelse(
  af_acc_direct >= 0,
  "Cases avoidable without anthropogenic warming",
  "Cases increased under no-anthropogenic-warming"
)]

# Region-level AF is calculated after aggregating city cases within each GCM.
region_gcm <- city_gcm[!is.na(region), .(
  total_cases = sum(total_cases, na.rm = TRUE),
  an_acc_direct = sum(an_acc_direct, na.rm = TRUE)
), by = .(tag, gcm, region)]

region_gcm[, af_acc_direct := (
  100 * an_acc_direct / total_cases
)]

national_forest_gcm <- national_gcm[, .(
  tag, gcm,
  region = "National overall",
  total_cases,
  an_acc_direct,
  af_acc_direct
)]

forest_gcm <- rbindlist(
  list(region_gcm, national_forest_gcm),
  fill = TRUE
)

set.seed(GCM_BOOT_SEED)

forest_summary <- forest_gcm[
  is.finite(af_acc_direct),
  {
    boot <- gcm_median_bootstrap(
      af_acc_direct
    )
    c(
      boot,
      list(cases = median(an_acc_direct))
    )
  },
  by = .(tag, region)
]

forest_summary[, `:=`(
  disease = factor(
    unname(disease_labels[tag]),
    levels = unname(disease_labels[disease_order])
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

# Regional case composition uses annual-average direct paired-GCM case
# differences. Division by the common number of study years does not alter
# the within-disease 100% composition, but makes all exported case quantities
# and annotations interpretable as cases per year.
region_cases <- region_gcm[, .(
  study_period_acc_cases = median(an_acc_direct),
  annual_acc_cases = median(an_acc_direct) / STUDY_YEARS,
  annual_gcm_low = min(an_acc_direct) / STUDY_YEARS,
  annual_gcm_high = max(an_acc_direct) / STUDY_YEARS
), by = .(tag, region)]

region_cases[, abs_total := sum(
  abs(annual_acc_cases), na.rm = TRUE
), by = tag]

region_cases[, abs_contribution := (
  100 * abs(annual_acc_cases) / abs_total
)]

region_cases[, contribution_label := fifelse(
  abs_contribution >= 5,
  paste0(
    fifelse(annual_acc_cases >= 0, "+", "\u2212"),
    number(abs_contribution, accuracy = 1),
    "%"
  ),
  ""
)]

region_cases[, `:=`(
  disease = factor(
    unname(disease_short[tag]),
    levels = unname(disease_short[disease_order])
  ),
  region = factor(region, levels = region_order),
  annual_cases_thousands = annual_acc_cases / 1000,
  composition_pct = abs_contribution
)]

fwrite(
  city,
  file.path(
    FIGURE4_DIR,
    "V01-Figure4-direct-city-map-data.csv"
  )
)

fwrite(
  forest_summary,
  file.path(
    FIGURE4_DIR,
    "V01-Figure4-direct-regional-forest-summary.csv"
  )
)

fwrite(
  region_cases,
  file.path(
    FIGURE4_DIR,
    "V01-Figure4-direct-regional-cases.csv"
  )
)
