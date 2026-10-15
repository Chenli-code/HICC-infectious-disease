# -*- coding: UTF-8 -*-
# ==============================================================================
# Figure 6 panel F2 | AF response to a uniform +1 deg C warming
#
# F2's subtitle is "AF response to uniform +1 deg C warming". That is not the
# observed delta AF divided by the observed delta T: dividing by each city's own
# warming leaves the GCM spread in delta T inside the answer, and no ratio-based
# reconstruction reproduces the published panel (all give 7/10 GCMs agreeing,
# against 10/10 in the figure). Applying the SAME +1 deg C everywhere removes
# the exposure contrast entirely, so what is left is the response contrast.
#
# For every city, GCM and year this computes
#     response_1c = AF(factual temperature + 1 deg C) - AF(factual temperature)
# using exactly the estimator in 08-gastrointestinal-gdp-sii.R, so F2 and the
# rest of Figure 6 share one exposure-response function and one weighting rule.
#
# Output: SII/output/figure6-objects/figure6-plus1c-response-by-city-gcm.csv
# Runtime is dominated by reading the ten climate-pair files.
# ==============================================================================

options(stringsAsFactors = FALSE, scipen = 999)

suppressPackageStartupMessages({
  library(data.table)
})

if (!dir.exists("SII") || !dir.exists("user_route")) {
  stop("Run from the project root (the folder containing SII/ and user_route/).")
}

SETUP_SCRIPT <- file.path("SII", "08-gastrointestinal-gdp-sii.R")
stopifnot(file.exists(SETUP_SCRIPT))

OUT_DIR <- file.path("SII", "output", "figure6-objects")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 1. Reuse the setup of script 08 without repeating its own AF loop.
#
# Everything up to "af_results <- vector(...)" defines the pooled fit, the case
# counts, the analysis calendar and the acc_* helpers. Evaluating only that part
# gives an identical estimator at a fraction of the cost.
# ------------------------------------------------------------------------------
message("Figure 6 F2 | loading the setup from 08-gastrointestinal-gdp-sii.R")

setup_expressions <- parse(SETUP_SCRIPT)
for (ii in seq_along(setup_expressions)) {
  expression_text <- paste(deparse(setup_expressions[[ii]]), collapse = " ")
  if (grepl("^af_results <- vector", expression_text)) break
  eval(setup_expressions[[ii]], envir = globalenv())
}

required <- c(
  "fit", "counts", "health_codes", "analysis_dates", "analysis_years",
  "gcm_names", "climate_dir", "annual_scenario_af_city",
  "acc_complete_climate_pair", "municipality_provinces"
)
absent <- required[!vapply(required, exists, logical(1L))]
if (length(absent)) {
  stop("Setup did not produce: ", paste(absent, collapse = ", "))
}

WARMING_INCREMENT <- 1.0

# ------------------------------------------------------------------------------
# 2. Per-GCM +1 deg C response
#
# annual_scenario_af_city() returns af_factual and af_counterfactual for the two
# temperature series it is handed. Passing (tf + 1) as the "factual" arm and tf
# as the "counterfactual" arm makes their difference the +1 deg C response.
# ------------------------------------------------------------------------------
response_results <- vector("list", length(gcm_names))

for (gcm_index in seq_along(gcm_names)) {
  current_gcm <- gcm_names[gcm_index]
  message(sprintf(
    "[%02d/%02d] +1 deg C response: %s",
    gcm_index, length(gcm_names), current_gcm
  ))

  climate_pair <- readRDS(
    file.path(climate_dir, paste0("climate_pair_", current_gcm, ".rds"))
  )
  climate <- acc_complete_climate_pair(
    pair = climate_pair,
    codes = health_codes,
    dates = analysis_dates
  )

  daily_data <- merge(climate, counts, by = c("code", "date"), all.x = TRUE)
  daily_data[is.na(n), n := 0]
  setorder(daily_data, code, date)

  # Two baselines, because the choice is not cosmetic. The factual series is
  # bias-corrected and therefore almost identical across GCMs (spread 0.07 deg
  # C), so a +1 deg C applied to it produces ten near-identical curves. The
  # hist-nat series keeps each model's own climate, so the same +1 deg C lands
  # on ten different baselines and the panel retains a visible GCM spread.
  from_factual <- daily_data[, {
    annual_scenario_af_city(
      factual_temperature = tf + WARMING_INCREMENT,
      counterfactual_temperature = tf,
      cases = n, dates = date, fit = fit
    )
  }, by = code]
  setnames(
    from_factual,
    c("af_factual", "af_counterfactual"),
    c("af_plus1c_fact", "af_base_fact")
  )

  from_natural <- daily_data[, {
    annual_scenario_af_city(
      factual_temperature = tn + WARMING_INCREMENT,
      counterfactual_temperature = tn,
      cases = n, dates = date, fit = fit
    )
  }, by = code]
  setnames(
    from_natural,
    c("af_factual", "af_counterfactual"),
    c("af_plus1c_nat", "af_base_nat")
  )

  current <- merge(
    from_factual,
    from_natural[, .(code, year, af_plus1c_nat, af_base_nat)],
    by = c("code", "year")
  )
  current[, gcm := current_gcm]
  response_results[[gcm_index]] <- current
  rm(from_factual, from_natural)

  rm(climate_pair, climate, daily_data, current)
  gc(verbose = FALSE)
}

response_raw <- rbindlist(response_results, use.names = TRUE)
response_raw[, response_1c := af_plus1c_fact - af_base_fact]
response_raw[, response_1c_nat := af_plus1c_nat - af_base_nat]

# ------------------------------------------------------------------------------
# 3. Same administrative harmonisation as script 08
# ------------------------------------------------------------------------------
response_raw[, code := sprintf("%06d", as.integer(code))]
response_raw[, code_merge := fcase(
  substr(code, 1L, 2L) %chin% municipality_provinces &
    substr(code, 5L, 6L) != "00",
  code,
  substr(code, 1L, 2L) %chin% municipality_provinces,
  NA_character_,
  default = paste0(substr(code, 1L, 4L), "00")
)]

response_raw <- response_raw[
  !is.na(code_merge) &
    is.finite(case_weight) & case_weight > 0 &
    is.finite(response_1c) & is.finite(response_1c_nat)
]

# Cases weight the merge of several raw units into one prefecture, matching
# af_gcm in script 08.
response_city <- response_raw[, .(
  case_weight = sum(case_weight),
  response_1c = weighted.mean(response_1c, case_weight),
  response_1c_nat = weighted.mean(response_1c_nat, case_weight)
), by = .(code_merge, year, gcm)]

setorder(response_city, gcm, code_merge, year)

out_file <- file.path(OUT_DIR, "figure6-plus1c-response-by-city-gcm.csv")
fwrite(response_city, out_file)

message("\nWritten: ", out_file)
message(sprintf(
  "cities %d | years %d-%d | GCMs %d | median response %.2f pp per +1 deg C",
  uniqueN(response_city$code_merge),
  min(response_city$year), max(response_city$year),
  uniqueN(response_city$gcm),
  median(response_city$response_1c)
))
