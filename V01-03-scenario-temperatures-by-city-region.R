###############################################################################
## V01 Supplementary Figure
## Regional classification and city temperature exposure under factual and
## natural-only counterfactual climates, 2008–2019
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

# ------------------------------- paths ---------------------------------------
FIG_DIR <- file.path(USER_ROOT, "figure1")

DATA_DIR <- file.path(USER_OUT, "figure1_climate", "V01")

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)

CITY_TEMP_GCM_FILE <- file.path(
  DATA_DIR, "V01-Supp-city-scenario-temperature-by-gcm.csv"
)

region_order <- c(
  "Northeast", "North", "East", "Central",
  "South", "Southwest", "Northwest"
)

province_region <- c(
  `11` = "North", `12` = "North", `13` = "North", `14` = "North", `15` = "North",
  `21` = "Northeast", `22` = "Northeast", `23` = "Northeast",
  `31` = "East", `32` = "East", `33` = "East", `34` = "East",
  `35` = "East", `36` = "East", `37` = "East",
  `41` = "Central", `42` = "Central", `43` = "Central",
  `44` = "South", `45` = "South", `46` = "South",
  `50` = "Southwest", `51` = "Southwest", `52` = "Southwest",
  `53` = "Southwest", `54` = "Southwest",
  `61` = "Northwest", `62` = "Northwest", `63` = "Northwest",
  `64` = "Northwest", `65` = "Northwest"
)

valid_provinces <- names(province_region)

# --------------------------- health-city cohort ------------------------------
# Supplementary city summaries are restricted to the study cities that have
# health outcome records. A city is retained once in the climate-exposure
# cohort even when it contributes to more than one disease class.
health_env <- new.env(parent = emptyenv())

stopifnot(exists("incidence", envir = health_env, inherits = FALSE))

health_inc <- as.data.table(health_env$incidence)

health_inc <- health_inc[
  year <= 2019 &
    code_p != "66" &
    ename != "COVID-19"
]

health_inc[
  ename %chin% c("influenza A H1N1", "avian influenza H1N1"),
  classification := "Respiratory diseases"
]

health_inc[, `:=`(
  code = user_pad_code(code_merge),
  date = as.Date(date)
)]

health_inc <- health_inc[
  date >= USER_START_DATE &
    date <= USER_END_DATE &
    classification %chin% USER_SPECS$classification
]

health_cities <- unique(
  health_inc[, .(
    code,
    province_code = substr(code, 1L, 2L)
  )]
)

health_cities <- health_cities[
  province_code %chin% valid_provinces
]

health_codes <- sort(unique(health_cities$code))

# ------------------------------- climate summary -----------------------------
# The paired files contain population-weighted daily city temperatures for
# factual (tf) and natural-only counterfactual (tn) conditions.
if (!file.exists(CITY_TEMP_GCM_FILE)) {
  message("Building city scenario-temperature summaries from 10 paired GCMs...")
  city_temp_list <- vector("list", length(USER_GCMS))

  for (ii in seq_along(USER_GCMS)) {
    gcm_i <- USER_GCMS[ii]
    message(sprintf("  [%02d/%02d] %s", ii, length(USER_GCMS), gcm_i))

    x <- as.data.table(readRDS(user_pair_file(gcm_i)))[, .(
      code = user_pad_code(code),
      tf = as.numeric(tf),
      tn = as.numeric(tn)
    )]
    x[, province_code := substr(code, 1L, 2L)]

    city_temp_list[[ii]] <- x[
      province_code %chin% valid_provinces,
      .(
        factual = mean(tf, na.rm = TRUE),
        histnat = mean(tn, na.rm = TRUE)
      ),
      by = .(code, province_code)
    ][, gcm := gcm_i]

    rm(x)
    invisible(gc())
  }

  city_temp_gcm <- rbindlist(city_temp_list, use.names = TRUE)
  fwrite(city_temp_gcm, CITY_TEMP_GCM_FILE)
} else {
  city_temp_gcm <- fread(CITY_TEMP_GCM_FILE)
}

city_temp_gcm[, code := user_pad_code(code)]

city_temp_gcm[, province_code := sprintf(
  "%02d", as.integer(province_code)
)]

# Keep one common city cohort: cities with health data in at least one of the
# three analysed disease classes. For a health city absent from the original
# climate-city file, use its province-GCM mean (and the national-GCM mean only
# if the province has no available city). This is the same pragmatic inclusion
# principle used by the health models and avoids restricting the cohort to the
# climate/health intersection.
available_temp <- city_temp_gcm[
  province_code %chin% valid_provinces,
  .(
    factual = mean(factual, na.rm = TRUE),
    histnat = mean(histnat, na.rm = TRUE)
  ),
  by = .(code, province_code, gcm)
]

province_temp <- available_temp[, .(
  factual_province = mean(factual, na.rm = TRUE),
  histnat_province = mean(histnat, na.rm = TRUE)
), by = .(province_code, gcm)]

national_temp <- available_temp[, .(
  factual_national = mean(factual, na.rm = TRUE),
  histnat_national = mean(histnat, na.rm = TRUE)
), by = gcm]

city_temp_gcm <- CJ(
  code = health_codes,
  gcm = USER_GCMS,
  unique = TRUE
)

city_temp_gcm[, province_code := substr(code, 1L, 2L)]

city_temp_gcm <- merge(
  city_temp_gcm,
  available_temp,
  by = c("code", "province_code", "gcm"),
  all.x = TRUE,
  sort = FALSE
)

city_temp_gcm <- merge(
  city_temp_gcm,
  province_temp,
  by = c("province_code", "gcm"),
  all.x = TRUE,
  sort = FALSE
)

city_temp_gcm <- merge(
  city_temp_gcm,
  national_temp,
  by = "gcm",
  all.x = TRUE,
  sort = FALSE
)

city_temp_gcm[, climate_source := fifelse(
  is.finite(factual) & is.finite(histnat),
  "direct_city",
  fifelse(
    is.finite(factual_province) & is.finite(histnat_province),
    "province_gcm_mean",
    "national_gcm_mean"
  )
)]

city_temp_gcm[, `:=`(
  factual = fcoalesce(
    factual,
    factual_province,
    factual_national
  ),
  histnat = fcoalesce(
    histnat,
    histnat_province,
    histnat_national
  )
)]

city_temp_gcm[, c(
  "factual_province",
  "histnat_province",
  "factual_national",
  "histnat_national"
) := NULL]

city_temp_gcm[, region := factor(
  unname(province_region[province_code]), levels = region_order
)]

city_temp <- city_temp_gcm[, .(
  factual = mean(factual, na.rm = TRUE),
  factual_low = min(factual, na.rm = TRUE),
  factual_high = max(factual, na.rm = TRUE),
  histnat = mean(histnat, na.rm = TRUE),
  histnat_low = min(histnat, na.rm = TRUE),
  histnat_high = max(histnat, na.rm = TRUE),
  gcm_count = .N,
  climate_source = if (
    all(climate_source == "direct_city")
  ) {
    "direct_city"
  } else if (
    any(climate_source == "province_gcm_mean")
  ) {
    "province_gcm_mean"
  } else {
    "national_gcm_mean"
  }
), by = .(code, province_code, region)]

region_gcm <- city_temp_gcm[, .(
  factual = mean(factual, na.rm = TRUE),
  histnat = mean(histnat, na.rm = TRUE),
  n_city = uniqueN(code)
), by = .(region, gcm)]

region_factual <- region_gcm[, .(
  mean_temperature = mean(factual),
  gcm_low = min(factual),
  gcm_high = max(factual),
  n_city = as.integer(round(mean(n_city)))
), by = region]

region_factual[, scenario := "Factual"]

region_histnat <- region_gcm[, .(
  mean_temperature = mean(histnat),
  gcm_low = min(histnat),
  gcm_high = max(histnat),
  n_city = as.integer(round(mean(n_city)))
), by = region]

region_histnat[, scenario := "Natural-only counterfactual"]

region_summary <- rbindlist(
  list(region_factual, region_histnat), use.names = TRUE
)

region_summary[, region := factor(region, levels = region_order)]

region_summary[, label := sprintf(
  "%.1f°C\n(%.1f–%.1f)", mean_temperature, gcm_low, gcm_high
)]

fwrite(
  city_temp,
  file.path(DATA_DIR, "V01-Supp-city-scenario-temperature-ensemble.csv")
)

fwrite(
  region_gcm,
  file.path(DATA_DIR, "V01-Supp-region-scenario-temperature-by-gcm.csv")
)

fwrite(
  region_summary,
  file.path(DATA_DIR, "V01-Supp-region-scenario-temperature-summary.csv")
)
