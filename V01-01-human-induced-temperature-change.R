###############################################################################
## V01 Figure 1 | Human-induced changes in temperature exposure across China
##
## A  Factual and natural-only counterfactual annual temperature
## B  Annual factual-minus-hist-nat temperature difference
## C  City-level mean anthropogenic temperature change
## D  Distribution of city-level anthropogenic warming across seven regions
##
## Important interpretation:
## - Temperatures in each city are already population-weighted over constituent
##   climate grids.
## - National summaries give equal weight to study cities; they are not weighted
##   again by city population.
## - Bands and error bars show the full range across 10 GCMs, not a 95% CI.
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
FIG_DIR <- file.path(USER_OUT, "figures", "V01")

DATA_DIR <- file.path(USER_OUT, "figure1_climate", "V01")

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)

ANNUAL_GCM_FILE <- file.path(USER_OUT, "figures",
                            "Figure1N-annual-by-gcm.csv")

CITY_GCM_FILE <- file.path(USER_OUT, "figures",
                          "Figure1N-city-by-gcm.csv")

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

# ------------------------------- summaries -----------------------------------
annual_gcm <- fread(ANNUAL_GCM_FILE)

city_gcm <- fread(CITY_GCM_FILE)

city_gcm[, code := user_pad_code(code)]

city_gcm[, province_code := sprintf("%02d", as.integer(province_code))]

city_gcm[, region := unname(province_region[province_code])]

annual <- annual_gcm[, .(
  factual = mean(factual, na.rm = TRUE),
  factual_low = min(factual, na.rm = TRUE),
  factual_high = max(factual, na.rm = TRUE),
  histnat = mean(histnat, na.rm = TRUE),
  histnat_low = min(histnat, na.rm = TRUE),
  histnat_high = max(histnat, na.rm = TRUE),
  change = mean(change, na.rm = TRUE),
  change_low = min(change, na.rm = TRUE),
  change_high = max(change, na.rm = TRUE)
), by = year]

period_gcm <- annual_gcm[, .(
  period_change = mean(change, na.rm = TRUE)
), by = gcm]

city <- city_gcm[province_code %chin% valid_provinces, .(
  change = mean(change, na.rm = TRUE),
  change_low = min(change, na.rm = TRUE),
  change_high = max(change, na.rm = TRUE),
  gcm_count = .N
), by = .(code, province_code, region)]

city[, region := factor(region, levels = region_order)]

region_gcm <- city_gcm[
  province_code %chin% valid_provinces & !is.na(region),
  .(region_change = mean(change, na.rm = TRUE), n_city = uniqueN(code)),
  by = .(region, gcm)
]

region_summary <- region_gcm[, .(
  mean_change = mean(region_change),
  gcm_low = min(region_change),
  gcm_high = max(region_change),
  n_city = as.integer(round(mean(n_city)))
), by = region]

region_summary[, region := factor(region, levels = region_order)]

region_summary[, label := sprintf(
  "%.2f°C\n(%.2f–%.2f)", mean_change, gcm_low, gcm_high
)]

fwrite(annual, file.path(DATA_DIR, "V01-Figure1-annual-ensemble.csv"))

fwrite(period_gcm, file.path(DATA_DIR, "V01-Figure1-period-change-by-gcm.csv"))

fwrite(city, file.path(DATA_DIR, "V01-Figure1-city-change-ensemble.csv"))

fwrite(region_gcm, file.path(DATA_DIR, "V01-Figure1-region-change-by-gcm.csv"))

fwrite(region_summary, file.path(DATA_DIR, "V01-Figure1-region-change-summary.csv"))
