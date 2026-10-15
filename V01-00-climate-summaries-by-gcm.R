###############################################################################
## Prepare Figure 1 climate summaries
##
## City temperatures are population-weighted within cities. National and
## regional summaries are equal-location means across study cities.
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(sf)
})

ROOT <- Sys.getenv(
  "ATTR_PROJECT_ROOT",
  getwd()
)

source(file.path(ROOT, "03.0-user-route-config.R"))

FIG_DIR <- file.path(USER_OUT, "figures")

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

province_names <- c(
  `11` = "Beijing", `12` = "Tianjin", `13` = "Hebei", `14` = "Shanxi",
  `15` = "Inner Mongolia", `21` = "Liaoning", `22` = "Jilin",
  `23` = "Heilongjiang", `31` = "Shanghai", `32` = "Jiangsu",
  `33` = "Zhejiang", `34` = "Anhui", `35` = "Fujian", `36` = "Jiangxi",
  `37` = "Shandong", `41` = "Henan", `42` = "Hubei", `43` = "Hunan",
  `44` = "Guangdong", `45` = "Guangxi", `46` = "Hainan", `50` = "Chongqing",
  `51` = "Sichuan", `52` = "Guizhou", `53` = "Yunnan", `54` = "Tibet",
  `61` = "Shaanxi", `62` = "Gansu", `63` = "Qinghai", `64` = "Ningxia",
  `65` = "Xinjiang"
)

valid_provinces <- names(province_names)

summarise_one_gcm <- function(gcm) {
  x <- as.data.table(readRDS(user_pair_file(gcm)))[,
    .(code = as.character(code), date = as.Date(date),
      tf = as.numeric(tf), tn = as.numeric(tn))
  ]
  x[, `:=`(
    province_code = substr(code, 1L, 2L),
    year = as.integer(format(date, "%Y")),
    warming = tf - tn
  )]
  x <- x[province_code %chin% valid_provinces]

  # Each city value is already population-weighted over grids within that city.
  # National/provincial summaries below are equal-location means, matching the
  # location-based presentation used in the reference attribution figure.
  annual <- x[, .(
    factual = mean(tf, na.rm = TRUE),
    histnat = mean(tn, na.rm = TRUE),
    warming = mean(warming, na.rm = TRUE)
  ), by = year][, gcm := gcm]

  city <- x[, .(
    warming = mean(warming, na.rm = TRUE)
  ), by = .(code, province_code)][, gcm := gcm]

  province <- city[, .(
    warming = mean(warming, na.rm = TRUE),
    locations = .N
  ), by = province_code][, gcm := gcm]

  list(annual = annual, city = city, province = province)
}

all_gcm <- lapply(USER_GCMS, summarise_one_gcm)

annual_gcm <- rbindlist(lapply(all_gcm, `[[`, "annual"))

city_gcm <- rbindlist(lapply(all_gcm, `[[`, "city"))

province_gcm <- rbindlist(lapply(all_gcm, `[[`, "province"))

annual <- annual_gcm[, .(
  factual = median(factual),
  factual_low = min(factual),
  factual_high = max(factual),
  histnat = median(histnat),
  histnat_low = min(histnat),
  histnat_high = max(histnat),
  warming = median(warming),
  warming_low = min(warming),
  warming_high = max(warming)
), by = year]

city <- city_gcm[, .(
  warming = median(warming),
  warming_low = min(warming),
  warming_high = max(warming),
  positive_gcm = sum(warming > 0),
  gcm_count = .N
), by = .(code, province_code)]

province <- province_gcm[, .(
  warming = median(warming),
  warming_low = min(warming),
  warming_high = max(warming),
  locations = as.integer(round(median(locations))),
  positive_gcm = sum(warming > 0),
  gcm_count = .N
), by = province_code]

province[, province := unname(province_names[province_code])]

setorder(province, warming)

province[, province := factor(province, levels = province)]

fwrite(annual_gcm, file.path(FIG_DIR, "Figure1A-annual-temperature-by-gcm.csv"))

fwrite(annual, file.path(FIG_DIR, "Figure1A-annual-temperature-ensemble.csv"))

fwrite(city, file.path(FIG_DIR, "Figure1C-city-warming-ensemble.csv"))

fwrite(province, file.path(FIG_DIR, "Figure1B-province-warming-ensemble.csv"))

# Canonical inputs used by the final V01 Figure 1 script.  Keep the descriptive
# files above for auditing, but expose a stable schema (`change`) so the
# preparation and plotting stages cannot drift apart.
annual_gcm_v01 <- copy(annual_gcm)

setnames(annual_gcm_v01, "warming", "change")

city_gcm_v01 <- copy(city_gcm)

setnames(city_gcm_v01, "warming", "change")

fwrite(
  annual_gcm_v01,
  file.path(FIG_DIR, "Figure1N-annual-by-gcm.csv")
)

fwrite(
  city_gcm_v01,
  file.path(FIG_DIR, "Figure1N-city-by-gcm.csv")
)
