# -*- coding: UTF-8 -*-
# ==============================================================================
# Figure 6 panel E | contrasts across China's two conventional dividing lines
#
# Panel E shows a break that follows a line. Two candidate lines are quantified
# here, using the SAME geometry that 10-figure6-prepare-objects.R draws, so the
# line in the figure and the split in the table cannot disagree.
#
#   Hu Huanyong line (Heihe-Tengchong)  the diagonal one. A straight segment
#       from (98.5 E, 25.0 N) to (127.5 E, 50.2 N). It is the classic
#       population divide: the northwest half holds most of the land and little
#       of the population. This is the boundary the bivariate map follows.
#
#   Qinling-Huaihe line                 the horizontal one, a polyline along the
#       Qinling ridge and the Huai River, 104-121 E. It is a climate boundary
#       (the January 0 C isotherm) and is kept as a secondary split.
#
# Method notes:
#   - Cities are assigned by polygon centroid. A straight line splits the whole
#     plane, so the Hu line classifies every city; for Qinling-Huaihe, cities
#     outside 104-121 E are reported separately because the line stops there.
#   - Every contrast is paired within GCM: the difference is the median of ten
#     within-model differences, not the difference of two medians.
#   - Intervals are 2.5th-97.5th GCM percentiles, an ensemble spread and not a
#     confidence interval.
#
# Output: SII/output/figure6-objects/figure6-dividing-line-*.csv
# ==============================================================================

options(stringsAsFactors = FALSE, scipen = 999)

suppressPackageStartupMessages({
  library(data.table)
  library(sf)
})

if (!dir.exists("SII") || !dir.exists("user_route")) {
  stop("Run from the project root (the folder containing SII/ and user_route/).")
}
sf_use_s2(FALSE)

OUT_DIR <- file.path("SII", "output", "figure6-objects")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

DISEASE <- "Gastrointestinal and enterovirus"
YEARS <- 2008:2019
LINE_WEST <- 104.0
LINE_EAST <- 121.0

source("00-paths.R")
CITY_GEO_DIR <- PATH_CITY_GEOJSON_DIR
OBS_FILE <- PATH_OBS_WEATHER
TEMP_CSV <- file.path(
  "user_route", "output", "figure1_climate", "V01",
  "V01-Supp-city-scenario-temperature-by-gcm.csv"
)

pad6 <- function(x) {
  x <- as.character(x)
  ifelse(nchar(x) < 6L, paste0(x, strrep("0", 6L - nchar(x))), x)
}
q025 <- function(x) quantile(x, 0.025, names = FALSE, na.rm = TRUE)
q975 <- function(x) quantile(x, 0.975, names = FALSE, na.rm = TRUE)

# Same geometry as the drawn lines.
HU_LINE <- data.table(
  lon = c(98.5, 127.5),
  lat = c(25.0,  50.2)
)
QINLING_HUAIHE <- data.table(
  lon = c(104.0, 107.0, 109.5, 111.5, 113.5, 115.5, 117.5, 119.5, 121.0),
  lat = c( 34.3,  33.7,  33.6,  33.3,  32.8,  32.7,  32.8,  33.2,  34.0)
)

# Sign of the cross product of (city - start) with the segment direction.
hu_side <- function(lon, lat) {
  dx <- HU_LINE$lon[2] - HU_LINE$lon[1]
  dy <- HU_LINE$lat[2] - HU_LINE$lat[1]
  cross <- (lon - HU_LINE$lon[1]) * dy - (lat - HU_LINE$lat[1]) * dx
  fifelse(cross > 0, "Southeast", "Northwest")
}

qh_side <- function(lon, lat) {
  inside <- lon >= 104.0 & lon <= 121.0
  line_lat <- approx(QINLING_HUAIHE$lon, QINLING_HUAIHE$lat, xout = lon,
                     rule = 2)$y
  fifelse(!inside, "Outside line extent",
          fifelse(lat > line_lat, "North", "South"))
}

# ------------------------------------------------------------------------------
# 1. City centroids and the north/south label
# ------------------------------------------------------------------------------
geo_files <- list.files(CITY_GEO_DIR, pattern = "[.]json$", full.names = TRUE)
city_geo <- do.call(rbind, lapply(geo_files, function(f) {
  st_read(f, quiet = TRUE)[, c("id", "name")]
}))
city_geo <- suppressWarnings(st_make_valid(st_transform(city_geo, 4326)))
centroids <- suppressWarnings(st_coordinates(st_point_on_surface(st_geometry(city_geo))))

geo <- data.table(
  code_merge = pad6(city_geo$id),
  city = as.character(city_geo$name),
  lon = centroids[, 1],
  lat = centroids[, 2]
)
geo[, hu := hu_side(lon, lat)]
geo[, qh := qh_side(lon, lat)]

# ------------------------------------------------------------------------------
# 2. Outcomes: delta AF, the two F panels, and January temperature as a check
# ------------------------------------------------------------------------------
cached <- readRDS(
  file.path("SII", "output", "quantitative-statistics", "analysis-inputs.rds")
)
a <- as.data.table(cached$analysis_data)[disease == DISEASE & year %in% YEARS]
a[, delta_af := af_factual - af_counterfactual]

city_gcm <- a[, .(
  delta_af = mean(delta_af),
  af_factual = mean(af_factual),
  af_natural = mean(af_counterfactual),
  case_weight = sum(case_weight),
  gdp_rank = socioeconomic_rank[1L]
), by = .(code_merge, gcm)]

temp_gcm <- fread(TEMP_CSV, colClasses = list(character = "code"))
temp_gcm[, code_merge := pad6(code)]
temp_gcm[, delta_t := factual - histnat]

response_file <- file.path(OUT_DIR, "figure6-plus1c-response-by-city-gcm.csv")
if (!file.exists(response_file)) {
  stop("Missing ", response_file, ". Run SII/09-figure6-plus1c-response.R first.")
}
response <- fread(response_file, colClasses = list(character = "code_merge"))
response <- response[, .(response_1c = mean(response_1c_nat)), by = .(code_merge, gcm)]

obs <- as.data.table(readRDS(OBS_FILE))
obs <- obs[, .(code_merge = pad6(id), date = as.Date(date), t = as.numeric(tmean))]
january <- obs[
  date >= as.Date("2008-01-01") & date <= as.Date("2019-12-31") &
    data.table::month(date) == 1L,
  .(january_t = mean(t, na.rm = TRUE)),
  by = code_merge
]

dt <- Reduce(
  function(x, y) merge(x, y, by = intersect(names(x), names(y))),
  list(city_gcm, temp_gcm[, .(code_merge, gcm, delta_t)], response,
       geo[, .(code_merge, city, lon, lat, hu, qh)], january)
)

# ------------------------------------------------------------------------------
# 3. Paired contrast across each line, weighting cities by their case totals
# ------------------------------------------------------------------------------
metrics <- c("delta_af", "delta_t", "response_1c", "af_factual", "af_natural")
units <- c(delta_af = "pp", delta_t = "deg C", response_1c = "pp per +1 deg C",
           af_factual = "%", af_natural = "%")

contrast_for <- function(side_col, a_side, b_side, line_name) {
  panel <- dt[get(side_col) %chin% c(a_side, b_side)]
  by_gcm <- panel[, lapply(.SD, function(v) weighted.mean(v, case_weight)),
                  .SDcols = metrics, by = c("gcm", side_col)]
  rbindlist(lapply(metrics, function(m) {
    w <- dcast(by_gcm, as.formula(paste("gcm ~", side_col)), value.var = m)
    w[, difference := get(a_side) - get(b_side)]
    data.table(
      line = line_name, metric = m, unit = units[[m]],
      side_a = a_side, side_b = b_side,
      value_a = median(w[[a_side]]), value_b = median(w[[b_side]]),
      difference = median(w$difference),
      diff_lower = q025(w$difference), diff_upper = q975(w$difference),
      n_gcm = nrow(w), n_a_higher = sum(w$difference > 0)
    )
  }))
}

contrast <- rbindlist(list(
  contrast_for("hu", "Northwest", "Southeast", "Hu Huanyong line"),
  contrast_for("qh", "North", "South", "Qinling-Huaihe line")
))
contrast[, label := sprintf(
  "%.2f vs %.2f; difference %+.2f [%.2f, %.2f] %s; %d/%d GCMs",
  value_a, value_b, difference, diff_lower, diff_upper, unit,
  n_a_higher, n_gcm
)]

# ------------------------------------------------------------------------------
# 4. Composition: is this just the GDP axis in disguise?
# ------------------------------------------------------------------------------
composition_for <- function(side_col, line_name) {
  out <- dt[, .(
    n_cities = uniqueN(code_merge),
    cases = sum(case_weight) / uniqueN(gcm),
    mean_gdp_rank = weighted.mean(gdp_rank, case_weight),
    mean_january_t = weighted.mean(january_t, case_weight),
    mean_latitude = weighted.mean(lat, case_weight),
    mean_longitude = weighted.mean(lon, case_weight)
  ), by = c(side_col)]
  setnames(out, side_col, "side")
  out[, line := line_name][]
}
composition <- rbindlist(list(
  composition_for("hu", "Hu Huanyong line"),
  composition_for("qh", "Qinling-Huaihe line")
))

city_table <- unique(
  dt[, .(code_merge, city, lon = round(lon, 2), lat = round(lat, 2),
         hu, qh, january_t = round(january_t, 2),
         gdp_rank = round(gdp_rank, 3))],
  by = "code_merge"
)[order(hu, -lat)]

fwrite(contrast, file.path(OUT_DIR, "figure6-dividing-line-contrast.csv"))
fwrite(composition, file.path(OUT_DIR, "figure6-dividing-line-composition.csv"))
fwrite(city_table, file.path(OUT_DIR, "figure6-dividing-line-city-side.csv"))

message("\n--- composition ---")
print(composition[, .(line, side, n_cities,
                      mean_gdp_rank = round(mean_gdp_rank, 3),
                      mean_january_t = round(mean_january_t, 2),
                      mean_latitude = round(mean_latitude, 1),
                      mean_longitude = round(mean_longitude, 1))])
message("\n--- paired contrasts ---")
for (ln in unique(contrast$line)) {
  message("")
  message(ln)
  print(contrast[line == ln, .(metric, label)])
}
message("")
message("Written to: ", OUT_DIR)
