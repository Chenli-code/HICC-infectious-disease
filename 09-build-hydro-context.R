# -*- coding: UTF-8 -*-
###############################################################################
## City-level hydrological context for the sanitation module
##
## A sanitation system is loaded hardest when water moves through it fastest:
## sustained wet-season rainfall and heavy-rain days drive combined-sewer
## overflow, surface runoff and the mobilisation of stored faecal material.
## This script builds those indicators from the observed city-day weather
## series already used by the main analysis.
##
##   Source: 24-wildfire/code/data/China_493_cities_total_fire_PM25_O3_
##           weather_2007_2019.rds (ERA5-derived city-day tmean, RH, rainfall)
##
##   warm_season_rainfall  mean daily rainfall in June-September (mm)
##   annual_rainfall       mean annual total rainfall (mm)
##   heavy_rain_days       days per year with rainfall of at least 25 mm
##   very_heavy_rain_days  days per year with rainfall of at least 50 mm
##   wet_day_fraction      fraction of days with rainfall of at least 0.1 mm
##   warm_season_humidity  mean relative humidity in June-September (%)
##
##   hydro_stress = mean of z(log1p(warm_season_rainfall)) and
##                  z(log1p(heavy_rain_days))
##
## Surface runoff and soil moisture would be the more direct quantities, but
## they require an ERA5-Land download that is not in the project; this is
## stated as a limitation rather than approximated.
###############################################################################

options(stringsAsFactors = FALSE, scipen = 999)

suppressPackageStartupMessages({
  library(data.table)
})

project_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
output_dir <- file.path(project_root, "Mechanism", "output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(project_root, "00-paths.R"))
observed_file <- PATH_OBS_WEATHER

start_date <- as.Date("2008-01-01")
end_date <- as.Date("2019-12-31")
warm_months <- 6:9
heavy_rain_threshold <- 25
very_heavy_rain_threshold <- 50

# The project pads administrative codes on the right, so that a four-digit
# prefecture code becomes the six-digit code_merge used everywhere else.
pad_code <- function(x) {
  x <- sub("[.]0+$", "", as.character(x))
  x <- sprintf("%-6s", x)
  gsub(" ", "0", x)
}

standardise <- function(x) {
  usable <- is.finite(x)
  result <- rep(NA_real_, length(x))
  result[usable] <- (x[usable] - mean(x[usable])) / sd(x[usable])
  result
}

message("Reading observed city-day weather ...")

observed <- as.data.table(readRDS(observed_file))
message("Columns: ", paste(names(observed), collapse = ", "))

observed <- observed[, .(
  code = pad_code(id),
  date = as.Date(date),
  rainfall = as.numeric(rainfall),
  relative_humidity = as.numeric(RH),
  tmean = as.numeric(tmean)
)]
observed <- observed[date >= start_date & date <= end_date]
observed[, `:=`(
  month = as.integer(format(date, "%m")),
  year = as.integer(format(date, "%Y"))
)]

n_years <- observed[, uniqueN(year)]

hydro <- observed[, .(
  warm_season_rainfall = mean(rainfall[month %in% warm_months], na.rm = TRUE),
  annual_rainfall = sum(rainfall, na.rm = TRUE) / n_years,
  heavy_rain_days = sum(rainfall >= heavy_rain_threshold, na.rm = TRUE) /
    n_years,
  very_heavy_rain_days = sum(rainfall >= very_heavy_rain_threshold,
                             na.rm = TRUE) / n_years,
  wet_day_fraction = mean(rainfall >= 0.1, na.rm = TRUE),
  rainfall_p95 = as.numeric(quantile(rainfall, 0.95, na.rm = TRUE)),
  warm_season_humidity = mean(
    relative_humidity[month %in% warm_months], na.rm = TRUE
  ),
  n_days = .N
), by = code]

hydro <- hydro[is.finite(warm_season_rainfall) & is.finite(heavy_rain_days)]

hydro[, `:=`(
  warm_rainfall_z = standardise(log1p(warm_season_rainfall)),
  heavy_rain_z = standardise(log1p(heavy_rain_days))
)]
hydro[, hydro_stress := standardise(rowMeans(
  cbind(warm_rainfall_z, heavy_rain_z), na.rm = TRUE
))]

setorder(hydro, code)
fwrite(hydro, file.path(output_dir, "city-hydrological-context.csv"))

message("Cities with hydrological indicators: ", nrow(hydro))
print(summary(hydro[, .(warm_season_rainfall, annual_rainfall,
                        heavy_rain_days, hydro_stress)]))
