###############################################################################
## Input paths
##
## Every script in this repository reads external data through the variables
## defined here. Run all scripts with the repository root as the working
## directory, or set ATTR_PROJECT_ROOT to the repository root.
##
## By default the input data are expected under <repository>/data with the
## layout documented in README.md. To keep the data elsewhere, set
## ATTR_DATA_DIR to that directory.
###############################################################################

ATTR_ROOT <- normalizePath(
  Sys.getenv("ATTR_PROJECT_ROOT", getwd()),
  winslash = "/",
  mustWork = TRUE
)
if (!file.exists(file.path(ATTR_ROOT, "00-paths.R"))) {
  stop(
    "Run from the repository root or set ATTR_PROJECT_ROOT to it: ",
    ATTR_ROOT
  )
}
Sys.setenv(ATTR_PROJECT_ROOT = ATTR_ROOT)

DATA_DIR <- Sys.getenv("ATTR_DATA_DIR", file.path(ATTR_ROOT, "data"))

## Notifiable infectious disease cases (restricted; see Data availability).
PATH_HEALTH <- file.path(DATA_DIR, "health", "DATA.incidence.id.RData")

## Child and adolescent population by city and year (restricted).
PATH_POP_TOTAL <- file.path(DATA_DIR, "population", "pop.total.RData")

## Observed city-level daily weather and PM2.5, population weighted.
PATH_OBS_WEATHER <- file.path(
  DATA_DIR, "weather",
  "China_493_cities_total_fire_PM25_O3_weather_2007_2019.rds"
)

## Bias-adjusted CMIP6 temperatures for 10 GCMs, aggregated to cities.
## The *_pop files are population weighted (primary analysis); the *_grid
## files are simple grid-cell means (sensitivity analysis only).
PATH_CLIMATE_DIR <- file.path(DATA_DIR, "climate")
PATH_GCM_FACTUAL <- file.path(
  PATH_CLIMATE_DIR, "all_gcms_all_years_factual_pop.rds"
)
PATH_GCM_HISTNAT <- file.path(
  PATH_CLIMATE_DIR, "all_gcms_all_years_hist_nat_pop.rds"
)

## Administrative boundaries.
PATH_CITY_MAP <- file.path(DATA_DIR, "boundaries", "city_map.shp")
PATH_PROVINCE_MAP <- file.path(DATA_DIR, "boundaries", "province.shp")
PATH_CITY_GEOJSON_DIR <- file.path(DATA_DIR, "boundaries", "geometryProvince")

## Socioeconomic rasters.
PATH_GDP_DIRS <- file.path(DATA_DIR, "socioeconomic", "gridded_gdp")
PATH_GDP_PER_CAPITA <- file.path(
  DATA_DIR, "socioeconomic", "rast_adm2_gdp_perCapita_1990_2024.tif"
)
PATH_HFP_DIR <- file.path(DATA_DIR, "socioeconomic", "human_footprint")
PATH_NTL_DIR <- file.path(DATA_DIR, "socioeconomic", "nighttime_light")

## Sectoral water-pollution discharge rasters (Zenodo 10.5281/zenodo.16930436),
## laid out as <dir>/<COD|NH3N>/<Urban|Rural|Industry>/chn_1km_pol_<sector>_<year>.tif
PATH_SANITATION_RASTER_DIR <- file.path(DATA_DIR, "sanitation", "water_pollution_raster")
