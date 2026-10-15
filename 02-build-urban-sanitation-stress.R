# -*- coding: UTF-8 -*-
###############################################################################
## Build city-level urban sanitation-system stress indicators
##
## Source data
##   High-resolution gridded dataset of sectoral water pollution discharges in
##   China, 2007-2022, 30 arc-second, Zenodo 10.5281/zenodo.16930436
##   (Yuan Ze and Ma Ting, 2025; CC BY 4.0). Cell values are annual discharges
##   to the environment in kg per cell per year, that is, after in-sector
##   treatment removal, downscaled from national and provincial pollution
##   source statistics with land use, night-time lights, population, livestock
##   and geocoded industrial enterprises.
##
## Sectors used
##   Urban    urban residential (municipal domestic sewage)
##   Rural    rural residential (largely uncollected domestic sewage)
##   Industry industrial point sources (non-domestic contrast)
##
## Analysis unit
##   code_merge, identical to the health and GDP analyses: six-digit district
##   codes inside the four municipalities, otherwise prefecture code + "00".
##
## Primary sanitation-stress metric
##   Urban residential COD discharged per resident per year. For a given
##   population, a larger residual discharge implies that a smaller share of
##   domestic sewage is intercepted, conveyed and treated by the municipal
##   system, i.e. a more heavily loaded urban sanitation system.
##
## DENOMINATOR WARNING
##   pop.total.RData is NOT total resident population. It is the study
##   denominator for children and adolescents, covers 2013-2021 only, and sums
##   to 328 million nationally in 2015. The columns ending in
##   "_per_child_adolescent" therefore divide an all-age discharge by a child
##   count whose share of the population itself varies with income, so they
##   must not be read as per-resident loads. The population-free density and
##   share columns are the ones used downstream.
###############################################################################

options(stringsAsFactors = FALSE, scipen = 999)

suppressPackageStartupMessages({
  library(data.table)
  library(sf)
  library(terra)
  library(exactextractr)
  library(stringr)
})

sf_use_s2(FALSE)

analysis_years <- 2008:2019
pollutants <- c("COD", "NH3N")
sectors <- c("Urban", "Rural", "Industry")
municipality_provinces <- c("11", "12", "31", "50")

project_root <- normalizePath(".", winslash = "/", mustWork = TRUE)

source(file.path(project_root, "00-paths.R"))

raster_root <- Sys.getenv("SANITATION_RASTER_ROOT", PATH_SANITATION_RASTER_DIR)

city_map_file <- PATH_CITY_MAP

population_file <- PATH_POP_TOTAL

output_dir <- file.path(
  project_root, "Mechanism", "data_external",
  "03_China_urban_water_pollution"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

output_file <- file.path(
  output_dir, "urban-sanitation-stress-city-year.csv"
)

# Six-digit administrative code, then the shared analysis-unit mapping.
pad_six <- function(x) {
  x <- sub("[.]0+$", "", as.character(x))
  x <- suppressWarnings(as.integer(x))
  ifelse(is.na(x), NA_character_, sprintf("%06d", x))
}

to_code_merge <- function(code) {
  province <- substr(code, 1L, 2L)
  fcase(
    is.na(code), NA_character_,
    province %chin% municipality_provinces & substr(code, 5L, 6L) != "00", code,
    province %chin% municipality_provinces, NA_character_,
    default = paste0(substr(code, 1L, 4L), "00")
  )
}

# ---------------------------------------------------------------------------
# 1. City polygons dissolved to code_merge
# ---------------------------------------------------------------------------

message("Reading city polygons ...")

city_map <- st_read(city_map_file, quiet = TRUE)
city_map$code_merge <- str_pad(
  sub("[.]0+$", "", as.character(city_map$code_merge)),
  width = 6L, side = "right", pad = "0"
)
city_map <- city_map[
  !is.na(city_map$code_merge) & nchar(city_map$code_merge) == 6L,
  "code_merge"
]
city_map <- st_make_valid(city_map)
city_map <- aggregate(
  city_map,
  by = list(code_merge = city_map$code_merge),
  FUN = function(x) x[1L],
  do_union = TRUE
)
city_map <- city_map[, "code_merge"]
message("City polygons: ", nrow(city_map))

city_area_km2 <- as.numeric(st_area(st_transform(city_map, 6933))) / 1e6

# ---------------------------------------------------------------------------
# 2. Zonal sums of gridded discharge and of loaded urban area
# ---------------------------------------------------------------------------

cell_area <- NULL
zonal_records <- list()

for (pollutant in pollutants) {
  for (sector in sectors) {
    for (year in analysis_years) {
      raster_file <- file.path(
        raster_root, pollutant, sector,
        sprintf("chn_1km_pol_%s_%d.tif", sector, year)
      )
      if (!file.exists(raster_file)) {
        stop("Missing raster: ", raster_file)
      }

      grid <- rast(raster_file)
      names(grid) <- "discharge"

      if (is.null(cell_area)) {
        cell_area <- cellSize(grid, unit = "km", mask = FALSE)
        names(cell_area) <- "cell_km2"
      }

      totals <- exact_extract(grid, city_map, fun = "sum", progress = FALSE)

      loaded_area <- exact_extract(
        c(grid, cell_area), city_map,
        fun = function(values, coverage_fraction) {
          loaded <- is.finite(values$discharge) & values$discharge > 0
          sum(values$cell_km2[loaded] * coverage_fraction[loaded], na.rm = TRUE)
        },
        progress = FALSE
      )

      zonal_records[[paste(pollutant, sector, year)]] <- data.table(
        code_merge = city_map$code_merge,
        year = year,
        pollutant = pollutant,
        sector = sector,
        discharge_kg = as.numeric(totals),
        loaded_km2 = as.numeric(loaded_area)
      )

      message(sprintf(
        "  %s / %-8s / %d  national sum = %.3e kg",
        pollutant, sector, year, sum(as.numeric(totals), na.rm = TRUE)
      ))

      rm(grid)
      invisible(gc(verbose = FALSE))
    }
  }
}

zonal <- rbindlist(zonal_records, use.names = TRUE)

# ---------------------------------------------------------------------------
# 3. Population denominator on the same analysis unit
# ---------------------------------------------------------------------------

message("Reading the child and adolescent denominator ...")

population_environment <- new.env(parent = emptyenv())
load(population_file, envir = population_environment)
population_raw <- as.data.table(population_environment[["pop.total"]])

population_raw[, code := pad_six(code)]
population_raw[, code_merge := to_code_merge(code)]
population_raw[, year := as.integer(year)]
population_raw[, pop := as.numeric(pop)]

population <- population_raw[
  year %in% analysis_years & !is.na(code_merge) &
    is.finite(pop) & pop > 0,
  .(child_adolescent_population = sum(pop, na.rm = TRUE)),
  by = .(code_merge, year)
]

message("Denominator units: ", uniqueN(population$code_merge))

# ---------------------------------------------------------------------------
# 4. City-year table with derived stress metrics
# ---------------------------------------------------------------------------

zonal[, variable := paste0(tolower(pollutant), "_", tolower(sector))]

discharge_wide <- dcast(
  zonal, code_merge + year ~ variable, value.var = "discharge_kg"
)

area_wide <- dcast(
  zonal[pollutant == "COD"],
  code_merge + year ~ paste0("loaded_km2_", tolower(sector)),
  value.var = "loaded_km2"
)

sanitation <- merge(discharge_wide, area_wide,
                    by = c("code_merge", "year"), all = TRUE)
sanitation <- merge(sanitation, population,
                    by = c("code_merge", "year"), all.x = TRUE)
sanitation <- merge(
  sanitation,
  data.table(code_merge = city_map$code_merge, city_area_km2 = city_area_km2),
  by = "code_merge", all.x = TRUE
)

sanitation[, `:=`(
  cod_domestic = cod_urban + cod_rural,
  nh3n_domestic = nh3n_urban + nh3n_rural
)]

per_capita_targets <- c(
  "cod_urban", "cod_rural", "cod_industry", "cod_domestic",
  "nh3n_urban", "nh3n_rural", "nh3n_industry", "nh3n_domestic"
)

for (variable in per_capita_targets) {
  sanitation[, (paste0(variable, "_per_child_adolescent")) := fifelse(
    is.finite(child_adolescent_population) & child_adolescent_population > 0,
    get(variable) / child_adolescent_population, NA_real_
  )]
}

sanitation[, `:=`(
  cod_urban_per_urban_km2 = fifelse(
    is.finite(loaded_km2_urban) & loaded_km2_urban > 0,
    cod_urban / loaded_km2_urban, NA_real_
  ),
  nh3n_urban_per_urban_km2 = fifelse(
    is.finite(loaded_km2_urban) & loaded_km2_urban > 0,
    nh3n_urban / loaded_km2_urban, NA_real_
  ),
  cod_rural_per_rural_km2 = fifelse(
    is.finite(loaded_km2_rural) & loaded_km2_rural > 0,
    cod_rural / loaded_km2_rural, NA_real_
  ),
  nh3n_rural_per_rural_km2 = fifelse(
    is.finite(loaded_km2_rural) & loaded_km2_rural > 0,
    nh3n_rural / loaded_km2_rural, NA_real_
  ),
  cod_domestic_per_loaded_km2 = fifelse(
    (loaded_km2_urban + loaded_km2_rural) > 0,
    cod_domestic / (loaded_km2_urban + loaded_km2_rural), NA_real_
  ),
  urban_land_share = fifelse(
    is.finite(city_area_km2) & city_area_km2 > 0,
    loaded_km2_urban / city_area_km2, NA_real_
  ),
  domestic_share_of_water_pollution = fifelse(
    (cod_domestic + cod_industry) > 0,
    cod_domestic / (cod_domestic + cod_industry), NA_real_
  ),
  rural_share_of_domestic = fifelse(
    cod_domestic > 0, cod_rural / cod_domestic, NA_real_
  )
)]

setorder(sanitation, code_merge, year)
fwrite(sanitation, output_file)

message("Written: ", output_file)
message("City-years: ", nrow(sanitation),
        " | cities: ", uniqueN(sanitation$code_merge),
        " | with a child denominator: ",
        uniqueN(sanitation[is.finite(child_adolescent_population)]$code_merge))
