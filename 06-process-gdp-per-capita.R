# -*- coding: UTF-8 -*-
# ==============================================================================
# Process gridded GDP per capita for the sanitation analysis
# Source: downscaled gridded GDP per capita (PPP), Kummu et al., Scientific
#         Data (2025); multiband raster with one layer per year
#         (rast_adm2_gdp_perCapita_1990_2024.tif, layers gdp_pc_YYYY).
# Output: one row per code_merge and one GDPpc column per year, 2008-2019.
# Aggregation: mean of the grid cells within each analysis unit, weighted by
#         the fraction of each cell covered by the polygon. This is a spatial
#         mean, not a population-weighted mean.
# Analysis units without their own polygon use the province polygon, as in
#         03-process-gridded-gdp.R.
# ==============================================================================

options(stringsAsFactors = FALSE, scipen = 999)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(sf)
  library(stringr)
  library(terra)
  library(tidyr)
})

sf_use_s2(FALSE)

analysis_years <- 2008:2019
source("00-paths.R")
gdp_pc_file <- PATH_GDP_PER_CAPITA

incidence_file <- PATH_HEALTH
city_map_file <- PATH_CITY_MAP
province_map_file <- PATH_PROVINCE_MAP

output_dir <- file.path("SII", "data", "gdp-per-capita")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Analysis-region codes ---------------------------------------------------------

load(incidence_file)
if (!exists("incidence") || !"code_merge" %in% names(incidence)) {
  stop("The incidence RData does not contain incidence$code_merge.")
}

code_infec <- sort(unique(str_pad(
  as.character(incidence$code_merge),
  width = 6L,
  side = "right",
  pad = "0"
)))
code_infec <- code_infec[!is.na(code_infec) & nchar(code_infec) == 6L]
rm(incidence)

# Administrative boundaries ----------------------------------------------------

normalise_province_code <- function(x) {
  raw_value <- trimws(as.character(x))
  numeric_value <- suppressWarnings(as.numeric(raw_value))
  standard_value <- ifelse(
    is.na(numeric_value),
    raw_value,
    format(numeric_value, scientific = FALSE, trim = TRUE)
  )
  standard_value <- trimws(sub("\\..*$", "", standard_value))
  standard_value <- ifelse(
    nchar(standard_value) <= 2L,
    str_pad(standard_value, 2L, side = "left", pad = "0"),
    str_pad(standard_value, 6L, side = "left", pad = "0")
  )
  ifelse(
    is.na(standard_value) | !nzchar(standard_value),
    NA_character_,
    substr(standard_value, 1L, 2L)
  )
}

city_map <- st_read(city_map_file, quiet = TRUE) %>%
  mutate(
    code_merge = str_pad(
      as.character(code_merge),
      width = 6L,
      side = "right",
      pad = "0"
    )
  ) %>%
  filter(code_merge %in% code_infec) %>%
  st_make_valid() %>%
  dplyr::select(code_merge, geometry) %>%
  group_by(code_merge) %>%
  summarise(.groups = "drop")

province_map <- st_read(province_map_file, quiet = TRUE)
province_type_field <- "类型"
province_code_field <- "省代码"
province_exclude_value <- "不统计"

province_map <- province_map %>%
  filter(.data[[province_type_field]] != province_exclude_value) %>%
  mutate(
    code_p = normalise_province_code(.data[[province_code_field]])
  ) %>%
  filter(!is.na(code_p), nzchar(code_p)) %>%
  dplyr::select(code_p, geometry) %>%
  st_make_valid() %>%
  group_by(code_p) %>%
  summarise(.groups = "drop")

code_crosswalk <- data.table(code_merge = code_infec)
code_crosswalk[, code_p := substr(code_merge, 1L, 2L)]
code_crosswalk[code_p == "66", code_p := "65"]
code_crosswalk[, exact_city_geometry := code_merge %chin% city_map$code_merge]
code_crosswalk[, source_code := fifelse(exact_city_geometry, code_merge, code_p)]

# Raster layers -----------------------------------------------------------------

gdp_pc <- rast(gdp_pc_file)
layer_names <- paste0("gdp_pc_", analysis_years)
if (!all(layer_names %in% names(gdp_pc))) {
  stop(
    "GDP per capita raster lacks layers: ",
    paste(setdiff(layer_names, names(gdp_pc)), collapse = ", ")
  )
}
gdp_pc <- gdp_pc[[layer_names]]

province_for_crop <- st_transform(province_map, crs(gdp_pc, proj = TRUE))
gdp_pc <- crop(gdp_pc, vect(province_for_crop), snap = "out")

# Polygon extraction ------------------------------------------------------------

extract_gdp_pc_mean <- function(raster_stack, polygons, id_column) {
  if (nrow(polygons) == 0L) return(data.table())
  polygons <- st_transform(polygons, crs(raster_stack, proj = TRUE))
  extracted <- terra::extract(
    raster_stack,
    vect(polygons),
    fun = mean,
    na.rm = TRUE,
    exact = TRUE
  )
  out <- as.data.table(extracted)[, -1L]
  setnames(out, paste0("GDPpc", analysis_years))
  out[, source_code := as.character(polygons[[id_column]])]
  out
}

exact_codes <- code_crosswalk[exact_city_geometry == TRUE, code_merge]
fallback_codes <- code_crosswalk[exact_city_geometry == FALSE, unique(source_code)]

values <- rbindlist(list(
  extract_gdp_pc_mean(
    gdp_pc, city_map[city_map$code_merge %in% exact_codes, ], "code_merge"
  ),
  extract_gdp_pc_mean(
    gdp_pc, province_map[province_map$code_p %in% fallback_codes, ], "code_p"
  )
), use.names = TRUE, fill = TRUE)

gdp_pc_wide <- merge(
  code_crosswalk[, .(code_merge, source_code)],
  values,
  by = "source_code",
  all.x = TRUE,
  sort = FALSE
)[, source_code := NULL]
setorder(gdp_pc_wide, code_merge)

final_csv <- file.path(output_dir, "gdp-per-capita-2008-2019.csv")
fwrite(gdp_pc_wide, final_csv, bom = TRUE, na = "")

message("GDP per capita processing completed. Final CSV: ", final_csv)
