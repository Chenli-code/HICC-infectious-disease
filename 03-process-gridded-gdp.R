# -*- coding: UTF-8 -*-
# ==============================================================================
# Process Wang & Sun annual gridded GDP for the SII analysis
# Output: one row per code_merge and one GDP column per year, 2008-2019.
# Aggregation: polygon-coverage-weighted sum; zero is a valid GDP value.
# Unit: 2005 PPP international dollars, as documented by the data repository.
# Population data are not used in this script.
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
gdp_dirs <- PATH_GDP_DIRS

incidence_file <- PATH_HEALTH
city_map_file <- PATH_CITY_MAP
province_map_file <- PATH_PROVINCE_MAP

output_dir <- file.path("SII", "data", "gridded-gdp")
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
province_type_field <- "\u7c7b\u578b"
province_code_field <- "\u7701\u4ee3\u7801"
province_exclude_value <- "\u4e0d\u7edf\u8ba1"

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
code_crosswalk[, `:=`(
  extraction_level = fifelse(
    exact_city_geometry,
    "city",
    "province_fallback"
  ),
  source_code = fifelse(exact_city_geometry, code_merge, code_p)
)]

# Annual raster manifest --------------------------------------------------------

gdp_files <- unlist(lapply(
  gdp_dirs,
  list.files,
  pattern = "^GDP[0-9]{4}\\.tif$",
  full.names = TRUE,
  ignore.case = TRUE
), use.names = FALSE)

gdp_manifest <- data.table(
  year = as.integer(str_extract(basename(gdp_files), "[0-9]{4}")),
  raster_file = gdp_files
)[year %in% analysis_years]
setorder(gdp_manifest, year)

if (!identical(gdp_manifest$year, analysis_years)) {
  stop("GDP files must contain exactly one TIFF for every year from 2008 to 2019.")
}

# Polygon extraction ------------------------------------------------------------

extract_gdp_sum <- function(gdp_raster, polygons, id_column) {
  if (nrow(polygons) == 0L) {
    return(data.table(source_code = character(), GDP = numeric()))
  }

  polygons <- st_transform(polygons, crs(gdp_raster, proj = TRUE))
  extracted <- terra::extract(
    gdp_raster,
    vect(polygons),
    fun = sum,
    na.rm = TRUE,
    exact = TRUE
  )

  data.table(
    source_code = as.character(polygons[[id_column]]),
    GDP = as.numeric(extracted[[2L]])
  )
}

extract_one_gdp_year <- function(year, raster_file) {
  message("Processing ", year, ": ", basename(raster_file))

  gdp_global <- rast(raster_file)
  if (nlyr(gdp_global) != 1L) {
    stop("GDP TIFF must contain one raster layer: ", raster_file)
  }

  province_for_crop <- st_transform(
    province_map,
    crs(gdp_global, proj = TRUE)
  )
  gdp_china <- crop(
    gdp_global,
    vect(province_for_crop),
    snap = "out"
  )

  exact_codes <- code_crosswalk[exact_city_geometry == TRUE, code_merge]
  city_subset <- city_map[city_map$code_merge %in% exact_codes, ]

  fallback_codes <- code_crosswalk[
    exact_city_geometry == FALSE,
    unique(source_code)
  ]
  province_subset <- province_map[province_map$code_p %in% fallback_codes, ]

  city_values <- extract_gdp_sum(gdp_china, city_subset, "code_merge")
  province_values <- extract_gdp_sum(
    gdp_china,
    province_subset,
    "code_p"
  )

  city_result <- merge(
    code_crosswalk[exact_city_geometry == TRUE],
    city_values,
    by = "source_code",
    all.x = TRUE,
    sort = FALSE
  )
  fallback_result <- merge(
    code_crosswalk[exact_city_geometry == FALSE],
    province_values,
    by = "source_code",
    all.x = TRUE,
    sort = FALSE
  )

  result <- rbindlist(
    list(city_result, fallback_result),
    use.names = TRUE,
    fill = TRUE
  )
  result[, year := as.integer(year)]
  setorder(result, code_merge)

  result[, .(code_merge, year, GDP)]
}

gdp_results <- lapply(seq_len(nrow(gdp_manifest)), function(i) {
  extract_one_gdp_year(
    year = gdp_manifest$year[i],
    raster_file = gdp_manifest$raster_file[i]
  )
})

gdp_long <- rbindlist(gdp_results, use.names = TRUE)
setorder(gdp_long, code_merge, year)

gdp_wide <- gdp_long %>%
  as.data.frame() %>%
  pivot_wider(
    names_from = year,
    values_from = GDP,
    names_prefix = "GDP"
  ) %>%
  arrange(code_merge)

final_csv <- file.path(output_dir, "gridded-gdp-2008-2019.csv")
fwrite(as.data.table(gdp_wide), final_csv, bom = TRUE, na = "")

message("GDP processing completed. Final CSV: ", final_csv)
