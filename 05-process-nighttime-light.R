# -*- coding: UTF-8 -*-
# ==============================================================================
# Script Name:  05-process-nighttime-light.R
# Description:  提取逐年 NPP-VIIRS-like 夜间灯光的行政区指标。
# Scope:        仅处理 22262545 中 SRUNet_NPP_VIIRS_V2_Like_YYYY.zip。
# Period:       2008-2019
# Raster:       WGS84，15 arc-second（约 500 m），Float32。
# Unit:         nW cm^-2 sr^-1（单位来自数据论文，TIFF unit tag 为空）。
# Metrics:      NTLMean = 面积加权平均辐亮度。
#               NTLTotal = 按边界覆盖比例求和的总灯光强度。
#               NTLLitFrac = 亮灯面积比例（像元值 > 0，范围 0-1）。
# Final output: 一个宽表 CSV，每行一个 code_merge。
# ==============================================================================

options(stringsAsFactors = FALSE, scipen = 999)

if (.Platform$OS.type == "windows" && identical(Sys.getlocale("LC_CTYPE"), "C")) {
  invisible(suppressWarnings(Sys.setlocale("LC_CTYPE", "Chinese")))
}

required_packages <- c(
  "archive", "data.table", "dplyr", "sf", "stringr", "terra", "tidyr"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1L), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Install required package(s): ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

library(archive)
library(data.table)
library(dplyr)
library(sf)
library(stringr)
library(terra)
library(tidyr)

sf_use_s2(FALSE)

# ==============================================================================
# 0. 路径与参数
# ==============================================================================

analysis_years <- 2008:2019
lit_threshold <- 0
ntl_unit <- "nW cm^-2 sr^-1"

source("00-paths.R")
ntl_dir <- PATH_NTL_DIR
incidence_file <- PATH_HEALTH
city_map_file <- PATH_CITY_MAP
province_map_file <- PATH_PROVINCE_MAP

output_dir <- file.path("SII", "data", "nighttime-light")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

message("==========================================================")
message("Nighttime-light processing started: ", Sys.time())
message("Years: ", min(analysis_years), "-", max(analysis_years))
message("==========================================================")

# ==============================================================================
# 1. 传染病分析代码
# ==============================================================================

message("Step 1: Loading incidence codes...")

load(incidence_file)
if (!exists("incidence") || !"code_merge" %in% names(incidence)) {
  stop("The incidence RData does not contain incidence$code_merge.", call. = FALSE)
}

code_infec <- sort(unique(str_pad(
  as.character(incidence$code_merge),
  width = 6L,
  side = "right",
  pad = "0"
)))
code_infec <- code_infec[!is.na(code_infec) & nchar(code_infec) == 6L]

rm(incidence)
gc(verbose = FALSE)

message("  Unique target codes: ", length(code_infec))

# ==============================================================================
# 2. 城市和省级边界
# ==============================================================================

message("Step 2: Loading city and province maps...")

normalise_province_code <- function(x) {
  raw_value <- trimws(as.character(x))
  numeric_value <- suppressWarnings(as.numeric(raw_value))
  standard_value <- ifelse(
    is.na(numeric_value),
    raw_value,
    format(numeric_value, scientific = FALSE, trim = TRUE)
  )
  standard_value <- trimws(sub("\\..*$", "", standard_value))

  short_code <- !is.na(standard_value) & nchar(standard_value) <= 2L
  standard_value[short_code] <- str_pad(
    standard_value[short_code], 2L, side = "left", pad = "0"
  )
  standard_value[!short_code & !is.na(standard_value)] <- str_pad(
    standard_value[!short_code & !is.na(standard_value)],
    6L,
    side = "left",
    pad = "0"
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
      as.character(code_merge), 6L, side = "right", pad = "0"
    )
  ) %>%
  filter(!is.na(code_merge), code_merge %in% code_infec) %>%
  st_make_valid() %>%
  dplyr::select(code_merge, geometry) %>%
  group_by(code_merge) %>%
  summarise(.groups = "drop")

province_map <- st_read(province_map_file, quiet = TRUE)
province_type_field <- "\u7c7b\u578b"
province_code_field <- "\u7701\u4ee3\u7801"
province_exclude_value <- "\u4e0d\u7edf\u8ba1"

if (!all(c(province_type_field, province_code_field) %in% names(province_map))) {
  stop("Province shapefile is missing its type and/or code field.", call. = FALSE)
}

province_map <- province_map %>%
  filter(.data[[province_type_field]] != province_exclude_value) %>%
  mutate(code_p = normalise_province_code(.data[[province_code_field]])) %>%
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
  extraction_level = fifelse(exact_city_geometry, "city", "province_fallback"),
  source_code = fifelse(exact_city_geometry, code_merge, code_p)
)]

unmatched_provinces <- code_crosswalk[
  exact_city_geometry == FALSE & !source_code %chin% province_map$code_p,
  unique(source_code)
]
if (length(unmatched_provinces) > 0L) {
  stop(
    "No province fallback geometry for province code(s): ",
    paste(unmatched_provinces, collapse = ", "),
    call. = FALSE
  )
}

message("  Exact city codes: ", code_crosswalk[extraction_level == "city", .N])
message(
  "  Province fallback codes: ",
  code_crosswalk[extraction_level == "province_fallback", .N]
)

# ==============================================================================
# 3. 仅索引统一主序列 SRUNet_NPP_VIIRS_V2_Like
# ==============================================================================

message("Step 3: Indexing nighttime-light archives...")

if (!dir.exists(ntl_dir)) {
  stop("Nighttime-light directory not found: ", ntl_dir, call. = FALSE)
}

ntl_zip_files <- list.files(
  ntl_dir,
  pattern = "^SRUNet_NPP_VIIRS_V2_Like_[0-9]{4}\\.zip$",
  full.names = TRUE,
  ignore.case = TRUE
)

ntl_manifest <- data.table(
  year = as.integer(str_extract(basename(ntl_zip_files), "[0-9]{4}")),
  zip_file = ntl_zip_files
)
ntl_manifest <- ntl_manifest[year %in% analysis_years]
setorder(ntl_manifest, year)

missing_years <- setdiff(analysis_years, ntl_manifest$year)
duplicate_years <- ntl_manifest[, .N, by = year][N != 1L]
if (length(missing_years) > 0L || nrow(duplicate_years) > 0L) {
  stop(
    "NTL manifest must contain exactly one main-series ZIP per year. ",
    "Missing years: ", paste(missing_years, collapse = ", "),
    "; duplicated years: ", paste(duplicate_years$year, collapse = ", "),
    call. = FALSE
  )
}

# ==============================================================================
# 4. 解压当年 LZMA ZIP
# ==============================================================================

# 这些 ZIP 使用 compression method 14 (LZMA)，不使用 utils::unzip()。
# 每年将 ZIP 复制到临时 ASCII 路径，解压当年约 11 GB TIFF，
# 提取完成后立即删除整个当年临时目录。
unpack_ntl_year <- function(year, zip_file, temp_year_dir) {
  expected_member <- paste0("SRUNet_NPP_VIIRS_V2_Like_", year, ".tif")
  local_zip <- file.path(temp_year_dir, paste0("ntl", year, ".zip"))

  if (!file.copy(zip_file, local_zip, overwrite = TRUE)) {
    stop("Could not copy NTL archive to temporary directory: ", zip_file)
  }

  archive_contents <- archive::archive(local_zip)
  matching_member <- archive_contents$path[
    tolower(basename(archive_contents$path)) == tolower(expected_member)
  ]
  if (length(matching_member) != 1L) {
    stop(
      "Expected exactly one ", expected_member,
      " in ", basename(zip_file),
      call. = FALSE
    )
  }

  archive::archive_extract(
    archive = local_zip,
    dir = temp_year_dir,
    files = matching_member
  )

  extracted_file <- file.path(temp_year_dir, matching_member)
  if (!file.exists(extracted_file)) {
    stop("NTL TIFF was not extracted: ", expected_member, call. = FALSE)
  }
  extracted_file
}

# ==============================================================================
# 5. 县域 NTL 统计提取
# ==============================================================================

build_ntl_extraction_stack <- function(ntl_raster) {
  # 15 arc-second 经纬网格的实际面积随纬度变化。用 cellSize()
  # 构建 km2 权重，使平均亮度和亮灯比例是真正的面积加权指标。
  cell_area_km2 <- cellSize(ntl_raster, unit = "km")
  valid_area_km2 <- ifel(is.na(ntl_raster), NA, cell_area_km2)
  radiance_area <- ntl_raster * cell_area_km2
  lit_area_km2 <- ifel(
    is.na(ntl_raster),
    NA,
    ifel(ntl_raster > lit_threshold, cell_area_km2, 0)
  )

  extraction_stack <- c(
    ntl_raster,
    radiance_area,
    valid_area_km2,
    lit_area_km2
  )
  names(extraction_stack) <- c(
    "ntl_total",
    "radiance_area",
    "valid_area_km2",
    "lit_area_km2"
  )

  extraction_stack
}

extract_ntl_statistics <- function(extraction_stack, polygons, id_column) {
  if (nrow(polygons) == 0L) {
    return(data.table(
      source_code = character(),
      ntl_mean = numeric(),
      ntl_total = numeric(),
      ntl_lit_fraction = numeric()
    ))
  }

  polygons_raster_crs <- st_transform(
    polygons,
    crs(extraction_stack, proj = TRUE)
  )

  extracted <- terra::extract(
    extraction_stack,
    vect(polygons_raster_crs),
    fun = sum,
    na.rm = TRUE,
    exact = TRUE
  )

  result <- data.table(
    source_code = as.character(polygons[[id_column]]),
    ntl_total = as.numeric(extracted$ntl_total),
    ntl_mean = as.numeric(extracted$radiance_area) /
      as.numeric(extracted$valid_area_km2),
    ntl_lit_fraction = as.numeric(extracted$lit_area_km2) /
      as.numeric(extracted$valid_area_km2)
  )

  numeric_columns <- c("ntl_total", "ntl_mean", "ntl_lit_fraction")
  for (column_name in numeric_columns) {
    set(
      result,
      i = which(!is.finite(result[[column_name]])),
      j = column_name,
      value = NA_real_
    )
  }
  result
}

# ==============================================================================
# 6. 处理单个年份
# ==============================================================================

extract_one_ntl_year <- function(year, zip_file) {
  message("  Processing ", year, ": ", basename(zip_file))

  temp_year_dir <- tempfile(pattern = paste0("ntl-", year, "-"))
  dir.create(temp_year_dir, recursive = TRUE, showWarnings = FALSE)

  temp_root <- normalizePath(tempdir(), winslash = "/", mustWork = TRUE)
  resolved_temp_year <- normalizePath(
    temp_year_dir,
    winslash = "/",
    mustWork = TRUE
  )
  if (!startsWith(resolved_temp_year, paste0(temp_root, "/"))) {
    stop("Refusing to use a temporary directory outside tempdir().")
  }
  on.exit(unlink(temp_year_dir, recursive = TRUE, force = TRUE), add = TRUE)

  tif_file <- unpack_ntl_year(year, zip_file, temp_year_dir)
  ntl_global <- rast(tif_file)

  inspected_extent <- c(
    -180.00208333335,
    180.00208621353,
    -65.00208445335,
    75.00208333335
  )
  actual_extent <- as.vector(ext(ntl_global))

  invalid_structure <-
    nlyr(ntl_global) != 1L ||
    names(ntl_global) != paste0("SRUNet_NPP_VIIRS_V2_Like_", year) ||
    nrow(ntl_global) != 33601L ||
    ncol(ntl_global) != 86401L ||
    any(abs(res(ntl_global) - c(0.0041666667, 0.0041666667)) > 1e-10) ||
    any(abs(actual_extent - inspected_extent) > 1e-5) ||
    as.character(crs(ntl_global, describe = TRUE)$code[1]) != "4326" ||
    datatype(ntl_global) != "FLT4S"

  if (invalid_structure) {
    stop(
      "NTL TIFF structure differs from the inspected product: ",
      basename(tif_file),
      call. = FALSE
    )
  }

  province_for_crop <- st_transform(province_map, crs(ntl_global, proj = TRUE))
  ntl_china <- crop(ntl_global, vect(province_for_crop), snap = "out")
  ntl_extraction_stack <- build_ntl_extraction_stack(ntl_china)

  exact_codes <- code_crosswalk[exact_city_geometry == TRUE, code_merge]
  city_subset <- city_map[city_map$code_merge %in% exact_codes, ]

  fallback_province_codes <- code_crosswalk[
    exact_city_geometry == FALSE,
    unique(source_code)
  ]
  province_subset <- province_map[
    province_map$code_p %in% fallback_province_codes,
  ]

  city_values <- extract_ntl_statistics(
    ntl_extraction_stack,
    city_subset,
    "code_merge"
  )
  province_values <- extract_ntl_statistics(
    ntl_extraction_stack,
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
  result[, `:=`(
    year = as.integer(year),
    ntl_unit = ntl_unit,
    lit_threshold = lit_threshold,
    source_file = basename(zip_file)
  )]
  setorder(result, code_merge)

  result[, .(
    code_merge,
    year,
    ntl_mean,
    ntl_total,
    ntl_lit_fraction,
    ntl_unit,
    lit_threshold,
    extraction_level,
    source_code,
    source_file
  )]
}

# ==============================================================================
# 7. 循环处理 2008-2019
# ==============================================================================

message("Step 4: Extracting nighttime-light indicators...")

ntl_results <- lapply(seq_len(nrow(ntl_manifest)), function(i) {
  extract_one_ntl_year(
    year = ntl_manifest$year[i],
    zip_file = ntl_manifest$zip_file[i]
  )
})

ntl_long <- rbindlist(ntl_results, use.names = TRUE)
setorder(ntl_long, code_merge, year)

expected_rows <- length(code_infec) * length(analysis_years)
if (nrow(ntl_long) != expected_rows) {
  stop(
    "Unexpected result size: expected ", expected_rows,
    " rows, received ", nrow(ntl_long),
    call. = FALSE
  )
}
if (anyDuplicated(ntl_long[, paste(code_merge, year)])) {
  stop("Duplicated code_merge-year rows in NTL result.", call. = FALSE)
}
if (any(ntl_long$ntl_mean < 0 | ntl_long$ntl_total < 0, na.rm = TRUE)) {
  stop("Negative NTL indicators were found after extraction.", call. = FALSE)
}
if (any(
  ntl_long$ntl_lit_fraction < 0 | ntl_long$ntl_lit_fraction > 1,
  na.rm = TRUE
)) {
  stop("Illuminated-area fractions outside 0-1 were found.", call. = FALSE)
}

# ==============================================================================
# 8. 生成一个最终宽表 CSV
# ==============================================================================

message("Step 5: Creating the final wide nighttime-light table...")

mean_wide <- ntl_long %>%
  as.data.frame() %>%
  dplyr::select(code_merge, year, ntl_mean) %>%
  pivot_wider(names_from = year, values_from = ntl_mean, names_prefix = "NTLMean")

total_wide <- ntl_long %>%
  as.data.frame() %>%
  dplyr::select(code_merge, year, ntl_total) %>%
  pivot_wider(names_from = year, values_from = ntl_total, names_prefix = "NTLTotal")

lit_fraction_wide <- ntl_long %>%
  as.data.frame() %>%
  dplyr::select(code_merge, year, ntl_lit_fraction) %>%
  pivot_wider(
    names_from = year,
    values_from = ntl_lit_fraction,
    names_prefix = "NTLLitFrac"
  )

ntl_wide <- mean_wide %>%
  left_join(total_wide, by = "code_merge") %>%
  left_join(lit_fraction_wide, by = "code_merge") %>%
  arrange(code_merge)

expected_columns <- c(
  "code_merge",
  paste0("NTLMean", analysis_years),
  paste0("NTLTotal", analysis_years),
  paste0("NTLLitFrac", analysis_years)
)
if (!identical(names(ntl_wide), expected_columns)) {
  stop(
    "Unexpected final NTL columns: ",
    paste(names(ntl_wide), collapse = ", "),
    call. = FALSE
  )
}
if (nrow(ntl_wide) != length(code_infec) || anyDuplicated(ntl_wide$code_merge)) {
  stop(
    "The final NTL table must contain exactly one row per code_merge.",
    call. = FALSE
  )
}

final_csv <- file.path(output_dir, "nighttime-light-2008-2019.csv")
fwrite(as.data.table(ntl_wide), final_csv, bom = TRUE, na = "")

ntl_qc <- ntl_long[, .(
  n_codes = uniqueN(code_merge),
  n_city_extractions = sum(extraction_level == "city"),
  n_province_fallbacks = sum(extraction_level == "province_fallback"),
  n_missing_mean = sum(is.na(ntl_mean)),
  median_ntl_mean = median(ntl_mean, na.rm = TRUE),
  median_ntl_total = median(ntl_total, na.rm = TRUE),
  median_lit_fraction = median(ntl_lit_fraction, na.rm = TRUE)
), by = year]
setorder(ntl_qc, year)

message("==========================================================")
message("Nighttime-light processing completed: ", Sys.time())
message("Final CSV: ", final_csv)
message("Rows (unique code_merge): ", nrow(ntl_wide))
message("Metrics: NTLMean, NTLTotal, NTLLitFrac for 2008-2019")
message("==========================================================")
print(ntl_qc)
