# -*- coding: UTF-8 -*-
# ==============================================================================
# Script Name:  04-process-human-footprint.R
# Description:  将逐年 Human Footprint (HFP) 栅格按行政区边界覆盖
#               比例计算面积加权平均值。
# Scope:        仅处理 Human Footprint；不处理 GDP、NTL、GHSL 和医院数据。
# Period:       2008-2019
# Raster:       hfpYYYY.tif，World Mollweide，1000 m，单波段 Float32。
# Value:        HFP 强度指数，代表年份的实测范围约为 0-50。
# Aggregation:  area-fraction weighted mean；0 为有效值，NaN/NA 为缺失。
# Final output: 一个宽表 CSV；每行一个 code_merge，每年一个 HFP 列。
# ==============================================================================

options(stringsAsFactors = FALSE, scipen = 999)

if (.Platform$OS.type == "windows" && identical(Sys.getlocale("LC_CTYPE"), "C")) {
  invisible(suppressWarnings(Sys.setlocale("LC_CTYPE", "Chinese")))
}

required_packages <- c("data.table", "dplyr", "sf", "stringr", "terra", "tidyr")
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
hfp_aggregation_method <- "area-fraction weighted mean"

source("00-paths.R")
hfp_dir <- PATH_HFP_DIR

incidence_file <- PATH_HEALTH
city_map_file <- PATH_CITY_MAP
province_map_file <- PATH_PROVINCE_MAP

output_dir <- file.path("SII", "data", "human-footprint")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

message("==========================================================")
message("Human Footprint processing started: ", Sys.time())
message("Years: ", min(analysis_years), "-", max(analysis_years))
message("==========================================================")

# ==============================================================================
# 1. 读取传染病分析代码
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
# 2. 读取城市和省级边界
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
    standard_value[short_code],
    width = 2L,
    side = "left",
    pad = "0"
  )
  standard_value[!short_code & !is.na(standard_value)] <- str_pad(
    standard_value[!short_code & !is.na(standard_value)],
    width = 6L,
    side = "left",
    pad = "0"
  )

  ifelse(
    is.na(standard_value) | !nzchar(standard_value),
    NA_character_,
    substr(standard_value, 1L, 2L)
  )
}

city_map_raw <- st_read(city_map_file, quiet = TRUE) %>%
  mutate(
    code_merge = str_pad(
      as.character(code_merge),
      width = 6L,
      side = "right",
      pad = "0"
    )
  ) %>%
  filter(!is.na(code_merge), code_merge %in% code_infec) %>%
  st_make_valid()

city_map <- city_map_raw %>%
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
# 3. 建立 2008-2019 ZIP 文件清单
# ==============================================================================

message("Step 3: Indexing Human Footprint archives...")

if (!dir.exists(hfp_dir)) {
  stop("Human Footprint directory not found: ", hfp_dir, call. = FALSE)
}

hfp_zip_files <- list.files(
  hfp_dir,
  pattern = "^hfp[0-9]{4}\\.zip$",
  full.names = TRUE,
  ignore.case = TRUE
)

hfp_manifest <- data.table(
  year = as.integer(str_extract(basename(hfp_zip_files), "[0-9]{4}")),
  zip_file = hfp_zip_files
)
hfp_manifest <- hfp_manifest[year %in% analysis_years]
setorder(hfp_manifest, year)

missing_years <- setdiff(analysis_years, hfp_manifest$year)
duplicate_years <- hfp_manifest[, .N, by = year][N != 1L]

if (length(missing_years) > 0L || nrow(duplicate_years) > 0L) {
  stop(
    "HFP archive manifest must contain exactly one ZIP per year. ",
    "Missing years: ", paste(missing_years, collapse = ", "),
    "; duplicated years: ", paste(duplicate_years$year, collapse = ", "),
    call. = FALSE
  )
}

# ==============================================================================
# 4. 年度 ZIP 解压函数
# ==============================================================================

# 为了避免 Windows/GDAL 对中文 ZIP 路径的兼容问题，每年先将 ZIP
# 复制到 R 临时目录，仅解压当年 TIFF。当年提取完成后立即删除临时文件。
unpack_hfp_year <- function(year, zip_file, temp_year_dir) {
  expected_member <- paste0("hfp", year, ".tif")
  local_zip <- file.path(temp_year_dir, basename(zip_file))

  if (!file.copy(zip_file, local_zip, overwrite = TRUE)) {
    stop("Could not copy HFP archive to temporary directory: ", zip_file)
  }

  zip_members <- utils::unzip(local_zip, list = TRUE)$Name
  matching_member <- zip_members[
    tolower(basename(zip_members)) == tolower(expected_member)
  ]

  if (length(matching_member) != 1L) {
    stop(
      "Expected exactly one ", expected_member,
      " in ", basename(zip_file),
      call. = FALSE
    )
  }

  utils::unzip(
    local_zip,
    files = matching_member,
    exdir = temp_year_dir,
    overwrite = TRUE
  )

  extracted_file <- file.path(temp_year_dir, matching_member)
  if (!file.exists(extracted_file)) {
    stop("HFP TIFF was not extracted: ", expected_member, call. = FALSE)
  }

  extracted_file
}

# ==============================================================================
# 5. HFP 面积加权平均值提取函数
# ==============================================================================

extract_hfp_mean <- function(hfp_raster, polygons, id_column) {
  if (nrow(polygons) == 0L) {
    return(data.table(source_code = character(), hfp_mean = numeric()))
  }

  polygons_raster_crs <- st_transform(polygons, crs(hfp_raster, proj = TRUE))
  extracted <- terra::extract(
    hfp_raster,
    vect(polygons_raster_crs),
    fun = mean,
    na.rm = TRUE,
    exact = TRUE
  )

  value_column <- setdiff(names(extracted), "ID")
  if (length(value_column) != 1L) {
    stop("Unexpected HFP extraction output columns.", call. = FALSE)
  }

  result <- data.table(
    source_code = as.character(polygons[[id_column]]),
    hfp_mean = as.numeric(extracted[[value_column]])
  )
  result[!is.finite(hfp_mean), hfp_mean := NA_real_]
  result
}

# ==============================================================================
# 6. 处理单个年份
# ==============================================================================

extract_one_hfp_year <- function(year, zip_file) {
  message("  Processing ", year, ": ", basename(zip_file))

  temp_year_dir <- tempfile(pattern = paste0("hfp-", year, "-"))
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

  on.exit(
    unlink(temp_year_dir, recursive = TRUE, force = TRUE),
    add = TRUE
  )

  tif_file <- unpack_hfp_year(year, zip_file, temp_year_dir)
  hfp_global <- rast(tif_file)

  hfp_extent <- as.vector(ext(hfp_global))
  hfp_crs_name <- crs(hfp_global, describe = TRUE)$name[1]

  invalid_structure <-
    nlyr(hfp_global) != 1L ||
    names(hfp_global) != paste0("hfp", year) ||
    nrow(hfp_global) != 16382L ||
    ncol(hfp_global) != 36081L ||
    any(abs(res(hfp_global) - c(1000, 1000)) > 1e-8) ||
    any(abs(hfp_extent - c(-18040094, 18040906, -7363043, 9018957)) > 1e-6) ||
    is.na(hfp_crs_name) ||
    !grepl("Mollweide", hfp_crs_name, ignore.case = TRUE) ||
    datatype(hfp_global) != "FLT4S"

  if (invalid_structure) {
    stop(
      "HFP TIFF structure differs from the inspected product: ",
      basename(tif_file),
      call. = FALSE
    )
  }

  province_for_crop <- st_transform(
    province_map,
    crs(hfp_global, proj = TRUE)
  )
  hfp_china <- crop(
    hfp_global,
    vect(province_for_crop),
    snap = "out"
  )

  exact_codes <- code_crosswalk[exact_city_geometry == TRUE, code_merge]
  city_subset <- city_map[city_map$code_merge %in% exact_codes, ]

  fallback_province_codes <- code_crosswalk[
    exact_city_geometry == FALSE,
    unique(source_code)
  ]
  province_subset <- province_map[
    province_map$code_p %in% fallback_province_codes,
  ]

  city_values <- extract_hfp_mean(hfp_china, city_subset, "code_merge")
  province_values <- extract_hfp_mean(hfp_china, province_subset, "code_p")

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
    aggregation_method = hfp_aggregation_method,
    source_file = basename(zip_file)
  )]
  setorder(result, code_merge)

  result[, .(
    code_merge,
    year,
    hfp_mean,
    aggregation_method,
    extraction_level,
    source_code,
    source_file
  )]
}

# ==============================================================================
# 7. 循环处理 2008-2019
# ==============================================================================

message("Step 4: Extracting Human Footprint for all analysis years...")

hfp_results <- lapply(seq_len(nrow(hfp_manifest)), function(i) {
  extract_one_hfp_year(
    year = hfp_manifest$year[i],
    zip_file = hfp_manifest$zip_file[i]
  )
})

hfp_long <- rbindlist(hfp_results, use.names = TRUE)
setorder(hfp_long, code_merge, year)

expected_rows <- length(code_infec) * length(analysis_years)
if (nrow(hfp_long) != expected_rows) {
  stop(
    "Unexpected result size: expected ", expected_rows,
    " rows, received ", nrow(hfp_long),
    call. = FALSE
  )
}
if (anyDuplicated(hfp_long[, paste(code_merge, year)])) {
  stop("Duplicated code_merge-year rows in HFP result.", call. = FALSE)
}
if (any(hfp_long$hfp_mean < 0 | hfp_long$hfp_mean > 50, na.rm = TRUE)) {
  stop("HFP means outside the inspected 0-50 range were found.", call. = FALSE)
}

# ==============================================================================
# 8. 生成唯一的最终 CSV：一行一个 code_merge，每年一列
# ==============================================================================

message("Step 5: Creating the final wide Human Footprint table...")

final_csv <- file.path(output_dir, "human-footprint-2008-2019.csv")

hfp_wide <- hfp_long %>%
  as.data.frame() %>%
  dplyr::select(code_merge, year, hfp_mean) %>%
  pivot_wider(
    names_from = year,
    values_from = hfp_mean,
    names_prefix = "HFP"
  ) %>%
  arrange(code_merge)

expected_columns <- c("code_merge", paste0("HFP", analysis_years))
if (!identical(names(hfp_wide), expected_columns)) {
  stop(
    "Unexpected final HFP columns: ",
    paste(names(hfp_wide), collapse = ", "),
    call. = FALSE
  )
}
if (nrow(hfp_wide) != length(code_infec) || anyDuplicated(hfp_wide$code_merge)) {
  stop(
    "The final HFP table must contain exactly one row per code_merge.",
    call. = FALSE
  )
}

fwrite(as.data.table(hfp_wide), final_csv, bom = TRUE, na = "")

hfp_qc <- hfp_long[, .(
  n_codes = uniqueN(code_merge),
  n_city_extractions = sum(extraction_level == "city"),
  n_province_fallbacks = sum(extraction_level == "province_fallback"),
  n_missing = sum(is.na(hfp_mean)),
  n_zero = sum(hfp_mean == 0, na.rm = TRUE),
  hfp_min = if (all(is.na(hfp_mean))) NA_real_ else min(hfp_mean, na.rm = TRUE),
  hfp_median = if (all(is.na(hfp_mean))) NA_real_ else median(hfp_mean, na.rm = TRUE),
  hfp_mean = if (all(is.na(hfp_mean))) NA_real_ else mean(hfp_mean, na.rm = TRUE),
  hfp_max = if (all(is.na(hfp_mean))) NA_real_ else max(hfp_mean, na.rm = TRUE)
), by = year]
setorder(hfp_qc, year)

message("==========================================================")
message("Human Footprint processing completed: ", Sys.time())
message("Final CSV: ", final_csv)
message("Rows (unique code_merge): ", nrow(hfp_wide))
message("HFP year columns: HFP", min(analysis_years), "-HFP", max(analysis_years))
message("==========================================================")
print(hfp_qc)
