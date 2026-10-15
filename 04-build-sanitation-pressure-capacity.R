# -*- coding: UTF-8 -*-
###############################################################################
## Sanitation-system stress as pressure minus capacity
##
## Stress = z(pressure) - z(capacity)
##
##   Pressure (higher = heavier residual load on the system)
##     P1 log urban residential COD per km2 of loaded urban land   [city]
##     P2 log rural residential COD per km2 of loaded rural land   [city]
##     P3 log domestic NH3-N per km2 of loaded domestic land       [city]
##
##   Capacity (higher = better served)
##     C1 logit of the sewered share of the domestic load, that is one minus
##        the rural-residential share                              [city]
##     C2 log treated wastewater volume per kg of residual domestic COD
##                                                                 [province]
##     C3 log reported treatment facilities per tonne of residual domestic COD
##                                                                 [province]
##
## NO POPULATION DENOMINATOR IS USED.
##   The project has no total resident population series. pop.total.RData is
##   the child and adolescent study denominator (2013-2021, 328 million
##   nationally in 2015), and the child share of the population itself varies
##   with income, so dividing an all-age discharge by it would build the
##   socioeconomic gradient into the indicator. Every term above is therefore a
##   ratio of two quantities that were measured directly: a discharge mass, a
##   loaded land area, a treated volume, or a facility count.
##
##   GDP per capita is unaffected: SII/06-process-gdp-per-capita.R reads a
##   dedicated per-capita GDP raster (rast_adm2_gdp_perCapita_1990_2024.tif)
##   and takes the within-city mean, so it needs no population input.
##
## Inputs are all pre-computed; nothing is re-modelled here.
##   Mechanism/data_external/03_China_urban_water_pollution/
##     urban-sanitation-stress-city-year.csv        (script 02)
##   Mechanism/data_external/05_China_WWTP_2006_2019/
##     province_year_WWTP_capacity_validation_2008_2019.csv
##   SII/data/gdp-per-capita/gdp-per-capita-2008-2019.csv
##   SII/data/gridded-gdp/gridded-gdp-2008-2019.csv  (robustness denominator)
##
## Caveat carried into the manuscript: C2 and C3 vary only between provinces,
## so the capacity term mixes a city-level coverage measure with a provincial
## treatment-intensity measure.
###############################################################################

options(stringsAsFactors = FALSE, scipen = 999)

suppressPackageStartupMessages({
  library(data.table)
})

project_root <- normalizePath(".", winslash = "/", mustWork = TRUE)

output_dir <- file.path(project_root, "Mechanism", "output")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

analysis_years <- 2008:2019

sanitation_file <- file.path(
  project_root, "Mechanism", "data_external",
  "03_China_urban_water_pollution",
  "urban-sanitation-stress-city-year.csv"
)

wwtp_file <- file.path(
  project_root, "Mechanism", "data_external", "05_China_WWTP_2006_2019",
  "province_year_WWTP_capacity_validation_2008_2019.csv"
)

gdp_per_capita_file <- file.path(
  project_root, "SII", "data", "gdp-per-capita",
  "gdp-per-capita-2008-2019.csv"
)

gridded_gdp_file <- file.path(
  project_root, "SII", "data", "gridded-gdp",
  "gridded-gdp-2008-2019.csv"
)

province_lookup <- data.table(
  province_name = c(
    "Beijing", "Tianjin", "Hebei", "Shanxi", "Inner Mongolia",
    "Liaoning", "Jilin", "Heilongjiang", "Shanghai", "Jiangsu",
    "Zhejiang", "Anhui", "Fujian", "Jiangxi", "Shandong",
    "Henan", "Hubei", "Hunan", "Guangdong", "Guangxi",
    "Hainan", "Chongqing", "Sichuan", "Guizhou", "Yunnan",
    "Tibet", "Shaanxi", "Gansu", "Qinghai", "Ningxia", "Xinjiang"
  ),
  province_code = c(
    "11", "12", "13", "14", "15",
    "21", "22", "23", "31", "32",
    "33", "34", "35", "36", "37",
    "41", "42", "43", "44", "45",
    "46", "50", "51", "52", "53",
    "54", "61", "62", "63", "64", "65"
  )
)

standardise <- function(x) {
  usable <- is.finite(x)
  result <- rep(NA_real_, length(x))
  result[usable] <- (x[usable] - mean(x[usable])) / sd(x[usable])
  result
}

logit <- function(p) {
  bounded <- pmin(pmax(p, 0.001), 0.999)
  log(bounded / (1 - bounded))
}

row_mean <- function(...) {
  rowMeans(cbind(...), na.rm = TRUE)
}

sanitation_panel <- fread(
  sanitation_file, colClasses = list(character = "code_merge")
)

sanitation_panel <- sanitation_panel[year %in% analysis_years]

sanitation_panel[, `:=`(
  nh3n_domestic_per_loaded_km2 = fifelse(
    (loaded_km2_urban + loaded_km2_rural) > 0,
    nh3n_domestic / (loaded_km2_urban + loaded_km2_rural), NA_real_
  ),
  rural_share_of_domestic_nh3n = fifelse(
    nh3n_domestic > 0, nh3n_rural / nh3n_domestic, NA_real_
  )
)]

city_indicators <- sanitation_panel[, .(
  cod_urban_per_urban_km2 = mean(cod_urban_per_urban_km2, na.rm = TRUE),
  cod_rural_per_rural_km2 = mean(cod_rural_per_rural_km2, na.rm = TRUE),
  cod_domestic_per_loaded_km2 = mean(cod_domestic_per_loaded_km2, na.rm = TRUE),
  nh3n_domestic_per_loaded_km2 = mean(nh3n_domestic_per_loaded_km2,
                                      na.rm = TRUE),
  rural_share_of_domestic = mean(rural_share_of_domestic, na.rm = TRUE),
  rural_share_of_domestic_nh3n = mean(rural_share_of_domestic_nh3n,
                                      na.rm = TRUE),
  cod_domestic = mean(cod_domestic, na.rm = TRUE),
  cod_industry = mean(cod_industry, na.rm = TRUE),
  cod_industry_per_industry_km2 = mean(
    fifelse(loaded_km2_industry > 0, cod_industry / loaded_km2_industry,
            NA_real_),
    na.rm = TRUE
  ),
  urban_land_share = mean(urban_land_share, na.rm = TRUE),
  child_adolescent_population = mean(child_adolescent_population, na.rm = TRUE)
), by = .(code = code_merge)]

city_indicators[, province_code := substr(code, 1L, 2L)]

city_indicators[, sewered_share_of_domestic := 1 - rural_share_of_domestic]

wwtp <- fread(wwtp_file, encoding = "UTF-8")

setnames(wwtp, names(wwtp), sub("^﻿", "", names(wwtp)))

wwtp <- merge(
  wwtp, province_lookup,
  by.x = "Province", by.y = "province_name", all.x = TRUE
)

if (anyNA(wwtp$province_code)) {
  stop("Unmatched province names: ",
       paste(unique(wwtp[is.na(province_code)]$Province), collapse = ", "))
}

wwtp_province <- wwtp[Year %in% analysis_years, .(
  treated_wastewater_m3 = 1e4 * mean(treated_wastewater_10k_m3, na.rm = TRUE),
  reported_wwtp_records = mean(reported_wwtp_records, na.rm = TRUE),
  enhanced_activated_sludge_share = mean(
    enhanced_activated_sludge_share, na.rm = TRUE
  ),
  wwtp_years = .N
), by = province_code]

province_load <- city_indicators[, .(
  province_domestic_cod = sum(cod_domestic, na.rm = TRUE)
), by = province_code]

capacity <- merge(wwtp_province, province_load, by = "province_code",
                  all.x = TRUE)

capacity[, `:=`(
  treated_volume_per_kg_cod = treated_wastewater_m3 / province_domestic_cod,
  wwtp_per_tonne_cod = 1e3 * reported_wwtp_records / province_domestic_cod
)]

index <- merge(city_indicators, capacity, by = "province_code", all.x = TRUE)

index[, `:=`(
  pressure_1 = standardise(log(cod_urban_per_urban_km2)),
  pressure_2 = standardise(log(cod_rural_per_rural_km2)),
  pressure_3 = standardise(log(nh3n_domestic_per_loaded_km2)),
  capacity_1 = standardise(logit(sewered_share_of_domestic)),
  capacity_2 = standardise(log(treated_volume_per_kg_cod)),
  capacity_3 = standardise(log(wwtp_per_tonne_cod))
)]

index[, `:=`(
  pressure_index = row_mean(pressure_1, pressure_2, pressure_3),
  capacity_index = row_mean(capacity_1, capacity_2, capacity_3)
)]

index[, sanitation_stress_raw := pressure_index - capacity_index]

index[, sanitation_stress := standardise(sanitation_stress_raw)]

# City-level-only variant, so that the provincial component can be shown not
# to drive the result.
index[, capacity_city_only := capacity_1]

index[, sanitation_stress_city_only := standardise(
  pressure_index - capacity_city_only
)]

gdp_per_capita_wide <- fread(gdp_per_capita_file,
                             colClasses = list(character = "code_merge"))

setnames(gdp_per_capita_wide, 1L, "code_merge")

gdp <- gdp_per_capita_wide[, .(
  code = code_merge,
  gdp_per_capita = rowMeans(.SD, na.rm = TRUE)
), .SDcols = paste0("GDPpc", analysis_years)]

# Total gridded GDP divided by gridded GDP per capita gives an approximate
# resident population. It is kept only as a robustness denominator and is not
# used in the index: the per-capita raster is a spatial mean, so cities with
# large rural areas get an inflated implied population.
gridded_gdp_wide <- fread(gridded_gdp_file,
                          colClasses = list(character = "code_merge"))

setnames(gridded_gdp_wide, 1L, "code_merge")

gdp_total <- gridded_gdp_wide[, .(
  code = code_merge,
  gdp_total = rowMeans(.SD, na.rm = TRUE)
), .SDcols = paste0("GDP", analysis_years)]

index <- merge(index, gdp, by = "code", all.x = TRUE)

index <- merge(index, gdp_total, by = "code", all.x = TRUE)

index[, log_gdp_per_capita := log(gdp_per_capita)]

index[, implied_population := fifelse(
  is.finite(gdp_total) & is.finite(gdp_per_capita) & gdp_per_capita > 0,
  gdp_total / gdp_per_capita, NA_real_
)]

index <- index[is.finite(sanitation_stress) & is.finite(log_gdp_per_capita)]

index[, gdp_tertile := cut(
  frank(log_gdp_per_capita, ties.method = "first"),
  breaks = quantile(
    frank(log_gdp_per_capita, ties.method = "first"),
    probs = c(0, 1 / 3, 2 / 3, 1)
  ),
  labels = c("Lower GDP", "Middle GDP", "Higher GDP"),
  include.lowest = TRUE
)]

component_correlations <- as.data.table(
  cor(index[, ..component_columns], use = "pairwise.complete.obs"),
  keep.rownames = "component"
)

stress_by_gdp <- index[, .(
  n_cities = .N,
  gdp_per_capita_median = median(gdp_per_capita),
  pressure_index_median = median(pressure_index),
  capacity_index_median = median(capacity_index),
  sanitation_stress_median = median(sanitation_stress),
  sanitation_stress_q25 = quantile(sanitation_stress, 0.25),
  sanitation_stress_q75 = quantile(sanitation_stress, 0.75),
  sewered_share_median = median(sewered_share_of_domestic, na.rm = TRUE)
), by = gdp_tertile][order(gdp_tertile)]

index_correlations <- data.table(
  quantity = c(
    "sanitation stress vs log GDP per capita",
    "pressure index vs log GDP per capita",
    "capacity index vs log GDP per capita",
    "city-only stress vs full stress",
    "sanitation stress vs unsewered share",
    "sanitation stress vs industrial COD density"
  ),
  spearman_rho = c(
    cor(index$sanitation_stress, index$log_gdp_per_capita, method = "spearman"),
    cor(index$pressure_index, index$log_gdp_per_capita, method = "spearman"),
    cor(index$capacity_index, index$log_gdp_per_capita, method = "spearman"),
    cor(index$sanitation_stress_city_only, index$sanitation_stress,
        method = "spearman"),
    cor(index$sanitation_stress, index$rural_share_of_domestic,
        method = "spearman"),
    cor(index$sanitation_stress, log1p(index$cod_industry_per_industry_km2),
        method = "spearman", use = "complete.obs")
  ),
  n_cities = nrow(index)
)

fwrite(index, file.path(output_dir, "sanitation-pressure-capacity-index.csv"))

fwrite(component_correlations,
       file.path(output_dir, "sanitation-index-component-correlations.csv"))

fwrite(stress_by_gdp,
       file.path(output_dir, "sanitation-index-by-gdp-tertile.csv"))

fwrite(index_correlations,
       file.path(output_dir, "sanitation-index-correlations.csv"))
