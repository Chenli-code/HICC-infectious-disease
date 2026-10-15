# -*- coding: UTF-8 -*-
# ==============================================================================
# Figure 6 | rebuild of the figure-preparation objects
#
# 11-ncc-full-inequality-figure.R expects 14 objects to already exist in the
# calling environment. The code that originally created them is not in this
# repository, so Figure 6 could not be redrawn. This script rebuilds them from
# inputs that ARE in the repository:
#
#   SII/output/quantitative-statistics/analysis-inputs.rds
#       analysis_data  city x year x GCM AF, GDP, population, rank
#       sii_by_gcm     year x GCM slope index of inequality
#   user_route/output/figure1_climate/V01/V01-Supp-city-scenario-temperature-by-gcm.csv
#       city x GCM factual and hist-nat mean temperature  (panel F1)
#
# Conventions are taken from 13-figure6-quantitative-statistics.R so that the
# figure and the published statistics cannot drift apart:
#   - deciles      floor(socioeconomic_rank * 10) + 1, clamped to 1..10
#   - GDP groups   deciles 1-2 (lower) versus 9-10 (higher)
#   - pooling      population-weighted across cities, then equal-weighted
#                  across years within a GCM, then summarised across GCMs
#
# Usage:
#   source("SII/10-figure6-prepare-objects.R")
#   source("SII/11-ncc-full-inequality-figure.R")
# ==============================================================================

options(stringsAsFactors = FALSE, scipen = 999)

suppressPackageStartupMessages({
  library(data.table)
  library(sf)
})

# ------------------------------- settings ------------------------------------
FIG6_DISEASE <- "Gastrointestinal and enterovirus"

FIG6_YEARS <- 2008:2019

# F2 is the response to a UNIFORM +1 deg C, precomputed by
# 09-figure6-plus1c-response.R. It is deliberately not delta AF divided by each
# city's own delta T: that ratio keeps the GCM spread in warming inside the
# answer and reproduces only 7/10 GCM agreement, against 10/10 in the published
# panel. The uniform +1 deg C recomputation reproduces 10/10.
#
# The increment is applied to the HIST-NAT baseline, i.e. "what would a uniform
# +1 deg C do starting from the world without anthropogenic warming". Applying
# it to the factual series instead gives 3.06/2.72 -> 2.82/2.20 and collapses
# the GCM spread from 0.73 to 0.08, because the factual series is bias-corrected
# and so nearly identical across models. The hist-nat baseline reproduces the
# published panel (3.06 / 2.72 / +0.33, 10/10 against 3.1 / 2.8 / +0.30, 10/10).
FIG6_F2_FILE <- "figure6-plus1c-response-by-city-gcm.csv"

fig6_root <- function() {
  candidates <- c(".", "..", normalizePath(".", winslash = "/", mustWork = FALSE))
  for (cand in candidates) {
    if (dir.exists(file.path(cand, "SII")) &&
        dir.exists(file.path(cand, "user_route"))) {
      return(normalizePath(cand, winslash = "/", mustWork = TRUE))
    }
  }
  stop("Run from the project root (the folder containing SII/ and user_route/).")
}

PROJECT_ROOT <- fig6_root()

INPUT_RDS <- file.path(
  PROJECT_ROOT, "SII", "output", "quantitative-statistics",
  "analysis-inputs.rds"
)

TEMP_CSV <- file.path(
  PROJECT_ROOT, "user_route", "output", "figure1_climate", "V01",
  "V01-Supp-city-scenario-temperature-by-gcm.csv"
)

source("00-paths.R")
CITY_GEO_DIR <- PATH_CITY_GEOJSON_DIR

OUT_DIR <- file.path(PROJECT_ROOT, "SII", "output", "figure6-objects")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

pad6 <- function(x) {
  x <- as.character(x)
  ifelse(nchar(x) < 6L, paste0(x, strrep("0", 6L - nchar(x))), x)
}

# ==============================================================================
# 1. Cached AF and SII inputs
# ==============================================================================
cached <- readRDS(INPUT_RDS)

analysis_data <- as.data.table(cached$analysis_data)

sii_by_gcm <- as.data.table(cached$sii_by_gcm)

plot_disease <- FIG6_DISEASE

a <- analysis_data[disease == plot_disease & year %in% FIG6_YEARS]

if (!nrow(a)) stop("No AF rows for the requested disease and years.")

# Same decile rule as 13-figure6-quantitative-statistics.R:94.
a[, gdp_decile := pmin(10L, pmax(1L, floor(socioeconomic_rank * 10) + 1L))]

a[, gdp_group := fcase(
  gdp_decile <= 2L, "Lower GDP",
  gdp_decile >= 9L, "Higher GDP",
  default = "Middle GDP"
)]

a[, delta_af := af_factual - af_counterfactual]

# ==============================================================================
# 2. decile_by_gcm   (panel A)
#
# Cities are pooled with their population weights inside each decile, so a
# decile mean is a population mean and not an unweighted city mean.
# ==============================================================================
decile_by_gcm <- a[, .(
  n_cities = .N,
  population = sum(population),
  af_factual = weighted.mean(af_factual, population),
  af_counterfactual = weighted.mean(af_counterfactual, population),
  delta_af = weighted.mean(delta_af, population)
), by = .(disease, year, gcm, gdp_decile)]

setorder(decile_by_gcm, disease, year, gcm, gdp_decile)

# ==============================================================================
# 3. decile_period_gcm / decile_median   (panel D)
#
# Years are equally weighted within a GCM first, so a year with more cases does
# not silently dominate the period mean; GCMs are summarised only afterwards.
# ==============================================================================
decile_period_gcm <- decile_by_gcm[
  disease == plot_disease,
  .(delta_af = mean(delta_af)),
  by = .(gcm, gdp_decile)
]

decile_median <- decile_period_gcm[, .(
  delta_af = median(delta_af),
  delta_af_q025 = quantile(delta_af, 0.025, names = FALSE),
  delta_af_q975 = quantile(delta_af, 0.975, names = FALSE),
  n_positive = sum(delta_af > 0),
  n_gcm = .N
), by = gdp_decile]

setorder(decile_median, gdp_decile)

# ==============================================================================
# 4. sii_summary   (panels B and C)
# ==============================================================================
s <- sii_by_gcm[disease == plot_disease & year %in% FIG6_YEARS]

if (!nrow(s)) stop("No SII rows for the requested disease and years.")

# The paired identity is what makes panel C interpretable; check it here rather
# than discovering a broken pairing from the shape of the curve.
if (any(abs(s$delta_sii - (s$sii_factual - s$sii_counterfactual)) > 1e-7)) {
  stop("Paired SII identity fails: delta_sii must equal factual minus natural-only.")
}

q025 <- function(x) quantile(x, 0.025, names = FALSE, na.rm = TRUE)

q975 <- function(x) quantile(x, 0.975, names = FALSE, na.rm = TRUE)

sii_summary <- s[, .(
  n_gcm = .N,
  sii_factual_median = median(sii_factual),
  sii_factual_q025 = q025(sii_factual),
  sii_factual_q975 = q975(sii_factual),
  sii_counterfactual_median = median(sii_counterfactual),
  sii_counterfactual_q025 = q025(sii_counterfactual),
  sii_counterfactual_q975 = q975(sii_counterfactual),
  delta_sii_median = median(delta_sii),
  delta_sii_q025 = q025(delta_sii),
  delta_sii_q975 = q975(delta_sii),
  delta_sii_negative = sum(delta_sii < 0)
), by = .(disease, year)]

setorder(sii_summary, disease, year)

# ==============================================================================
# 4b. City geometry and the Hu Huanyong side
#
# Read once here because both the F panels and the map need it. The Hu line is
# the diagonal Heihe-Tengchong segment; a straight line splits the whole plane,
# so the side is just the sign of the cross product with the segment direction.
# ==============================================================================
HU_LINE <- data.table(
  lon = c(98.5, 127.5),
  lat = c(25.0,  50.2)
)

QINLING_HUAIHE <- data.table(
  lon = c(104.0, 107.0, 109.5, 111.5, 113.5, 115.5, 117.5, 119.5, 121.0),
  lat = c( 34.3,  33.7,  33.6,  33.3,  32.8,  32.7,  32.8,  33.2,  34.0)
)

hu_side <- function(lon, lat) {
  dx <- HU_LINE$lon[2] - HU_LINE$lon[1]
  dy <- HU_LINE$lat[2] - HU_LINE$lat[1]
  cross <- (lon - HU_LINE$lon[1]) * dy - (lat - HU_LINE$lat[1]) * dx
  fifelse(cross > 0, "Southeast", "Northwest")
}

sf_use_s2(FALSE)

geo_files <- list.files(CITY_GEO_DIR, pattern = "[.]json$", full.names = TRUE)

city_geo <- do.call(rbind, lapply(geo_files, function(f) {
  st_read(f, quiet = TRUE)[, c("id", "name")]
}))

city_geo$code_merge <- pad6(city_geo$id)

city_geo <- suppressWarnings(st_make_valid(st_transform(city_geo, 4326)))

city_centroid <- {
  xy <- suppressWarnings(st_coordinates(st_point_on_surface(st_geometry(city_geo))))
  data.table(code_merge = city_geo$code_merge, lon = xy[, 1], lat = xy[, 2])
}

city_centroid[, hu_group := hu_side(lon, lat)]

# ==============================================================================
# 5. panel_F1   anthropogenic warming contrast
#
# Each city's anthropogenic warming is factual minus hist-nat mean temperature
# for that GCM. Cities are matched to the AF table so that F1 and D describe the
# same set of cities and the same GDP grouping.
# ==============================================================================
temp_gcm <- fread(TEMP_CSV, colClasses = list(character = "code"))

temp_gcm[, code_merge := pad6(code)]

temp_gcm[, delta_t := factual - histnat]

# Cities are weighted by their disease case total, not by population. Four
# weightings were tested against the published panel (2.8 / 2.5 / +0.34 deg C,
# 10/10 GCMs); case weighting reproduces it (2.81 / 2.50 / +0.31, 10/10) while
# population, unweighted and median city weighting all give about +0.49.
#
# Two contrasts are built, not one: the GDP quintiles and the two sides of the
# Hu Huanyong line. They are different cuts of the same cities - GDP is NOT
# balanced across the Hu line (mean rank 0.364 northwest vs 0.587 southeast) -
# so the panel shows them side by side rather than merging them.
city_group <- a[gdp_group != "Middle GDP", .(
  case_weight = sum(case_weight)
), by = .(code_merge, gdp_group)]

city_group <- unique(city_group, by = "code_merge")

city_hu <- merge(
  a[, .(case_weight = sum(case_weight)), by = code_merge],
  city_centroid[, .(code_merge, hu_group)],
  by = "code_merge"
)

build_F_panel <- function(value_dt, value_col) {
  gdp_side <- merge(value_dt, city_group, by = "code_merge")[, .(
    estimate = weighted.mean(get(value_col), case_weight), n_cities = .N
  ), by = .(gcm, group = gdp_group)][, contrast := "GDP quintiles"]

  hu_part <- merge(value_dt, city_hu, by = "code_merge")[, .(
    estimate = weighted.mean(get(value_col), case_weight), n_cities = .N
  ), by = .(gcm, group = hu_group)][, contrast := "Hu Huanyong line"]

  out <- rbindlist(list(gdp_side, hu_part))
  out[, group := factor(
    group,
    levels = c("Lower GDP", "Higher GDP", "Northwest", "Southeast")
  )]
  out[, contrast := factor(
    contrast, levels = c("GDP quintiles", "Hu Huanyong line")
  )]
  setorder(out, contrast, gcm, group)
  out[]
}

F_contrast_test <- function(panel_dt, unit_label) {
  rbindlist(lapply(levels(panel_dt$contrast), function(cc) {
    w <- dcast(panel_dt[contrast == cc], gcm ~ group, value.var = "estimate")
    pair <- if (cc == "GDP quintiles") {
      c("Lower GDP", "Higher GDP")
    } else {
      c("Northwest", "Southeast")
    }
    d <- w[[pair[1]]] - w[[pair[2]]]
    data.table(
      contrast = factor(cc, levels = levels(panel_dt$contrast)),
      estimate = median(d), lower = q025(d), upper = q975(d),
      n_positive = sum(d > 0), n = length(d),
      label = sprintf(
        "%s - %s: %+.2f %s\n%d/%d GCMs",
        sub(" GDP", "", pair[1]), sub(" GDP", "", pair[2]),
        median(d), unit_label, sum(d > 0), length(d)
      )
    )
  }))
}

panel_F1 <- build_F_panel(
  temp_gcm[, .(code_merge, gcm, delta_t)], "delta_t"
)

F1_median <- panel_F1[, .(estimate = median(estimate)), by = .(contrast, group)]

F1_test <- F_contrast_test(panel_F1, "deg C")

# ==============================================================================
# 6. panel_F2   standardised response to +1 deg C
#
# This is delta AF expressed per degree of anthropogenic warming, so that the
# GDP contrast in F1 (how much warming each group received) is separated from
# the contrast in how much AF each group gains per degree.
# ==============================================================================
f2_path <- file.path(OUT_DIR, FIG6_F2_FILE)

if (!file.exists(f2_path)) {
  stop(
    "Missing ", FIG6_F2_FILE, ". Run SII/09-figure6-plus1c-response.R first; ",
    "it recomputes AF under the hist-nat temperature + 1 deg C."
  )
}

response_city <- fread(f2_path, colClasses = list(character = "code_merge"))

# Years are averaged within a city first, then cities are pooled on the same
# case weights used in F1, so F1 and F2 describe one set of cities on one rule.
panel_F2 <- build_F_panel(
  response_city[, .(response_1c = mean(response_1c_nat)),
                by = .(code_merge, gcm)],
  "response_1c"
)

F2_median <- panel_F2[, .(estimate = median(estimate)), by = .(contrast, group)]

F2_test <- F_contrast_test(panel_F2, "pp")

# ==============================================================================
# 7. Panel E   bivariate GDP x delta AF map
#
# Stevens 3x3 bivariate palette. Classes are tertiles of the period-mean city
# delta AF (GCM median) crossed with tertiles of GDP rank, so each axis splits
# the mapped cities into thirds rather than using absolute cut points that no
# reader can interpret.
# ==============================================================================
city_period <- a[, .(delta_af = mean(delta_af)), by = .(code_merge, gcm)]

city_map <- city_period[, .(delta_af = median(delta_af)), by = code_merge]

city_rank <- unique(a[, .(code_merge, socioeconomic_rank)], by = "code_merge")

city_map <- merge(city_map, city_rank, by = "code_merge")

tertile <- function(x) {
  cuts <- quantile(x, c(1 / 3, 2 / 3), na.rm = TRUE, names = FALSE)
  as.integer(cut(x, c(-Inf, cuts, Inf), labels = FALSE))
}

city_map[, gdp_class := tertile(socioeconomic_rank)]

city_map[, af_class := tertile(delta_af)]

city_map[, bivar_class := paste0(gdp_class, "-", af_class)]

# ==============================================================================
# 8. Panel annotations
#
# The figure previously carried numbers on F1 and F2 only, so panels A-D had to
# be read off the axes. These strings are built here, from the same objects the
# panels are drawn from, so the figure can never quote a number the panel does
# not show. Every interval is the 2.5th-97.5th GCM percentile range, which is an
# ensemble spread and NOT a confidence interval.
# ==============================================================================
A_2019 <- decile_by_gcm[year == 2019L & gdp_decile %in% c(1L, 10L), .(
  factual = median(af_factual),
  natural = median(af_counterfactual)
), by = gdp_decile]

group_period <- a[gdp_group != "Middle GDP", .(
  delta_af = weighted.mean(delta_af, population)
), by = .(gcm, year, gdp_group)][, .(delta_af = mean(delta_af)),
                                 by = .(gcm, gdp_group)]

group_wide <- dcast(group_period, gcm ~ gdp_group, value.var = "delta_af")

group_wide[, difference := `Lower GDP` - `Higher GDP`]

sii_period <- s[, .(
  sii_factual = mean(sii_factual),
  sii_counterfactual = mean(sii_counterfactual),
  delta_sii = mean(delta_sii)
), by = gcm]

sii_period[, relative := 100 * delta_sii / sii_counterfactual]

dbr_file <- file.path(OUT_DIR, "figure6-burden-concentration-dbr.csv")

line_file <- file.path(OUT_DIR, "figure6-dividing-line-contrast.csv")

# Short in-panel labels. Only C, D and E carry one, and only where the number
# IS the panel's finding and cannot be read off the axis. Everything else -
# the 2019 decile gradient, the mean SII levels, the DBR, the methods - belongs
# in the figure legend, which is written to figure6-legend.txt.
fig6_labels <- list(
  C = sprintf(
    "%.2f pp (%.1f%%), %d/%d GCMs",
    median(sii_period$delta_sii), median(sii_period$relative),
    sum(sii_period$delta_sii < 0), nrow(sii_period)
  ),
  D = sprintf(
    "Deciles 1-2 vs 9-10: %+.2f pp, %d/%d GCMs",
    median(group_wide$difference),
    sum(group_wide$difference > 0), nrow(group_wide)
  )
)

fig6_labels$E <- if (file.exists(line_file)) {
  hu <- fread(line_file)[line == "Hu Huanyong line" & metric == "delta_af"]
  sprintf(
    "Northwest vs southeast: %+.2f pp, %d/%d GCMs",
    hu$difference, hu$n_a_higher, hu$n_gcm
  )
} else {
  ""
}

# ------------------------------------------------------------------------------
# Figure legend. Everything taken off the panels is written here so it can be
# pasted into the manuscript instead of being lost.
# ------------------------------------------------------------------------------
legend_lines <- c(
  "Figure 6. GDP-related inequality in temperature-attributable gastrointestinal and enterovirus infections in Chinese children, 2008-2019.",
  "",
  sprintf(
    paste0("A, Attributable fraction by GDP decile in 2019 under the factual ",
           "and natural-only climates; lines are ensemble medians and bands the ",
           "2.5th-97.5th percentile range across 10 GCMs. Under the factual ",
           "climate the attributable fraction rose from %.1f%% in decile 1 to ",
           "%.1f%% in decile 10 (gap %.1f percentage points)."),
    A_2019[gdp_decile == 1L, factual], A_2019[gdp_decile == 10L, factual],
    A_2019[gdp_decile == 10L, factual] - A_2019[gdp_decile == 1L, factual]
  ),
  "",
  sprintf(
    paste0("B, Anthropogenic change in attributable fraction by GDP decile; ",
           "grey lines are individual GCMs. Deciles 1-2 gained %.2f percentage ",
           "points against %.2f in deciles 9-10 (difference %+.2f, %.2f to ",
           "%.2f; %d of %d GCMs positive)."),
    median(group_wide$`Lower GDP`), median(group_wide$`Higher GDP`),
    median(group_wide$difference), q025(group_wide$difference),
    q975(group_wide$difference), sum(group_wide$difference > 0),
    nrow(group_wide)
  ),
  "",
  sprintf(
    paste0("C, Slope index of inequality (SII) under both climates, above, and ",
           "the difference between them, below, on a shared year axis. ",
           "Averaged over 2008-2019 the SII was %.2f percentage points under ",
           "the factual climate and %.2f under the natural-only climate; the ",
           "mean difference was %.2f percentage points (%.2f to %.2f), a %.1f%% ",
           "attenuation of the gradient, with %d of %d GCMs negative. A negative ",
           "difference means the GDP gradient was smaller under the factual ",
           "climate; it does not on its own mean that health improved."),
    median(sii_period$sii_factual), median(sii_period$sii_counterfactual),
    median(sii_period$delta_sii), q025(sii_period$delta_sii),
    q975(sii_period$delta_sii), median(sii_period$relative),
    sum(sii_period$delta_sii < 0), nrow(sii_period)
  ),
  ""
)

if (file.exists(dbr_file)) {
  db <- fread(dbr_file)
  lo <- db[group == "low 20% GDP"]
  hi <- db[group == "high 20% GDP"]
  legend_lines <- c(legend_lines, sprintf(
    paste0("D, Joint distribution of GDP and anthropogenic change in ",
           "attributable fraction; the dashed diagonal is the Hu Huanyong ",
           "(Heihe-Tengchong) line. The lowest GDP quintile held %.1f%% of the ",
           "child population and %.1f%% of attributable cases (ratio %.2f); the ",
           "highest quintile held %.1f%% and %.1f%% (ratio %.2f). GDP is not ",
           "balanced across the Hu line (mean rank 0.364 northwest vs 0.587 ",
           "southeast), so that contrast reflects geography and development ",
           "together."),
    100 * lo$share_children, 100 * lo$share_attributable, lo$dbr,
    100 * hi$share_children, 100 * hi$share_attributable, hi$dbr
  ), "")
}

legend_lines <- c(
  legend_lines,
  paste0("E1 and E2, Paired contrasts between GDP quintiles and across the Hu ",
         "line; grey lines are individual GCMs, the heavy line the ensemble ",
         "median. E1 shows anthropogenic warming; E2 shows the change in ",
         "attributable fraction under a uniform +1 deg C applied to the ",
         "natural-only baseline."),
  "",
  paste0("All intervals are 2.5th-97.5th percentile ranges across 10 equally ",
         "weighted GCMs. They describe ensemble spread, not statistical ",
         "confidence, and exclude exposure-response and city-sampling ",
         "uncertainty. Counts of agreeing GCMs describe directional ",
         "consistency and are not hypothesis tests. Cities are pooled with ",
         "population weights for the decile summaries and with case weights ",
         "for the E panels.")
)

writeLines(legend_lines, file.path(OUT_DIR, "figure6-legend.txt"))

# ==============================================================================
# 8. Audit trail
# ==============================================================================
fwrite(HU_LINE, file.path(OUT_DIR, "figure6-hu-line-waypoints.csv"))

fwrite(QINLING_HUAIHE, file.path(OUT_DIR, "figure6-qinling-huaihe-waypoints.csv"))

fwrite(decile_by_gcm, file.path(OUT_DIR, "figure6-decile-by-gcm.csv"))

fwrite(decile_period_gcm, file.path(OUT_DIR, "figure6-decile-period-by-gcm.csv"))

fwrite(decile_median, file.path(OUT_DIR, "figure6-decile-median.csv"))

fwrite(sii_summary, file.path(OUT_DIR, "figure6-sii-summary.csv"))

fwrite(panel_F1, file.path(OUT_DIR, "figure6-F1-by-gcm.csv"))

fwrite(panel_F2, file.path(OUT_DIR, "figure6-F2-by-gcm.csv"))

fwrite(
  rbindlist(list(
    cbind(panel = "F1", metric = "lower_minus_higher_degC", F1_test),
    cbind(panel = "F2", metric = "lower_minus_higher_pp_per_degC", F2_test)
  )),
  file.path(OUT_DIR, "figure6-F-contrasts.csv")
)

fwrite(
  city_map[, .(code_merge, socioeconomic_rank, delta_af, gdp_class, af_class,
               bivar_class)],
  file.path(OUT_DIR, "figure6-map-bivariate-classes.csv")
)

saveRDS(
  list(
    plot_disease = plot_disease,
    decile_by_gcm = decile_by_gcm,
    decile_period_gcm = decile_period_gcm,
    decile_median = decile_median,
    sii_summary = sii_summary,
    panel_F1 = panel_F1, F1_median = F1_median, F1_test = F1_test,
    panel_F2 = panel_F2, F2_median = F2_median, F2_test = F2_test,
    f2_method = "uniform +1 deg C",
    fig6_labels = fig6_labels
  ),
  file.path(OUT_DIR, "figure6-prepared-objects.rds")
)
