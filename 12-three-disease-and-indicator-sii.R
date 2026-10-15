# -*- coding: UTF-8 -*-
# ==============================================================================
# Three-category GDP-related inequality preview
#
# This script leaves the existing single-category figure unchanged. It applies
# the same annual AF and SII workflow to the three major disease categories and
# creates a 3-category x 3-question comparison figure.
# ==============================================================================

options(stringsAsFactors = FALSE, scipen = 999)

suppressPackageStartupMessages({
  library(data.table)
})

project_root <- "."

source_script <- file.path(
  project_root,
  "SII",
  "08-gastrointestinal-gdp-sii.R"
)

disease_specs <- data.table(
  disease = c(
    "Respiratory infectious diseases",
    "Gastrointestinal and enterovirus",
    "Vector-borne diseases"
  ),
  short_label = c(
    "Respiratory",
    "Gastrointestinal / enterovirus",
    "Vector-borne"
  ),
  tag = c(
    "respiratory",
    "gastrointestinal",
    "vectorborne"
  ),
  data_file = c(
    "resp_fit_allcities.rds",
    "gastrointestinal_fit_allcities.rds",
    "vectorborne_fit_allcities.rds"
  )
)

source_text_original <- readLines(
  source_script,
  encoding = "UTF-8",
  warn = FALSE
)

all_sii <- vector("list", nrow(disease_specs))

all_sii_by_gcm <- vector("list", nrow(disease_specs))

all_analysis <- vector("list", nrow(disease_specs))

for (i in seq_len(nrow(disease_specs))) {
  current_spec <- disease_specs[i]
  message(
    "Calculating annual inequality: ",
    current_spec$short_label
  )

  source_text <- source_text_original

  source_text <- sub(
    'project_root <- normalizePath(".", winslash = "/", mustWork = TRUE)',
    'project_root <- "."',
    source_text,
    fixed = TRUE
  )

  source_text <- sub(
    'disease_name <- "Gastrointestinal and enterovirus"',
    paste0('disease_name <- "', current_spec$disease, '"'),
    source_text,
    fixed = TRUE
  )

  source_text <- gsub(
    "erf_pooled_gastrointestinal.rds",
    paste0("erf_pooled_", current_spec$tag, ".rds"),
    source_text,
    fixed = TRUE
  )

  source_text <- gsub(
    "gastrointestinal_fit_allcities.rds",
    current_spec$data_file,
    source_text,
    fixed = TRUE
  )

  analysis_environment <- new.env(parent = globalenv())

  eval(
    parse(text = source_text),
    envir = analysis_environment
  )

  all_sii[[i]] <- copy(analysis_environment$sii_summary)
  all_sii_by_gcm[[i]] <- copy(analysis_environment$sii_by_gcm)
  all_analysis[[i]] <- copy(analysis_environment$analysis_data)

  rm(analysis_environment)
  gc(verbose = FALSE)
}

sii_all <- rbindlist(all_sii, use.names = TRUE, fill = TRUE)

sii_by_gcm_all <- rbindlist(
  all_sii_by_gcm,
  use.names = TRUE,
  fill = TRUE
)

analysis_all <- rbindlist(all_analysis, use.names = TRUE, fill = TRUE)

label_lookup <- setNames(
  disease_specs$short_label,
  disease_specs$disease
)

sii_all[, disease_label := unname(label_lookup[disease])]

sii_by_gcm_all[, disease_label := unname(label_lookup[disease])]

analysis_all[, disease_label := unname(label_lookup[disease])]

disease_levels <- disease_specs$short_label

sii_all[, disease_label := factor(disease_label, levels = disease_levels)]

sii_by_gcm_all[, disease_label := factor(
  disease_label,
  levels = disease_levels
)]

analysis_all[, disease_label := factor(
  disease_label,
  levels = disease_levels
)]

# Population-ranked GDP deciles for the 2019 gradient panels.
analysis_all[, gdp_decile := pmin(
  10L,
  pmax(1L, floor(socioeconomic_rank * 10) + 1L)
)]

# A concentration/Lorenz curve is not appropriate here because AF can be
# negative. Instead, the first row shows the observed GDP-decile pattern and
# the population-weighted regression line whose slope is the SII.
gradient_2019_gcm <- analysis_all[
  year == 2019L &
    is.finite(GDP) &
    is.finite(population) &
    population > 0 &
    is.finite(af_factual) &
    is.finite(af_counterfactual),
  .(
    af_factual = weighted.mean(
      af_factual,
      w = population,
      na.rm = TRUE
    ),
    af_counterfactual = weighted.mean(
      af_counterfactual,
      w = population,
      na.rm = TRUE
    )
  ),
  by = .(
    disease_label,
    gcm,
    gdp_decile
  )
]

gradient_2019 <- melt(
  gradient_2019_gcm,
  id.vars = c(
    "disease_label",
    "gcm",
    "gdp_decile"
  ),
  measure.vars = c(
    "af_counterfactual",
    "af_factual"
  ),
  variable.name = "scenario",
  value.name = "af"
)

gradient_2019[, scenario := fifelse(
  scenario == "af_counterfactual",
  "Natural-only",
  "Factual"
)]

gradient_2019 <- gradient_2019[, .(
  estimate = median(af, na.rm = TRUE)
), by = .(
  disease_label,
  scenario,
  gdp_decile
)]

gradient_2019[, socioeconomic_rank :=
  (gdp_decile - 0.5) / 10]

gradient_fit_2019_gcm <- analysis_all[
  year == 2019L &
    is.finite(socioeconomic_rank) &
    is.finite(population) &
    population > 0 &
    is.finite(af_factual) &
    is.finite(af_counterfactual),
  {
    factual_model <- lm(
      af_factual ~ socioeconomic_rank,
      data = .SD,
      weights = population
    )
    natural_model <- lm(
      af_counterfactual ~ socioeconomic_rank,
      data = .SD,
      weights = population
    )

    data.table(
      scenario = c("Natural-only", "Factual"),
      intercept = c(
        unname(coef(natural_model)["(Intercept)"]),
        unname(coef(factual_model)["(Intercept)"])
      ),
      SII = c(
        unname(coef(natural_model)["socioeconomic_rank"]),
        unname(coef(factual_model)["socioeconomic_rank"])
      )
    )
  },
  by = .(
    disease_label,
    gcm
  )
]

gradient_fit_2019 <- gradient_fit_2019_gcm[, .(
  intercept = median(intercept, na.rm = TRUE),
  SII = median(SII, na.rm = TRUE)
), by = .(
  disease_label,
  scenario
)]

gradient_2019[, scenario := factor(
  scenario,
  levels = c("Natural-only", "Factual")
)]

gradient_fit_2019[, scenario := factor(
  scenario,
  levels = c("Natural-only", "Factual")
)]

gradient_fit_2019[, `:=`(
  label = sprintf("%s SII: %+.1f pp", scenario, SII),
  label_vjust = fifelse(scenario == "Natural-only", 1.3, 2.8)
)]

ntl_wide <- fread(
  file.path(
    project_root,
    "SII",
    "data",
    "nighttime-light",
    "nighttime-light-2008-2019.csv"
  ),
  colClasses = list(character = "code_merge")
)

hfp_wide <- fread(
  file.path(
    project_root,
    "SII",
    "data",
    "human-footprint",
    "human-footprint-2008-2019.csv"
  ),
  colClasses = list(character = "code_merge")
)

ntl_long <- melt(
  ntl_wide,
  id.vars = "code_merge",
  measure.vars = paste0("NTLMean", 2008:2019),
  variable.name = "year_column",
  value.name = "indicator_value"
)

ntl_long[, `:=`(
  year = as.integer(sub("NTLMean", "", year_column)),
  indicator = "Night-time lights"
)]

ntl_long <- ntl_long[, .(
  code_merge,
  year,
  indicator,
  indicator_value = as.numeric(indicator_value)
)]

# The original annual HFP product ends in 2018. The 2019 proxy is deliberately
# excluded so that the robustness result is based only on observed HFP years.
hfp_long <- melt(
  hfp_wide,
  id.vars = "code_merge",
  measure.vars = paste0("HFP", 2008:2018),
  variable.name = "year_column",
  value.name = "indicator_value"
)

hfp_long[, `:=`(
  year = as.integer(sub("HFP", "", year_column)),
  indicator = "Human footprint"
)]

hfp_long <- hfp_long[, .(
  code_merge,
  year,
  indicator,
  indicator_value = as.numeric(indicator_value)
)]

af_indicator_base <- analysis_all[, .(
  disease_label,
  code_merge,
  year,
  gcm,
  population,
  af_factual,
  af_counterfactual,
  GDP
)]

indicator_data <- rbindlist(list(
  af_indicator_base[, .(
    disease_label,
    code_merge,
    year,
    gcm,
    population,
    af_factual,
    af_counterfactual,
    indicator = "GDP",
    indicator_value = as.numeric(GDP)
  )],
  merge(
    af_indicator_base[, .(
      disease_label,
      code_merge,
      year,
      gcm,
      population,
      af_factual,
      af_counterfactual
    )],
    ntl_long,
    by = c("code_merge", "year"),
    all = FALSE
  ),
  merge(
    af_indicator_base[, .(
      disease_label,
      code_merge,
      year,
      gcm,
      population,
      af_factual,
      af_counterfactual
    )],
    hfp_long,
    by = c("code_merge", "year"),
    all = FALSE
  )
), use.names = TRUE)

indicator_data <- indicator_data[
  is.finite(indicator_value) &
    is.finite(population) &
    population > 0 &
    is.finite(af_factual) &
    is.finite(af_counterfactual)
]

indicator_rank_base <- unique(indicator_data[, .(
  disease_label,
  code_merge,
  year,
  indicator,
  indicator_value,
  population
)])

indicator_rank_groups <- indicator_rank_base[, .(
  group_population = sum(population, na.rm = TRUE)
), by = .(
  disease_label,
  year,
  indicator,
  indicator_value
)]

setorder(
  indicator_rank_groups,
  disease_label,
  year,
  indicator,
  indicator_value
)

indicator_rank_groups[, cumulative_population_before :=
  data.table::shift(
    cumsum(group_population),
    fill = 0
  ),
by = .(
  disease_label,
  year,
  indicator
)]

indicator_rank_groups[, contextual_rank :=
  (
    cumulative_population_before +
      0.5 * group_population
  ) /
  sum(group_population),
by = .(
  disease_label,
  year,
  indicator
)]

indicator_rank_table <- merge(
  indicator_rank_base,
  indicator_rank_groups[, .(
    disease_label,
    year,
    indicator,
    indicator_value,
    contextual_rank
  )],
  by = c(
    "disease_label",
    "year",
    "indicator",
    "indicator_value"
  ),
  all.x = TRUE
)

indicator_data <- merge(
  indicator_data,
  indicator_rank_table[, .(
    disease_label,
    code_merge,
    year,
    indicator,
    contextual_rank
  )],
  by = c(
    "disease_label",
    "code_merge",
    "year",
    "indicator"
  ),
  all = FALSE
)

indicator_sii_year_gcm <- indicator_data[, {
  current_data <- .SD[
    is.finite(contextual_rank) &
      is.finite(population) &
      population > 0
  ]

  if (
    nrow(current_data) < 20L ||
      uniqueN(current_data$contextual_rank) < 2L
  ) {
    data.table(
      sii_factual = NA_real_,
      sii_natural = NA_real_,
      delta_sii = NA_real_
    )
  } else {
    factual_model <- lm(
      af_factual ~ contextual_rank,
      data = current_data,
      weights = population
    )

    natural_model <- lm(
      af_counterfactual ~ contextual_rank,
      data = current_data,
      weights = population
    )

    factual_sii <- unname(coef(factual_model)["contextual_rank"])
    natural_sii <- unname(coef(natural_model)["contextual_rank"])

    data.table(
      sii_factual = factual_sii,
      sii_natural = natural_sii,
      delta_sii = factual_sii - natural_sii
    )
  }
}, by = .(
  disease_label,
  indicator,
  year,
  gcm
)]

# First average the annual SII within each GCM, then summarise across GCMs.
indicator_sii_period_gcm <- indicator_sii_year_gcm[, .(
  sii_factual = mean(sii_factual, na.rm = TRUE),
  sii_natural = mean(sii_natural, na.rm = TRUE),
  delta_sii = mean(delta_sii, na.rm = TRUE),
  first_year = min(year),
  last_year = max(year)
), by = .(
  disease_label,
  indicator,
  gcm
)]

indicator_sii_summary <- indicator_sii_period_gcm[, .(
  factual_median = median(sii_factual, na.rm = TRUE),
  factual_lower = quantile(sii_factual, 0.025, na.rm = TRUE),
  factual_upper = quantile(sii_factual, 0.975, na.rm = TRUE),
  natural_median = median(sii_natural, na.rm = TRUE),
  natural_lower = quantile(sii_natural, 0.025, na.rm = TRUE),
  natural_upper = quantile(sii_natural, 0.975, na.rm = TRUE),
  delta_median = median(delta_sii, na.rm = TRUE),
  first_year = min(first_year),
  last_year = max(last_year)
), by = .(
  disease_label,
  indicator
)]

indicator_levels <- c(
  "GDP",
  "Night-time lights",
  "Human footprint"
)

indicator_sii_summary[, indicator := factor(
  indicator,
  levels = indicator_levels
)]

indicator_sii_summary[, row_label := paste(
  disease_label,
  indicator,
  sep = "  |  "
)]

indicator_row_levels <- unlist(lapply(
  disease_levels,
  function(disease_i) {
    paste(disease_i, indicator_levels, sep = "  |  ")
  }
))

indicator_sii_summary[, row_id :=
  length(indicator_row_levels) -
  match(row_label, indicator_row_levels) + 1L]

output_dir <- file.path(
  project_root,
  "SII",
  "output",
  "gdp-inequality-figure"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


# ==============================================================================
# Save the numerical results (Supplementary Figs. S7 and S8)
# ==============================================================================

out_dir <- file.path("SII", "output", "three-disease-inequality")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
fwrite(sii_all, file.path(out_dir, "sii-by-year-summary.csv"))
fwrite(sii_by_gcm_all, file.path(out_dir, "sii-by-year-gcm.csv"))
fwrite(gradient_2019, file.path(out_dir, "gdp-decile-af-2019.csv"))
fwrite(gradient_fit_2019, file.path(out_dir, "sii-fit-2019.csv"))
fwrite(indicator_sii_year_gcm, file.path(out_dir, "indicator-sii-by-year-gcm.csv"))
fwrite(indicator_sii_period_gcm, file.path(out_dir, "indicator-sii-period-by-gcm.csv"))
fwrite(indicator_sii_summary, file.path(out_dir, "indicator-sii-period-summary.csv"))
