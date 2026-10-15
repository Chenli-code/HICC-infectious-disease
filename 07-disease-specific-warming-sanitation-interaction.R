# -*- coding: UTF-8 -*-
###############################################################################
## Supplementary: warming x sanitation interaction, disease by disease
##
## The city-level anthropogenic attributable fraction already exists for the
## three disease classes but not for individual diseases, so it is computed
## here from the stored city BLUP coefficients with the shared two-stage
## attribution helpers. No exposure-response model is refitted.
##
## For every disease with a prepared two-stage fit:
##   1. city x GCM anthropogenic attributable fraction, using the paired
##      factual versus natural-only contrast;
##   2. the GCM-median city value is regressed on
##        delta AF ~ delta T * sanitation stress + GDP + city climate;
##   3. the interaction term is shown as a forest plot by transmission route.
##
## Set MECHANISM_DISEASE_REUSE=false to recompute the attribution cache.
###############################################################################

options(stringsAsFactors = FALSE, scipen = 999)

suppressPackageStartupMessages({
  library(data.table)
  library(dlnm)
})

project_root <- normalizePath(".", winslash = "/", mustWork = TRUE)

Sys.setenv(ATTR_PROJECT_ROOT = project_root)

source(file.path(project_root, "R", "attribution-core.R"))

source(file.path(project_root, "two_stage_analysis", "R", "two-stage-core.R"))

output_dir <- file.path(project_root, "Mechanism", "output")

cache_dir <- file.path(output_dir, "disease_city_attribution")

dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

analysis_dates <- seq.Date(as.Date("2008-01-01"), as.Date("2019-12-31"),
                           by = "day")

gcm_names <- c(
  "ACCESS-CM2", "ACCESS-ESM1-5", "BCC-CSM2-MR", "CanESM5", "FGOALS-g3",
  "GFDL-ESM4", "ipsl-cm6a-lr", "MIROC6", "MRI-ESM2-0", "NorESM2-LM"
)

minimum_cities <- 60L

reuse <- tolower(Sys.getenv("MECHANISM_DISEASE_REUSE", "true")) == "true"

climate_dir <- file.path(project_root, "full_analysis", "cache",
                         "climate_pairs")

disease_fit_dir <- file.path(
  project_root, "two_stage_analysis", "output", "individual_diseases", "fits"
)

index_file <- file.path(output_dir, "sanitation-pressure-capacity-index.csv")

city_temperature_file <- file.path(
  project_root, "user_route", "output", "figure1_climate", "V01",
  "V01-Supp-city-scenario-temperature-by-gcm.csv"
)

standardise <- function(x) {
  usable <- is.finite(x)
  result <- rep(NA_real_, length(x))
  result[usable] <- (x[usable] - mean(x[usable])) / sd(x[usable])
  result
}

disease_files <- list.files(
  disease_fit_dir, pattern = "^prepared-mixmeta-BLUP-", full.names = TRUE
)

disease_fits <- lapply(disease_files, readRDS)

names(disease_fits) <- vapply(disease_fits, function(fit) fit$disease,
                              character(1))

all_codes <- sort(unique(unlist(lapply(
  disease_fits, function(fit) names(fit$coefficient)
))))

attribution_for_gcm <- function(gcm) {
  cache_path <- file.path(cache_dir, paste0("disease-city-", gcm, ".rds"))
  if (reuse && file.exists(cache_path)) {
    message("Reusing ", basename(cache_path))
    return(readRDS(cache_path))
  }

  message("Attribution: ", gcm)
  pair <- readRDS(file.path(climate_dir, paste0("climate_pair_", gcm, ".rds")))
  climate_all <- acc_complete_climate_pair(
    pair = pair, codes = all_codes, dates = analysis_dates
  )
  rm(pair)
  invisible(gc(verbose = FALSE))

  disease_rows <- lapply(names(disease_fits), function(disease) {
    fit <- disease_fits[[disease]]
    codes <- names(fit$coefficient)
    counts <- as.data.table(fit$dat)[, .(code, date, n)]
    analysis <- merge(
      climate_all[code %chin% codes], counts,
      by = c("code", "date"), all.x = TRUE
    )
    analysis[is.na(n), n := 0]
    setorder(analysis, code, date)
    setkey(analysis, code)

    city_rows <- lapply(codes, function(city) {
      city_data <- analysis[.(city)]
      daily <- ts_city_daily_attribution(
        factual_temperature = city_data$tf,
        histnat_temperature = city_data$tn,
        cases = city_data$n,
        dates = city_data$date,
        coefficient = fit$coefficient[[city]],
        argvar = fit$argvar,
        lag = fit$lag,
        reference = fit$reference[[city]],
        tag = fit$tag
      )
      summary_row <- ts_summarise_city_daily(
        daily, sum(city_data$n, na.rm = TRUE)
      )
      summary_row[, code := city]
      summary_row[]
    })

    result <- rbindlist(city_rows, fill = TRUE)
    result[, `:=`(disease = disease, class_tag = fit$tag)]
    result[]
  })

  result <- rbindlist(disease_rows, fill = TRUE)
  result[, gcm := gcm]
  saveRDS(result, cache_path, compress = TRUE)
  rm(climate_all)
  invisible(gc(verbose = FALSE))
  result
}

attribution <- rbindlist(lapply(gcm_names, attribution_for_gcm), fill = TRUE)

if (!"af_acc_direct" %chin% names(attribution)) {
  attribution[, af_acc_direct := fifelse(
    total_cases > 0, 100 * an_acc_direct / total_cases, NA_real_
  )]
}

fwrite(attribution,
       file.path(output_dir, "disease-city-attribution-by-gcm.csv"))

city_temperature <- fread(city_temperature_file,
                          colClasses = list(character = "code"))

city_temperature[, delta_t := factual - histnat]

index <- fread(index_file, colClasses = list(character = c("code",
                                                           "province_code")))

city_meta <- rbindlist(lapply(names(disease_fits), function(disease) {
  meta <- as.data.table(disease_fits[[disease]]$city_meta)[
    , .(code, mean_temp, temp_iqr)
  ]
  meta[, disease := disease]
  meta[]
}))

panel <- merge(
  attribution[, .(disease, class_tag, code, gcm, total_cases, af_acc_direct)],
  city_temperature[, .(code, gcm, delta_t)],
  by = c("code", "gcm"), all.x = TRUE
)

city_median <- panel[, .(
  total_cases = first(total_cases),
  delta_af = median(af_acc_direct, na.rm = TRUE),
  delta_t = median(delta_t, na.rm = TRUE)
), by = .(disease, class_tag, code)]

city_median <- merge(city_median, city_meta, by = c("disease", "code"),
                     all.x = TRUE)

city_median <- merge(
  city_median,
  index[, .(code, sanitation_stress, pressure_index, log_gdp_per_capita)],
  by = "code", all.x = TRUE
)

city_median <- city_median[
  is.finite(delta_af) & is.finite(delta_t) & is.finite(sanitation_stress) &
    is.finite(log_gdp_per_capita) & is.finite(total_cases) & total_cases > 0
]

estimates <- city_median[, {
  current <- copy(.SD)
  if (uniqueN(current$code) < minimum_cities) {
    data.table()
  } else {
    current[, `:=`(
      delta_t_z = standardise(delta_t),
      sanitation_z = standardise(sanitation_stress),
      gdp_z = standardise(log_gdp_per_capita),
      mean_temp_z = standardise(mean_temp),
      temp_iqr_z = standardise(temp_iqr)
    )]
    model <- lm(
      delta_af ~ delta_t_z * sanitation_z + gdp_z + mean_temp_z + temp_iqr_z,
      data = current, weights = total_cases
    )
    full_model <- lm(
      delta_af ~ delta_t_z * (sanitation_z + gdp_z + mean_temp_z + temp_iqr_z),
      data = current, weights = total_cases
    )
    coefficients <- summary(model)$coefficients
    confidence <- confint(model)
    full_coefficients <- summary(full_model)$coefficients
    data.table(
      n_cities = nrow(current),
      total_cases = sum(current$total_cases),
      mean_delta_af = weighted.mean(current$delta_af, current$total_cases),
      estimate = coefficients["delta_t_z:sanitation_z", "Estimate"],
      lower = confidence["delta_t_z:sanitation_z", 1],
      upper = confidence["delta_t_z:sanitation_z", 2],
      p_value = coefficients["delta_t_z:sanitation_z", "Pr(>|t|)"],
      full_estimate = full_coefficients["delta_t_z:sanitation_z", "Estimate"],
      full_p_value = full_coefficients["delta_t_z:sanitation_z", "Pr(>|t|)"]
    )
  }
}, by = disease]

pathway_groups <- data.table(
  disease = c(
    "Infectious diarrhoea", "Dysentery", "Hepatitis A",
    "Typhoid and paratyphoid",
    "Hand, foot, and mouth disease", "Acute haemorrhagic conjunctivitis",
    "Seasonal influenza", "Mumps", "Rubella", "Scarlet fever", "Tuberculosis"
  ),
  pathway = c(
    rep("Faecal-oral and water-related", 4),
    rep("Other enterovirus and contact", 2),
    rep("Respiratory and person-to-person", 5)
  )
)

estimates <- merge(estimates, pathway_groups, by = "disease", all.x = TRUE)

estimates[is.na(pathway), pathway := "Other"]

setorder(estimates, pathway, -estimate)

fwrite(city_median,
       file.path(output_dir, "disease-city-interaction-data.csv"))

fwrite(estimates,
       file.path(output_dir, "disease-warming-sanitation-interaction.csv"))
