# -*- coding: UTF-8 -*-
###############################################################################
## Supplementary: is the sanitation signal specific to faecal-oral infections?
##
## For every disease and disease class that already has a fitted two-stage
## model, the city exposure-response function is evaluated on one common
## national warm-season temperature distribution, warmed by 1 degree Celsius.
## The resulting standardised response is regressed on sanitation stress with
## GDP and city climate held constant, and the sanitation coefficient is shown
## as a forest plot grouped by transmission route.
##
## Nothing is refitted: the city BLUP coefficients are read from
##   two_stage_analysis/output/individual_diseases/fits/prepared-mixmeta-BLUP-*
##   full_analysis/output/{gastrointestinal,resp,vectorborne}_fit_allcities.rds
###############################################################################

options(stringsAsFactors = FALSE, scipen = 999)

suppressPackageStartupMessages({
  library(data.table)
  library(dlnm)
})

project_root <- normalizePath(".", winslash = "/", mustWork = TRUE)

output_dir <- file.path(project_root, "Mechanism", "output")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

warm_months <- 6:9

minimum_cities <- 60L

index_file <- file.path(output_dir, "sanitation-pressure-capacity-index.csv")

class_fit_files <- c(
  gastrointestinal = file.path(project_root, "full_analysis", "output",
                               "gastrointestinal_fit_allcities.rds"),
  respiratory = file.path(project_root, "full_analysis", "output",
                          "resp_fit_allcities.rds"),
  vectorborne = file.path(project_root, "full_analysis", "output",
                          "vectorborne_fit_allcities.rds")
)

disease_fit_dir <- file.path(
  project_root, "two_stage_analysis", "output", "individual_diseases", "fits"
)

pathway_groups <- data.table(
  label = c(
    "Infectious diarrhoea", "Dysentery", "Hepatitis A",
    "Typhoid and paratyphoid",
    "Hand, foot, and mouth disease", "Acute haemorrhagic conjunctivitis",
    "Seasonal influenza", "Mumps", "Rubella", "Scarlet fever", "Tuberculosis",
    "Gastrointestinal and enterovirus (class)",
    "Respiratory infections (class)", "Vector-borne diseases (class)"
  ),
  pathway = c(
    rep("Faecal-oral and water-related", 4),
    rep("Other enterovirus and contact", 2),
    rep("Respiratory and person-to-person", 5),
    "Disease class", "Disease class", "Disease class"
  )
)

standardise <- function(x) {
  usable <- is.finite(x)
  result <- rep(NA_real_, length(x))
  result[usable] <- (x[usable] - mean(x[usable])) / sd(x[usable])
  result
}

reference_fit <- readRDS(class_fit_files[["gastrointestinal"]])

reference_daily <- as.data.table(reference_fit$dat)[, .(date, tmean)]

reference_daily[, month := as.integer(format(date, "%m"))]

reference_grid <- as.numeric(quantile(
  reference_daily[month %in% warm_months & is.finite(tmean), tmean],
  probs = seq(0.005, 0.995, by = 0.01), na.rm = TRUE, names = FALSE
))

standardised_response <- function(coefficients, argvar) {
  difference <-
    do.call(onebasis, c(list(x = reference_grid + 1), argvar)) -
    do.call(onebasis, c(list(x = reference_grid), argvar))
  delta_eta <- drop(difference %*% as.numeric(coefficients))
  100 * mean(1 - exp(-delta_eta))
}

index <- fread(index_file, colClasses = list(character = c("code",
                                                           "province_code")))

index <- index[, .(code, sanitation_stress, pressure_index,
                   capacity_index, log_gdp_per_capita,
                   cod_industry_per_industry_km2)]

collect_responses <- function(coefficient_list, argvar, city_meta, label) {
  codes <- names(coefficient_list)
  responses <- data.table(
    code = codes,
    response = vapply(
      codes,
      function(city) standardised_response(coefficient_list[[city]], argvar),
      numeric(1)
    )
  )
  meta <- as.data.table(city_meta)[, .(code, total_cases, mean_temp, temp_iqr)]
  merged <- merge(responses, meta, by = "code")
  merged[, label := label]
  merged
}

response_records <- list()

for (class_tag in names(class_fit_files)) {
  fit <- readRDS(class_fit_files[[class_tag]])
  label <- switch(
    class_tag,
    gastrointestinal = "Gastrointestinal and enterovirus (class)",
    respiratory = "Respiratory infections (class)",
    vectorborne = "Vector-borne diseases (class)"
  )
  coefficient_list <- lapply(fit$blups, function(entry) entry$blup)
  response_records[[label]] <- collect_responses(
    coefficient_list, fit$argvar, fit$city_meta, label
  )
  message("  ", label, ": ", nrow(response_records[[label]]), " cities")
  rm(fit)
  invisible(gc(verbose = FALSE))
}

disease_files <- list.files(disease_fit_dir, pattern = "^prepared-mixmeta-BLUP-",
                            full.names = TRUE)

for (path in disease_files) {
  fit <- readRDS(path)
  label <- fit$disease
  response_records[[label]] <- collect_responses(
    fit$coefficient, fit$argvar, fit$city_meta, label
  )
  message("  ", label, ": ", nrow(response_records[[label]]), " cities")
  rm(fit)
  invisible(gc(verbose = FALSE))
}

responses <- rbindlist(response_records, use.names = TRUE)

responses <- merge(responses, index, by = "code")

estimates <- responses[, {
  current <- .SD[is.finite(response) & is.finite(sanitation_stress) &
                   is.finite(log_gdp_per_capita) & is.finite(total_cases) &
                   total_cases > 0]
  if (nrow(current) < minimum_cities) {
    data.table()
  } else {
    current[, `:=`(
      sanitation_z = standardise(sanitation_stress),
      pressure_z = standardise(pressure_index),
      gdp_z = standardise(log_gdp_per_capita),
      mean_temp_z = standardise(mean_temp),
      temp_iqr_z = standardise(temp_iqr)
    )]
    model <- lm(response ~ sanitation_z + gdp_z + mean_temp_z + temp_iqr_z,
                data = current, weights = total_cases)
    pressure_model <- lm(
      response ~ pressure_z + gdp_z + mean_temp_z + temp_iqr_z,
      data = current, weights = total_cases
    )
    coefficients <- summary(model)$coefficients
    confidence <- confint(model)
    pressure_coefficients <- summary(pressure_model)$coefficients
    data.table(
      n_cities = nrow(current),
      total_cases = sum(current$total_cases),
      mean_response = weighted.mean(current$response, current$total_cases),
      estimate = coefficients["sanitation_z", "Estimate"],
      lower = confidence["sanitation_z", 1],
      upper = confidence["sanitation_z", 2],
      p_value = coefficients["sanitation_z", "Pr(>|t|)"],
      pressure_estimate = pressure_coefficients["pressure_z", "Estimate"],
      pressure_p_value = pressure_coefficients["pressure_z", "Pr(>|t|)"]
    )
  }
}, by = label]

estimates <- merge(estimates, pathway_groups, by = "label", all.x = TRUE)

estimates[is.na(pathway), pathway := "Other"]

# The response scale differs between diseases, so a relative version is also
# reported: the sanitation coefficient divided by the mean response.
estimates[, relative_estimate := fifelse(
  abs(mean_response) > 0.1, estimate / mean_response, NA_real_
)]

setorder(estimates, pathway, -estimate)

fwrite(responses, file.path(output_dir, "disease-standardised-responses.csv"))

fwrite(estimates,
       file.path(output_dir, "disease-sanitation-response-estimates.csv"))
