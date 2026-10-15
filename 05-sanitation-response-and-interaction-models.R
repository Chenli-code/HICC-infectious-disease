# -*- coding: UTF-8 -*-
###############################################################################
## Figure 7b-d: sanitation-system stress and the anthropogenic burden
##
## Everything here reuses estimates that already exist in the project. No
## exposure-response model is refitted.
##
##   Outcome 1  city anthropogenic attributable fraction
##              two_stage_analysis/output/major_classes/two-stage-city-by-gcm.csv
##              (af_acc_direct, city BLUP two-stage, per GCM)
##   Outcome 2  city response to a uniform +1 degree Celsius, computed from the
##              stored city BLUP coefficients; two versions
##                own-climate     the city's own daily series is warmed
##                standardised    a common national warm-season temperature
##                                distribution is warmed, so that between-city
##                                differences reflect the response function
##                                only
##   Exposure   city anthropogenic warming, factual minus natural-only
##              user_route/output/figure1_climate/V01/
##              V01-Supp-city-scenario-temperature-by-gcm.csv
##   Context    Mechanism/output/sanitation-pressure-capacity-index.csv
##              (script 04), GDP per capita, city climate metadata
##
## Panels
##   7b  delta AF ~ delta T * sanitation stress + GDP + climate
##   7c  response to +1 C ~ sanitation stress
##   7d  GDP coefficient on the response, crude / + climate / + sanitation
###############################################################################

options(stringsAsFactors = FALSE, scipen = 999)

suppressPackageStartupMessages({
  library(data.table)
  library(dlnm)
  library(mixmeta)
  library(tsModel)
})

set.seed(20260904L)

project_root <- normalizePath(".", winslash = "/", mustWork = TRUE)

output_dir <- file.path(project_root, "Mechanism", "output")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

warm_months <- 6:9

disease_tag <- "gastrointestinal"

city_af_file <- file.path(
  project_root, "two_stage_analysis", "output", "major_classes",
  "two-stage-city-by-gcm.csv"
)

city_temperature_file <- file.path(
  project_root, "user_route", "output", "figure1_climate", "V01",
  "V01-Supp-city-scenario-temperature-by-gcm.csv"
)

index_file <- file.path(output_dir, "sanitation-pressure-capacity-index.csv")

fit_file <- file.path(
  project_root, "full_analysis", "output",
  "gastrointestinal_fit_allcities.rds"
)

standardise <- function(x) {
  usable <- is.finite(x)
  result <- rep(NA_real_, length(x))
  result[usable] <- (x[usable] - mean(x[usable])) / sd(x[usable])
  result
}

city_af <- fread(city_af_file, colClasses = list(character = "code"))

city_af <- city_af[tag == disease_tag, .(
  code, gcm, total_cases,
  delta_af = af_acc_direct,
  delta_cases = an_acc_direct
)]

city_temperature <- fread(
  city_temperature_file, colClasses = list(character = "code")
)

city_temperature[, delta_t := factual - histnat]

city_panel <- merge(
  city_af,
  city_temperature[, .(code, gcm, factual, histnat, delta_t)],
  by = c("code", "gcm"), all.x = TRUE
)

city_median <- city_panel[, .(
  total_cases = first(total_cases),
  delta_af = median(delta_af, na.rm = TRUE),
  delta_cases = median(delta_cases, na.rm = TRUE),
  delta_t = median(delta_t, na.rm = TRUE),
  n_gcm = uniqueN(gcm)
), by = code]

fit <- readRDS(fit_file)

argvar <- fit$argvar

lag_days <- as.integer(fit$config$lag)

blups <- fit$blups

daily <- as.data.table(fit$dat)[, .(code, date, tmean, n)]

daily[, month := as.integer(format(date, "%m"))]

setkey(daily, code)

forward_weights <- function(cases, lag) {
  as.numeric(rowMeans(as.matrix(Lag(cases, -seq.int(0L, as.integer(lag))))))
}

# A common national warm-season temperature distribution. Using a fixed grid
# with national day weights removes each city's own climate from the
# comparison, leaving only its response function.
national_warm <- daily[month %in% warm_months & is.finite(tmean), tmean]

reference_grid <- quantile(
  national_warm, probs = seq(0.005, 0.995, by = 0.01), na.rm = TRUE,
  names = FALSE
)

response_for_city <- function(city) {
  coefficients <- as.numeric(blups[[city]]$blup)
  city_daily <- daily[.(city)]

  weights <- forward_weights(city_daily$n, lag_days)
  basis_observed <- do.call(onebasis, c(list(x = city_daily$tmean), argvar))
  basis_warmed <- do.call(onebasis, c(list(x = city_daily$tmean + 1), argvar))
  difference <- basis_warmed - basis_observed
  valid <- is.finite(weights) & rowSums(!is.finite(difference)) == 0L
  own_delta_eta <- drop(difference[valid, , drop = FALSE] %*% coefficients)
  own_attributable <- weights[valid] * (1 - exp(-own_delta_eta))
  warm <- city_daily$month[valid] %in% warm_months

  standard_difference <-
    do.call(onebasis, c(list(x = reference_grid + 1), argvar)) -
    do.call(onebasis, c(list(x = reference_grid), argvar))
  standard_delta_eta <- drop(standard_difference %*% coefficients)

  data.table(
    code = city,
    response_own_climate = 100 * sum(own_attributable) / sum(weights[valid]),
    response_own_climate_warm = if (any(warm)) {
      100 * sum(own_attributable[warm]) / sum(weights[valid][warm])
    } else {
      NA_real_
    },
    response_standardised = 100 * mean(1 - exp(-standard_delta_eta))
  )
}

responses <- rbindlist(lapply(names(blups), response_for_city))

index <- fread(index_file, colClasses = list(character = c("code",
                                                           "province_code")))

city_meta <- as.data.table(fit$city_meta)[, .(code, mean_temp, temp_iqr)]

city_level <- Reduce(
  function(x, y) merge(x, y, by = "code", all.x = TRUE),
  list(city_median, responses, city_meta,
       index[, .(code, province_code, sanitation_stress,
                 sanitation_stress_city_only, pressure_index, capacity_index,
                 gdp_per_capita, log_gdp_per_capita, gdp_tertile,
                 rural_share_of_domestic,
                 cod_industry_per_industry_km2)])
)

city_level <- city_level[
  is.finite(delta_af) & is.finite(delta_t) & is.finite(sanitation_stress) &
    is.finite(log_gdp_per_capita) & is.finite(total_cases) & total_cases > 0
]

city_level[, `:=`(
  delta_t_z = standardise(delta_t),
  sanitation_z = standardise(sanitation_stress),
  sanitation_city_only_z = standardise(sanitation_stress_city_only),
  gdp_z = standardise(log_gdp_per_capita),
  mean_temp_z = standardise(mean_temp),
  temp_iqr_z = standardise(temp_iqr),
  industry_z = standardise(log1p(cod_industry_per_industry_km2))
)]

city_level[, stress_tertile := cut(
  frank(sanitation_stress, ties.method = "first"),
  breaks = quantile(frank(sanitation_stress, ties.method = "first"),
                    probs = c(0, 1 / 3, 2 / 3, 1)),
  labels = c("Lower stress", "Middle stress", "Higher stress"),
  include.lowest = TRUE
)]

# fread returns the tertile label as plain text, so the display order has to be
# restored explicitly.
city_level[, gdp_tertile := factor(
  gdp_tertile, levels = c("Lower GDP", "Middle GDP", "Higher GDP")
)]

interaction_formula <- delta_af ~ delta_t_z * sanitation_z + gdp_z +
  mean_temp_z + temp_iqr_z

interaction_model <- lm(interaction_formula, data = city_level,
                        weights = total_cases)

interaction_variants <- list(
  main = interaction_model,
  city_only_stress = lm(
    delta_af ~ delta_t_z * sanitation_city_only_z + gdp_z + mean_temp_z +
      temp_iqr_z,
    data = city_level, weights = total_cases
  ),
  industrial_control = lm(
    delta_af ~ delta_t_z * industry_z + gdp_z + mean_temp_z + temp_iqr_z,
    data = city_level, weights = total_cases
  ),
  unweighted = lm(interaction_formula, data = city_level),
  no_gdp = lm(
    delta_af ~ delta_t_z * sanitation_z + mean_temp_z + temp_iqr_z,
    data = city_level, weights = total_cases
  ),
  # The slope of the exposure-response function varies with local climate and
  # with GDP, so those interactions with warming must be in the model before
  # a sanitation interaction can be claimed.
  all_interactions = lm(
    delta_af ~ delta_t_z * (sanitation_z + gdp_z + mean_temp_z + temp_iqr_z),
    data = city_level, weights = total_cases
  ),
  # Sanitation against the industrial negative control, mutually adjusted.
  sanitation_vs_industry = lm(
    delta_af ~ delta_t_z * (sanitation_z + industry_z) + gdp_z + mean_temp_z +
      temp_iqr_z,
    data = city_level, weights = total_cases
  )
)

extract_terms <- function(model, label) {
  estimates <- summary(model)$coefficients
  confidence <- confint(model)
  terms <- setdiff(rownames(estimates), "(Intercept)")
  data.table(
    model = label,
    term = terms,
    estimate = estimates[terms, "Estimate"],
    std_error = estimates[terms, "Std. Error"],
    lower = confidence[terms, 1],
    upper = confidence[terms, 2],
    p_value = estimates[terms, "Pr(>|t|)"],
    r_squared = summary(model)$r.squared,
    n = nobs(model)
  )
}

interaction_table <- rbindlist(lapply(
  names(interaction_variants),
  function(label) extract_terms(interaction_variants[[label]], label)
))

# The same model inside each GCM, to show that the interaction does not hinge
# on the ensemble median.
by_gcm <- merge(
  city_panel,
  city_level[, .(code, sanitation_z, gdp_z, mean_temp_z, temp_iqr_z)],
  by = "code"
)

by_gcm[, delta_t_z := standardise(delta_t), by = gcm]

interaction_by_gcm <- by_gcm[, {
  model <- lm(delta_af ~ delta_t_z * sanitation_z + gdp_z + mean_temp_z +
                temp_iqr_z, data = .SD, weights = total_cases)
  estimates <- summary(model)$coefficients
  data.table(
    interaction_estimate = estimates["delta_t_z:sanitation_z", "Estimate"],
    interaction_p = estimates["delta_t_z:sanitation_z", "Pr(>|t|)"],
    n_cities = nobs(model)
  )
}, by = gcm]

response_outcomes <- c(
  own_climate_warm = "response_own_climate_warm",
  standardised = "response_standardised"
)

response_specifications <- c(
  crude = "~ gdp_z",
  plus_climate = "~ gdp_z + mean_temp_z + temp_iqr_z",
  plus_sanitation = "~ gdp_z + mean_temp_z + temp_iqr_z + sanitation_z"
)

response_table <- rbindlist(lapply(
  names(response_outcomes),
  function(outcome_name) {
    outcome_variable <- response_outcomes[[outcome_name]]
    outcome_data <- city_level[is.finite(get(outcome_variable))]
    rbindlist(lapply(names(response_specifications), function(spec_name) {
      model <- lm(
        as.formula(paste(outcome_variable, response_specifications[[spec_name]])),
        data = outcome_data, weights = total_cases
      )
      terms <- extract_terms(model, spec_name)
      terms[, outcome := outcome_name]
      terms[]
    }))
  }
))

sanitation_response_specifications <- c(
  "sanitation adjusted" = "~ sanitation_z + mean_temp_z + temp_iqr_z + gdp_z",
  "plus industrial control" = paste(
    "~ sanitation_z + industry_z + mean_temp_z + temp_iqr_z + gdp_z"
  ),
  "city-level capacity only" = paste(
    "~ sanitation_city_only_z + mean_temp_z + temp_iqr_z + gdp_z"
  ),
  "pressure and capacity separately" = paste(
    "~ pressure_index + capacity_index + mean_temp_z + temp_iqr_z + gdp_z"
  )
)

sanitation_response_models <- rbindlist(lapply(
  names(response_outcomes),
  function(outcome_name) {
    outcome_variable <- response_outcomes[[outcome_name]]
    outcome_data <- city_level[is.finite(get(outcome_variable))]
    rbindlist(lapply(
      names(sanitation_response_specifications),
      function(spec_name) {
        weighted_model <- lm(
          as.formula(paste(outcome_variable,
                           sanitation_response_specifications[[spec_name]])),
          data = outcome_data, weights = total_cases
        )
        unweighted_model <- lm(
          as.formula(paste(outcome_variable,
                           sanitation_response_specifications[[spec_name]])),
          data = outcome_data
        )
        rbind(
          extract_terms(weighted_model, spec_name)[, `:=`(
            outcome = outcome_name, weighting = "case-weighted"
          )],
          extract_terms(unweighted_model, spec_name)[, `:=`(
            outcome = outcome_name, weighting = "unweighted"
          )]
        )
      }
    ))
  }
))

fwrite(city_level, file.path(output_dir, "figure7-city-level-data.csv"))

# Rank correlation between sanitation-system stress and log GDP per capita
# across the cities in the burden models (reported in the Supplementary Text).
fwrite(
  data.table(
    n_cities = nrow(city_level),
    spearman_rho = cor(
      city_level$sanitation_stress, city_level$log_gdp_per_capita,
      method = "spearman"
    )
  ),
  file.path(output_dir, "sanitation-gdp-spearman.csv")
)

fwrite(interaction_table,
       file.path(output_dir, "figure7-interaction-models.csv"))

fwrite(interaction_by_gcm,
       file.path(output_dir, "figure7-interaction-by-gcm.csv"))

fwrite(response_table, file.path(output_dir, "figure7-response-models.csv"))

fwrite(sanitation_response_models,
       file.path(output_dir, "figure7-response-sanitation-terms.csv"))
