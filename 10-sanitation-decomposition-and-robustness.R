# -*- coding: UTF-8 -*-
###############################################################################
## Figure 7, panels d and e, plus the supplementary robustness set
##
##   d  pressure versus capacity: which side of the sanitation balance tracks
##      the response to equivalent warming
##   e  hydrological context: does sanitation stress matter more where the
##      system is also hydraulically loaded
##
##   S1 non-linearity of the sanitation-response relationship (natural spline)
##   S2 seasonal specificity: warm-season versus cold-season reference
##      distribution
##   S3 spatial robustness: region fixed effects and leave-one-region-out
##
## Inputs are all already computed:
##   Mechanism/output/figure7-city-level-data.csv        (script 05)
##   Mechanism/output/city-hydrological-context.csv      (script 09)
##   two_stage_analysis/output/major_classes/two-stage-city-by-gcm.csv (region)
##   full_analysis/output/gastrointestinal_fit_allcities.rds (BLUPs, for S2)
###############################################################################

options(stringsAsFactors = FALSE, scipen = 999)

suppressPackageStartupMessages({
  library(data.table)
  library(dlnm)
  library(splines)
})

set.seed(20260904L)

project_root <- normalizePath(".", winslash = "/", mustWork = TRUE)

output_dir <- file.path(project_root, "Mechanism", "output")


city_file <- file.path(output_dir, "figure7-city-level-data.csv")

hydro_file <- file.path(output_dir, "city-hydrological-context.csv")

region_file <- file.path(
  project_root, "two_stage_analysis", "output", "major_classes",
  "two-stage-city-by-gcm.csv"
)

fit_file <- file.path(
  project_root, "full_analysis", "output",
  "gastrointestinal_fit_allcities.rds"
)

warm_months <- 6:9

cold_months <- c(11L, 12L, 1L, 2L)

standardise <- function(x) {
  usable <- is.finite(x)
  result <- rep(NA_real_, length(x))
  result[usable] <- (x[usable] - mean(x[usable])) / sd(x[usable])
  result
}

extract_terms <- function(model, label, keep = NULL) {
  estimates <- summary(model)$coefficients
  confidence <- confint(model)
  terms <- setdiff(rownames(estimates), "(Intercept)")
  if (!is.null(keep)) {
    terms <- intersect(terms, keep)
  }
  data.table(
    model = label,
    term = terms,
    estimate = estimates[terms, "Estimate"],
    std_error = estimates[terms, "Std. Error"],
    lower = confidence[terms, 1],
    upper = confidence[terms, 2],
    p_value = estimates[terms, "Pr(>|t|)"],
    n = nobs(model)
  )
}

partial_r_squared <- function(full_model, term, data, weights_column) {
  local_data <- as.data.frame(data)
  local_data$analysis_weight <- local_data[[weights_column]]
  reduced_formula <- update(formula(full_model), paste(". ~ . -", term))
  reduced_model <- lm(reduced_formula, data = local_data,
                      weights = analysis_weight)
  sse_full <- sum(weighted.residuals(full_model)^2)
  sse_reduced <- sum(weighted.residuals(reduced_model)^2)
  (sse_reduced - sse_full) / sse_reduced
}

city_level <- fread(city_file,
                    colClasses = list(character = c("code", "province_code")))

hydro <- fread(hydro_file, colClasses = list(character = "code"))

region_lookup <- unique(
  fread(region_file, colClasses = list(character = "code"))[
    tag == "gastrointestinal", .(code, region)
  ]
)

city_level <- merge(
  city_level,
  hydro[, .(code, warm_season_rainfall, annual_rainfall, heavy_rain_days,
            wet_day_fraction, warm_season_humidity, hydro_stress)],
  by = "code", all.x = TRUE
)

city_level <- merge(city_level, region_lookup, by = "code", all.x = TRUE)

analysis <- city_level[
  is.finite(response_standardised) & is.finite(sanitation_z) &
    is.finite(total_cases) & total_cases > 0
]

analysis[, `:=`(
  pressure_z = standardise(pressure_index),
  capacity_z = standardise(capacity_index),
  hydro_z = standardise(hydro_stress)
)]

covariates <- "gdp_z + mean_temp_z + temp_iqr_z"

decomposition_models <- list(
  "Domestic pollution pressure" = paste(
    "response_standardised ~ pressure_z +", covariates
  ),
  "Sanitation capacity" = paste(
    "response_standardised ~ capacity_z +", covariates
  ),
  "Pressure and capacity together" = paste(
    "response_standardised ~ pressure_z + capacity_z +", covariates
  ),
  "Composite stress" = paste(
    "response_standardised ~ sanitation_z +", covariates
  )
)

decomposition <- rbindlist(lapply(
  names(decomposition_models),
  function(label) {
    fitted_model <- lm(as.formula(decomposition_models[[label]]),
                       data = analysis, weights = total_cases)
    keep <- intersect(c("pressure_z", "capacity_z", "sanitation_z"),
                      names(coef(fitted_model)))
    terms <- extract_terms(fitted_model, label, keep = keep)
    partial_values <- vapply(
      terms$term,
      function(current) partial_r_squared(fitted_model, current, analysis,
                                          "total_cases"),
      numeric(1)
    )
    terms[, partial_r2 := partial_values]
    terms[, model_r2 := summary(fitted_model)$r.squared]
    terms[]
  }
))

decomposition[, term_label := fcase(
  term == "pressure_z", "Domestic pollution pressure",
  term == "capacity_z", "Sanitation capacity",
  term == "sanitation_z", "Composite stress",
  default = term
)]

fwrite(decomposition, file.path(output_dir,
                                "figure7-panelD-decomposition.csv"))

hydro_analysis <- analysis[is.finite(hydro_z)]

hydro_analysis[, `:=`(
  sanitation_group = fifelse(
    sanitation_z > median(sanitation_z), "Higher sanitation stress",
    "Lower sanitation stress"
  ),
  hydro_group = fifelse(
    hydro_z > median(hydro_z), "Wetter", "Drier"
  )
)]

hydro_analysis[, combined_group := factor(
  paste(sanitation_group, hydro_group, sep = " / "),
  levels = c(
    "Lower sanitation stress / Drier", "Lower sanitation stress / Wetter",
    "Higher sanitation stress / Drier", "Higher sanitation stress / Wetter"
  )
)]

# One model answers both questions, so that the plotted contrast and the
# reported test cannot disagree: the four cells are marginal predictions from
# the continuous interaction model, evaluated at each cell's own centroid.
interaction_model <- lm(
  response_standardised ~ sanitation_z * hydro_z + gdp_z + mean_temp_z +
    temp_iqr_z,
  data = hydro_analysis, weights = total_cases
)

interaction_terms <- extract_terms(
  interaction_model, "continuous interaction (plotted model)",
  keep = c("sanitation_z", "hydro_z", "sanitation_z:hydro_z")
)

cell_centroids <- hydro_analysis[, .(
  sanitation_z = mean(sanitation_z),
  hydro_z = mean(hydro_z),
  n_cities = .N,
  total_cases = sum(total_cases)
), by = combined_group]

cell_centroids[, `:=`(gdp_z = 0, mean_temp_z = 0, temp_iqr_z = 0)]

group_predictions <- as.data.table(predict(
  interaction_model, newdata = cell_centroids, interval = "confidence"
))

group_summary <- cbind(
  cell_centroids[, .(combined_group, n_cities, total_cases,
                     sanitation_centroid = sanitation_z,
                     hydro_centroid = hydro_z)],
  group_predictions
)

setnames(group_summary, c("fit", "lwr", "upr"),
         c("adjusted_response", "lower", "upper"))

setorder(group_summary, combined_group)

# The dichotomised model is kept only as a documented sensitivity. It gives a
# much larger difference-in-differences than the continuous test, so the
# region jackknife below shows where that comes from.
categorical_model <- lm(
  response_standardised ~ sanitation_group * hydro_group + gdp_z +
    mean_temp_z + temp_iqr_z,
  data = hydro_analysis, weights = total_cases
)

categorical_term <- grep(":", names(coef(categorical_model)), value = TRUE)[1L]

categorical_terms <- extract_terms(
  categorical_model, "dichotomised difference-in-differences",
  keep = categorical_term
)

categorical_jackknife <- rbindlist(lapply(
  sort(unique(hydro_analysis[!is.na(region) & nzchar(region)]$region)),
  function(dropped) {
    subset_data <- hydro_analysis[region != dropped]
    model <- lm(
      response_standardised ~ sanitation_group * hydro_group + gdp_z +
        mean_temp_z + temp_iqr_z,
      data = subset_data, weights = total_cases
    )
    current_term <- grep(":", names(coef(model)), value = TRUE)[1L]
    terms <- extract_terms(model, paste("excluding", dropped),
                           keep = current_term)
    terms[, dropped_region := dropped]
    terms[]
  }
))

stratified_terms <- rbindlist(lapply(c("Drier", "Wetter"), function(group) {
  subset_data <- hydro_analysis[hydro_group == group]
  model <- lm(
    response_standardised ~ sanitation_z + gdp_z + mean_temp_z + temp_iqr_z,
    data = subset_data, weights = total_cases
  )
  extract_terms(model, paste("sanitation slope,", group),
                keep = "sanitation_z")
}))

fwrite(group_summary, file.path(output_dir, "figure7-panelE-hydro-groups.csv"))

fwrite(rbindlist(list(interaction_terms, categorical_terms, stratified_terms,
                      categorical_jackknife), fill = TRUE),
       file.path(output_dir, "figure7-panelE-hydro-models.csv"))

linear_model <- lm(
  response_standardised ~ sanitation_z + gdp_z + mean_temp_z + temp_iqr_z,
  data = analysis, weights = total_cases
)

spline_model <- lm(
  response_standardised ~ ns(sanitation_z, df = 3) + gdp_z + mean_temp_z +
    temp_iqr_z,
  data = analysis, weights = total_cases
)

non_linearity_test <- anova(linear_model, spline_model)

spline_grid <- data.table(
  sanitation_z = seq(quantile(analysis$sanitation_z, 0.02),
                     quantile(analysis$sanitation_z, 0.98), length.out = 120),
  gdp_z = 0, mean_temp_z = 0, temp_iqr_z = 0
)

spline_predictions <- as.data.table(predict(
  spline_model, newdata = spline_grid, interval = "confidence"
))

spline_curve <- cbind(spline_grid[, .(sanitation_z)], spline_predictions)

setnames(spline_curve, c("fit", "lwr", "upr"),
         c("response", "lower", "upper"))

non_linearity <- data.table(
  comparison = "linear versus natural spline with 3 df",
  f_statistic = non_linearity_test$F[2],
  df = non_linearity_test$Df[2],
  p_value = non_linearity_test$`Pr(>F)`[2],
  linear_r2 = summary(linear_model)$r.squared,
  spline_r2 = summary(spline_model)$r.squared
)

fwrite(spline_curve, file.path(output_dir, "figureS-sanitation-spline.csv"))

fwrite(non_linearity,
       file.path(output_dir, "figureS-sanitation-nonlinearity-test.csv"))

fit <- readRDS(fit_file)

argvar <- fit$argvar

blups <- fit$blups

daily <- as.data.table(fit$dat)[, .(code, date, tmean)]

daily[, month := as.integer(format(date, "%m"))]

season_grid <- function(months) {
  as.numeric(quantile(
    daily[month %in% months & is.finite(tmean), tmean],
    probs = seq(0.005, 0.995, by = 0.01), na.rm = TRUE, names = FALSE
  ))
}

warm_grid <- season_grid(warm_months)

cold_grid <- season_grid(cold_months)

season_response <- function(grid, coefficients) {
  difference <- do.call(onebasis, c(list(x = grid + 1), argvar)) -
    do.call(onebasis, c(list(x = grid), argvar))
  100 * mean(1 - exp(-drop(difference %*% coefficients)))
}

season_responses <- rbindlist(lapply(names(blups), function(city) {
  coefficients <- as.numeric(blups[[city]]$blup)
  data.table(
    code = city,
    response_warm_grid = season_response(warm_grid, coefficients),
    response_cold_grid = season_response(cold_grid, coefficients)
  )
}))

season_data <- merge(analysis, season_responses, by = "code", all.x = TRUE)

seasonal <- rbindlist(lapply(
  c(warm = "response_warm_grid", cold = "response_cold_grid"),
  function(outcome) {
    model <- lm(
      as.formula(paste(outcome,
                       "~ sanitation_z + gdp_z + mean_temp_z + temp_iqr_z")),
      data = season_data[is.finite(get(outcome))], weights = total_cases
    )
    terms <- extract_terms(model, outcome, keep = "sanitation_z")
    terms[, mean_response := weighted.mean(
      season_data[[outcome]], season_data$total_cases, na.rm = TRUE
    )]
    terms[]
  }
), idcol = "season")

fwrite(seasonal, file.path(output_dir, "figureS-sanitation-seasonal.csv"))

region_data <- analysis[!is.na(region) & nzchar(region)]

region_models <- rbindlist(list(
  extract_terms(
    lm(response_standardised ~ sanitation_z + gdp_z + mean_temp_z +
         temp_iqr_z, data = region_data, weights = total_cases),
    "no spatial adjustment", keep = "sanitation_z"
  ),
  extract_terms(
    lm(response_standardised ~ sanitation_z + gdp_z + mean_temp_z +
         temp_iqr_z + factor(region), data = region_data,
       weights = total_cases),
    "region fixed effects", keep = "sanitation_z"
  ),
  extract_terms(
    lm(response_standardised ~ sanitation_z + gdp_z + mean_temp_z +
         temp_iqr_z + factor(province_code), data = region_data,
       weights = total_cases),
    "province fixed effects", keep = "sanitation_z"
  )
))

leave_one_out <- rbindlist(lapply(
  sort(unique(region_data$region)),
  function(dropped) {
    subset_data <- region_data[region != dropped]
    model <- lm(
      response_standardised ~ sanitation_z + gdp_z + mean_temp_z + temp_iqr_z,
      data = subset_data, weights = total_cases
    )
    terms <- extract_terms(model, paste("excluding", dropped),
                           keep = "sanitation_z")
    terms[, dropped_region := dropped]
    terms[]
  }
))

spatial <- rbindlist(list(region_models, leave_one_out), fill = TRUE)

fwrite(spatial, file.path(output_dir, "figureS-sanitation-spatial.csv"))
