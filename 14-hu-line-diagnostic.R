# -*- coding: UTF-8 -*-
# Compare the Figure 7 sanitation sample on each side of the Hu line.
# Pass an ASCII-safe project root on Windows if R cannot decode the Chinese path.

suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[[1]] else "."
out <- if (length(args) > 1L) args[[2]] else file.path(root, "Mechanism", "output")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

city <- fread(
  file.path(root, "Mechanism", "output", "figure7-city-level-data.csv"),
  colClasses = list(character = c("code", "province_code"))
)
side <- fread(
  file.path(root, "SII", "output", "figure6-objects", "figure6-dividing-line-city-side.csv"),
  colClasses = list(character = "code_merge")
)
setnames(side, "code_merge", "code")
side <- unique(side[, .(code, hu)])
joined <- merge(city, side, by = "code")
stopifnot(nrow(joined) > 0L)

figure6_response <- fread(
  file.path(root, "SII", "output", "figure6-objects",
            "figure6-plus1c-response-by-city-gcm.csv"),
  colClasses = list(character = c("code_merge", "gcm"))
)[, .(
  figure6_response_nat = mean(response_1c_nat, na.rm = TRUE),
  figure6_response_fact = mean(response_1c, na.rm = TRUE)
), by = code_merge]
setnames(figure6_response, "code_merge", "code")
joined <- merge(joined, figure6_response, by = "code", all.x = TRUE)

safe_weighted_mean <- function(value, weight) {
  ok <- is.finite(value) & is.finite(weight) & weight > 0
  if (!any(ok)) return(NA_real_)
  weighted.mean(value[ok], weight[ok])
}

by_side <- joined[, .(
  cities = .N,
  cases = sum(total_cases, na.rm = TRUE),
  median_gdp_per_capita = median(gdp_per_capita, na.rm = TRUE),
  median_delta_af = median(delta_af, na.rm = TRUE),
  case_weighted_delta_af = safe_weighted_mean(delta_af, total_cases),
  median_response_plus1c = median(response_standardised, na.rm = TRUE),
  case_weighted_response_plus1c = safe_weighted_mean(response_standardised, total_cases),
  median_figure6_response_nat = median(figure6_response_nat, na.rm = TRUE),
  case_weighted_figure6_response_nat = safe_weighted_mean(figure6_response_nat, total_cases),
  median_figure6_response_fact = median(figure6_response_fact, na.rm = TRUE),
  median_sanitation_z = median(sanitation_z, na.rm = TRUE),
  case_weighted_sanitation_z = safe_weighted_mean(sanitation_z, total_cases),
  median_pressure = median(pressure_index, na.rm = TRUE),
  median_capacity = median(capacity_index, na.rm = TRUE)
), by = hu]

by_side[, `:=`(
  figure7_cities = nrow(city),
  hu_classified_cities = nrow(side),
  joined_cities = nrow(joined)
)]

correlations <- rbindlist(lapply(
  c("All matched cities", "Northwest", "Southeast"),
  function(g) {
    d <- if (g == "All matched cities") joined else joined[hu == g]
    data.table(
      group = g,
      n = nrow(d),
      rho_gdp_stress = cor(d$log_gdp_per_capita, d$sanitation_z,
                           method = "spearman", use = "complete.obs"),
      rho_stress_response = cor(d$sanitation_z, d$response_standardised,
                                method = "spearman", use = "complete.obs"),
      rho_stress_figure6_nat = cor(d$sanitation_z, d$figure6_response_nat,
                                   method = "spearman", use = "complete.obs")
    )
  }
))
cross_estimator_correlation <- joined[, .(
  n = sum(is.finite(response_standardised) & is.finite(figure6_response_nat)),
  pearson = cor(response_standardised, figure6_response_nat,
                use = "complete.obs"),
  spearman = cor(response_standardised, figure6_response_nat,
                 method = "spearman", use = "complete.obs")
)]

model_data <- joined[
  is.finite(response_standardised) & is.finite(sanitation_z) &
    is.finite(gdp_z) & is.finite(mean_temp_z) & is.finite(temp_iqr_z) &
    is.finite(total_cases) & total_cases > 0
]
base_model <- lm(
  response_standardised ~ sanitation_z + gdp_z + mean_temp_z + temp_iqr_z,
  data = model_data, weights = total_cases
)
side_model <- lm(
  response_standardised ~ sanitation_z + gdp_z + mean_temp_z + temp_iqr_z + hu,
  data = model_data, weights = total_cases
)
interaction_model <- lm(
  response_standardised ~ sanitation_z * hu + gdp_z + mean_temp_z + temp_iqr_z,
  data = model_data, weights = total_cases
)

comparable_models <- rbindlist(lapply(c(
  "response_standardised", "response_own_climate",
  "figure6_response_nat", "figure6_response_fact"
), function(outcome) {
  d <- model_data[is.finite(get(outcome))]
  m <- lm(
    as.formula(paste(outcome,
      "~ sanitation_z + gdp_z + mean_temp_z + temp_iqr_z")),
    data = d, weights = total_cases
  )
  ct <- summary(m)$coefficients["sanitation_z", ]
  ci <- confint(m, "sanitation_z")
  data.table(
    outcome = outcome,
    n = nobs(m),
    estimate = unname(ct[["Estimate"]]),
    lower = unname(ci[[1]]),
    upper = unname(ci[[2]]),
    p_value = unname(ct[["Pr(>|t|)"]])
  )
}))

model_rows <- rbindlist(lapply(list(
  "Matched cities; original adjustment" = base_model,
  "Add Hu-line side" = side_model,
  "Allow side-specific slope" = interaction_model
), function(m) {
  ct <- summary(m)$coefficients
  ci <- confint(m)
  terms <- intersect(rownames(ct), c("sanitation_z", "sanitation_z:huSoutheast", "huSoutheast"))
  data.table(
    model = deparse(formula(m), width.cutoff = 500),
    term = terms,
    estimate = ct[terms, "Estimate"],
    lower = ci[terms, 1],
    upper = ci[terms, 2],
    p_value = ct[terms, "Pr(>|t|)"],
    n = nobs(m)
  )
}), use.names = TRUE)

fwrite(by_side, file.path(out, "figure7-hu-line-by-side.csv"))
fwrite(correlations, file.path(out, "figure7-hu-line-correlations.csv"))
fwrite(cross_estimator_correlation, file.path(out, "figure7-cross-estimator-correlation.csv"))
fwrite(model_rows, file.path(out, "figure7-hu-line-models.csv"))
fwrite(comparable_models, file.path(out, "figure7-comparable-response-models.csv"))
print(by_side)
print(correlations)
print(cross_estimator_correlation)
print(model_rows[, .(term, estimate, lower, upper, p_value, n)])
print(comparable_models)
