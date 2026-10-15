###############################################################################
## Core functions for city-specific DLNM -> mixmeta -> BLUP attribution
###############################################################################

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("Package 'data.table' is required.")
}
if (!requireNamespace("dlnm", quietly = TRUE)) {
  stop("Package 'dlnm' is required.")
}
if (!requireNamespace("mixmeta", quietly = TRUE)) {
  stop("Package 'mixmeta' is required.")
}

ts_basis <- function(x, argvar) {
  do.call(dlnm::onebasis, c(list(x = x), argvar))
}

ts_centered_basis <- function(x, reference, argvar) {
  basis <- ts_basis(x, argvar)
  centre <- as.numeric(ts_basis(reference, argvar))
  sweep(basis, 2L, centre, "-")
}

ts_find_minimum_risk_temperature <- function(
    temperature,
    coefficient,
    argvar
) {
  observed <- temperature[is.finite(temperature)]
  if (!length(observed)) {
    return(NA_real_)
  }
  limits <- as.numeric(stats::quantile(
    observed, c(0.01, 0.99), na.rm = TRUE, names = FALSE
  ))
  grid <- seq(limits[1L], limits[2L], length.out = 400L)
  linear_predictor <- drop(ts_basis(grid, argvar) %*% coefficient)
  grid[which.min(linear_predictor)]
}

ts_curve_is_plausible <- function(
    coefficient,
    temperature,
    reference,
    argvar,
    lower = 0.10,
    upper = 10
) {
  temperature <- temperature[is.finite(temperature)]
  if (!length(temperature) || any(!is.finite(coefficient))) {
    return(FALSE)
  }
  limits <- as.numeric(stats::quantile(
    temperature, c(0.01, 0.99),
    na.rm = TRUE, names = FALSE
  ))
  grid <- seq(limits[1L], limits[2L], length.out = 80L)
  log_rr <- drop(
    ts_centered_basis(grid, reference, argvar) %*% coefficient
  )
  rr <- exp(log_rr)
  all(is.finite(rr) & rr >= lower & rr <= upper)
}

ts_refit_robust_meta <- function(fit) {
  metadata <- data.table::as.data.table(fit$city_meta)
  metadata[, code := as.character(code)]
  fit_data <- data.table::as.data.table(fit$dat)
  fit_data[, code := as.character(code)]
  temperature_by_city <- split(fit_data$tmean, fit_data$code)
  minimum_city_cases <- if (!is.null(fit$minimum_city_cases)) {
    as.numeric(fit$minimum_city_cases)
  } else {
    TS_MIN_CITY_CASES
  }
  candidate_codes <- intersect(names(fit$raw_coef), metadata$code)
  candidate_codes <- intersect(
    candidate_codes,
    metadata[total_cases >= minimum_city_cases, code]
  )
  keep <- vapply(candidate_codes, function(code_i) {
    ts_curve_is_plausible(
      coefficient = fit$raw_coef[[code_i]],
      temperature = temperature_by_city[[code_i]],
      reference = fit$cen_common,
      argvar = fit$argvar
    )
  }, logical(1))
  retained_codes <- candidate_codes[keep]
  if (length(retained_codes) < 8L) {
    stop(
      "Fewer than eight city models passed the prespecified ",
      "RR-range quality control."
    )
  }

  meta_data <- metadata[match(retained_codes, code)]
  scaling <- list(
    mean_temp_mean = mean(meta_data$mean_temp, na.rm = TRUE),
    mean_temp_sd = stats::sd(meta_data$mean_temp, na.rm = TRUE),
    temp_iqr_mean = mean(meta_data$temp_iqr, na.rm = TRUE),
    temp_iqr_sd = stats::sd(meta_data$temp_iqr, na.rm = TRUE)
  )
  meta_data[, `:=`(
    mean_temp_z = (
      mean_temp - scaling$mean_temp_mean
    ) / scaling$mean_temp_sd,
    temp_iqr_z = (
      temp_iqr - scaling$temp_iqr_mean
    ) / scaling$temp_iqr_sd
  )]
  meta_data <- as.data.frame(meta_data)
  rownames(meta_data) <- retained_codes
  response <- do.call(rbind, fit$raw_coef[retained_codes])
  rownames(response) <- retained_codes
  covariance <- fit$raw_vcov[retained_codes]
  meta_model <- if (length(retained_codes) < 20L) {
    mixmeta::mixmeta(
      response ~ 1,
      covariance,
      data = meta_data,
      method = "reml",
      control = list(maxiter = 1000)
    )
  } else {
    tryCatch(
      mixmeta::mixmeta(
        response ~ mean_temp_z + temp_iqr_z,
        covariance,
        data = meta_data,
        method = "reml",
        control = list(
          showiter = FALSE,
          igls.inititer = 10,
          maxiter = 1000
        )
      ),
      error = function(e) {
        mixmeta::mixmeta(
          response ~ 1,
          covariance,
          data = meta_data,
          method = "reml",
          control = list(maxiter = 1000)
        )
      }
    )
  }
  blups <- mixmeta::blup(meta_model, vcov = TRUE)
  if (length(blups) == length(retained_codes)) {
    names(blups) <- retained_codes
  } else if (is.null(names(blups))) {
    stop("Second-stage BLUP output could not be mapped to city codes.")
  }
  list(
    model = meta_model,
    blups = blups,
    retained_codes = retained_codes,
    rejected_codes = setdiff(candidate_codes, retained_codes),
    scaling = scaling,
    temperature_by_city = temperature_by_city
  )
}

ts_meta_predicted_coefficients <- function(
    fit,
    codes,
    robust_meta
) {
  codes <- as.character(codes)
  response_dimension <- length(fit$raw_coef[[1L]])
  beta <- stats::coef(robust_meta$model)
  beta_matrix <- matrix(beta, nrow = response_dimension)
  term_count <- ncol(beta_matrix)

  metadata <- data.table::as.data.table(fit$city_meta)
  metadata[, code := as.character(code)]

  if (term_count > 1L) {
    metadata[, mean_temp_z := (
      mean_temp - robust_meta$scaling$mean_temp_mean
    ) / robust_meta$scaling$mean_temp_sd]
    metadata[, temp_iqr_z := (
      temp_iqr - robust_meta$scaling$temp_iqr_mean
    ) / robust_meta$scaling$temp_iqr_sd]
    # Avoid extrapolating second-stage covariates beyond the range represented
    # by retained cities.
    metadata[, `:=`(
      mean_temp_z = pmax(-2.5, pmin(2.5, mean_temp_z)),
      temp_iqr_z = pmax(-2.5, pmin(2.5, temp_iqr_z))
    )]
  }

  requested <- metadata[match(codes, code)]
  if (anyNA(requested$code)) {
    stop("One or more requested cities are absent from city metadata.")
  }
  design <- if (term_count == 1L) {
    matrix(1, nrow = nrow(requested), ncol = 1L)
  } else if (term_count == 3L) {
    cbind(1, requested$mean_temp_z, requested$temp_iqr_z)
  } else {
    stop("Unsupported second-stage meta-regression design.")
  }

  predicted <- design %*% t(beta_matrix)
  rownames(predicted) <- requested$code
  predicted
}

ts_city_coefficient_map <- function(fit) {
  selected_codes <- as.character(fit$selected_codes)
  robust_meta <- ts_refit_robust_meta(fit)
  predicted <- ts_meta_predicted_coefficients(
    fit, selected_codes, robust_meta
  )
  coefficients <- setNames(
    lapply(seq_along(selected_codes), function(i) {
      as.numeric(predicted[i, ])
    }),
    selected_codes
  )
  source <- setNames(
    rep("robust second-stage prediction", length(selected_codes)),
    selected_codes
  )
  response_dimension <- nrow(matrix(
    stats::coef(robust_meta$model),
    nrow = length(fit$raw_coef[[1L]])
  ))
  intercept_coefficient <- as.numeric(matrix(
    stats::coef(robust_meta$model),
    nrow = response_dimension
  )[, 1L])

  # A second-stage prediction can still extrapolate for a city with unusual
  # climate covariates. Replace such predictions with the meta-regression
  # intercept before considering city BLUPs.
  for (code_i in selected_codes) {
    if (!ts_curve_is_plausible(
      coefficient = coefficients[[code_i]],
      temperature = robust_meta$temperature_by_city[[code_i]],
      reference = fit$cen_common,
      argvar = fit$argvar
    )) {
      coefficients[[code_i]] <- intercept_coefficient
      source[[code_i]] <- "robust second-stage intercept"
    }
  }

  successful_codes <- intersect(
    names(robust_meta$blups), selected_codes
  )
  for (code_i in successful_codes) {
    candidate <- as.numeric(robust_meta$blups[[code_i]]$blup)
    if (ts_curve_is_plausible(
      coefficient = candidate,
      temperature = robust_meta$temperature_by_city[[code_i]],
      reference = fit$cen_common,
      argvar = fit$argvar
    )) {
      coefficients[[code_i]] <- candidate
      source[[code_i]] <- "quality-controlled city BLUP"
    }
  }
  qc <- data.table::data.table(
    code = selected_codes,
    retained_first_stage = selected_codes %in%
      robust_meta$retained_codes,
    coefficient_source = unname(source)
  )
  list(
    coefficients = coefficients,
    source = source,
    robust_meta = robust_meta,
    qc = qc
  )
}

ts_prepare_major_fit <- function(fit, tag) {
  city_map <- ts_city_coefficient_map(fit)
  dat <- data.table::as.data.table(fit$dat)
  dat[, `:=`(
    code = as.character(code),
    date = as.Date(date),
    n = as.numeric(n)
  )]
  reference <- if (tag == "respiratory") {
    NA_real_
  } else {
    as.numeric(stats::quantile(
      dat$tmean,
      TS_NONRESP_REFERENCE_PROBABILITY,
      na.rm = TRUE
    ))
  }

  city_reference <- setNames(
    rep(reference, length(city_map$coefficients)),
    names(city_map$coefficients)
  )
  if (tag == "respiratory") {
    for (code_i in names(city_map$coefficients)) {
      city_reference[[code_i]] <- ts_find_minimum_risk_temperature(
        dat[code == code_i, tmean],
        city_map$coefficients[[code_i]],
        fit$argvar
      )
    }
  }

  list(
    tag = tag,
    dat = dat,
    lag = as.integer(fit$config$lag),
    argvar = fit$argvar,
    coefficient = city_map$coefficients,
    coefficient_source = city_map$source,
    reference = city_reference,
    meta_model = city_map$robust_meta$model,
    meta_scaling = city_map$robust_meta$scaling,
    coefficient_qc = city_map$qc,
    raw_coef = fit$raw_coef,
    raw_vcov = fit$raw_vcov,
    city_meta = fit$city_meta,
    fit_status = fit$fit_status
  )
}

ts_city_daily_attribution <- function(
    factual_temperature,
    histnat_temperature,
    cases,
    dates,
    coefficient,
    argvar,
    lag,
    reference,
    tag
) {
  forward_cases <- acc_forward_weights(cases, lag)
  factual_basis <- ts_basis(factual_temperature, argvar)
  histnat_basis <- ts_basis(histnat_temperature, argvar)
  basis_difference <- factual_basis - histnat_basis

  valid <- is.finite(forward_cases) &
    rowSums(!is.finite(basis_difference)) == 0L
  total_cases <- sum(cases, na.rm = TRUE)
  if (!any(valid) || total_cases <= 0) {
    return(data.table::data.table())
  }
  forward_total <- sum(forward_cases[valid], na.rm = TRUE)
  weights <- forward_cases[valid] * total_cases / forward_total

  factual_temperature <- factual_temperature[valid]
  histnat_temperature <- histnat_temperature[valid]
  dates <- as.Date(dates[valid])
  factual_basis <- factual_basis[valid, , drop = FALSE]
  histnat_basis <- histnat_basis[valid, , drop = FALSE]
  basis_difference <- basis_difference[valid, , drop = FALSE]

  direct_eta <- drop(basis_difference %*% coefficient)
  direct_fraction <- 1 - exp(-direct_eta)

  centre_basis <- as.numeric(ts_basis(reference, argvar))
  factual_eta <- drop(
    sweep(factual_basis, 2L, centre_basis, "-") %*% coefficient
  )
  histnat_eta <- drop(
    sweep(histnat_basis, 2L, centre_basis, "-") %*% coefficient
  )
  factual_attributable_fraction <- 1 - exp(-factual_eta)
  histnat_attributable_fraction <- 1 - exp(-histnat_eta)

  if (tag == "respiratory") {
    relevant_factual <- factual_attributable_fraction
    relevant_histnat <- histnat_attributable_fraction
  } else {
    factual_high_basis <- ts_basis(
      pmax(factual_temperature, reference), argvar
    )
    histnat_high_basis <- ts_basis(
      pmax(histnat_temperature, reference), argvar
    )
    relevant_factual <- 1 - exp(-drop(
      sweep(factual_high_basis, 2L, centre_basis, "-") %*%
        coefficient
    ))
    relevant_histnat <- 1 - exp(-drop(
      sweep(histnat_high_basis, 2L, centre_basis, "-") %*%
        coefficient
    ))
  }

  data.table::data.table(
    date = dates,
    factual_temperature = factual_temperature,
    histnat_temperature = histnat_temperature,
    paired_midpoint = (
      factual_temperature + histnat_temperature
    ) / 2,
    case_weight = weights,
    direct_cases = direct_fraction * weights,
    factual_temperature_cases = (
      factual_attributable_fraction * weights
    ),
    histnat_temperature_cases = (
      histnat_attributable_fraction * weights
    ),
    relevant_factual_cases = relevant_factual * weights,
    relevant_histnat_cases = relevant_histnat * weights
  )
}

ts_summarise_city_daily <- function(daily, total_cases) {
  if (!nrow(daily) || total_cases <= 0) {
    return(data.table::data.table(
      total_cases = total_cases,
      an_acc_direct = 0,
      af_acc_direct = 0,
      factual_temperature_cases = 0,
      histnat_temperature_cases = 0,
      factual_temperature_af = 0,
      histnat_temperature_af = 0,
      relevant_factual_cases = 0,
      relevant_histnat_cases = 0,
      disease_relevant_ac = NA_real_
    ))
  }
  direct_cases <- sum(daily$direct_cases, na.rm = TRUE)
  factual_cases <- sum(
    daily$factual_temperature_cases, na.rm = TRUE
  )
  histnat_cases <- sum(
    daily$histnat_temperature_cases, na.rm = TRUE
  )
  relevant_factual <- sum(
    daily$relevant_factual_cases, na.rm = TRUE
  )
  relevant_histnat <- sum(
    daily$relevant_histnat_cases, na.rm = TRUE
  )
  data.table::data.table(
    total_cases = total_cases,
    an_acc_direct = direct_cases,
    af_acc_direct = 100 * direct_cases / total_cases,
    factual_temperature_cases = factual_cases,
    histnat_temperature_cases = histnat_cases,
    factual_temperature_af = 100 * factual_cases / total_cases,
    histnat_temperature_af = 100 * histnat_cases / total_cases,
    relevant_factual_cases = relevant_factual,
    relevant_histnat_cases = relevant_histnat,
    disease_relevant_ac = if (
      is.finite(relevant_factual) &&
        abs(relevant_factual) > .Machine$double.eps
    ) {
      100 * (relevant_factual - relevant_histnat) /
        relevant_factual
    } else {
      NA_real_
    }
  )
}
