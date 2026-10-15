###############################################################################
## Shared attribution functions
##
## Final estimand
##   counterfactual risk ratio = RR(hist-nat) / RR(factual)
##   direct ACC cases          = sum[y * (1 - counterfactual risk ratio)]
##   direct ACC AF             = direct ACC cases / sum(y)
##
## The health-response reference temperature cancels from this paired
## factual-versus-hist-nat contrast. All final figures and tables must use
## these functions instead of subtracting two reference-based burdens.
###############################################################################

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("Package 'data.table' is required.")
}
if (!requireNamespace("dlnm", quietly = TRUE)) {
  stop("Package 'dlnm' is required.")
}
if (!requireNamespace("tsModel", quietly = TRUE)) {
  stop("Package 'tsModel' is required.")
}

ACC_FORMULA_VERSION <- "direct-paired-v2"
ACC_BOOTSTRAP_REPETITIONS <- 10000L
ACC_BOOTSTRAP_SEED <- 20260724L

acc_file_md5 <- function(path) {
  if (!file.exists(path)) {
    return(NA_character_)
  }
  unname(tools::md5sum(path))
}

acc_cache_metadata <- function(
    fit_file,
    climate_file,
    baseline = "observed_daily"
) {
  list(
    formula_version = ACC_FORMULA_VERSION,
    fit_md5 = acc_file_md5(fit_file),
    climate_md5 = acc_file_md5(climate_file),
    baseline = as.character(baseline)
  )
}

acc_cache_is_current <- function(path, expected_metadata) {
  if (!file.exists(path)) {
    return(FALSE)
  }
  cached <- tryCatch(readRDS(path), error = function(e) NULL)
  if (is.null(cached) || is.null(cached$metadata)) {
    return(FALSE)
  }
  identical(cached$metadata, expected_metadata)
}

acc_complete_climate_pair <- function(pair, codes, dates) {
  pair <- data.table::as.data.table(pair)[, .(
    code = as.character(code),
    date = as.Date(date),
    tf = as.numeric(tf),
    tn = as.numeric(tn)
  )]
  pair <- pair[code %in% as.character(codes)]
  base <- data.table::CJ(
    code = as.character(codes),
    date = as.Date(dates),
    unique = TRUE
  )
  base[, province_code := substr(code, 1L, 2L)]
  pair[, province_code := substr(code, 1L, 2L)]

  province_daily <- pair[, .(
    tf_province = mean(tf, na.rm = TRUE),
    tn_province = mean(tn, na.rm = TRUE)
  ), by = .(province_code, date)]
  national_daily <- pair[, .(
    tf_national = mean(tf, na.rm = TRUE),
    tn_national = mean(tn, na.rm = TRUE)
  ), by = date]

  out <- merge(
    base,
    pair[, .(code, date, tf, tn)],
    by = c("code", "date"),
    all.x = TRUE
  )
  out <- merge(
    out,
    province_daily,
    by = c("province_code", "date"),
    all.x = TRUE
  )
  out <- merge(out, national_daily, by = "date", all.x = TRUE)
  out[, `:=`(
    tf = data.table::fcoalesce(tf, tf_province, tf_national),
    tn = data.table::fcoalesce(tn, tn_province, tn_national)
  )]
  out[, c(
    "tf_province", "tn_province",
    "tf_national", "tn_national"
  ) := NULL]
  data.table::setorder(out, code, date)

  if (any(!is.finite(out$tf)) || any(!is.finite(out$tn))) {
    stop("Temperature completion failed for one or more city-days.")
  }
  out
}

acc_prepare_case_weights <- function(
    cases,
    dates,
    baseline = c("observed_daily", "day_of_year")
) {
  baseline <- match.arg(baseline)
  cases <- data.table::as.data.table(cases)[, .(
    code = as.character(code),
    date = as.Date(date),
    n = as.numeric(n)
  )]
  grid <- data.table::CJ(
    code = sort(unique(cases$code)),
    date = as.Date(dates),
    unique = TRUE
  )

  if (baseline == "observed_daily") {
    out <- merge(grid, cases, by = c("code", "date"), all.x = TRUE)
    out[is.na(n), n := 0]
  } else {
    cases[, mmdd := format(date, "%m-%d")]
    climatology <- cases[, .(
      n = mean(n, na.rm = TRUE)
    ), by = .(code, mmdd)]
    grid[, mmdd := format(date, "%m-%d")]
    out <- merge(grid, climatology, by = c("code", "mmdd"), all.x = TRUE)
    out[is.na(n), n := 0]
    out[, mmdd := NULL]

    # Keep the total baseline equal to the observed total within every city.
    totals <- cases[, .(observed_total = sum(n, na.rm = TRUE)), by = code]
    out <- merge(out, totals, by = "code", all.x = TRUE)
    out[, climatology_total := sum(n, na.rm = TRUE), by = code]
    out[
      is.finite(climatology_total) & climatology_total > 0,
      n := n * observed_total / climatology_total
    ]
    out[, c("observed_total", "climatology_total") := NULL]
  }

  data.table::setorder(out, code, date)
  out
}

acc_forward_weights <- function(cases, lag) {
  lag_days <- seq.int(0L, as.integer(lag))
  as.numeric(rowMeans(
    as.matrix(tsModel::Lag(cases, -lag_days))
  ))
}

acc_one_basis <- function(temperature, argvar) {
  do.call(
    dlnm::onebasis,
    c(list(x = temperature), argvar)
  )
}

acc_direct_daily <- function(
    factual_temperature,
    histnat_temperature,
    cases,
    fit,
    dates = seq_along(cases)
) {
  forward_cases <- acc_forward_weights(cases, fit$lag)
  x_factual <- acc_one_basis(factual_temperature, fit$argvar)
  x_histnat <- acc_one_basis(histnat_temperature, fit$argvar)
  x_delta <- x_factual - x_histnat

  valid <- is.finite(forward_cases) &
    rowSums(!is.finite(x_delta)) == 0L
  total_cases <- sum(cases, na.rm = TRUE)
  if (!any(valid) || total_cases <= 0) {
    return(data.table::data.table())
  }

  # Scale the forward weights once and use exactly the same weights for point
  # estimates and all reported derived quantities.
  forward_total <- sum(forward_cases[valid], na.rm = TRUE)
  if (!is.finite(forward_total) || forward_total <= 0) {
    stop("The forward case-weight total is invalid.")
  }
  weights <- forward_cases[valid] * total_cases / forward_total
  delta_eta <- drop(
    x_delta[valid, , drop = FALSE] %*% fit$reduced_coef
  )
  daily_direct_af <- 1 - exp(-delta_eta)
  data.table::data.table(
    date = dates[valid],
    case_weight = weights,
    counterfactual_risk_ratio = exp(-delta_eta),
    daily_direct_af = daily_direct_af,
    direct_cases = daily_direct_af * weights
  )
}

acc_direct_city <- function(
    factual_temperature,
    histnat_temperature,
    cases,
    fit
) {
  total_cases <- sum(cases, na.rm = TRUE)
  daily <- acc_direct_daily(
    factual_temperature = factual_temperature,
    histnat_temperature = histnat_temperature,
    cases = cases,
    fit = fit
  )
  if (!nrow(daily) || total_cases <= 0) {
    return(data.table::data.table(
      total_cases = total_cases,
      included_forward_cases = 0,
      an_acc_direct = 0,
      counterfactual_cases = total_cases,
      af_acc_direct = 0
    ))
  }
  an_acc_direct <- sum(daily$direct_cases)
  counterfactual_cases <- total_cases - an_acc_direct
  af_acc_direct <- 100 * an_acc_direct / total_cases

  data.table::data.table(
    total_cases = total_cases,
    included_forward_cases = sum(daily$case_weight),
    an_acc_direct = an_acc_direct,
    counterfactual_cases = counterfactual_cases,
    af_acc_direct = af_acc_direct
  )
}

acc_check_identity <- function(
    total_cases,
    counterfactual_cases,
    an_acc_direct,
    af_acc_direct,
    tolerance = 1e-8
) {
  scale <- pmax(1, abs(total_cases), abs(an_acc_direct))
  case_error <- abs(
    total_cases - counterfactual_cases - an_acc_direct
  ) / scale
  af_error <- abs(
    af_acc_direct - 100 * an_acc_direct / total_cases
  )
  ok <- is.finite(case_error) & is.finite(af_error) &
    case_error <= tolerance & af_error <= tolerance
  if (any(!ok)) {
    stop("Direct ACC identity check failed.")
  }
  invisible(TRUE)
}

acc_bootstrap_indices <- function(
    gcm_count,
    repetitions = ACC_BOOTSTRAP_REPETITIONS,
    seed = ACC_BOOTSTRAP_SEED
) {
  set.seed(seed)
  matrix(
    sample.int(
      gcm_count,
      size = gcm_count * repetitions,
      replace = TRUE
    ),
    nrow = repetitions,
    ncol = gcm_count
  )
}

acc_bootstrap_median <- function(x, index_matrix) {
  x <- as.numeric(x)
  if (length(x) != ncol(index_matrix) || any(!is.finite(x))) {
    stop("Bootstrap input must contain one finite estimate per GCM.")
  }
  simulated <- apply(index_matrix, 1L, function(index) {
    median(x[index], na.rm = TRUE)
  })
  c(
    estimate = median(x),
    lower = unname(stats::quantile(simulated, 0.025)),
    upper = unname(stats::quantile(simulated, 0.975)),
    gcm_min = min(x),
    gcm_max = max(x)
  )
}

acc_summarise_gcm_metrics <- function(
    data,
    grouping_columns,
    metric_columns,
    gcm_levels,
    repetitions = ACC_BOOTSTRAP_REPETITIONS,
    seed = ACC_BOOTSTRAP_SEED
) {
  data <- data.table::as.data.table(data)
  grouping_columns <- as.character(grouping_columns)
  metric_columns <- as.character(metric_columns)
  gcm_levels <- as.character(gcm_levels)
  index_matrix <- acc_bootstrap_indices(
    length(gcm_levels),
    repetitions = repetitions,
    seed = seed
  )

  data[
    ,
    {
      group_data <- .SD[match(gcm_levels, gcm)]
      if (any(is.na(group_data$gcm))) {
        stop("A summary group does not contain all prespecified GCMs.")
      }
      values <- lapply(metric_columns, function(metric) {
        acc_bootstrap_median(group_data[[metric]], index_matrix)
      })
      names(values) <- metric_columns
      flattened <- unlist(lapply(names(values), function(metric) {
        stats::setNames(
          values[[metric]],
          paste(metric, names(values[[metric]]), sep = "_")
        )
      }))
      as.list(flattened)
    },
    by = grouping_columns,
    .SDcols = c("gcm", metric_columns)
  ]
}
