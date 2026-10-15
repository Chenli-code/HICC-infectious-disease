###############################################################################
## Shared engines for disease-relevant attributable contribution
##
## Respiratory diseases:
##   primary AN = all temperature-attributable AN relative to the fitted MMT.
##
## Gastrointestinal/enterovirus and vector-borne diseases:
##   primary AN = positive temperature-attributable AN accumulated during each
##   city's four hottest consecutive calendar months.
##
## This file defines functions only. It intentionally performs no analysis and
## writes no output, so it can be sourced safely by Figures 3-5 and the
## warm-season attribution script.
###############################################################################

reference_ac_city <- function(
    factual_temperature,
    histnat_temperature,
    cases,
    fit
) {
  reference <- as.numeric(fit$cen)
  forward_cases <- acc_forward_weights(cases, fit$lag)
  basis_factual <- acc_one_basis(
    factual_temperature,
    fit$argvar
  )
  basis_histnat <- acc_one_basis(
    histnat_temperature,
    fit$argvar
  )
  basis_reference <- as.numeric(acc_one_basis(
    reference,
    fit$argvar
  ))

  valid <- is.finite(forward_cases) &
    is.finite(factual_temperature) &
    is.finite(histnat_temperature) &
    rowSums(!is.finite(basis_factual)) == 0L &
    rowSums(!is.finite(basis_histnat)) == 0L

  total_cases <- sum(cases, na.rm = TRUE)
  empty_result <- function() {
    data.table::data.table(
      total_cases = total_cases,
      included_forward_cases = 0,
      factual_days_above_ref = 0L,
      histnat_days_above_ref = 0L,
      factual_days_below_ref = 0L,
      histnat_days_below_ref = 0L,
      an_factual_above_ref = 0,
      an_histnat_above_ref = 0,
      an_factual_below_ref = 0,
      an_histnat_below_ref = 0,
      an_acc_above_ref = 0,
      ac_above_ref = NA_real_,
      an_acc_below_ref = 0,
      cold_avoided_fraction = NA_real_,
      af_factual_above_ref_all_cases = 0,
      af_acc_above_ref_all_cases = 0,
      af_acc_below_ref_all_cases = 0
    )
  }
  if (!any(valid) || total_cases <= 0) {
    return(empty_result())
  }

  forward_total <- sum(forward_cases[valid], na.rm = TRUE)
  if (!is.finite(forward_total) || forward_total <= 0) {
    return(empty_result())
  }
  weights <- forward_cases[valid] * total_cases / forward_total

  centred_factual <- sweep(
    basis_factual[valid, , drop = FALSE],
    2L,
    basis_reference,
    "-"
  )
  centred_histnat <- sweep(
    basis_histnat[valid, , drop = FALSE],
    2L,
    basis_reference,
    "-"
  )
  eta_factual <- drop(
    centred_factual %*% fit$reduced_coef
  )
  eta_histnat <- drop(
    centred_histnat %*% fit$reduced_coef
  )

  factual_valid <- factual_temperature[valid]
  histnat_valid <- histnat_temperature[valid]
  factual_above <- factual_valid > reference
  histnat_above <- histnat_valid > reference
  factual_below <- factual_valid < reference
  histnat_below <- histnat_valid < reference

  daily_an_factual <- weights * (1 - exp(-eta_factual))
  daily_an_histnat <- weights * (1 - exp(-eta_histnat))

  an_factual_above <- sum(
    daily_an_factual * factual_above,
    na.rm = TRUE
  )
  an_histnat_above <- sum(
    daily_an_histnat * histnat_above,
    na.rm = TRUE
  )
  an_factual_below <- sum(
    daily_an_factual * factual_below,
    na.rm = TRUE
  )
  an_histnat_below <- sum(
    daily_an_histnat * histnat_below,
    na.rm = TRUE
  )
  an_acc_above <- an_factual_above - an_histnat_above
  an_acc_below <- an_factual_below - an_histnat_below

  data.table::data.table(
    total_cases = total_cases,
    included_forward_cases = sum(weights),
    factual_days_above_ref = sum(factual_above),
    histnat_days_above_ref = sum(histnat_above),
    factual_days_below_ref = sum(factual_below),
    histnat_days_below_ref = sum(histnat_below),
    an_factual_above_ref = an_factual_above,
    an_histnat_above_ref = an_histnat_above,
    an_factual_below_ref = an_factual_below,
    an_histnat_below_ref = an_histnat_below,
    an_acc_above_ref = an_acc_above,
    ac_above_ref = if (
      abs(an_factual_above) > sqrt(.Machine$double.eps)
    ) {
      100 * an_acc_above / an_factual_above
    } else {
      NA_real_
    },
    an_acc_below_ref = an_acc_below,
    cold_avoided_fraction = if (
      an_factual_below > 0 && an_histnat_below > 0
    ) {
      100 * (
        an_histnat_below - an_factual_below
      ) / an_histnat_below
    } else {
      NA_real_
    },
    af_factual_above_ref_all_cases =
      100 * an_factual_above / total_cases,
    af_acc_above_ref_all_cases =
      100 * an_acc_above / total_cases,
    af_acc_below_ref_all_cases =
      100 * an_acc_below / total_cases
  )
}

reference_ac_select_primary <- function(
    result,
    class_tag
) {
  result <- data.table::as.data.table(result)
  respiratory <- identical(
    as.character(class_tag),
    "respiratory"
  )
  result[, `:=`(
    reference_burden_definition = if (respiratory) {
      "All temperature-attributable AN relative to MMT"
    } else {
      "AN above the fitted numerical reference"
    },
    an_factual_reference = if (respiratory) {
      an_factual_above_ref + an_factual_below_ref
    } else {
      an_factual_above_ref
    },
    an_histnat_reference = if (respiratory) {
      an_histnat_above_ref + an_histnat_below_ref
    } else {
      an_histnat_above_ref
    }
  )]
  result[, an_acc_reference :=
    an_factual_reference - an_histnat_reference]
  result[, ac_reference := data.table::fifelse(
    is.finite(an_factual_reference) &
      abs(an_factual_reference) >
        sqrt(.Machine$double.eps),
    100 * an_acc_reference / an_factual_reference,
    NA_real_
  )]
  result
}

warm4_ac_city <- function(
    factual_temperature,
    histnat_temperature,
    cases,
    dates,
    warm_month_numbers,
    fit
) {
  total_cases <- sum(cases, na.rm = TRUE)
  forward_cases <- acc_forward_weights(cases, fit$lag)
  basis_factual <- acc_one_basis(
    factual_temperature,
    fit$argvar
  )
  basis_histnat <- acc_one_basis(
    histnat_temperature,
    fit$argvar
  )
  basis_reference <- as.numeric(acc_one_basis(
    as.numeric(fit$cen),
    fit$argvar
  ))
  date_month <- as.integer(format(
    as.Date(dates),
    "%m"
  ))

  valid <- is.finite(forward_cases) &
    is.finite(factual_temperature) &
    is.finite(histnat_temperature) &
    rowSums(!is.finite(basis_factual)) == 0L &
    rowSums(!is.finite(basis_histnat)) == 0L

  empty_result <- function() {
    data.table::data.table(
      total_cases_warm4 = total_cases,
      included_forward_cases_warm4 = 0,
      warm4_case_weight = 0,
      warm4_city_days = 0L,
      an_factual_warm4 = 0,
      an_histnat_warm4 = 0,
      an_acc_warm4 = 0,
      ac_warm4 = NA_real_,
      an_factual_warm4_signed = 0,
      an_histnat_warm4_signed = 0
    )
  }
  if (!any(valid) || total_cases <= 0) {
    return(empty_result())
  }

  forward_total <- sum(
    forward_cases[valid],
    na.rm = TRUE
  )
  if (!is.finite(forward_total) || forward_total <= 0) {
    return(empty_result())
  }

  weights <- forward_cases[valid] *
    total_cases / forward_total
  warm <- date_month[valid] %in%
    as.integer(warm_month_numbers)
  warm_case_weight <- sum(
    weights[warm],
    na.rm = TRUE
  )
  if (!any(warm) || warm_case_weight <= 0) {
    return(empty_result())
  }

  centred_factual <- sweep(
    basis_factual[valid, , drop = FALSE],
    2L,
    basis_reference,
    "-"
  )
  centred_histnat <- sweep(
    basis_histnat[valid, , drop = FALSE],
    2L,
    basis_reference,
    "-"
  )
  eta_factual <- drop(
    centred_factual %*% fit$reduced_coef
  )
  eta_histnat <- drop(
    centred_histnat %*% fit$reduced_coef
  )
  daily_an_factual <- weights * (
    1 - exp(-eta_factual)
  )
  daily_an_histnat <- weights * (
    1 - exp(-eta_histnat)
  )
  an_factual_warm4 <- sum(
    pmax(daily_an_factual[warm], 0),
    na.rm = TRUE
  )
  an_histnat_warm4 <- sum(
    pmax(daily_an_histnat[warm], 0),
    na.rm = TRUE
  )
  an_acc_warm4 <-
    an_factual_warm4 - an_histnat_warm4

  data.table::data.table(
    total_cases_warm4 = total_cases,
    included_forward_cases_warm4 =
      sum(weights),
    warm4_case_weight = warm_case_weight,
    warm4_city_days = sum(warm),
    an_factual_warm4 = an_factual_warm4,
    an_histnat_warm4 = an_histnat_warm4,
    an_acc_warm4 = an_acc_warm4,
    ac_warm4 = if (
      an_factual_warm4 >
        sqrt(.Machine$double.eps)
    ) {
      100 * an_acc_warm4 /
        an_factual_warm4
    } else {
      NA_real_
    },
    an_factual_warm4_signed = sum(
      daily_an_factual[warm],
      na.rm = TRUE
    ),
    an_histnat_warm4_signed = sum(
      daily_an_histnat[warm],
      na.rm = TRUE
    )
  )
}

select_disease_relevant_ac <- function(
    reference_result,
    warm4_result,
    class_tag
) {
  reference_result <- data.table::as.data.table(
    reference_result
  )
  warm4_result <- data.table::as.data.table(
    warm4_result
  )
  if (nrow(reference_result) != nrow(warm4_result)) {
    stop(
      "Reference and warm-season results must have the same number of rows."
    )
  }

  respiratory <- identical(
    as.character(class_tag),
    "respiratory"
  )
  result <- data.table::copy(reference_result)
  result[, `:=`(
    primary_burden_definition = if (respiratory) {
      "All temperature-attributable AN relative to MMT"
    } else {
      paste0(
        "Positive temperature-attributable AN during each city's ",
        "four hottest consecutive months"
      )
    },
    an_factual_primary = if (respiratory) {
      an_factual_above_ref + an_factual_below_ref
    } else {
      warm4_result$an_factual_warm4
    },
    an_histnat_primary = if (respiratory) {
      an_histnat_above_ref + an_histnat_below_ref
    } else {
      warm4_result$an_histnat_warm4
    }
  )]
  result[, an_acc_primary :=
    an_factual_primary - an_histnat_primary]
  result[, ac_primary := data.table::fifelse(
    is.finite(an_factual_primary) &
      abs(an_factual_primary) >
        sqrt(.Machine$double.eps),
    100 * an_acc_primary / an_factual_primary,
    NA_real_
  )]
  result
}

reference_ac_bootstrap_summary <- function(
    factual,
    histnat,
    total_cases,
    bootstrap_indices
) {
  factual <- as.numeric(factual)
  histnat <- as.numeric(histnat)
  total_cases <- as.numeric(total_cases)
  gcm_count <- ncol(bootstrap_indices)
  if (
    length(factual) != gcm_count ||
      length(histnat) != gcm_count ||
      length(total_cases) != gcm_count ||
      any(!is.finite(factual)) ||
      any(!is.finite(histnat)) ||
      any(!is.finite(total_cases))
  ) {
    stop("One finite paired reference-burden estimate is required per GCM.")
  }

  derive <- function(index) {
    factual_i <- mean(factual[index])
    histnat_i <- mean(histnat[index])
    total_i <- mean(total_cases[index])
    acc_i <- factual_i - histnat_i
    ac_i <- if (
      abs(factual_i) > sqrt(.Machine$double.eps)
    ) {
      100 * acc_i / factual_i
    } else {
      NA_real_
    }
    c(
      factual = factual_i,
      histnat = histnat_i,
      acc = acc_i,
      ac = ac_i,
      factual_af_all_cases = 100 * factual_i / total_i,
      acc_af_all_cases = 100 * acc_i / total_i
    )
  }

  point <- derive(seq_along(factual))
  simulated <- t(apply(
    bootstrap_indices,
    1L,
    derive
  ))
  interval <- apply(
    simulated,
    2L,
    stats::quantile,
    probs = c(0.025, 0.975),
    na.rm = TRUE,
    names = FALSE
  )
  list(
    point = point,
    lower = interval[1L, ],
    upper = interval[2L, ],
    simulated = simulated
  )
}
