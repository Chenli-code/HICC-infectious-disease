###############################################################################
## Sensitivity analyses for anthropogenic climate-attributable burden
##
## This script does not overwrite the main analysis. It evaluates:
##   S1. Observed daily case weights versus day-of-year climatological weights.
##   S2. Alternative disease-specific lag windows.
##   S3. Population-weighted versus simple grid-weighted GCM temperatures.
##
## Reference-temperature sensitivity is intentionally omitted because the
## direct paired contrast is exactly reference-invariant.
##
## Primary direct estimand, for city i, day t, disease d and GCM g:
##   AF_ACC(t) = 1 - RR_hist-nat(t) / RR_factual(t)
##             = 1 - exp{eta_hist-nat(t) - eta_factual(t)}
##
## Run:
##   source("V01-10-run-sensitivity-analysis.R")
##
## Optional environment switches:
##   SENS_PREPARE_CLIMATE = true/false  (default true)
##   SENS_FIT_MODELS      = true/false  (default true)
##   SENS_RUN_ATTRIBUTION = true/false  (default true)
##   SENS_REBUILD_GRID    = true/false  (default false)
##   SENS_REFIT_MODELS    = true/false  (default false)
##   SENS_RECALCULATE     = true/false  (default false)
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(lubridate)
  library(stringr)
  library(dlnm)
  library(splines)
  library(gnm)
})

invisible(Sys.setlocale("LC_TIME", "C"))

ROOT <- Sys.getenv(
  "ATTR_PROJECT_ROOT",
  getwd()
)
source(file.path(ROOT, "03.0-user-route-config.R"))
source(USER_CORE)

SENSITIVITY_CACHE_VERSION <- paste0(
  ACC_FORMULA_VERSION,
  "-per-fit-city-set-v2"
)

SENS_DIR <- file.path(ROOT, "Sensitivity")
CACHE_DIR <- file.path(SENS_DIR, "cache")
GRID_PAIR_DIR <- file.path(CACHE_DIR, "climate_pairs_grid")
FIT_DIR <- file.path(CACHE_DIR, "lag_models")
ATTR_CACHE_DIR <- file.path(CACHE_DIR, "attribution")
RESULT_DIR <- file.path(SENS_DIR, "results")

invisible(lapply(
  c(SENS_DIR, CACHE_DIR, GRID_PAIR_DIR, FIT_DIR, ATTR_CACHE_DIR, RESULT_DIR),
  dir.create, recursive = TRUE, showWarnings = FALSE
))

CLIMATE_ROOT <- PATH_CLIMATE_DIR

CLIMATE_FILES <- list(
  population = list(
    factual = file.path(
      CLIMATE_ROOT, "all_gcms_all_years_factual_pop.rds"
    ),
    histnat = file.path(
      CLIMATE_ROOT, "all_gcms_all_years_hist_nat_pop.rds"
    )
  ),
  grid = list(
    factual = file.path(
      CLIMATE_ROOT, "all_gcms_all_years_factual_grid.rds"
    ),
    histnat = file.path(
      CLIMATE_ROOT, "all_gcms_all_years_hist_nat_grid.rds"
    )
  )
)

POP_PAIR_DIR <- USER_CACHE

flag <- function(name, default = FALSE) {
  value <- Sys.getenv(name, if (default) "true" else "false")
  tolower(value) %chin% c("true", "1", "yes", "y")
}

DO_PREPARE <- flag("SENS_PREPARE_CLIMATE", TRUE)
DO_FIT <- flag("SENS_FIT_MODELS", TRUE)
DO_ATTRIBUTION <- flag("SENS_RUN_ATTRIBUTION", TRUE)
REBUILD_GRID <- flag("SENS_REBUILD_GRID", FALSE)
REFIT_MODELS <- flag("SENS_REFIT_MODELS", FALSE)
RECALCULATE <- flag("SENS_RECALCULATE", FALSE)

safe_name <- function(x) gsub("[^A-Za-z0-9_-]", "_", x)
pad_code <- function(x) {
  str_pad(as.character(x), width = 6L, side = "right", pad = "0")
}

pair_file <- function(pair_dir, gcm) {
  file.path(
    pair_dir,
    paste0("climate_pair_", safe_name(gcm), ".rds")
  )
}

# Main lags reproduce the current primary models. Alternative lags are fixed
# before analysis and are not selected according to the resulting effect size.
LAG_SPECS <- data.table(
  classification = rep(
    c(
      "Respiratory diseases",
      "Gastrointestinal and enterovirus",
      "Vectorborne"
    ),
    each = 3L
  ),
  tag = rep(
    c("respiratory", "gastrointestinal", "vectorborne"),
    each = 3L
  ),
  lag_set = rep(c("Shorter lag", "Main lag", "Longer lag"), 3L),
  lag = c(
    3L, 7L, 14L,
    2L, 3L, 7L,
    7L, 14L, 21L
  ),
  lagnk = c(
    2L, 2L, 2L,
    1L, 2L, 2L,
    2L, 2L, 2L
  ),
  add_rain = rep(c(FALSE, FALSE, TRUE), each = 3L),
  reference = rep(
    c(
      "mrt",
      USER_NONRESP_REFERENCE_TYPE,
      USER_NONRESP_REFERENCE_TYPE
    ),
    each = 3L
  )
)

setkey(LAG_SPECS, tag, lag_set)

fwrite(
  LAG_SPECS,
  file.path(RESULT_DIR, "V01-Sensitivity-analysis-specifications.csv")
)

# ------------------------- climate-pair preparation -------------------------

prepare_climate_pairs <- function(
    weighting = c("population", "grid"),
    rebuild = FALSE
) {
  weighting <- match.arg(weighting)
  pair_dir <- if (weighting == "population") {
    POP_PAIR_DIR
  } else {
    GRID_PAIR_DIR
  }

  targets <- vapply(
    USER_GCMS,
    function(g) pair_file(pair_dir, g),
    character(1)
  )
  if (all(file.exists(targets)) && !rebuild) {
    message(
      "Reusing ", weighting, "-weighted climate pairs: ", pair_dir
    )
    return(pair_dir)
  }

  dir.create(pair_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- CLIMATE_FILES[[weighting]]
  stopifnot(file.exists(paths$factual), file.exists(paths$histnat))

  temporary_factual <- function(gcm) {
    file.path(
      pair_dir,
      paste0("_factual_", safe_name(gcm), ".rds")
    )
  }

  message("Loading full ", weighting, "-weighted factual temperatures...")
  fact <- as.data.table(readRDS(paths$factual))
  required <- c(
    "CityID", "date", "year", "tas_corrected_weighted", "gcm"
  )
  stopifnot(all(required %in% names(fact)))
  fact <- fact[year %in% 2008:2019, .(
    code = pad_code(CityID),
    date = as.Date(date),
    gcm = as.character(gcm),
    tf = as.numeric(tas_corrected_weighted)
  )]
  if (fact[, anyDuplicated(.SD), .SDcols = c("code", "date", "gcm")] > 0L) {
    stop("Duplicate city-date-GCM keys in factual ", weighting, " data.")
  }

  gcms <- sort(unique(fact$gcm))
  if (!setequal(gcms, USER_GCMS)) {
    stop("Factual GCM names do not match USER_GCMS for ", weighting, ".")
  }
  for (g in USER_GCMS) {
    saveRDS(
      fact[gcm == g, .(code, date, tf)],
      temporary_factual(g),
      compress = TRUE
    )
  }
  rm(fact)
  invisible(gc())

  message("Loading full ", weighting, "-weighted hist-nat temperatures...")
  nat <- as.data.table(readRDS(paths$histnat))
  stopifnot(all(required %in% names(nat)))
  nat <- nat[year %in% 2008:2019, .(
    code = pad_code(CityID),
    date = as.Date(date),
    gcm = as.character(gcm),
    tn = as.numeric(tas_corrected_weighted)
  )]
  if (nat[, anyDuplicated(.SD), .SDcols = c("code", "date", "gcm")] > 0L) {
    stop("Duplicate city-date-GCM keys in hist-nat ", weighting, " data.")
  }
  if (!setequal(unique(nat$gcm), USER_GCMS)) {
    stop("Hist-nat GCM names do not match USER_GCMS for ", weighting, ".")
  }

  qc <- vector("list", length(USER_GCMS))
  for (ii in seq_along(USER_GCMS)) {
    g <- USER_GCMS[ii]
    target <- pair_file(pair_dir, g)
    message(
      sprintf(
        "Pairing %s temperatures [%02d/%02d]: %s",
        weighting, ii, length(USER_GCMS), g
      )
    )
    f <- as.data.table(readRDS(temporary_factual(g)))
    n <- nat[gcm == g, .(code, date, tn)]
    pair <- merge(f, n, by = c("code", "date"), all = FALSE)
    pair <- pair[is.finite(tf) & is.finite(tn)]
    setorder(pair, code, date)
    saveRDS(pair, target, compress = TRUE)

    qc[[ii]] <- pair[, .(
      weighting = weighting,
      gcm = g,
      cities = uniqueN(code),
      paired_rows = .N,
      first_date = min(date),
      last_date = max(date),
      mean_signal = mean(tf - tn),
      median_signal = median(tf - tn)
    )]
    rm(f, n, pair)
    invisible(gc())
  }
  rm(nat)
  invisible(gc())

  temporary <- vapply(USER_GCMS, temporary_factual, character(1))
  unlink(temporary[file.exists(temporary)])
  fwrite(
    rbindlist(qc),
    file.path(
      RESULT_DIR,
      paste0("V01-climate-pair-QC-", weighting, ".csv")
    )
  )
  pair_dir
}

# ----------------------------- ERF estimation -------------------------------

main_fit_file <- function(tag) {
  file.path(USER_OUT, paste0("erf_pooled_", tag, ".rds"))
}

sensitivity_fit_file <- function(tag, lag_set, lag) {
  file.path(
    FIT_DIR,
    paste0(
      "erf_", tag, "_",
      safe_name(tolower(lag_set)), "_lag", lag, ".rds"
    )
  )
}

load_model_context <- function() {
  message("Loading health data for lag sensitivity models...")
  health_env <- new.env(parent = emptyenv())
  load(USER_PATH_HEALTH, envir = health_env)
  if (!exists("incidence", envir = health_env, inherits = FALSE)) {
    stop("Object 'incidence' was not found in USER_PATH_HEALTH.")
  }
  inc <- as.data.table(get("incidence", envir = health_env))
  rm(health_env)
  inc <- inc[year <= 2019 & code_p != "66" & ename != "COVID-19"]
  inc[
    ename %in% c("influenza A H1N1", "avian influenza H1N1"),
    classification := "Respiratory diseases"
  ]
  inc[, `:=`(
    code = pad_code(code_merge),
    date = as.Date(date)
  )]
  inc <- inc[
    date >= USER_START_DATE & date <= USER_END_DATE
  ]

  message("Loading observed weather for lag sensitivity models...")
  obs <- as.data.table(readRDS(USER_PATH_OBS))
  obs[, `:=`(
    code = pad_code(id),
    date = as.Date(date)
  )]
  obs <- obs[
    date >= USER_START_DATE & date <= USER_END_DATE,
    .(
      code, date,
      tmean = as.numeric(tmean),
      RH = as.numeric(RH),
      rainfall = as.numeric(rainfall),
      total_PM25 = as.numeric(total_PM25)
    )
  ]
  setorder(obs, code, date)

  vars <- c("tmean", "RH", "rainfall", "total_PM25")
  obs[, pro := substr(code, 1L, 2L)]
  province <- obs[
    ,
    lapply(.SD, function(z) mean(z, na.rm = TRUE)),
    by = .(pro, date),
    .SDcols = vars
  ]
  setnames(province, vars, paste0(vars, "_prov"))
  national <- obs[
    ,
    lapply(.SD, function(z) mean(z, na.rm = TRUE)),
    by = date,
    .SDcols = vars
  ]
  setnames(national, vars, paste0(vars, "_nat"))

  list(
    inc = inc,
    obs = obs,
    province = province,
    national = national,
    variables = vars,
    nonresp_reference_temperature = as.numeric(quantile(
      obs$tmean,
      USER_NONRESP_REFERENCE_PROBABILITY,
      na.rm = TRUE
    ))
  )
}

complete_observed_climate <- function(codes, context) {
  panel <- CJ(
    code = sort(unique(codes)),
    date = USER_DATES
  )
  panel[, pro := substr(code, 1L, 2L)]
  panel <- merge(
    panel,
    context$obs[
      ,
      c("code", "date", context$variables),
      with = FALSE
    ],
    by = c("code", "date"),
    all.x = TRUE
  )
  panel <- merge(
    panel,
    context$province,
    by = c("pro", "date"),
    all.x = TRUE
  )
  panel <- merge(
    panel,
    context$national,
    by = "date",
    all.x = TRUE
  )
  for (v in context$variables) {
    panel[, (v) := fcoalesce(
      get(v),
      get(paste0(v, "_prov")),
      get(paste0(v, "_nat"))
    )]
  }
  panel[
    ,
    c(
      "pro",
      paste0(context$variables, "_prov"),
      paste0(context$variables, "_nat")
    ) := NULL
  ]
  setorder(panel, code, date)
  panel
}

one_basis <- function(x, argvar) {
  do.call(onebasis, c(list(x = x), argvar))
}

lp_reduced <- function(x, coefficient, argvar, reference) {
  bx <- one_basis(x, argvar)
  b0 <- as.numeric(one_basis(reference, argvar))
  drop(sweep(bx, 2L, b0, "-") %*% coefficient)
}

fit_lag_model <- function(spec, context) {
  output <- sensitivity_fit_file(
    spec$tag, spec$lag_set, spec$lag
  )
  if (file.exists(output) && !REFIT_MODELS) {
    message("Reusing lag-sensitivity ERF: ", basename(output))
    return(output)
  }

  message(
    "Fitting ", spec$classification,
    " | ", spec$lag_set,
    " | lag 0-", spec$lag
  )
  counts <- context$inc[
    classification == spec$classification,
    .(n = .N),
    by = .(code, date)
  ]
  health_codes <- sort(unique(counts$code))
  dat <- complete_observed_climate(health_codes, context)
  dat <- merge(
    dat, counts,
    by = c("code", "date"),
    all.x = TRUE
  )
  dat[is.na(n), n := 0L]
  setorder(dat, code, date)
  dat[, `:=`(
    year = year(date),
    month = month(date),
    dow = wday(date),
    lag_PM25 = frollmean(
      total_PM25,
      n = spec$lag + 1L,
      align = "right",
      na.rm = TRUE
    )
  ), by = code]
  dat[, stratum := interaction(
    code, year, month, dow, drop = TRUE
  )]

  argvar <- list(
    fun = "ns",
    df = 3,
    knots = as.numeric(
      quantile(
        dat$tmean,
        c(0.10, 0.75, 0.90),
        na.rm = TRUE
      )
    ),
    Bound = range(dat$tmean, na.rm = TRUE)
  )
  arglag <- list(
    fun = "ns",
    knots = logknots(spec$lag, spec$lagnk),
    df = 3
  )

  cb_t <- crossbasis(
    dat$tmean,
    lag = spec$lag,
    argvar = argvar,
    arglag = arglag,
    group = dat$code
  )
  cb_h <- crossbasis(
    dat$RH,
    lag = spec$lag,
    argvar = list(fun = "ns", df = 3),
    arglag = arglag,
    group = dat$code
  )

  if (isTRUE(spec$add_rain)) {
    cb_rain <- crossbasis(
      dat$rainfall,
      lag = spec$lag,
      argvar = list(fun = "ns", df = 3),
      arglag = arglag,
      group = dat$code
    )
    model <- gnm(
      n ~ cb_t + cb_h + cb_rain + lag_PM25,
      family = quasipoisson(),
      eliminate = factor(stratum),
      data = dat,
      na.action = na.exclude
    )
  } else {
    model <- gnm(
      n ~ cb_t + cb_h + lag_PM25,
      family = quasipoisson(),
      eliminate = factor(stratum),
      data = dat,
      na.action = na.exclude
    )
  }

  initial_reference <- median(dat$tmean, na.rm = TRUE)
  reduced_initial <- crossreduce(
    cb_t, model,
    type = "overall",
    cen = initial_reference
  )

  if (spec$reference == USER_NONRESP_REFERENCE_TYPE) {
    reference <- context$nonresp_reference_temperature
  } else {
    temperature_grid <- seq(
      floor(quantile(dat$tmean, 0.01, na.rm = TRUE)),
      ceiling(quantile(dat$tmean, 0.99, na.rm = TRUE)),
      by = 0.1
    )
    reference <- temperature_grid[
      which.min(
        lp_reduced(
          temperature_grid,
          coef(reduced_initial),
          argvar,
          initial_reference
        )
      )
    ]
  }

  reduced <- crossreduce(
    cb_t, model,
    type = "overall",
    cen = reference
  )

  model_coefficient <- coef(model)
  temperature_index <- grep("^cb_t", names(model_coefficient))
  if (length(temperature_index) != ncol(cb_t)) {
    stop(
      "Could not identify all temperature cross-basis coefficients for ",
      spec$tag, " lag ", spec$lag, "."
    )
  }

  fit <- list(
    disease = spec$classification,
    tag = spec$tag,
    lag_set = spec$lag_set,
    lag = as.integer(spec$lag),
    lagnk = as.integer(spec$lagnk),
    add_rain = spec$add_rain,
    reference = spec$reference,
    cen = reference,
    crossover = dat[, .(code, date, n)],
    argvar = argvar,
    arglag = arglag,
    reduced_coef = coef(reduced),
    reduced_vcov = reduced$vcov,
    temp_coef = model_coefficient[temperature_index],
    temp_vcov = vcov(model)[
      temperature_index,
      temperature_index,
      drop = FALSE
    ],
    coverage = data.table(
      classification = spec$classification,
      tag = spec$tag,
      lag_set = spec$lag_set,
      lag = spec$lag,
      lagnk = spec$lagnk,
      reference = spec$reference,
      reference_temperature = reference,
      health_cities = length(health_codes),
      total_cases = sum(dat$n),
      observed_climate_cities = sum(
        health_codes %in% unique(context$obs$code)
      ),
      climate_imputed_cities = sum(
        !health_codes %in% unique(context$obs$code)
      )
    )
  )
  saveRDS(fit, output, compress = TRUE)
  fwrite(
    fit$coverage,
    file.path(
      RESULT_DIR,
      paste0(
        "V01-lag-model-coverage-", spec$tag,
        "-", safe_name(tolower(spec$lag_set)), ".csv"
      )
    )
  )

  rm(
    counts, dat, cb_t, cb_h, model,
    reduced_initial, reduced, fit
  )
  if (exists("cb_rain")) {
    rm(cb_rain)
  }
  invisible(gc())
  output
}

fit_file_for_spec <- function(spec) {
  if (spec$lag_set == "Main lag") {
    path <- main_fit_file(spec$tag)
    if (!file.exists(path)) {
      stop(
        "Main fit is missing: ", path,
        ". Run 03.1-fit-user-route-pooled-erf.R first."
      )
    }
    return(path)
  }
  sensitivity_fit_file(spec$tag, spec$lag_set, spec$lag)
}

# ---------------------------- attribution engine ----------------------------

evaluate_fit_for_gcm <- function(
    fit,
    climate_split,
    gcm,
    weighting,
    lag_set,
    baseline
) {
  cases <- as.data.table(fit$crossover)[, .(
    code = as.character(code),
    date = as.Date(date),
    n = as.numeric(n)
  )]
  case_weights <- acc_prepare_case_weights(
    cases = cases,
    dates = USER_DATES,
    baseline = baseline
  )
  cases_split <- split(case_weights, by = "code", keep.by = FALSE)
  codes <- sort(unique(case_weights$code))

  city_rows <- vector("list", length(codes))

  for (cc in seq_along(codes)) {
    code_i <- codes[cc]
    climate_i <- climate_split[[code_i]]
    cases_i <- cases_split[[code_i]]
    if (is.null(climate_i) || is.null(cases_i)) {
      next
    }

    dat <- merge(
      climate_i,
      cases_i,
      by = "date",
      all.x = TRUE
    )
    dat[is.na(n), n := 0]
    setorder(dat, date)
    direct <- acc_direct_city(
      factual_temperature = dat$tf,
      histnat_temperature = dat$tn,
      cases = dat$n,
      fit = fit
    )
    city_rows[[cc]] <- data.table(
      code = code_i,
      province_code = substr(code_i, 1L, 2L),
      classification = fit$disease,
      tag = fit$tag,
      weighting = weighting,
      lag_set = lag_set,
      lag = as.integer(fit$lag),
      baseline = baseline,
      gcm = gcm,
      total_cases = direct$total_cases,
      included_forward_cases = direct$included_forward_cases,
      direct_an = direct$an_acc_direct,
      counterfactual_cases = direct$counterfactual_cases,
      direct_af = direct$af_acc_direct,
      formula_version = ACC_FORMULA_VERSION
    )
  }

  city <- rbindlist(city_rows, fill = TRUE)
  if (!nrow(city)) {
    stop(
      "No valid attribution rows for ",
      fit$tag, ", ", weighting, ", ", lag_set, ", ", gcm
    )
  }
  city[!is.finite(direct_af), direct_af := NA_real_]

  national <- city[, .(
    total_cases = sum(total_cases),
    included_forward_cases = sum(included_forward_cases),
    direct_an = sum(direct_an),
    counterfactual_cases = sum(counterfactual_cases)
  ), by = .(
    classification, tag, weighting,
    lag_set, lag, baseline, gcm,
    formula_version
  )]
  national[, direct_af := 100 * direct_an / total_cases]
  acc_check_identity(
    national$total_cases,
    national$counterfactual_cases,
    national$direct_an,
    national$direct_af
  )

  list(city = city, national = national)
}

attribution_cache_file <- function(weighting, baseline, spec, gcm) {
  file.path(
    ATTR_CACHE_DIR,
    paste0(
      SENSITIVITY_CACHE_VERSION, "__",
      weighting, "__", baseline, "__", spec$tag, "__",
      safe_name(tolower(spec$lag_set)), "__",
      safe_name(gcm), ".rds"
    )
  )
}

run_attribution <- function(
    population_pair_dir,
    grid_pair_dir
) {
  # Population-weighted temperatures are evaluated under all lag settings.
  # Grid-weighted temperatures are evaluated only under the primary lag.
  analysis_specs <- rbindlist(list(
    copy(LAG_SPECS)[, `:=`(
      weighting = "population",
      baseline = "observed_daily"
    )],
    copy(LAG_SPECS[lag_set == "Main lag"])[, `:=`(
      weighting = "grid",
      baseline = "observed_daily"
    )],
    copy(LAG_SPECS[lag_set == "Main lag"])[, `:=`(
      weighting = "population",
      baseline = "day_of_year"
    )]
  ))

  city_outputs <- list()
  national_outputs <- list()
  output_index <- 0L

  for (weighting_i in c("population", "grid")) {
    pair_dir_i <- if (weighting_i == "population") {
      population_pair_dir
    } else {
      grid_pair_dir
    }
    specs_i <- analysis_specs[weighting == weighting_i]

    for (gg in seq_along(USER_GCMS)) {
      gcm_i <- USER_GCMS[gg]
      message(
        sprintf(
          "Attribution: %s weighting | GCM %02d/%02d | %s",
          weighting_i, gg, length(USER_GCMS), gcm_i
        )
      )

      pair <- readRDS(pair_file(pair_dir_i, gcm_i))

      for (ss in seq_len(nrow(specs_i))) {
        spec <- specs_i[ss]
        cache_i <- attribution_cache_file(
          weighting_i, spec$baseline, spec, gcm_i
        )
        fit_file_i <- fit_file_for_spec(spec)
        cache_metadata <- list(
          sensitivity_cache_version = SENSITIVITY_CACHE_VERSION,
          formula_version = ACC_FORMULA_VERSION,
          weighting = weighting_i,
          baseline = spec$baseline,
          lag_set = spec$lag_set,
          lag = as.integer(spec$lag),
          fit_md5 = acc_file_md5(fit_file_i),
          climate_md5 = acc_file_md5(
            pair_file(pair_dir_i, gcm_i)
          )
        )

        if (
          !RECALCULATE &&
            acc_cache_is_current(cache_i, cache_metadata)
        ) {
          result <- readRDS(cache_i)$result
        } else {
          fit <- readRDS(fit_file_i)
          fit_codes <- sort(unique(as.character(
            fit$crossover$code
          )))
          complete_pair <- acc_complete_climate_pair(
            pair = pair,
            codes = fit_codes,
            dates = USER_DATES
          )
          climate_split <- split(
            complete_pair,
            by = "code",
            keep.by = FALSE
          )
          result <- evaluate_fit_for_gcm(
            fit = fit,
            climate_split = climate_split,
            gcm = gcm_i,
            weighting = weighting_i,
            lag_set = spec$lag_set,
            baseline = spec$baseline
          )
          saveRDS(
            list(
              metadata = cache_metadata,
              result = result
            ),
            cache_i,
            compress = TRUE
          )
          rm(fit, fit_codes, complete_pair, climate_split)
          invisible(gc())
        }

        output_index <- output_index + 1L
        city_outputs[[output_index]] <- result$city
        national_outputs[[output_index]] <- result$national
        rm(result)
      }
      rm(pair)
      invisible(gc())
    }
  }

  city <- rbindlist(city_outputs, fill = TRUE)
  national <- rbindlist(national_outputs, fill = TRUE)
  setorder(
    city,
    weighting, baseline, tag, lag_set, gcm, code
  )
  setorder(
    national,
    weighting, baseline, tag, lag_set, gcm
  )

  fwrite(
    city,
    file.path(
      RESULT_DIR,
      "V01-sensitivity-attribution-city-by-gcm.csv"
    )
  )
  fwrite(
    national,
    file.path(
      RESULT_DIR,
      "V01-sensitivity-attribution-national-by-gcm.csv"
    )
  )

  invisible(list(city = city, national = national))
}

# ------------------------------- execution ----------------------------------

population_pair_dir <- POP_PAIR_DIR
grid_pair_dir <- GRID_PAIR_DIR

if (DO_PREPARE) {
  population_pair_dir <- prepare_climate_pairs(
    "population",
    rebuild = FALSE
  )
  grid_pair_dir <- prepare_climate_pairs(
    "grid",
    rebuild = REBUILD_GRID
  )
}

for (path_i in unlist(CLIMATE_FILES, recursive = TRUE)) {
  if (!file.exists(path_i)) {
    stop("Missing climate source file: ", path_i)
  }
}

if (DO_FIT) {
  alternative_specs <- LAG_SPECS[lag_set != "Main lag"]
  context <- load_model_context()
  for (ii in seq_len(nrow(alternative_specs))) {
    fit_lag_model(alternative_specs[ii], context)
  }
  rm(context)
  invisible(gc())
}

if (DO_ATTRIBUTION) {
  result <- run_attribution(
    population_pair_dir = population_pair_dir,
    grid_pair_dir = grid_pair_dir
  )
  message(
    "Completed ", nrow(result$national),
    " disease–GCM sensitivity estimates."
  )
}

message(
  "Sensitivity-analysis results written to: ",
  normalizePath(RESULT_DIR, winslash = "/", mustWork = FALSE)
)
