###############################################################################
## Figure 3 | Fit 24 disease-specific national pooled ERFs
##
## Strategy is identical to the major-class route used in Figure 2:
##   - time-stratified national pooled case-crossover model;
##   - city-specific lag reset through crossbasis(group = code);
##   - temperature ns(df=3; knots at P10/P75/P90);
##   - distributed lag from the published disease-specific lag table;
##   - humidity DLNM and lagged PM2.5 adjustment;
##   - additional rainfall DLNM for vector-borne diseases;
##   - respiratory reference = minimum-risk temperature;
##   - gastrointestinal/vector-borne curve reference = national observed P50;
##   - their disease-relevant attribution window is defined separately as each
##     city's four hottest consecutive months.
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(lubridate)
  library(stringr)
  library(dlnm)
  library(splines)
  library(gnm)
})
Sys.setlocale("LC_TIME", "C")
ROOT <- Sys.getenv(
  "ATTR_PROJECT_ROOT",
  getwd()
)
source(file.path(ROOT, "10.0-Figure3-disease-specification.R"))

REUSE <- tolower(Sys.getenv("FIG3_REUSE_ERF", "true")) == "true"
FILTER <- trimws(Sys.getenv("FIG3_DISEASE_FILTER", ""))
STATUS_FILE <- file.path(FIG3_DIR, "10.1-disease-erf-fit-status.csv")
CURVE_FILE <- file.path(FIG3_DIR, "10.1-disease-erf-curves.csv")
COVERAGE_FILE <- file.path(FIG3_DIR, "10.1-disease-erf-coverage.csv")

specs <- copy(FIG3_DISEASE_SPECS)
if (nzchar(FILTER)) {
  wanted <- trimws(strsplit(FILTER, "\\|")[[1]])
  specs <- specs[disease %chin% wanted | disease_tag %chin% wanted]
  if (!nrow(specs)) stop("FIG3_DISEASE_FILTER did not match any disease.")
}

message("Loading incidence data once for disease-specific ERFs...")
load(USER_PATH_HEALTH)
stopifnot(exists("incidence"))
inc <- as.data.table(incidence)
rm(incidence); invisible(gc())
inc <- inc[
  year <= 2019 & code_p != "66" & ename != "COVID-19" &
    as.Date(date) >= USER_START_DATE & as.Date(date) <= USER_END_DATE
]
inc[, `:=`(
  code = user_pad_code(code_merge),
  date = as.Date(date),
  ename = as.character(ename),
  classification = as.character(classification)
)]
inc <- inc[ename %chin% specs$disease, .(code, date, ename, classification)]

message("Loading observed weather once...")
obs <- as.data.table(readRDS(USER_PATH_OBS))
obs[, `:=`(code = user_pad_code(id), date = as.Date(date))]
obs <- obs[
  date >= USER_START_DATE & date <= USER_END_DATE,
  .(code, date, tmean = as.numeric(tmean), RH = as.numeric(RH),
    rainfall = as.numeric(rainfall), total_PM25 = as.numeric(total_PM25))
]
setorder(obs, code, date)
weather_vars <- c("tmean", "RH", "rainfall", "total_PM25")
obs[, pro := substr(code, 1, 2)]
prov <- obs[, lapply(.SD, function(z) mean(z, na.rm = TRUE)),
            by = .(pro, date), .SDcols = weather_vars]
setnames(prov, weather_vars, paste0(weather_vars, "_prov"))
nat <- obs[, lapply(.SD, function(z) mean(z, na.rm = TRUE)),
           by = date, .SDcols = weather_vars]
setnames(nat, weather_vars, paste0(weather_vars, "_nat"))
global_nonresp_reference <- as.numeric(quantile(
  obs$tmean,
  USER_NONRESP_REFERENCE_PROBABILITY,
  na.rm = TRUE
))
complete_observed_weather <- function(codes) {
  panel <- CJ(code = sort(unique(codes)), date = USER_DATES)
  panel[, pro := substr(code, 1, 2)]
  panel <- merge(
    panel, obs[, c("code", "date", weather_vars), with = FALSE],
    by = c("code", "date"), all.x = TRUE
  )
  panel <- merge(panel, prov, by = c("pro", "date"), all.x = TRUE)
  panel <- merge(panel, nat, by = "date", all.x = TRUE)
  for (v in weather_vars) {
    panel[, (v) := fcoalesce(
      get(v), get(paste0(v, "_prov")), get(paste0(v, "_nat"))
    )]
  }
  panel[, c("pro", paste0(weather_vars, "_prov"),
            paste0(weather_vars, "_nat")) := NULL]
  setorder(panel, code, date)
  panel
}

one_basis <- function(x, argvar) do.call(onebasis, c(list(x = x), argvar))
lp_reduced <- function(x, b, argvar, cen) {
  bx <- one_basis(x, argvar)
  bc <- as.numeric(one_basis(cen, argvar))
  drop(sweep(bx, 2, bc, "-") %*% b)
}

read_existing <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(readRDS(path), error = function(e) NULL)
}

status_rows <- if (file.exists(STATUS_FILE)) fread(STATUS_FILE) else data.table(
  disease = character(), disease_tag = character(), classification = character(),
  lag = integer(), lagnk = integer(), lag_source = character(),
  status = character(), message = character(), fit_file = character()
)
curve_rows <- list()
coverage_rows <- list()

for (ii in seq_len(nrow(specs))) {
  s <- specs[ii]
  fit_file <- figure3_fit_file(s$disease_tag)
  old <- if (REUSE) read_existing(fit_file) else NULL
  if (
    !is.null(old) &&
      identical(as.integer(old$lag), as.integer(s$lag)) &&
      identical(as.integer(old$lagnk), as.integer(s$lagnk)) &&
      identical(
        tolower(as.character(old$reference)),
        tolower(as.character(s$reference))
      )
  ) {
    message(sprintf("[%d/%d] Reusing ERF: %s", ii, nrow(specs), s$disease))
    status_rows <- rbind(
      status_rows[disease != s$disease],
      data.table(
        disease = s$disease, disease_tag = s$disease_tag,
        classification = s$classification, lag = s$lag, lagnk = s$lagnk,
        lag_source = s$lag_source, status = "reused", message = "",
        fit_file = fit_file
      ), fill = TRUE
    )
    if (!is.null(old$curve)) curve_rows[[length(curve_rows) + 1L]] <- old$curve
    if (!is.null(old$coverage)) coverage_rows[[length(coverage_rows) + 1L]] <- old$coverage
    fwrite(status_rows, STATUS_FILE)
    next
  }

  message(sprintf("[%d/%d] Fitting disease ERF: %s | lag 0-%d",
                  ii, nrow(specs), s$disease, s$lag))
  counts <- inc[ename == s$disease, .(n = .N), by = .(code, date)]
  total_cases <- sum(counts$n)
  health_codes <- sort(unique(counts$code))
  warning_text <- character()

  result <- tryCatch(
    withCallingHandlers({
      if (!total_cases || !length(health_codes)) stop("No eligible cases.")
      dat <- complete_observed_weather(health_codes)
      dat <- merge(dat, counts, by = c("code", "date"), all.x = TRUE)
      dat[is.na(n), n := 0L]
      setorder(dat, code, date)
      dat[, `:=`(
        year = year(date), month = month(date), dow = wday(date),
        lag_PM25 = frollmean(
          total_PM25, n = s$lag + 1L, align = "right", na.rm = TRUE
        )
      ), by = code]
      dat[, stratum := interaction(code, year, month, dow, drop = TRUE)]

      argvar <- list(
        fun = "ns", df = 3,
        knots = as.numeric(quantile(
          dat$tmean, c(0.10, 0.75, 0.90), na.rm = TRUE
        )),
        Bound = range(dat$tmean, na.rm = TRUE)
      )
      arglag <- list(
        fun = "ns", knots = logknots(s$lag, s$lagnk), df = 3
      )
      cb_t <- crossbasis(
        dat$tmean, lag = s$lag, argvar = argvar, arglag = arglag,
        group = dat$code
      )
      cb_h <- crossbasis(
        dat$RH, lag = s$lag,
        argvar = list(fun = "ns", df = 3), arglag = arglag,
        group = dat$code
      )

      if (isTRUE(s$add_rain)) {
        cb_rain <- crossbasis(
          dat$rainfall, lag = s$lag,
          argvar = list(fun = "ns", df = 3), arglag = arglag,
          group = dat$code
        )
        model <- gnm(
          n ~ cb_t + cb_h + cb_rain + lag_PM25,
          family = quasipoisson(), eliminate = factor(stratum),
          data = dat, na.action = na.exclude
        )
      } else {
        model <- gnm(
          n ~ cb_t + cb_h + lag_PM25,
          family = quasipoisson(), eliminate = factor(stratum),
          data = dat, na.action = na.exclude
        )
      }

      initial_cen <- median(dat$tmean, na.rm = TRUE)
      red0 <- crossreduce(cb_t, model, type = "overall", cen = initial_cen)
      if (s$reference == USER_NONRESP_REFERENCE_TYPE) {
        cen <- global_nonresp_reference
      } else {
        grid_ref <- seq(
          floor(quantile(dat$tmean, 0.01, na.rm = TRUE)),
          ceiling(quantile(dat$tmean, 0.99, na.rm = TRUE)), by = 0.1
        )
        lp0 <- lp_reduced(grid_ref, coef(red0), argvar, initial_cen)
        if (!any(is.finite(lp0))) stop("Minimum-risk temperature was not estimable.")
        cen <- grid_ref[which.min(replace(lp0, !is.finite(lp0), Inf))]
      }
      red <- crossreduce(cb_t, model, type = "overall", cen = cen)

      model_coef <- coef(model)
      temp_index <- grep("^cb_t", names(model_coef))
      if (length(temp_index) != ncol(cb_t)) {
        stop("Could not identify all temperature cross-basis coefficients.")
      }
      temp_coef <- model_coef[temp_index]
      temp_vcov <- vcov(model)[temp_index, temp_index, drop = FALSE]
      if (any(!is.finite(temp_coef)) || any(!is.finite(temp_vcov))) {
        stop("Temperature coefficients/covariance are non-finite.")
      }

      basis_template <- crossbasis(
        rep(cen, max(45L, s$lag + 5L)), lag = s$lag,
        argvar = argvar, arglag = arglag
      )
      grid_plot <- seq(
        floor(quantile(dat$tmean, 0.01, na.rm = TRUE)),
        ceiling(quantile(dat$tmean, 0.99, na.rm = TRUE)), length.out = 300
      )
      bx <- one_basis(grid_plot, argvar)
      bc <- as.numeric(one_basis(cen, argvar))
      Xred <- sweep(bx, 2, bc, "-")
      eta <- drop(Xred %*% coef(red))
      se <- sqrt(pmax(0, rowSums((Xred %*% red$vcov) * Xred)))
      curve <- data.table(
        disease = s$disease, disease_tag = s$disease_tag,
        classification = s$classification, lag = s$lag, cen = cen,
        temperature = grid_plot,
        rr = exp(eta), rr_low = exp(eta - 1.96 * se),
        rr_high = exp(eta + 1.96 * se)
      )
      coverage <- data.table(
        disease = s$disease, disease_tag = s$disease_tag,
        classification = s$classification, lag = s$lag, lagnk = s$lagnk,
        lag_source = s$lag_source, reference = s$reference, cen = cen,
        health_cities = length(health_codes), total_cases = total_cases,
        nonzero_case_days = nrow(counts),
        observed_climate_cities = sum(health_codes %in% unique(obs$code)),
        climate_imputed_cities = sum(!health_codes %in% unique(obs$code)),
        temp_condition_number = kappa(temp_vcov)
      )
      fit <- list(
        disease = s$disease, disease_tag = s$disease_tag,
        classification = s$classification, class_tag = s$class_tag,
        lag = as.integer(s$lag), lagnk = as.integer(s$lagnk),
        lag_source = s$lag_source, add_rain = isTRUE(s$add_rain),
        reference = s$reference, cen = cen,
        basis_template = basis_template,
        temp_coef = temp_coef, temp_vcov = temp_vcov,
        argvar = argvar, arglag = arglag,
        reduced_coef = coef(red), reduced_vcov = red$vcov,
        sparse_counts = counts, health_codes = health_codes,
        total_cases = total_cases, coverage = coverage, curve = curve,
        warnings = warning_text
      )
      saveRDS(fit, fit_file, compress = TRUE)
      list(curve = curve, coverage = coverage)
    }, warning = function(w) {
      warning_text <<- unique(c(warning_text, conditionMessage(w)))
      invokeRestart("muffleWarning")
    }),
    error = function(e) e
  )

  if (inherits(result, "error")) {
    message("  FAILED: ", conditionMessage(result))
    status_row <- data.table(
      disease = s$disease, disease_tag = s$disease_tag,
      classification = s$classification, lag = s$lag, lagnk = s$lagnk,
      lag_source = s$lag_source, status = "failed",
      message = conditionMessage(result), fit_file = fit_file
    )
  } else {
    curve_rows[[length(curve_rows) + 1L]] <- result$curve
    coverage_rows[[length(coverage_rows) + 1L]] <- result$coverage
    status_row <- data.table(
      disease = s$disease, disease_tag = s$disease_tag,
      classification = s$classification, lag = s$lag, lagnk = s$lagnk,
      lag_source = s$lag_source, status = "fitted",
      message = paste(warning_text, collapse = " | "), fit_file = fit_file
    )
  }
  status_rows <- rbind(status_rows[disease != s$disease], status_row, fill = TRUE)
  fwrite(status_rows, STATUS_FILE)
  if (length(curve_rows)) fwrite(rbindlist(curve_rows, fill = TRUE), CURVE_FILE)
  if (length(coverage_rows)) fwrite(rbindlist(coverage_rows, fill = TRUE), COVERAGE_FILE)
  rm(counts, result); invisible(gc())
}

status_rows[, order_id__ := match(disease, FIG3_DISEASE_SPECS$disease)]
setorder(status_rows, order_id__)
status_rows[, order_id__ := NULL]
fwrite(status_rows, STATUS_FILE)
print(status_rows[, .(disease, lag, status, message)])
