###############################################################################
## 03.1 Fit one national pooled case-crossover/DLNM ERF per major class
##
## This follows the user's original route:
##   - one national ERF for each major disease class;
##   - all health cities are retained;
##   - missing observed-climate cities use province-day means, then national-day
##     means as a final fallback;
##   - the lag basis is reset at every city boundary with group = code.
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
source(file.path(ROOT, "03.0-user-route-config.R"))

REUSE <- tolower(Sys.getenv("USER_ROUTE_REUSE_ERF", "true")) == "true"
LAGNK <- 2L

message("Loading health data once for all three major classes...")
load(USER_PATH_HEALTH)
stopifnot(exists("incidence"))
inc <- as.data.table(incidence)
rm(incidence); invisible(gc())
inc <- inc[year <= 2019 & code_p != "66" & ename != "COVID-19"]
inc[ename %in% c("influenza A H1N1", "avian influenza H1N1"),
    classification := "Respiratory diseases"]
inc[, `:=`(code = user_pad_code(code_merge), date = as.Date(date))]
inc <- inc[date >= USER_START_DATE & date <= USER_END_DATE]

message("Loading observed climate once...")
obs <- as.data.table(readRDS(USER_PATH_OBS))
obs[, `:=`(code = user_pad_code(id), date = as.Date(date))]
obs <- obs[date >= USER_START_DATE & date <= USER_END_DATE,
           .(code, date,
             tmean = as.numeric(tmean), RH = as.numeric(RH),
             rainfall = as.numeric(rainfall),
             total_PM25 = as.numeric(total_PM25))]
setorder(obs, code, date)
vars <- c("tmean", "RH", "rainfall", "total_PM25")
obs[, pro := substr(code, 1, 2)]
prov <- obs[, lapply(.SD, function(z) mean(z, na.rm = TRUE)),
            by = .(pro, date), .SDcols = vars]
setnames(prov, vars, paste0(vars, "_prov"))
nat <- obs[, lapply(.SD, function(z) mean(z, na.rm = TRUE)),
           by = date, .SDcols = vars]
setnames(nat, vars, paste0(vars, "_nat"))
global_nonresp_reference <- as.numeric(quantile(
  obs$tmean,
  USER_NONRESP_REFERENCE_PROBABILITY,
  na.rm = TRUE
))
complete_observed_climate <- function(codes) {
  panel <- CJ(code = sort(unique(codes)), date = USER_DATES)
  panel[, pro := substr(code, 1, 2)]
  panel <- merge(panel, obs[, c("code", "date", vars), with = FALSE],
                 by = c("code", "date"), all.x = TRUE)
  panel <- merge(panel, prov, by = c("pro", "date"), all.x = TRUE)
  panel <- merge(panel, nat, by = "date", all.x = TRUE)
  for (v in vars) {
    panel[, (v) := fcoalesce(
      get(v), get(paste0(v, "_prov")), get(paste0(v, "_nat"))
    )]
  }
  panel[, c("pro", paste0(vars, "_prov"), paste0(vars, "_nat")) := NULL]
  setorder(panel, code, date)
  panel
}

one_basis <- function(x, argvar) do.call(onebasis, c(list(x = x), argvar))
lp_reduced <- function(x, b, argvar, cen) {
  bx <- one_basis(x, argvar)
  bc <- as.numeric(one_basis(cen, argvar))
  drop(sweep(bx, 2, bc, "-") %*% b)
}

fit_status <- vector("list", nrow(USER_SPECS))
for (i in seq_len(nrow(USER_SPECS))) {
  s <- USER_SPECS[i, ]
  fit_file <- file.path(USER_OUT, paste0("erf_pooled_", s$tag, ".rds"))
  if (REUSE && file.exists(fit_file)) {
    existing_fit <- tryCatch(
      readRDS(fit_file),
      error = function(e) NULL
    )
    reference_current <- !is.null(existing_fit) &&
      identical(
        tolower(as.character(existing_fit$reference)),
        tolower(as.character(s$reference))
      )
    if (reference_current) {
      message("Reusing national ERF: ", s$classification)
      fit_status[[i]] <- data.table(
        classification = s$classification, tag = s$tag,
        status = "reused", fit_file = fit_file
      )
      next
    }
    message(
      "Existing ERF uses an outdated reference; refitting: ",
      s$classification
    )
  }

  message("Preparing all health cities: ", s$classification)
  counts <- inc[classification == s$classification,
                .(n = .N), by = .(code, date)]
  health_codes <- sort(unique(counts$code))
  dat <- complete_observed_climate(health_codes)
  dat <- merge(dat, counts, by = c("code", "date"), all.x = TRUE)
  dat[is.na(n), n := 0L]
  setorder(dat, code, date)
  dat[, `:=`(
    year = year(date), month = month(date), dow = wday(date),
    lag_PM25 = frollmean(total_PM25, n = s$lag + 1L,
                         align = "right", na.rm = TRUE)
  ), by = code]
  dat[, stratum := interaction(code, year, month, dow, drop = TRUE)]

  argvar <- list(
    fun = "ns", df = 3,
    knots = as.numeric(quantile(dat$tmean, c(0.10, 0.75, 0.90), na.rm = TRUE)),
    Bound = range(dat$tmean, na.rm = TRUE)
  )
  arglag <- list(fun = "ns", knots = logknots(s$lag, LAGNK), df = 3)
  cb_t <- crossbasis(dat$tmean, lag = s$lag, argvar = argvar,
                     arglag = arglag, group = dat$code)
  cb_h <- crossbasis(dat$RH, lag = s$lag,
                     argvar = list(fun = "ns", df = 3),
                     arglag = arglag, group = dat$code)

  message("Fitting national pooled ERF: ", s$classification,
          " | cities=", length(health_codes), " | cases=", sum(dat$n))
  if (isTRUE(s$add_rain)) {
    cb_rain <- crossbasis(dat$rainfall, lag = s$lag,
                          argvar = list(fun = "ns", df = 3),
                          arglag = arglag, group = dat$code)
    model <- gnm(n ~ cb_t + cb_h + cb_rain + lag_PM25,
                 family = quasipoisson(), eliminate = factor(stratum),
                 data = dat, na.action = na.exclude)
  } else {
    model <- gnm(n ~ cb_t + cb_h + lag_PM25,
                 family = quasipoisson(), eliminate = factor(stratum),
                 data = dat, na.action = na.exclude)
  }

  initial_cen <- median(dat$tmean, na.rm = TRUE)
  red0 <- crossreduce(cb_t, model, type = "overall", cen = initial_cen)
  if (s$reference == USER_NONRESP_REFERENCE_TYPE) {
    cen <- global_nonresp_reference
  } else {
    grid <- seq(
      floor(quantile(dat$tmean, 0.01, na.rm = TRUE)),
      ceiling(quantile(dat$tmean, 0.99, na.rm = TRUE)), by = 0.1
    )
    cen <- grid[which.min(lp_reduced(grid, coef(red0), argvar, initial_cen))]
  }
  red <- crossreduce(cb_t, model, type = "overall", cen = cen)

  model_coef <- coef(model)
  temp_index <- grep("^cb_t", names(model_coef))
  if (length(temp_index) != ncol(cb_t)) {
    stop("Could not identify all temperature cross-basis coefficients in the pooled model.")
  }
  temp_coef <- model_coef[temp_index]
  temp_vcov <- vcov(model)[temp_index, temp_index, drop = FALSE]
  basis_template <- crossbasis(
    rep(cen, max(30L, s$lag + 5L)), lag = s$lag,
    argvar = argvar, arglag = arglag
  )

  crossover <- dat[, .(code, date, n)]
  coverage <- data.table(
    classification = s$classification,
    tag = s$tag,
    lag = s$lag,
    reference = s$reference,
    cen = cen,
    health_cities = length(health_codes),
    total_cases = sum(dat$n),
    observed_climate_cities = sum(health_codes %in% unique(obs$code)),
    climate_imputed_cities = sum(!health_codes %in% unique(obs$code))
  )
  fit <- list(
    disease = s$classification, tag = s$tag,
    lag = s$lag, add_rain = s$add_rain, reference = s$reference,
    basis_template = basis_template,
    temp_coef = temp_coef, temp_vcov = temp_vcov,
    cen = cen, crossover = crossover,
    argvar = argvar, arglag = arglag,
    reduced_coef = coef(red), reduced_vcov = red$vcov,
    coverage = coverage
  )
  saveRDS(fit, fit_file, compress = TRUE)
  fwrite(coverage, file.path(USER_OUT, paste0("erf_pooled_", s$tag, "_coverage.csv")))

  grid_plot <- seq(
    floor(quantile(dat$tmean, 0.01, na.rm = TRUE)),
    ceiling(quantile(dat$tmean, 0.99, na.rm = TRUE)), length.out = 300
  )
  curve <- data.table(
    classification = s$classification,
    temp = grid_plot,
    rr = exp(lp_reduced(grid_plot, coef(red), argvar, cen))
  )
  fwrite(curve, file.path(USER_OUT, paste0("erf_pooled_", s$tag, "_curve.csv")))

  fit_status[[i]] <- data.table(
    classification = s$classification, tag = s$tag,
    status = "fitted", fit_file = fit_file
  )
  rm(dat, counts, cb_t, cb_h, model, red0, red, fit, crossover, curve)
  if (exists("cb_rain")) rm(cb_rain)
  invisible(gc())
}

status <- rbindlist(fit_status, fill = TRUE)
fwrite(status, file.path(USER_OUT, "03.1-pooled-erf-fit-status.csv"))
print(status)
