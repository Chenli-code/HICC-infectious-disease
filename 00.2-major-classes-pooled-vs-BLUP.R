###############################################################################
## 00.2 Major disease classes: pooled national curve vs city-specific BLUP
##
## Purpose
##   Fit a pooled national response curve and city-specific BLUP curves for
##   model diagnostics and the regional heterogeneity panel in Figure 2.
##
## If the optional one-GCM attribution diagnostic is enabled, it uses the same
## direct paired contrast as the final analysis:
##   AF_ACC = 1 - RR_hist-nat / RR_factual.
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(lubridate)
  library(stringr)
  library(dlnm)
  library(splines)
  library(gnm)
  library(mixmeta)
})

invisible(Sys.setlocale("LC_TIME", "C"))

## ---------------------------------------------------------------------------
## 0. Configuration: override with environment variables if needed
## ---------------------------------------------------------------------------
ROOT <- Sys.getenv(
  "ATTR_PROJECT_ROOT",
  getwd()
)

OUT  <- Sys.getenv(
  "ATTR_OUTPUT_DIR",
  file.path(ROOT, "tests/resp_strategies/output")
)

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

source(file.path(ROOT, "00-paths.R"))

PATH_OBS    <- PATH_OBS_WEATHER

DIS          <- Sys.getenv("ATTR_TEST_CLASS", "Respiratory diseases")

CLASS_SPEC <- switch(
  DIS,
  "Respiratory diseases" = list(tag = "respiratory", lag = 7L,
                                add_rain = FALSE, reference = "mrt"),
  "Gastrointestinal and enterovirus" = list(tag = "gastrointestinal", lag = 3L,
                                             add_rain = FALSE, reference = "p50"),
  "Vectorborne" = list(tag = "vectorborne", lag = 14L,
                        add_rain = TRUE, reference = "p50"),
  stop("Unsupported ATTR_TEST_CLASS: ", DIS)
)

CLASS_TAG    <- CLASS_SPEC$tag

START_DATE   <- as.Date("2008-01-01")

END_DATE     <- as.Date("2019-12-31")

LAG          <- as.integer(Sys.getenv("ATTR_TEST_LAG", as.character(CLASS_SPEC$lag)))

ADD_RAIN     <- CLASS_SPEC$add_rain

REFERENCE_MODE <- CLASS_SPEC$reference

LAGNK        <- 2L

MIN_CASES    <- as.integer(Sys.getenv("RESP_TEST_MIN_CASES", "100"))

MAX_CITIES   <- as.integer(Sys.getenv("RESP_TEST_MAX_CITIES", "36"))

FIT_SUFFIX   <- if (MAX_CITIES <= 0L) "allcities" else sprintf("%dcities", MAX_CITIES)

TEST_GCM     <- Sys.getenv("RESP_TEST_GCM", "ACCESS-CM2")

SEED         <- as.integer(Sys.getenv("RESP_TEST_SEED", "20240601"))

RUN_ATTRIB   <- tolower(Sys.getenv("RESP_TEST_RUN_ATTRIB", "false")) == "true"

REUSE_FITS   <- tolower(Sys.getenv("RESP_TEST_REUSE_FITS", "true")) == "true"

MIN_FIT_SUCCESS <- as.numeric(Sys.getenv("ATTR_TEST_MIN_FIT_SUCCESS", "0.80"))

RR_GATE_LOW <- as.numeric(Sys.getenv("ATTR_TEST_RR_GATE_LOW", "0.10"))

RR_GATE_HIGH <- as.numeric(Sys.getenv("ATTR_TEST_RR_GATE_HIGH", "10"))

FIT_FILE     <- if (DIS == "Respiratory diseases") {
  file.path(OUT, sprintf("resp_fit_%s.rds", FIT_SUFFIX))
} else {
  file.path(OUT, sprintf("%s_fit_%s.rds", CLASS_TAG, FIT_SUFFIX))
}

log_file <- file.path(OUT, sprintf("%s_run.log", CLASS_TAG))

zz <- file(log_file, open = "wt")

on.exit({
  try(sink(type = "message"), silent = TRUE)
  try(sink(type = "output"), silent = TRUE)
  try(close(zz), silent = TRUE)
}, add = TRUE)

cat(DIS, "pooled-vs-BLUP pilot\n")

cat("Started:", format(Sys.time()), "\n")

cat("MAX_CITIES:", MAX_CITIES, "| MIN_CASES:", MIN_CASES,
    "| LAG:", LAG, "| ADD_RAIN:", ADD_RAIN,
    "| REFERENCE:", REFERENCE_MODE,
    "| GCM:", TEST_GCM, "| RUN_ATTRIB:", RUN_ATTRIB, "\n\n")

## ---------------------------------------------------------------------------
## Helpers
## ---------------------------------------------------------------------------
pad_code <- function(x) str_pad(as.character(x), width = 6, side = "right", pad = "0")

is_posdef <- function(x) {
  if (any(!is.finite(x))) return(FALSE)
  isTRUE(tryCatch({ chol(x); TRUE }, error = function(e) FALSE))
}

one_basis <- function(x, argvar) {
  do.call(onebasis, c(list(x = x), argvar))
}

lp_from_coef <- function(x, coef, argvar, cen) {
  bx <- one_basis(x, argvar)
  bc <- as.numeric(one_basis(cen, argvar))
  drop(sweep(bx, 2, bc, "-") %*% coef)
}

rr_from_coef <- function(x, coef, argvar, cen) {
  exp(lp_from_coef(x, coef, argvar, cen))
}

find_mrt <- function(x_obs, coef, argvar, probs = seq(0.05, 0.99, by = 0.01)) {
  grid <- unique(as.numeric(quantile(x_obs, probs, na.rm = TRUE)))
  grid[which.min(lp_from_coef(grid, coef, argvar, cen = median(x_obs, na.rm = TRUE)))]
}

select_pilot_cities <- function(meta, n, seed) {
  meta <- copy(meta[total_cases >= MIN_CASES & is.finite(mean_temp) & is.finite(temp_iqr)])
  if (n <= 0L) return(meta$code)
  if (nrow(meta) <= n) return(meta$code)
  meta[, case_group := cut(log1p(total_cases),
                           breaks = unique(quantile(log1p(total_cases), seq(0, 1, length.out = 4), na.rm = TRUE)),
                           include.lowest = TRUE, labels = FALSE)]
  meta[, temp_group := cut(mean_temp,
                           breaks = unique(quantile(mean_temp, seq(0, 1, length.out = 5), na.rm = TRUE)),
                           include.lowest = TRUE, labels = FALSE)]
  meta[is.na(case_group), case_group := 1L]
  meta[is.na(temp_group), temp_group := 1L]
  set.seed(seed)
  meta[, rand := runif(.N)]
  setorder(meta, case_group, temp_group, rand)
  strata <- unique(meta[, .(case_group, temp_group)])
  per_stratum <- max(1L, ceiling(n / nrow(strata)))
  picked <- meta[, head(.SD, per_stratum), by = .(case_group, temp_group)]$code
  if (length(picked) < n) {
    remaining <- setdiff(meta$code, picked)
    picked <- c(picked, head(remaining, n - length(picked)))
  }
  head(picked, n)
}

## ---------------------------------------------------------------------------
## 1. Prepare observed respiratory city-day panel
## ---------------------------------------------------------------------------
if (REUSE_FITS && file.exists(FIT_FILE)) {
  cat("Reusing fitted pilot:", FIT_FILE, "\n")
  fit <- readRDS(FIT_FILE)
  dat <- fit$dat
  city_meta <- fit$city_meta
  selected_codes <- fit$selected_codes
  argvar <- fit$argvar
  arglag <- fit$arglag
  red_pool <- fit$red_pool
  meta_model <- fit$meta_model
  blups <- fit$blups
  raw_coef <- fit$raw_coef
  raw_vcov <- fit$raw_vcov
  cen_common <- fit$cen_common
  city_cen <- fit$city_cen
  fit_status <- fit$fit_status
} else {
  cat("Loading health data...\n")
  load(PATH_HEALTH)
  stopifnot(exists("incidence"))
  inc <- as.data.table(incidence)
  rm(incidence); invisible(gc())
  inc <- inc[classification == DIS & year <= 2019 & code_p != "66" & ename != "COVID-19"]
  inc[, code := pad_code(code_merge)]
  inc[, date := as.Date(date)]
  counts <- inc[, .(n = .N), by = .(code, date)]
  rm(inc); invisible(gc())

  cat("Loading observed climate...\n")
  obs <- as.data.table(readRDS(PATH_OBS))
  obs[, `:=`(code = pad_code(id), date = as.Date(date))]
  obs <- obs[date >= START_DATE & date <= END_DATE,
             .(code, date, tmean = as.numeric(tmean), RH = as.numeric(RH),
               rainfall = as.numeric(rainfall), total_PM25 = as.numeric(total_PM25))]
  setorder(obs, code, date)
  dat_all <- merge(obs, counts, by = c("code", "date"), all.x = TRUE)
  dat_all[is.na(n), n := 0L]
  dat_all[, `:=`(year = year(date), month = month(date), dow = wday(date),
                 mmdd = format(date, "%m-%d"))]
  dat_all[, lag_PM25 := frollmean(total_PM25, n = LAG + 1L,
                                  align = "right", na.rm = TRUE), by = code]

  city_meta_all <- dat_all[, .(total_cases = sum(n, na.rm = TRUE),
                               mean_temp = mean(tmean, na.rm = TRUE),
                               temp_iqr = IQR(tmean, na.rm = TRUE),
                               n_days = .N), by = code]
  selected_codes <- select_pilot_cities(city_meta_all, MAX_CITIES, SEED)
  dat <- dat_all[code %in% selected_codes]
  city_meta <- city_meta_all[code %in% selected_codes]
  setorder(dat, code, date)
  rm(dat_all, obs, counts); invisible(gc())
  fwrite(city_meta, file.path(OUT, sprintf("%s_selected_city_metadata.csv", CLASS_TAG)))
  cat("Selected cities:", length(selected_codes),
      "| cases:", sum(dat$n), "\n")

  ## Same absolute-temperature basis in both strategies: isolates coefficient strategy.
  argvar <- list(fun = "ns", df = 3,
                 knots = as.numeric(quantile(dat$tmean, c(0.10, 0.75, 0.90), na.rm = TRUE)),
                 Bound = range(dat$tmean, na.rm = TRUE))
  arglag <- list(fun = "ns", knots = logknots(LAG, LAGNK), df = 3)
  arglag_h <- list(fun = "ns", knots = logknots(LAG, LAGNK), df = 3)

  ## -------------------------------------------------------------------------
  ## 2A. Pooled national curve
  ## -------------------------------------------------------------------------
  cat("Fitting pooled case-crossover model...\n")
  dat[, stratum_pool := interaction(code, year, month, dow, drop = TRUE)]
  cb_t_pool <- crossbasis(dat$tmean, lag = LAG, argvar = argvar,
                          arglag = arglag, group = dat$code)
  cb_h_pool <- crossbasis(dat$RH, lag = LAG,
                          argvar = list(fun = "ns", df = 3),
                          arglag = arglag_h, group = dat$code)
  if (ADD_RAIN) {
    cb_rain_pool <- crossbasis(dat$rainfall, lag = LAG,
                               argvar = list(fun = "ns", df = 3),
                               arglag = arglag_h, group = dat$code)
    model_pool <- gnm(n ~ cb_t_pool + cb_h_pool + cb_rain_pool + lag_PM25,
                      family = quasipoisson(), eliminate = factor(stratum_pool),
                      data = dat, na.action = na.exclude)
  } else {
    model_pool <- gnm(n ~ cb_t_pool + cb_h_pool + lag_PM25,
                      family = quasipoisson(), eliminate = factor(stratum_pool),
                      data = dat, na.action = na.exclude)
  }
  red_pool0 <- crossreduce(cb_t_pool, model_pool, type = "overall",
                           cen = median(dat$tmean, na.rm = TRUE))
  cen_common <- if (REFERENCE_MODE == "p50") {
    median(dat$tmean, na.rm = TRUE)
  } else {
    find_mrt(dat$tmean, coef(red_pool0), argvar)
  }
  red_pool <- crossreduce(cb_t_pool, model_pool, type = "overall", cen = cen_common)
  cat("Common reference temperature:", round(cen_common, 2), "C\n")
  rm(cb_t_pool, cb_h_pool, model_pool, red_pool0); invisible(gc())
  if (exists("cb_rain_pool")) rm(cb_rain_pool)

  ## -------------------------------------------------------------------------
  ## 2B. City-specific first stage
  ## -------------------------------------------------------------------------
  cat("Fitting city-specific first-stage models...\n")
  raw_coef <- raw_vcov <- list()
  fit_rows <- list()
  for (j in seq_along(selected_codes)) {
    cc <- selected_codes[j]
    d <- copy(dat[code == cc])
    d[, stratum_city := interaction(year, month, dow, drop = TRUE)]
    result <- tryCatch({
      cb_t <- crossbasis(d$tmean, lag = LAG, argvar = argvar, arglag = arglag)
      cb_h <- crossbasis(d$RH, lag = LAG,
                         argvar = list(fun = "ns", df = 3), arglag = arglag_h)
      if (ADD_RAIN) {
        cb_rain <- crossbasis(d$rainfall, lag = LAG,
                              argvar = list(fun = "ns", df = 3), arglag = arglag_h)
        m <- gnm(n ~ cb_t + cb_h + cb_rain + lag_PM25,
                 family = quasipoisson(), eliminate = factor(stratum_city),
                 data = d, na.action = na.exclude)
      } else {
        m <- gnm(n ~ cb_t + cb_h + lag_PM25,
                 family = quasipoisson(), eliminate = factor(stratum_city),
                 data = d, na.action = na.exclude)
      }
      red <- crossreduce(cb_t, m, type = "overall", cen = cen_common)
      v <- vcov(red); b <- coef(red)
      if (any(!is.finite(b)) || !is_posdef(v)) stop("non-finite coefficient or non-PD vcov")
      list(ok = TRUE, coef = b, vcov = v, message = NA_character_)
    }, error = function(e) list(ok = FALSE, message = conditionMessage(e)))
    if (isTRUE(result$ok)) {
      raw_coef[[cc]] <- result$coef
      raw_vcov[[cc]] <- result$vcov
    }
    fit_rows[[j]] <- data.table(code = cc, success = isTRUE(result$ok),
                                message = result$message)
    if (j %% 6 == 0 || j == length(selected_codes))
      cat(sprintf("  city %d/%d | successful %d\n", j, length(selected_codes), length(raw_coef)))
  }
  fit_status <- rbindlist(fit_rows, fill = TRUE)
  fwrite(fit_status, file.path(OUT, sprintf("%s_city_fit_status.csv", CLASS_TAG)))
  if (length(raw_coef) < 8) stop("Fewer than 8 city models succeeded; BLUP pilot is not estimable.")

  ## -------------------------------------------------------------------------
  ## 3. Multivariate meta-regression and city BLUPs
  ## -------------------------------------------------------------------------
  good_codes <- names(raw_coef)
  ymat <- do.call(rbind, raw_coef)
  rownames(ymat) <- good_codes
  Slist <- raw_vcov[good_codes]
  meta_dat <- city_meta[match(good_codes, code)]
  meta_dat[, `:=`(mean_temp_z = as.numeric(scale(mean_temp)),
                  temp_iqr_z = as.numeric(scale(temp_iqr)))]
  meta_model <- tryCatch(
    mixmeta(ymat ~ mean_temp_z + temp_iqr_z, Slist,
            data = meta_dat, method = "reml",
            control = list(showiter = FALSE, igls.inititer = 10)),
    error = function(e) {
      message("Meta-regression with climate predictors failed; using intercept-only REML: ", e$message)
      mixmeta(ymat ~ 1, Slist, data = meta_dat, method = "reml")
    }
  )
  blups <- blup(meta_model, vcov = TRUE)
  if (is.null(names(blups))) names(blups) <- good_codes
  if (!all(good_codes %in% names(blups))) names(blups) <- good_codes

  city_cen <- setNames(rep(NA_real_, length(good_codes)), good_codes)
  for (cc in good_codes) {
    b <- blups[[cc]]$blup
    city_cen[cc] <- find_mrt(dat[code == cc, tmean], b, argvar)
  }

  fit <- list(dat = dat, city_meta = city_meta, selected_codes = selected_codes,
              argvar = argvar, arglag = arglag, red_pool = red_pool,
              meta_model = meta_model, blups = blups,
              raw_coef = raw_coef, raw_vcov = raw_vcov,
              cen_common = cen_common, city_cen = city_cen,
              fit_status = fit_status,
              config = list(dis = DIS, lag = LAG, min_cases = MIN_CASES,
                            max_cities = MAX_CITIES, seed = SEED))
  saveRDS(fit, FIT_FILE)
  capture.output(summary(meta_model),
                 file = file.path(OUT, sprintf("%s_mixmeta_summary.txt", CLASS_TAG)))
  cat("Saved fitted pilot:", FIT_FILE, "\n")
}

## ---------------------------------------------------------------------------
## 4. Curve and shrinkage diagnostics
## ---------------------------------------------------------------------------
good_codes <- intersect(names(raw_coef), names(blups))

temp_grid <- seq(as.numeric(quantile(dat$tmean, 0.01, na.rm = TRUE)),
                 as.numeric(quantile(dat$tmean, 0.99, na.rm = TRUE)), length.out = 200)

curve_rows <- list(data.table(code = "Pooled", strategy = "Pooled national",
                              temp = temp_grid,
                              rr = rr_from_coef(temp_grid, coef(red_pool), argvar, cen_common)))

for (cc in good_codes) {
  curve_rows[[length(curve_rows) + 1L]] <- data.table(
    code = cc, strategy = "City BLUP", temp = temp_grid,
    rr = rr_from_coef(temp_grid, blups[[cc]]$blup, argvar, cen_common))
}

curve_dt <- rbindlist(curve_rows)

fwrite(curve_dt, file.path(OUT, sprintf("%s_curve_comparison.csv", CLASS_TAG)))

rr_diag <- rbindlist(lapply(good_codes, function(cc) {
  t95 <- quantile(dat[code == cc, tmean], 0.95, na.rm = TRUE)
  data.table(code = cc, temp_p95 = as.numeric(t95),
             rr_raw_p95 = rr_from_coef(t95, raw_coef[[cc]], argvar, cen_common),
             rr_blup_p95 = rr_from_coef(t95, blups[[cc]]$blup, argvar, cen_common),
             city_mrt = city_cen[[cc]])
}))

rr_diag <- merge(rr_diag, city_meta, by = "code", all.x = TRUE)

fwrite(rr_diag, file.path(OUT, sprintf("%s_blup_shrinkage_diagnostics.csv", CLASS_TAG)))

fit_success_fraction <- mean(fit_status$success, na.rm = TRUE)

rr_flag <- !is.finite(rr_diag$rr_blup_p95) |
  rr_diag$rr_blup_p95 < RR_GATE_LOW | rr_diag$rr_blup_p95 > RR_GATE_HIGH

rr_flag_fraction <- mean(rr_flag)

model_qc <- data.table(
  classification = DIS,
  attempted_cities = nrow(fit_status),
  successful_cities = sum(fit_status$success, na.rm = TRUE),
  fit_success_fraction = fit_success_fraction,
  p95_blup_rr_min = min(rr_diag$rr_blup_p95, na.rm = TRUE),
  p95_blup_rr_max = max(rr_diag$rr_blup_p95, na.rm = TRUE),
  p95_blup_rr_flag_fraction = rr_flag_fraction,
  passed = fit_success_fraction >= MIN_FIT_SUCCESS && rr_flag_fraction <= 0.10
)

fwrite(model_qc, file.path(OUT, sprintf("%s_model_qc.csv", CLASS_TAG)))

cat("\nFinished:", format(Sys.time()), "\n")

cat("Outputs:", OUT, "\n")
