###############################################################################
## Annual direct paired ACC attribution, 2008-2019
##
## Annual estimates use the same direct factual-versus-hist-nat contrast as
## the main analysis. Intervals resample only the 10 GCM annual point estimates.
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(dlnm)
})

ROOT <- Sys.getenv(
  "ATTR_PROJECT_ROOT",
  getwd()
)
source(file.path(ROOT, "03.0-user-route-config.R"))
source(USER_CORE)

BY_GCM_FILE <- file.path(
  USER_OUT, "04.1-annual-major-class-by-gcm.csv"
)
SUMMARY_FILE <- file.path(
  USER_OUT, "04.1-annual-major-class-components-summary.csv"
)
TREND_GCM_FILE <- file.path(
  USER_OUT, "04.1-annual-major-class-trend-by-gcm.csv"
)
TREND_SUMMARY_FILE <- file.path(
  USER_OUT, "04.1-annual-major-class-trend-summary.csv"
)

fit_files <- setNames(
  file.path(USER_OUT, paste0("erf_pooled_", USER_SPECS$tag, ".rds")),
  USER_SPECS$tag
)
stopifnot(
  all(file.exists(fit_files)),
  all(file.exists(vapply(USER_GCMS, user_pair_file, character(1))))
)

run_annual_gcm <- function(gcm) {
  message("Annual direct ACC: ", gcm)
  pair <- readRDS(user_pair_file(gcm))
  rows <- vector("list", length(fit_files))

  for (ii in seq_along(fit_files)) {
    tag <- names(fit_files)[ii]
    fit <- readRDS(fit_files[[ii]])
    cases <- as.data.table(fit$crossover)[, .(
      code = as.character(code),
      date = as.Date(date),
      n = as.numeric(n)
    )]
    climate <- acc_complete_climate_pair(
      pair = pair,
      codes = sort(unique(cases$code)),
      dates = USER_DATES
    )
    dat <- merge(climate, cases, by = c("code", "date"), all.x = TRUE)
    dat[is.na(n), n := 0]
    setorder(dat, code, date)

    daily <- dat[, {
      result <- acc_direct_daily(
        factual_temperature = tf,
        histnat_temperature = tn,
        cases = n,
        fit = fit,
        dates = date
      )
      result
    }, by = code]
    daily[, year := as.integer(format(as.Date(date), "%Y"))]
    rows[[ii]] <- daily[, .(
      annual_cases = sum(case_weight),
      an_acc = sum(direct_cases)
    ), by = year][, `:=`(
      tag = tag,
      classification = unname(USER_DISEASE_LABELS[tag]),
      gcm = gcm,
      component = "overall",
      formula_version = ACC_FORMULA_VERSION
    )]
    rows[[ii]][, af_acc := 100 * an_acc / annual_cases]
  }
  rbindlist(rows, use.names = TRUE, fill = TRUE)
}

annual_gcm <- rbindlist(
  lapply(USER_GCMS, run_annual_gcm),
  use.names = TRUE,
  fill = TRUE
)
annual_gcm[, gcm_order__ := match(gcm, USER_GCMS)]
setorder(annual_gcm, tag, year, gcm_order__)
annual_gcm[, gcm_order__ := NULL]

bootstrap_indices <- acc_bootstrap_indices(
  length(USER_GCMS),
  repetitions = ACC_BOOTSTRAP_REPETITIONS,
  seed = ACC_BOOTSTRAP_SEED
)

annual_summary <- annual_gcm[, {
  group_data <- .SD[match(USER_GCMS, gcm)]
  af <- acc_bootstrap_median(group_data$af_acc, bootstrap_indices)
  cases <- acc_bootstrap_median(group_data$an_acc, bootstrap_indices)
  list(
    annual_cases = median(group_data$annual_cases),
    af_acc_point = unname(af["estimate"]),
    af_acc_low = unname(af["lower"]),
    af_acc_high = unname(af["upper"]),
    an_acc_point = unname(cases["estimate"]),
    an_acc_low = unname(cases["lower"]),
    an_acc_high = unname(cases["upper"]),
    gcm_count = length(USER_GCMS)
  )
}, by = .(
  tag, classification, component, formula_version, year
)]

trend_gcm <- annual_gcm[, {
  model <- lm(af_acc ~ year)
  list(
    slope_pp_per_year = unname(coef(model)["year"])
  )
}, by = .(tag, classification, component, formula_version, gcm)]

trend_summary <- trend_gcm[, {
  group_data <- .SD[match(USER_GCMS, gcm)]
  slope <- acc_bootstrap_median(
    group_data$slope_pp_per_year,
    bootstrap_indices
  )
  list(
    slope_pp_per_year = unname(slope["estimate"]),
    slope_low = unname(slope["lower"]),
    slope_high = unname(slope["upper"]),
    gcm_count = length(USER_GCMS)
  )
}, by = .(tag, classification, component, formula_version)]

fwrite(annual_gcm, BY_GCM_FILE)
fwrite(annual_summary, SUMMARY_FILE)
fwrite(trend_gcm, TREND_GCM_FILE)
fwrite(trend_summary, TREND_SUMMARY_FILE)

message("Annual direct ACC results written to: ", USER_OUT)

