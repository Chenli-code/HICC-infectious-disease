###############################################################################
## Sensitivity analyses: ensemble summaries (Supplementary Fig. S10)
##
## Summarises the GCM-specific sensitivity estimates written by
## V01-10-run-sensitivity-analysis.R.
##
## A  Observed daily versus day-of-year climatological case weights
## B  Alternative lag windows
## C  Population-weighted versus simple grid-weighted temperature
## D  Primary specification (GCM-specific values are in the input file)
##
## All summaries use the direct paired ACC estimand. Point estimates are the
## median of the 10 GCM estimates; intervals are percentiles of 10,000
## bootstrap resamples of the 10 GCM point estimates.
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
})

ROOT <- Sys.getenv("ATTR_PROJECT_ROOT", getwd())
source(file.path(ROOT, "03.0-user-route-config.R"))
source(USER_CORE)

RESULT_DIR <- file.path(ROOT, "Sensitivity", "results")
RESULT_FILE <- file.path(
  RESULT_DIR, "V01-sensitivity-attribution-national-by-gcm.csv"
)
stopifnot(file.exists(RESULT_FILE))

dat <- fread(RESULT_FILE)
required_columns <- c(
  "tag", "classification", "weighting", "baseline",
  "lag_set", "lag", "gcm", "direct_af"
)
if (any(!required_columns %in% names(dat))) {
  stop(
    "Sensitivity results are outdated. Rerun ",
    "V01-10-run-sensitivity-analysis.R."
  )
}

bootstrap_indices <- acc_bootstrap_indices(
  length(USER_GCMS),
  repetitions = 10000L,
  seed = 20260724L
)

summarise_direct_af <- function(data, grouping_columns) {
  data[
    ,
    {
      group_data <- .SD[match(USER_GCMS, gcm)]
      result <- acc_bootstrap_median(
        group_data$direct_af,
        bootstrap_indices
      )
      list(
        estimate = unname(result["estimate"]),
        lower = unname(result["lower"]),
        upper = unname(result["upper"])
      )
    },
    by = grouping_columns,
    .SDcols = c("gcm", "direct_af")
  ]
}

# A. Case-weight baseline (population weighting, main lag).
a_summary <- summarise_direct_af(
  dat[weighting == "population" & lag_set == "Main lag"],
  c("tag", "baseline")
)
a_summary[, `:=`(analysis = "A_case_weights", specification = baseline)]

# B. Lag window (population weighting, observed daily case weights).
b_summary <- summarise_direct_af(
  dat[weighting == "population" & baseline == "observed_daily"],
  c("tag", "lag_set", "lag")
)
b_summary[, `:=`(
  analysis = "B_lag_window",
  specification = paste0(lag_set, " (0-", lag, ")")
)]

# C. Exposure weighting (main lag, observed daily case weights).
c_summary <- summarise_direct_af(
  dat[lag_set == "Main lag" & baseline == "observed_daily"],
  c("tag", "weighting")
)
c_summary[, `:=`(analysis = "C_exposure_weighting", specification = weighting)]

# D. Primary specification.
d_summary <- summarise_direct_af(
  dat[
    weighting == "population" &
      baseline == "observed_daily" &
      lag_set == "Main lag"
  ],
  "tag"
)
d_summary[, `:=`(analysis = "D_primary", specification = "primary")]

summary_all <- rbindlist(
  list(a_summary, b_summary, c_summary, d_summary),
  use.names = TRUE, fill = TRUE
)
setcolorder(summary_all, c(
  "analysis", "tag", "specification", "estimate", "lower", "upper"
))
setorder(summary_all, analysis, tag, specification)

fwrite(
  summary_all,
  file.path(RESULT_DIR, "V01-sensitivity-attribution-national-summary.csv")
)
print(summary_all)
