###############################################################################
## Figure 3 | Infectious disease burden attributable to human-induced climate
## change
##
## Row 1  A observed factual cases versus cases predicted under the natural-only
##          counterfactual climate;
##        B direct paired ACC AF among all observed cases
## Row 2  C-E direct ACC AF contributions across temperature
## Row 3  F temperature-range decomposition; G direct ACC cases
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(dlnm)
  library(scales)
  library(parallel)
})

ROOT <- Sys.getenv(
  "ATTR_PROJECT_ROOT",
  getwd()
)

source(file.path(ROOT, "03.0-user-route-config.R"))

source(USER_CORE)

# ------------------------------- paths ---------------------------------------
FIGURE3_DIR <- file.path(USER_ROOT, "Figure3")

dir.create(FIGURE3_DIR, recursive = TRUE, showWarnings = FALSE)

DIRECT_GCM_FILE <- file.path(
  USER_OUT,
  "direct_counterfactual_af",
  "V01-direct-ACC-AF-national-by-gcm.csv"
)

REFERENCE_AC_SUMMARY_FILE <- file.path(
  USER_OUT,
  "warm4_attribution",
  "V01-warm4-temperature-burden-national-summary.csv"
)

BIN_GCM_FILE <- file.path(
  FIGURE3_DIR,
  "V01-Figure3-direct-temperature-range-contributions-by-gcm.csv"
)

CUM_GCM_FILE <- file.path(
  FIGURE3_DIR,
  "V01-Figure3-direct-cumulative-burden-by-gcm.csv"
)

DETAIL_GCM_FILE <- file.path(
  FIGURE3_DIR,
  "V01-Figure3-direct-one-degree-contributions-by-gcm.csv"
)

CHECK_FILE <- file.path(
  FIGURE3_DIR,
  "V01-Figure3-direct-daily-decomposition-check.csv"
)

stopifnot(
  file.exists(DIRECT_GCM_FILE),
  file.exists(REFERENCE_AC_SUMMARY_FILE)
)

REUSE_DECOMP <- tolower(
  Sys.getenv("FIGURE3_REUSE_DECOMP", "true")
) == "true"

NCORES <- as.integer(Sys.getenv("FIGURE3_NCORES", "3"))

NCORES <- max(1L, min(NCORES, 3L))

disease_order <- c("respiratory", "gastrointestinal", "vectorborne")

disease_labels <- c(
  respiratory = "Respiratory infectious diseases",
  gastrointestinal = "Gastrointestinal and enterovirus diseases",
  vectorborne = "Vector-borne diseases"
)

disease_short <- c(
  respiratory = "Respiratory",
  gastrointestinal = "Gastrointestinal /\nenterovirus",
  vectorborne = "Vector-borne"
)

temperature_order <- c(
  "<P2.5", "P2.5–P25", "P25–P75", "P75–P97.5", ">P97.5"
)

signed_number <- function(x, accuracy = 0.1) {
  paste0(
    ifelse(x > 0, "+", ""),
    number(x, accuracy = accuracy)
  )
}

GCM_BOOT_N <- 10000L

GCM_BOOT_SEED <- 20260723L

DECOMP_VERSION <- paste0(ACC_FORMULA_VERSION, "-daily-v3")

gcm_bootstrap_indices <- acc_bootstrap_indices(
  length(USER_GCMS),
  repetitions = GCM_BOOT_N,
  seed = GCM_BOOT_SEED
)

gcm_median_bootstrap <- function(
  x
) {
  result <- acc_bootstrap_median(x, gcm_bootstrap_indices)
  list(
    estimate = unname(result["estimate"]),
    low = unname(result["lower"]),
    high = unname(result["upper"]),
    gcm_min = unname(result["gcm_min"]),
    gcm_max = unname(result["gcm_max"]),
    gcm_count = length(x)
  )
}

run_decomposition_task <- function(task) {
  suppressPackageStartupMessages({
    library(data.table)
    library(dlnm)
  })

  tag_i <- task$tag
  gcm_i <- task$gcm
  fit <- readRDS(file.path(
    USER_OUT, paste0("erf_pooled_", tag_i, ".rds")
  ))
  cases <- as.data.table(fit$crossover)
  cases[, `:=`(
    code = as.character(code),
    date = as.Date(date),
    n = as.numeric(n)
  )]
  case_weights <- acc_prepare_case_weights(
    cases,
    dates = USER_DATES,
    baseline = "observed_daily"
  )
  codes <- sort(unique(cases$code))
  pair <- readRDS(user_pair_file(gcm_i))
  clim <- acc_complete_climate_pair(
    pair = pair,
    codes = codes,
    dates = USER_DATES
  )
  dat <- merge(
    clim, case_weights,
    by = c("code", "date"), all.x = TRUE
  )
  dat[is.na(n), n := 0]
  setorder(dat, code, date)
  daily <- dat[, {
    contribution <- acc_direct_daily(
      factual_temperature = tf,
      histnat_temperature = tn,
      cases = n,
      fit = fit,
      dates = date
    )
    if (!nrow(contribution)) {
      NULL
    } else {
      index <- match(as.Date(contribution$date), as.Date(date))
      data.table(
        tf = tf[index],
        tn = tn[index],
        acc_daily = contribution$direct_cases
      )
    }
  }, by = code]

  cutpoints <- as.numeric(quantile(
    daily$tf, probs = c(0.025, 0.25, 0.75, 0.975),
    na.rm = TRUE, names = FALSE
  ))
  daily[, temperature_range := cut(
    tf,
    breaks = c(-Inf, cutpoints, Inf),
    labels = temperature_order,
    right = FALSE, include.lowest = TRUE
  )]
  bin <- daily[, .(
    acc_cases = sum(acc_daily, na.rm = TRUE)
  ), by = temperature_range]
  bin[, `:=`(
    tag = tag_i,
    gcm = gcm_i,
    formula_version = DECOMP_VERSION,
    total_acc_cases = sum(daily$acc_daily, na.rm = TRUE)
  )]

  existing_row <- fread(DIRECT_GCM_FILE)[
    tag == tag_i & gcm == gcm_i
  ][1L]
  included_cases_i <- existing_row$total_cases

  daily[, paired_midpoint := (tf + tn) / 2]
  daily[, temperature_bin := floor(paired_midpoint) + 0.5]
  detail <- daily[, .(
    acc_cases = sum(acc_daily, na.rm = TRUE)
  ), by = temperature_bin]
  detail[, `:=`(
    tag = tag_i,
    gcm = gcm_i,
    formula_version = DECOMP_VERSION,
    contribution_af_pp = 100 * acc_cases / included_cases_i
  )]

  setorder(daily, tf)
  daily[, percentile := 100 * seq_len(.N) / .N]
  daily[, percentile_bin := pmin(
    100L, pmax(1L, as.integer(ceiling(percentile)))
  )]
  cum <- daily[, .(
    acc_cases = sum(acc_daily, na.rm = TRUE),
    temperature = max(tf, na.rm = TRUE)
  ), by = percentile_bin]
  cum <- merge(
    data.table(percentile_bin = 1:100),
    cum, by = "percentile_bin", all.x = TRUE
  )
  cum[is.na(acc_cases), acc_cases := 0]
  if (anyNA(cum$temperature)) {
    ok <- which(is.finite(cum$temperature))
    cum[, temperature := approx(
      x = percentile_bin[ok],
      y = temperature[ok],
      xout = percentile_bin,
      rule = 2
    )$y]
  }
  setorder(cum, percentile_bin)
  cum[, cumulative_cases := cumsum(acc_cases)]
  cum <- rbindlist(list(
    data.table(
      percentile_bin = 0L, acc_cases = 0,
      cumulative_cases = 0,
      temperature = min(daily$tf, na.rm = TRUE)
    ),
    cum
  ), use.names = TRUE, fill = TRUE)
  total_acc <- sum(daily$acc_daily, na.rm = TRUE)
  cum[, `:=`(
    tag = tag_i,
    gcm = gcm_i,
    formula_version = DECOMP_VERSION,
    cumulative_percent = 100 * cumulative_cases /
      abs(total_acc),
    total_acc_cases = total_acc
  )]

  existing <- existing_row$an_acc_direct
  check <- data.table(
    tag = tag_i, gcm = gcm_i,
    formula_version = DECOMP_VERSION,
    daily_sum = total_acc,
    existing_total = existing,
    absolute_difference = total_acc - existing,
    relative_difference = (total_acc - existing) /
      max(1, abs(existing))
  )
  list(
    bin = bin, cumulative = cum,
    detail = detail, check = check
  )
}

CUM_HAS_TEMPERATURE <- file.exists(CUM_GCM_FILE) &&
  "temperature" %in% names(fread(CUM_GCM_FILE, nrows = 0))

DETAIL_EXISTS <- file.exists(DETAIL_GCM_FILE)

DECOMP_IS_CURRENT <- FALSE

if (all(file.exists(c(
  BIN_GCM_FILE, CUM_GCM_FILE, DETAIL_GCM_FILE, CHECK_FILE
)))) {
  decomp_headers <- lapply(
    c(BIN_GCM_FILE, CUM_GCM_FILE, DETAIL_GCM_FILE, CHECK_FILE),
    function(path) fread(path, nrows = 1L)
  )
  DECOMP_IS_CURRENT <- all(vapply(
    decomp_headers,
    function(x) {
      "formula_version" %in% names(x) &&
        identical(
          unique(as.character(x$formula_version)),
          DECOMP_VERSION
        )
    },
    logical(1)
  ))
}

if (
  !REUSE_DECOMP ||
    !all(file.exists(c(BIN_GCM_FILE, CUM_GCM_FILE, CHECK_FILE))) ||
    !CUM_HAS_TEMPERATURE ||
    !DETAIL_EXISTS ||
    !DECOMP_IS_CURRENT
) {
  tasks <- lapply(disease_order, function(tag_i) {
    lapply(USER_GCMS, function(gcm_i) {
      list(tag = tag_i, gcm = gcm_i)
    })
  })
  tasks <- unlist(tasks, recursive = FALSE)

  message(
    "Computing additive daily decomposition | tasks=",
    length(tasks), " | workers=", NCORES
  )
  if (NCORES == 1L) {
    result <- lapply(tasks, run_decomposition_task)
  } else {
    cl <- makeCluster(NCORES)
    clusterEvalQ(cl, {
      suppressPackageStartupMessages(library(data.table))
      suppressPackageStartupMessages(library(dlnm))
      NULL
    })
    clusterExport(
      cl,
      c(
        "USER_OUT", "USER_CACHE", "USER_DATES", "DIRECT_GCM_FILE",
        "ACC_FORMULA_VERSION", "DECOMP_VERSION",
        "temperature_order",
        "user_safe_gcm", "user_pair_file",
        "acc_complete_climate_pair",
        "acc_prepare_case_weights",
        "acc_forward_weights", "acc_one_basis",
        "acc_direct_daily",
        "run_decomposition_task"
      ),
      envir = environment()
    )
    result <- tryCatch(
      parLapplyLB(cl, tasks, run_decomposition_task),
      finally = stopCluster(cl)
    )
  }

  bin_gcm <- rbindlist(lapply(result, `[[`, "bin"))
  cum_gcm <- rbindlist(lapply(result, `[[`, "cumulative"))
  detail_gcm <- rbindlist(lapply(result, `[[`, "detail"))
  check <- rbindlist(lapply(result, `[[`, "check"))
  fwrite(bin_gcm, BIN_GCM_FILE)
  fwrite(cum_gcm, CUM_GCM_FILE)
  fwrite(detail_gcm, DETAIL_GCM_FILE)
  fwrite(check, CHECK_FILE)
} else {
  bin_gcm <- fread(BIN_GCM_FILE)
  cum_gcm <- fread(CUM_GCM_FILE)
  detail_gcm <- fread(DETAIL_GCM_FILE)
  check <- fread(CHECK_FILE)
}

if (max(abs(check$relative_difference), na.rm = TRUE) > 1e-8) {
  stop("Daily decomposition does not reproduce existing ACC totals.")
}

# ------------------------------- analysis data -------------------------------
direct_gcm <- fread(DIRECT_GCM_FILE)

reference_ac_summary <- fread(
  REFERENCE_AC_SUMMARY_FILE
)

required_relevant_columns <- c(
  "tag", "ac_primary",
  "ac_primary_low", "ac_primary_high",
  "primary_burden_definition"
)

if (
  any(
    !required_relevant_columns %in%
      names(reference_ac_summary)
  )
) {
  stop(
    "Warm-season attribution summary is outdated. Rerun ",
    "V01-13-hottest-four-months-AF.R."
  )
}

required_direct_columns <- c(
  "tag", "classification", "gcm", "total_cases",
  "counterfactual_cases", "an_acc_direct", "af_acc_direct"
)

if (any(!required_direct_columns %in% names(direct_gcm))) {
  stop(
    "The direct-ACC file is outdated. Rerun ",
    "V01-07-calculate-direct-counterfactual-AF.R."
  )
}

set.seed(GCM_BOOT_SEED + 1L)

ci_cases <- direct_gcm[, {
  boot_factual <- gcm_median_bootstrap(total_cases)
  boot_counterfactual <- gcm_median_bootstrap(counterfactual_cases)
  list(
    factual_cases = boot_factual$estimate,
    factual_cases_low = boot_factual$low,
    factual_cases_high = boot_factual$high,
    counterfactual_cases = boot_counterfactual$estimate,
    counterfactual_cases_low = boot_counterfactual$low,
    counterfactual_cases_high = boot_counterfactual$high
  )
}, by = .(classification, tag)]

set.seed(GCM_BOOT_SEED + 2L)

direct_ci <- direct_gcm[, {
  boot <- gcm_median_bootstrap(
    an_acc_direct
  )
  total_i <- median(total_cases, na.rm = TRUE)
  list(
    total_cases = total_i,
    an_acc = boot$estimate,
    an_acc_low = boot$low,
    an_acc_high = boot$high,
    an_acc_gcm_min = boot$gcm_min,
    an_acc_gcm_max = boot$gcm_max,
    af_acc = 100 * boot$estimate / total_i,
    af_acc_low = 100 * boot$low / total_i,
    af_acc_high = 100 * boot$high / total_i,
    af_acc_gcm_min = 100 * boot$gcm_min / total_i,
    af_acc_gcm_max = 100 * boot$gcm_max / total_i,
    gcm_count = boot$gcm_count
  )
}, by = tag]

ci <- merge(
  ci_cases,
  direct_ci,
  by = "tag",
  all.x = TRUE,
  sort = FALSE
)

ci[, `:=`(
  disease = factor(
    unname(disease_short[tag]),
    levels = rev(unname(disease_short[disease_order]))
  ),
  disease_long = factor(
    unname(disease_labels[tag]),
    levels = rev(unname(disease_labels[disease_order]))
  )
)]

bin_gcm[, temperature_range := factor(
  temperature_range, levels = temperature_order
)]

set.seed(GCM_BOOT_SEED + 3L)

bin_point <- bin_gcm[, {
  boot <- gcm_median_bootstrap(acc_cases)
  list(
    acc_cases = boot$estimate,
    ci_low = boot$low,
    ci_high = boot$high,
    gcm_min = boot$gcm_min,
    gcm_max = boot$gcm_max
  )
}, by = .(tag, temperature_range)]

bin_point[, net_cases := sum(acc_cases), by = tag]

bin_point[, absolute_contribution_percent := (
  100 * abs(acc_cases) /
    sum(abs(acc_cases), na.rm = TRUE)
), by = tag]

bin_point[, disease := factor(
  unname(disease_short[tag]),
  levels = rev(unname(disease_short[disease_order]))
)]

bin_point[, segment_label := fifelse(
  absolute_contribution_percent >= 5,
  paste0(
    number(
      absolute_contribution_percent,
      accuracy = 1
    ),
    "%"
  ),
  ""
)]

cum_gcm[, disease := factor(
  unname(disease_labels[tag]),
  levels = unname(disease_labels[disease_order])
)]

cum_summary <- cum_gcm[, .(
  temperature = median(temperature, na.rm = TRUE),
  cumulative_cases = median(
    cumulative_cases / 1000, na.rm = TRUE
  )
), by = .(tag, disease, percentile_bin)]

detail_grid <- rbindlist(lapply(disease_order, function(tag_i) {
  bins_i <- seq(
    min(
      detail_gcm[tag == tag_i, temperature_bin],
      na.rm = TRUE
    ),
    max(
      detail_gcm[tag == tag_i, temperature_bin],
      na.rm = TRUE
    ),
    by = 1
  )
  CJ(
    tag = tag_i,
    gcm = USER_GCMS,
    temperature_bin = bins_i,
    unique = TRUE
  )
}))

detail_full <- merge(
  detail_grid,
  detail_gcm[, .(
    tag, gcm, temperature_bin, contribution_af_pp
  )],
  by = c("tag", "gcm", "temperature_bin"),
  all.x = TRUE
)

detail_full[is.na(contribution_af_pp), contribution_af_pp := 0]

set.seed(GCM_BOOT_SEED + 4L)

detail_summary <- detail_full[, {
  boot <- gcm_median_bootstrap(
    contribution_af_pp
  )
  list(
    estimate = boot$estimate,
    ci_low = boot$low,
    ci_high = boot$high,
    gcm_min = boot$gcm_min,
    gcm_max = boot$gcm_max
  )
}, by = .(tag, temperature_bin)]

detail_summary[, direction := fifelse(
  estimate >= 0, "Increase", "Reduction"
)]

# -------------------------------- panel B ------------------------------------
b_direct <- ci[, .(
  tag,
  disease,
  estimate = af_acc,
  low = af_acc_low,
  high = af_acc_high,
  metric = "Direct ACC AF\namong all observed cases"
)]

b_reference <- reference_ac_summary[, .(
  tag,
  estimate = ac_primary,
  low = ac_primary_low,
  high = ac_primary_high,
  metric = paste0(
    "Disease-relevant AC\n",
    "(MMT total or warm-season burden)"
  )
)]

b_reference[, disease := factor(
  unname(disease_short[tag]),
  levels = rev(unname(disease_short[disease_order]))
)]

b_plot <- rbindlist(
  list(b_direct, b_reference),
  use.names = TRUE,
  fill = TRUE
)

b_plot[, metric := factor(
  metric,
  levels = c(
    "Direct ACC AF\namong all observed cases",
    paste0(
      "Disease-relevant AC\n",
      "(MMT total or warm-season burden)"
    )
  )
)]

b_plot[, y_plot := as.numeric(disease)]

b_plot[, value_label := paste0(
  signed_number(estimate, 0.1), "% (",
  signed_number(low, 0.1), " to ",
  signed_number(high, 0.1), "%)"
)]

bin_point[, disease_c := factor(
  as.character(disease),
  levels = unname(disease_short[disease_order])
)]

fwrite(
  bin_point,
  file.path(
    FIGURE3_DIR,
    "V01-Figure3-temperature-range-contributions-summary.csv"
  )
)

fwrite(
  cum_summary,
  file.path(
    FIGURE3_DIR,
    "V01-Figure3-cumulative-burden-summary.csv"
  )
)

fwrite(
  detail_summary,
  file.path(
    FIGURE3_DIR,
    "V01-Figure3-one-degree-contributions-summary.csv"
  )
)

fwrite(
  b_plot,
  file.path(
    FIGURE3_DIR,
    "V01-Figure3-headline-relative-burden-change.csv"
  )
)
