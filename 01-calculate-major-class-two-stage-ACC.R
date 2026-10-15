###############################################################################
## Calculate major-class anthropogenic attribution using city BLUP curves
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(dlnm)
})

ROOT <- Sys.getenv(
  "ATTR_PROJECT_ROOT",
  getwd()
)
source(file.path(
  ROOT, "two_stage_analysis/00-two-stage-config.R"
))
source(file.path(
  TS_ROOT, "R/two-stage-core.R"
))

REUSE <- tolower(Sys.getenv("TS_REUSE_MAJOR", "true")) == "true"
FORMULA_VERSION <- paste0(
  "city-dlnm-mixmeta-qc300-blup-direct-paired-p80-v6-",
  TS_BASELINE
)

prepared_fits <- lapply(TS_TAGS, function(tag_i) {
  ts_prepare_major_fit(
    readRDS(TS_MAJOR_FIT_FILES[[tag_i]]),
    tag_i
  )
})
names(prepared_fits) <- TS_TAGS
saveRDS(
  prepared_fits,
  file.path(TS_MAJOR_OUTPUT, "two-stage-prepared-major-fits.rds"),
  compress = TRUE
)

calculate_gcm <- function(gcm_i) {
  cache_path <- ts_major_cache_file(gcm_i)
  metadata <- list(
    formula_version = FORMULA_VERSION,
    gcm = gcm_i,
    climate_md5 = acc_file_md5(user_pair_file(gcm_i)),
    fit_md5 = vapply(
      TS_MAJOR_FIT_FILES, acc_file_md5, character(1)
    )
  )
  if (
    REUSE &&
      acc_cache_is_current(cache_path, metadata)
  ) {
    message("Reusing two-stage cache: ", gcm_i)
    return(readRDS(cache_path)$result)
  }

  message("Two-stage attribution: ", gcm_i)
  pair <- readRDS(user_pair_file(gcm_i))
  city_outputs <- list()
  bin_outputs <- list()
  annual_outputs <- list()
  output_index <- 0L

  for (tag_i in TS_TAGS) {
    fit <- prepared_fits[[tag_i]]
    case_weights <- acc_prepare_case_weights(
      fit$dat[, .(code, date, n)],
      dates = TS_DATES,
      baseline = TS_BASELINE
    )
    codes <- intersect(
      sort(unique(case_weights$code)),
      names(fit$coefficient)
    )
    climate <- acc_complete_climate_pair(
      pair = pair,
      codes = codes,
      dates = TS_DATES
    )
    analysis <- merge(
      climate,
      case_weights[code %in% codes],
      by = c("code", "date"),
      all.x = TRUE
    )
    analysis[is.na(n), n := 0]
    data.table::setorder(analysis, code, date)

    daily_rows <- vector("list", length(codes))
    city_rows <- vector("list", length(codes))
    for (j in seq_along(codes)) {
      code_i <- codes[j]
      city <- analysis[code == code_i]
      daily <- ts_city_daily_attribution(
        factual_temperature = city$tf,
        histnat_temperature = city$tn,
        cases = city$n,
        dates = city$date,
        coefficient = fit$coefficient[[code_i]],
        argvar = fit$argvar,
        lag = fit$lag,
        reference = fit$reference[[code_i]],
        tag = tag_i
      )
      if (nrow(daily)) {
        daily[, code := code_i]
        daily_rows[[j]] <- daily
      }
      city_summary <- ts_summarise_city_daily(
        daily, sum(city$n, na.rm = TRUE)
      )
      city_summary[, `:=`(
        code = code_i,
        reference_temperature = fit$reference[[code_i]],
        coefficient_source = fit$coefficient_source[[code_i]]
      )]
      city_rows[[j]] <- city_summary
    }

    daily <- rbindlist(daily_rows, use.names = TRUE, fill = TRUE)
    city <- rbindlist(city_rows, use.names = TRUE, fill = TRUE)
    city[, `:=`(
      tag = tag_i,
      classification = unname(TS_LABELS[[tag_i]]),
      gcm = gcm_i,
      province_code = substr(code, 1L, 2L),
      region = unname(USER_REGION_BY_PROVINCE[
        substr(code, 1L, 2L)
      ]),
      lag = fit$lag,
      formula_version = FORMULA_VERSION
    )]
    output_index <- output_index + 1L
    city_outputs[[output_index]] <- city

    cutpoints <- as.numeric(stats::quantile(
      daily$factual_temperature,
      c(0.025, 0.25, 0.75, 0.975),
      na.rm = TRUE,
      names = FALSE
    ))
    daily[, temperature_range := cut(
      factual_temperature,
      breaks = c(-Inf, cutpoints, Inf),
      labels = c(
        "Extreme cold (<P2.5)",
        "Moderate cold (P2.5-P25)",
        "Moderate (P25-P75)",
        "Moderate heat (P75-P97.5)",
        "Extreme heat (>P97.5)"
      ),
      right = FALSE,
      include.lowest = TRUE
    )]
    daily[, temperature_bin := floor(paired_midpoint) + 0.5]
    bin_outputs[[output_index]] <- daily[, .(
      acc_cases = sum(direct_cases, na.rm = TRUE),
      contribution_af_pp = 100 * sum(
        direct_cases, na.rm = TRUE
      ) / sum(city$total_cases, na.rm = TRUE)
    ), by = .(temperature_range, temperature_bin)][
      , `:=`(tag = tag_i, gcm = gcm_i)
    ]
    daily[, year := as.integer(format(date, "%Y"))]
    annual_outputs[[output_index]] <- daily[, .(
      an_acc_direct = sum(direct_cases, na.rm = TRUE),
      case_weight = sum(case_weight, na.rm = TRUE)
    ), by = year][
      , `:=`(
        af_acc_direct = 100 * an_acc_direct / case_weight,
        tag = tag_i,
        gcm = gcm_i
      )
    ]
  }

  result <- list(
    city = rbindlist(city_outputs, use.names = TRUE, fill = TRUE),
    bins = rbindlist(bin_outputs, use.names = TRUE, fill = TRUE),
    annual = rbindlist(
      annual_outputs, use.names = TRUE, fill = TRUE
    )
  )
  saveRDS(
    list(metadata = metadata, result = result),
    cache_path,
    compress = TRUE
  )
  result
}

all_results <- lapply(TS_GCMS, calculate_gcm)
city_gcm <- rbindlist(lapply(all_results, `[[`, "city"))
bin_gcm <- rbindlist(lapply(all_results, `[[`, "bins"))
annual_gcm <- rbindlist(lapply(all_results, `[[`, "annual"))

national_gcm <- city_gcm[, .(
  total_cases = sum(total_cases, na.rm = TRUE),
  an_acc_direct = sum(an_acc_direct, na.rm = TRUE),
  factual_temperature_cases = sum(
    factual_temperature_cases, na.rm = TRUE
  ),
  histnat_temperature_cases = sum(
    histnat_temperature_cases, na.rm = TRUE
  ),
  relevant_factual_cases = sum(
    relevant_factual_cases, na.rm = TRUE
  ),
  relevant_histnat_cases = sum(
    relevant_histnat_cases, na.rm = TRUE
  )
), by = .(tag, classification, gcm, lag, formula_version)]
national_gcm[, `:=`(
  af_acc_direct = 100 * an_acc_direct / total_cases,
  factual_temperature_af = (
    100 * factual_temperature_cases / total_cases
  ),
  histnat_temperature_af = (
    100 * histnat_temperature_cases / total_cases
  ),
  disease_relevant_ac = (
    100 * (
      relevant_factual_cases - relevant_histnat_cases
    ) / relevant_factual_cases
  )
)]

region_gcm <- city_gcm[
  !is.na(region),
  .(
    total_cases = sum(total_cases, na.rm = TRUE),
    an_acc_direct = sum(an_acc_direct, na.rm = TRUE),
    relevant_factual_cases = sum(
      relevant_factual_cases, na.rm = TRUE
    ),
    relevant_histnat_cases = sum(
      relevant_histnat_cases, na.rm = TRUE
    )
  ),
  by = .(tag, classification, gcm, region)
]
region_gcm[, `:=`(
  af_acc_direct = 100 * an_acc_direct / total_cases,
  disease_relevant_ac = 100 * (
    relevant_factual_cases - relevant_histnat_cases
  ) / relevant_factual_cases
)]

summary_gcm <- acc_summarise_gcm_metrics(
  national_gcm,
  grouping_columns = c("tag", "classification"),
  metric_columns = c(
    "an_acc_direct", "af_acc_direct",
    "factual_temperature_cases",
    "histnat_temperature_cases",
    "factual_temperature_af", "histnat_temperature_af",
    "disease_relevant_ac"
  ),
  gcm_levels = TS_GCMS,
  repetitions = TS_BOOTSTRAP_REPETITIONS,
  seed = TS_BOOTSTRAP_SEED
)

fwrite(
  city_gcm,
  file.path(TS_MAJOR_OUTPUT, "two-stage-city-by-gcm.csv")
)
fwrite(
  national_gcm,
  file.path(TS_MAJOR_OUTPUT, "two-stage-national-by-gcm.csv")
)
fwrite(
  region_gcm,
  file.path(TS_MAJOR_OUTPUT, "two-stage-region-by-gcm.csv")
)
fwrite(
  bin_gcm,
  file.path(TS_MAJOR_OUTPUT, "two-stage-temperature-bins-by-gcm.csv")
)
fwrite(
  annual_gcm,
  file.path(TS_MAJOR_OUTPUT, "two-stage-annual-by-gcm.csv")
)
fwrite(
  summary_gcm,
  file.path(TS_MAJOR_OUTPUT, "two-stage-national-summary.csv")
)

coefficient_audit <- rbindlist(lapply(TS_TAGS, function(tag_i) {
  fit <- prepared_fits[[tag_i]]
  data.table(
    tag = tag_i,
    code = names(fit$coefficient),
    coefficient_source = unname(fit$coefficient_source),
    reference_temperature = unname(fit$reference)
  )
}))
fwrite(
  coefficient_audit,
  file.path(TS_MAJOR_OUTPUT, "two-stage-coefficient-audit.csv")
)

print(summary_gcm[, .(
  classification,
  direct_af = af_acc_direct_estimate,
  direct_low = af_acc_direct_lower,
  direct_high = af_acc_direct_upper,
  relevant_ac = disease_relevant_ac_estimate,
  relevant_low = disease_relevant_ac_lower,
  relevant_high = disease_relevant_ac_upper
)])
