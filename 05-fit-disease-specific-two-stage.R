###############################################################################
## Disease-specific first-stage city DLNMs and second-stage mixmeta/BLUP models
##
## Each disease is fitted independently.  The lag window follows the published
## disease-specific temperature analysis.  Cities are eligible at >=300 cases;
## if this yields fewer than eight cities, the threshold is relaxed to >=100.
## Diseases with fewer than eight eligible or quality-controlled city models
## are reported as not estimable and are not forced into Figure 5.
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(lubridate)
  library(dlnm)
  library(gnm)
  library(mixmeta)
})
invisible(Sys.setlocale("LC_TIME", "C"))

ROOT <- Sys.getenv(
  "ATTR_PROJECT_ROOT",
  getwd()
)
source(file.path(ROOT, "two_stage_analysis/00-two-stage-config.R"))
source(file.path(TS_ROOT, "R/two-stage-core.R"))
source(file.path(ROOT, "10.0-Figure3-disease-specification.R"))

FIT_DIR <- file.path(TS_DISEASE_OUTPUT, "fits")
dir.create(FIT_DIR, recursive = TRUE, showWarnings = FALSE)
STATUS_FILE <- file.path(
  TS_DISEASE_OUTPUT, "two-stage-disease-fit-status.csv"
)
REUSE <- tolower(Sys.getenv("TS_REUSE_DISEASE_FITS", "true")) == "true"
FILTER <- trimws(Sys.getenv("TS_DISEASE_FILTER", ""))

specs <- copy(FIG3_DISEASE_SPECS)
if (nzchar(FILTER)) {
  wanted <- trimws(strsplit(FILTER, "\\|")[[1L]])
  specs <- specs[
    disease %chin% wanted | disease_tag %chin% wanted |
      class_tag %chin% wanted
  ]
  if (!nrow(specs)) stop("TS_DISEASE_FILTER did not match a disease.")
}

raw_fit_file <- function(tag) {
  file.path(FIT_DIR, paste0("raw-city-DLNM-", tag, ".rds"))
}
prepared_fit_file <- function(tag) {
  file.path(FIT_DIR, paste0("prepared-mixmeta-BLUP-", tag, ".rds"))
}
is_posdef <- function(x) {
  if (any(!is.finite(x))) return(FALSE)
  isTRUE(tryCatch({
    chol(x)
    TRUE
  }, error = function(e) FALSE))
}

message("Loading observed weather for city-specific disease models...")
weather <- as.data.table(readRDS(USER_PATH_OBS))
weather[, `:=`(
  code = user_pad_code(id),
  date = as.Date(date)
)]
weather <- weather[
  date >= USER_START_DATE & date <= USER_END_DATE,
  .(
    code, date,
    tmean = as.numeric(tmean),
    RH = as.numeric(RH),
    rainfall = as.numeric(rainfall),
    total_PM25 = as.numeric(total_PM25)
  )
]
setorder(weather, code, date)
weather_vars <- c("tmean", "RH", "rainfall", "total_PM25")
weather[, province_code := substr(code, 1L, 2L)]
province_weather <- weather[, lapply(.SD, mean, na.rm = TRUE),
  by = .(province_code, date), .SDcols = weather_vars
]
setnames(
  province_weather, weather_vars,
  paste0(weather_vars, "_province")
)
national_weather <- weather[, lapply(.SD, mean, na.rm = TRUE),
  by = date, .SDcols = weather_vars
]
setnames(
  national_weather, weather_vars,
  paste0(weather_vars, "_national")
)

complete_weather <- function(codes) {
  panel <- CJ(
    code = sort(unique(as.character(codes))),
    date = TS_DATES
  )
  panel[, province_code := substr(code, 1L, 2L)]
  panel <- merge(
    panel,
    weather[, c("code", "date", weather_vars), with = FALSE],
    by = c("code", "date"),
    all.x = TRUE
  )
  panel <- merge(
    panel,
    province_weather,
    by = c("province_code", "date"),
    all.x = TRUE
  )
  panel <- merge(
    panel,
    national_weather,
    by = "date",
    all.x = TRUE
  )
  for (variable in weather_vars) {
    panel[, (variable) := fcoalesce(
      get(variable),
      get(paste0(variable, "_province")),
      get(paste0(variable, "_national"))
    )]
  }
  panel[, c(
    "province_code",
    paste0(weather_vars, "_province"),
    paste0(weather_vars, "_national")
  ) := NULL]
  setorder(panel, code, date)
  panel
}

status_rows <- vector("list", nrow(specs))

for (ii in seq_len(nrow(specs))) {
  spec <- specs[ii]
  pooled_path <- figure3_fit_file(spec$disease_tag)
  raw_path <- raw_fit_file(spec$disease_tag)
  prepared_path <- prepared_fit_file(spec$disease_tag)
  message(sprintf(
    "[%02d/%02d] %s", ii, nrow(specs), spec$disease
  ))

  if (REUSE && file.exists(prepared_path)) {
    prepared <- tryCatch(readRDS(prepared_path), error = function(e) NULL)
    if (
      !is.null(prepared) &&
        identical(as.integer(prepared$lag), as.integer(spec$lag))
    ) {
      status_rows[[ii]] <- data.table(
        disease = spec$disease,
        disease_tag = spec$disease_tag,
        classification = spec$classification,
        class_tag = spec$class_tag,
        lag = spec$lag,
        minimum_city_cases = prepared$minimum_city_cases,
        eligible_cities = length(prepared$coefficient),
        successful_city_models = sum(
          prepared$fit_status$success, na.rm = TRUE
        ),
        retained_first_stage = sum(
          prepared$coefficient_qc$retained_first_stage,
          na.rm = TRUE
        ),
        status = "reused",
        message = "",
        prepared_fit = prepared_path
      )
      fwrite(rbindlist(status_rows, fill = TRUE), STATUS_FILE)
      next
    }
  }

  result <- tryCatch({
    if (!file.exists(pooled_path)) {
      stop("No disease-specific pooled model/count file.")
    }
    pooled_fit <- readRDS(pooled_path)
    counts <- as.data.table(pooled_fit$sparse_counts)[, .(
      code = as.character(code),
      date = as.Date(date),
      n = as.numeric(n)
    )]
    city_totals <- counts[, .(total_cases = sum(n)), by = code]
    minimum_city_cases <- if (
      city_totals[total_cases >= 300, .N] >= 8L
    ) {
      300L
    } else {
      100L
    }
    selected_codes <- city_totals[
      total_cases >= minimum_city_cases, code
    ]
    if (length(selected_codes) < 8L) {
      stop("Fewer than eight cities have at least 100 cases.")
    }

    dat <- complete_weather(selected_codes)
    dat <- merge(
      dat,
      counts[code %chin% selected_codes],
      by = c("code", "date"),
      all.x = TRUE
    )
    dat[is.na(n), n := 0]
    setorder(dat, code, date)
    dat[, `:=`(
      year = year(date),
      month = month(date),
      dow = wday(date)
    )]
    dat[, lag_PM25 := frollmean(
      total_PM25,
      n = as.integer(spec$lag) + 1L,
      align = "right",
      na.rm = TRUE
    ), by = code]

    city_meta <- dat[, .(
      total_cases = sum(n, na.rm = TRUE),
      mean_temp = mean(tmean, na.rm = TRUE),
      temp_iqr = IQR(tmean, na.rm = TRUE),
      n_days = .N
    ), by = code]

    argvar <- list(
      fun = "ns",
      df = 3,
      knots = as.numeric(quantile(
        dat$tmean, c(0.10, 0.75, 0.90), na.rm = TRUE
      )),
      Bound = range(dat$tmean, na.rm = TRUE)
    )
    arglag <- list(
      fun = "ns",
      knots = logknots(
        as.integer(spec$lag), as.integer(spec$lagnk)
      ),
      df = 3
    )
    cen_common <- if (
      spec$class_tag == "respiratory"
    ) {
      as.numeric(pooled_fit$cen)
    } else {
      as.numeric(quantile(
        dat$tmean,
        TS_NONRESP_REFERENCE_PROBABILITY,
        na.rm = TRUE
      ))
    }

    raw_coef <- list()
    raw_vcov <- list()
    fit_rows <- vector("list", length(selected_codes))
    for (jj in seq_along(selected_codes)) {
      code_i <- selected_codes[jj]
      city_data <- copy(dat[code == code_i])
      city_data[, stratum_city := interaction(
        year, month, dow, drop = TRUE
      )]
      city_result <- tryCatch({
        cb_temperature <- crossbasis(
          city_data$tmean,
          lag = as.integer(spec$lag),
          argvar = argvar,
          arglag = arglag
        )
        cb_humidity <- crossbasis(
          city_data$RH,
          lag = as.integer(spec$lag),
          argvar = list(fun = "ns", df = 3),
          arglag = arglag
        )
        if (isTRUE(spec$add_rain)) {
          cb_rainfall <- crossbasis(
            city_data$rainfall,
            lag = as.integer(spec$lag),
            argvar = list(fun = "ns", df = 3),
            arglag = arglag
          )
          model <- gnm(
            n ~ cb_temperature + cb_humidity +
              cb_rainfall + lag_PM25,
            family = quasipoisson(),
            eliminate = factor(stratum_city),
            data = city_data,
            na.action = na.exclude
          )
        } else {
          model <- gnm(
            n ~ cb_temperature + cb_humidity + lag_PM25,
            family = quasipoisson(),
            eliminate = factor(stratum_city),
            data = city_data,
            na.action = na.exclude
          )
        }
        reduced <- crossreduce(
          cb_temperature,
          model,
          type = "overall",
          cen = cen_common
        )
        coefficient <- coef(reduced)
        covariance <- vcov(reduced)
        if (
          any(!is.finite(coefficient)) ||
            !is_posdef(covariance)
        ) {
          stop("Non-finite coefficient or non-positive-definite vcov.")
        }
        list(
          success = TRUE,
          coefficient = coefficient,
          covariance = covariance,
          message = ""
        )
      }, error = function(e) {
        list(
          success = FALSE,
          coefficient = NULL,
          covariance = NULL,
          message = conditionMessage(e)
        )
      })
      if (isTRUE(city_result$success)) {
        raw_coef[[code_i]] <- city_result$coefficient
        raw_vcov[[code_i]] <- city_result$covariance
      }
      fit_rows[[jj]] <- data.table(
        code = code_i,
        success = isTRUE(city_result$success),
        message = city_result$message
      )
      if (jj %% 25L == 0L || jj == length(selected_codes)) {
        message(sprintf(
          "  city %d/%d; successful %d",
          jj, length(selected_codes), length(raw_coef)
        ))
      }
    }
    fit_status <- rbindlist(fit_rows, fill = TRUE)
    if (length(raw_coef) < 8L) {
      stop("Fewer than eight city-specific models succeeded.")
    }

    raw_fit <- list(
      disease = spec$disease,
      disease_tag = spec$disease_tag,
      classification = spec$classification,
      class_tag = spec$class_tag,
      dat = dat,
      city_meta = city_meta,
      selected_codes = selected_codes,
      argvar = argvar,
      arglag = arglag,
      raw_coef = raw_coef,
      raw_vcov = raw_vcov,
      cen_common = cen_common,
      fit_status = fit_status,
      minimum_city_cases = minimum_city_cases,
      config = list(
        lag = as.integer(spec$lag),
        lagnk = as.integer(spec$lagnk),
        add_rain = isTRUE(spec$add_rain)
      )
    )
    saveRDS(raw_fit, raw_path, compress = TRUE)
    prepared <- ts_prepare_major_fit(raw_fit, spec$class_tag)
    prepared$disease <- spec$disease
    prepared$disease_tag <- spec$disease_tag
    prepared$classification <- spec$classification
    prepared$minimum_city_cases <- minimum_city_cases
    saveRDS(prepared, prepared_path, compress = TRUE)

    data.table(
      disease = spec$disease,
      disease_tag = spec$disease_tag,
      classification = spec$classification,
      class_tag = spec$class_tag,
      lag = spec$lag,
      minimum_city_cases = minimum_city_cases,
      eligible_cities = length(selected_codes),
      successful_city_models = sum(fit_status$success),
      retained_first_stage = sum(
        prepared$coefficient_qc$retained_first_stage
      ),
      status = "fitted",
      message = "",
      prepared_fit = prepared_path
    )
  }, error = function(e) {
    data.table(
      disease = spec$disease,
      disease_tag = spec$disease_tag,
      classification = spec$classification,
      class_tag = spec$class_tag,
      lag = spec$lag,
      minimum_city_cases = NA_integer_,
      eligible_cities = NA_integer_,
      successful_city_models = NA_integer_,
      retained_first_stage = NA_integer_,
      status = "not estimable",
      message = conditionMessage(e),
      prepared_fit = prepared_path
    )
  })
  status_rows[[ii]] <- result
  fwrite(rbindlist(status_rows, fill = TRUE), STATUS_FILE)
  invisible(gc())
}

status <- rbindlist(status_rows, fill = TRUE)
status[, disease_order := match(
  disease, FIG3_DISEASE_SPECS$disease
)]
setorder(status, disease_order)
status[, disease_order := NULL]
fwrite(status, STATUS_FILE)
print(status)
