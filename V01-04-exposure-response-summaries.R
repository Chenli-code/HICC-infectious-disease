###############################################################################
## Figure 2 | Contrasting temperature-infectious disease associations
##
## A-C  National pooled cumulative DLNM exposure-response curves using the
##      original model reference-temperature definitions and coefficients.
## D    Heatmap of RR (95% CI) at eight observed temperature percentiles.
## E    National P5 and P95 RR forest plot.
## F    Region-specific P95 RR forest plot from random-effects meta-analysis
##      of the existing city-specific DLNM coefficients.
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(dlnm)
  library(scales)
  library(mixmeta)
})

ROOT <- Sys.getenv(
  "ATTR_PROJECT_ROOT",
  getwd()
)

source(file.path(ROOT, "03.0-user-route-config.R"))

# ------------------------------- paths ---------------------------------------
FIGURE2_DIR <- file.path(ROOT, "Figure2")

dir.create(FIGURE2_DIR, recursive = TRUE, showWarnings = FALSE)

FIT_FILES <- setNames(
  file.path(USER_OUT, paste0("erf_pooled_", USER_SPECS$tag, ".rds")),
  USER_SPECS$tag
)

stopifnot(all(file.exists(FIT_FILES)), file.exists(USER_PATH_OBS))

disease_order <- c("respiratory", "gastrointestinal", "vectorborne")

disease_labels <- c(
  respiratory = "Respiratory infectious diseases",
  gastrointestinal = "Gastrointestinal and enterovirus diseases",
  vectorborne = "Vector-borne infectious diseases"
)

disease_short <- c(
  respiratory = "Respiratory",
  gastrointestinal = "Gastrointestinal and\nenterovirus",
  vectorborne = "Vector-borne"
)

percentile_order <- c(
  "P5", "P10", "P25", "P50", "P75", "P90", "P95", "P99"
)

region_order <- c(
  "Northeast", "North", "East", "Central",
  "South", "Southwest", "Northwest"
)

province_region <- c(
  `11` = "North", `12` = "North", `13` = "North",
  `14` = "North", `15` = "North",
  `21` = "Northeast", `22` = "Northeast", `23` = "Northeast",
  `31` = "East", `32` = "East", `33` = "East", `34` = "East",
  `35` = "East", `36` = "East", `37` = "East",
  `41` = "Central", `42` = "Central", `43` = "Central",
  `44` = "South", `45` = "South", `46` = "South",
  `50` = "Southwest", `51` = "Southwest", `52` = "Southwest",
  `53` = "Southwest", `54` = "Southwest",
  `61` = "Northwest", `62` = "Northwest", `63` = "Northwest",
  `64` = "Northwest", `65` = "Northwest"
)

REGIONAL_FIT_FILES <- c(
  respiratory = file.path(
    USER_ROOT, "full_analysis/output/resp_fit_allcities.rds"
  ),
  gastrointestinal = file.path(
    USER_ROOT, "full_analysis/output/gastrointestinal_fit_allcities.rds"
  ),
  vectorborne = file.path(
    USER_ROOT, "full_analysis/output/vectorborne_fit_allcities.rds"
  )
)

stopifnot(all(file.exists(REGIONAL_FIT_FILES)))

# ------------------------------- exposure data -------------------------------
fits <- lapply(FIT_FILES, readRDS)

fit_reference_current <- vapply(
  names(fits),
  function(tag_i) {
    expected_i <- USER_SPECS$reference[
      match(tag_i, USER_SPECS$tag)
    ]
    identical(
      tolower(as.character(fits[[tag_i]]$reference)),
      tolower(as.character(expected_i))
    )
  },
  logical(1)
)

if (!all(fit_reference_current)) {
  stop(
    "One or more pooled ERFs use an outdated reference. Rerun ",
    "03.1-fit-user-route-pooled-erf.R before Figure 2."
  )
}

# Reconstruct the observed temperature panel using the same completion rule as
# the pooled ERF models: city observation, province-day mean, then national-day
# mean. This makes the temperature percentiles match the model input.
obs <- as.data.table(readRDS(USER_PATH_OBS))

obs[, `:=`(
  code = user_pad_code(id),
  date = as.Date(date),
  tmean = as.numeric(tmean)
)]

obs <- obs[
  date >= USER_START_DATE & date <= USER_END_DATE,
  .(code, date, tmean)
]

obs[, province_code := substr(code, 1L, 2L)]

province_day <- obs[, .(
  tmean_province = mean(tmean, na.rm = TRUE)
), by = .(province_code, date)]

national_day <- obs[, .(
  tmean_national = mean(tmean, na.rm = TRUE)
), by = date]

health_codes <- lapply(fits, function(fit) {
  unique(as.character(fit$crossover$code))
})

all_health_codes <- sort(unique(unlist(health_codes)))

temperature_panel <- CJ(
  code = all_health_codes,
  date = USER_DATES,
  unique = TRUE
)

temperature_panel[, province_code := substr(code, 1L, 2L)]

temperature_panel <- merge(
  temperature_panel,
  obs[, .(code, date, tmean)],
  by = c("code", "date"), all.x = TRUE, sort = FALSE
)

temperature_panel <- merge(
  temperature_panel,
  province_day,
  by = c("province_code", "date"), all.x = TRUE, sort = FALSE
)

temperature_panel <- merge(
  temperature_panel,
  national_day,
  by = "date", all.x = TRUE, sort = FALSE
)

temperature_panel[, tmean := fcoalesce(
  tmean, tmean_province, tmean_national
)]

temperature_panel[, c("tmean_province", "tmean_national") := NULL]

stopifnot(all(is.finite(temperature_panel$tmean)))

national_percentiles <- as.numeric(quantile(
  temperature_panel$tmean,
  probs = c(0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99),
  na.rm = TRUE, names = FALSE
))

names(national_percentiles) <- c(
  "P1", "P5", "P10", "P25", "P50", "P75", "P90", "P95", "P99"
)

# A common absolute-temperature axis is restricted to observed P1-P99 support
# and to the intersection of the fitted spline boundaries.
bounds <- do.call(rbind, lapply(fits, function(fit) fit$argvar$Bound))

common_x_limits <- c(
  max(floor(national_percentiles["P1"]), max(bounds[, 1])),
  min(ceiling(national_percentiles["P99"]), min(bounds[, 2]))
)

if (common_x_limits[1] >= common_x_limits[2]) {
  stop("No common exposure range across the three fitted temperature bases.")
}

temperature_grid <- seq(
  common_x_limits[1], common_x_limits[2], length.out = 600
)

fwrite(
  data.table(
    percentile = names(national_percentiles),
    probability = c(
      0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99
    ),
    temperature = national_percentiles
  ),
  file.path(FIGURE2_DIR, "V01-Figure2-temperature-percentiles.csv")
)

# ------------------------------- curve estimates -----------------------------
one_basis <- function(x, argvar) {
  do.call(onebasis, c(list(x = x), argvar))
}

curve_rows <- vector("list", length(disease_order))

forest_rows <- vector("list", length(disease_order))

reference_rows <- vector("list", length(disease_order))

for (ii in seq_along(disease_order)) {
  tag_i <- disease_order[ii]
  fit_i <- fits[[tag_i]]

  basis_grid <- one_basis(temperature_grid, fit_i$argvar)
  reference_i <- as.numeric(fit_i$cen)
  reference_type_i <- if (tolower(fit_i$reference) == "mrt") {
    "MMT"
  } else if (grepl(
    "^p[0-9]+$",
    tolower(fit_i$reference)
  )) {
    toupper(as.character(fit_i$reference))
  } else {
    toupper(as.character(fit_i$reference))
  }

  basis_reference <- as.numeric(
    one_basis(reference_i, fit_i$argvar)
  )
  contrast_grid <- sweep(
    basis_grid, 2L, basis_reference, "-"
  )
  lp <- drop(contrast_grid %*% fit_i$reduced_coef)
  se <- sqrt(pmax(
    0,
    rowSums(
      (contrast_grid %*% fit_i$reduced_vcov) * contrast_grid
    )
  ))

  curve_rows[[ii]] <- data.table(
    tag = tag_i,
    disease = unname(disease_labels[tag_i]),
    temperature = temperature_grid,
    rr = exp(lp),
    rr_low = exp(lp - qnorm(0.975) * se),
    rr_high = exp(lp + qnorm(0.975) * se),
    reference_temperature = reference_i,
    reference_type = reference_type_i
  )

  selected_temps <- national_percentiles[percentile_order]
  selected_basis <- one_basis(
    as.numeric(selected_temps), fit_i$argvar
  )
  selected_contrast <- sweep(
    selected_basis, 2L, basis_reference, "-"
  )
  selected_lp <- drop(
    selected_contrast %*% fit_i$reduced_coef
  )
  selected_se <- sqrt(pmax(
    0,
    rowSums(
      (selected_contrast %*% fit_i$reduced_vcov) *
        selected_contrast
    )
  ))

  forest_rows[[ii]] <- data.table(
    tag = tag_i,
    disease = unname(disease_short[tag_i]),
    percentile = names(selected_temps),
    temperature = as.numeric(selected_temps),
    reference_temperature = reference_i,
    reference_type = reference_type_i,
    rr = exp(selected_lp),
    rr_low = exp(selected_lp - qnorm(0.975) * selected_se),
    rr_high = exp(selected_lp + qnorm(0.975) * selected_se)
  )

  reference_rows[[ii]] <- data.table(
    tag = tag_i,
    disease = unname(disease_labels[tag_i]),
    reference_temperature = reference_i,
    reference_type = reference_type_i,
    original_reference_setting = fit_i$reference,
    common_temperature_min = common_x_limits[1],
    common_temperature_max = common_x_limits[2]
  )
}

curve_dt <- rbindlist(curve_rows)

forest_dt <- rbindlist(forest_rows)

reference_dt <- rbindlist(reference_rows)

forest_dt[, percentile := factor(
  percentile, levels = percentile_order
)]

forest_dt[, disease := factor(
  disease,
  levels = rev(unname(disease_short[disease_order]))
)]

forest_dt[, estimate_label := sprintf(
  "%.2f (%.2f–%.2f)", rr, rr_low, rr_high
)]

regional_rows <- list()

row_id <- 0L

for (tag_i in disease_order) {
  source_i <- readRDS(REGIONAL_FIT_FILES[tag_i])
  source_i$dat$code <- as.character(source_i$dat$code)
  raw_codes <- intersect(
    names(source_i$raw_coef),
    names(source_i$raw_vcov)
  )
  code_region <- unname(
    province_region[substr(raw_codes, 1L, 2L)]
  )

  for (region_i in region_order) {
    region_codes <- raw_codes[
      !is.na(code_region) & code_region == region_i
    ]
    if (!length(region_codes)) next

    ymat <- do.call(rbind, source_i$raw_coef[region_codes])
    slist <- source_i$raw_vcov[region_codes]
    if (length(region_codes) == 1L) {
      regional_coef <- as.numeric(source_i$raw_coef[[region_codes]])
      regional_vcov <- as.matrix(source_i$raw_vcov[[region_codes]])
      regional_method <- "single city"
    } else {
      convergence_warning <- NULL
      regional_meta <- withCallingHandlers(
        tryCatch(
          mixmeta(
            ymat ~ 1, slist, method = "reml",
            control = list(maxiter = 1000)
          ),
          error = function(e) NULL
        ),
        warning = function(w) {
          if (grepl(
            "convergence not reached",
            conditionMessage(w), fixed = TRUE
          )) {
            convergence_warning <<- conditionMessage(w)
            invokeRestart("muffleWarning")
          }
        }
      )
      if (
        is.null(regional_meta) ||
          !is.null(convergence_warning)
      ) {
        regional_meta <- mixmeta(
          ymat ~ 1, slist, method = "fixed"
        )
      }
      regional_coef <- as.numeric(coef(regional_meta))
      regional_vcov <- as.matrix(vcov(regional_meta))
      regional_method <- regional_meta$method
    }

    region_temp <- source_i$dat[
      source_i$dat$code %in% region_codes,
      "tmean"
    ][[1L]]
    region_p95 <- as.numeric(quantile(
      region_temp, 0.95, na.rm = TRUE
    ))

    region_reference <- as.numeric(fits[[tag_i]]$cen)
    reference_type <- if (
      tolower(fits[[tag_i]]$reference) == "mrt"
    ) {
      "National MMT"
    } else {
      paste0(
        "National ",
        toupper(as.character(fits[[tag_i]]$reference))
      )
    }

    basis_p95 <- one_basis(region_p95, source_i$argvar)
    basis_ref <- as.numeric(
      one_basis(region_reference, source_i$argvar)
    )
    contrast <- sweep(basis_p95, 2L, basis_ref, "-")
    regional_lp <- drop(contrast %*% regional_coef)
    regional_se <- sqrt(pmax(
      0,
      drop((contrast %*% regional_vcov) %*% t(contrast))
    ))

    row_id <- row_id + 1L
    regional_rows[[row_id]] <- data.table(
      tag = tag_i,
      disease = unname(disease_short[tag_i]),
      region = region_i,
      percentile = "P95",
      temperature = region_p95,
      reference_temperature = region_reference,
      reference_type = reference_type,
      n_city_models = length(region_codes),
      meta_method = regional_method,
      rr = exp(regional_lp),
      rr_low = exp(regional_lp - qnorm(0.975) * regional_se),
      rr_high = exp(regional_lp + qnorm(0.975) * regional_se)
    )
  }
}

regional_dt <- rbindlist(regional_rows, fill = TRUE)

regional_dt[, `:=`(
  disease = factor(
    disease, levels = unname(disease_short[disease_order])
  ),
  region = factor(region, levels = rev(region_order)),
  estimate_label = sprintf(
    "%.2f (%.2f–%.2f)", rr, rr_low, rr_high
  )
)]

regional_dt[, estimate_label_f := sprintf(
  "%.2f (%.2f–%.2f)", rr, rr_low, rr_high
)]

regional_dt[, label_x_f := max(rr_high, na.rm = TRUE) * 1.07,
            by = disease]

fwrite(
  curve_dt,
  file.path(FIGURE2_DIR, "V01-Figure2-response-curves.csv")
)

fwrite(
  forest_dt,
  file.path(FIGURE2_DIR, "V01-Figure2-selected-percentile-RR.csv")
)

fwrite(
  reference_dt,
  file.path(FIGURE2_DIR, "V01-Figure2-model-reference-temperatures.csv")
)

fwrite(
  regional_dt,
  file.path(FIGURE2_DIR, "V01-Figure2-regional-P95-RR.csv")
)
