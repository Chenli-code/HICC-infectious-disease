# -*- coding: UTF-8 -*-
# ==============================================================================
# 胃肠道疾病组：GDP相关SII
#
# 行政编码：
#   直辖市保留六位区县编码；其他省份使用前四位 + 00。
# 结局：
#   factual AF（%）与 counterfactual AF（%）。
# 分析单位：
#   每个年份 × 每个GCM分别计算SII。
# 人口权重：
#   固定使用2015年人口；无法匹配的地区不纳入分析。
# 输出对象（仅保留在内存中，不保存文件）：
#   af_gcm_raw、af_gcm、af_summary、analysis_data、
#   city_number_check、sii_by_gcm、sii_summary。
# ==============================================================================

options(stringsAsFactors = FALSE, scipen = 999)

suppressPackageStartupMessages({
  library(data.table)
  library(dlnm)
  library(tsModel)
})

# ==============================================================================
# 0. 参数
# ==============================================================================

analysis_years <- 2008:2019
analysis_dates <- seq.Date(
  as.Date("2008-01-01"),
  as.Date("2019-12-31"),
  by = "day"
)

population_weight_year <- 2015L
minimum_cities <- 20L
disease_name <- "Gastrointestinal and enterovirus"

gcm_names <- c(
  "ACCESS-CM2",
  "ACCESS-ESM1-5",
  "BCC-CSM2-MR",
  "CanESM5",
  "FGOALS-g3",
  "GFDL-ESM4",
  "ipsl-cm6a-lr",
  "MIROC6",
  "MRI-ESM2-0",
  "NorESM2-LM"
)

municipality_provinces <- c("11", "12", "31", "50")

# ==============================================================================
# 1. 文件路径
# ==============================================================================

project_root <- normalizePath(".", winslash = "/", mustWork = TRUE)

core_file <- file.path(
  project_root,
  "R",
  "attribution-core.R"
)

pooled_fit_file <- file.path(
  project_root,
  "user_route",
  "output",
  "erf_pooled_gastrointestinal.rds"
)

gastrointestinal_data_file <- file.path(
  project_root,
  "full_analysis",
  "output",
  "gastrointestinal_fit_allcities.rds"
)

climate_dir <- file.path(
  project_root,
  "full_analysis",
  "cache",
  "climate_pairs"
)

gdp_file <- file.path(
  project_root,
  "SII",
  "data",
  "gridded-gdp",
  "gridded-gdp-2008-2019.csv"
)

source(file.path(project_root, "00-paths.R"))
population_file_candidates <- PATH_POP_TOTAL

population_file <- population_file_candidates[
  file.exists(population_file_candidates)
][1L]

if (length(population_file) == 0L || is.na(population_file)) {
  stop("找不到人口文件 pop.total.RData。")
}

# ==============================================================================
# 2. 读取模型和胃肠道病例数据
# ==============================================================================

source(core_file)

# 合并暴露-反应模型。
fit <- readRDS(pooled_fit_file)

# 每个地区每日胃肠道疾病病例数据。
gastrointestinal_data <- readRDS(gastrointestinal_data_file)

health_codes <- sprintf(
  "%06d",
  as.integer(gastrointestinal_data$selected_codes)
)

counts <- as.data.table(gastrointestinal_data$dat)[, .(
  n = sum(as.numeric(n), na.rm = TRUE)
), by = .(
  code = sprintf("%06d", as.integer(code)),
  date = as.Date(date)
)]

counts <- counts[
  code %in% health_codes &
    date %in% analysis_dates
]

setorder(counts, code, date)

# ==============================================================================
# 3. 工具函数
# ==============================================================================

median_or_na <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else median(x)
}

extract_sii <- function(model) {
  estimate <- unname(coef(model)["socioeconomic_rank"])
  standard_error <- summary(model)$coefficients[
    "socioeconomic_rank",
    "Std. Error"
  ]
  confidence_interval <- confint(
    model,
    "socioeconomic_rank",
    level = 0.95
  )

  c(
    estimate = estimate,
    standard_error = standard_error,
    lower = confidence_interval[1L],
    upper = confidence_interval[2L]
  )
}

annual_scenario_af_city <- function(
    factual_temperature,
    counterfactual_temperature,
    cases,
    dates,
    fit
) {
  forward_cases <- acc_forward_weights(cases, fit$lag)

  basis_factual <- acc_one_basis(
    factual_temperature,
    fit$argvar
  )
  basis_counterfactual <- acc_one_basis(
    counterfactual_temperature,
    fit$argvar
  )
  basis_reference <- as.numeric(acc_one_basis(
    as.numeric(fit$cen),
    fit$argvar
  ))

  valid <- is.finite(forward_cases) &
    is.finite(factual_temperature) &
    is.finite(counterfactual_temperature) &
    rowSums(!is.finite(basis_factual)) == 0L &
    rowSums(!is.finite(basis_counterfactual)) == 0L

  if (!any(valid)) {
    return(data.table(
      year = analysis_years,
      case_weight = NA_real_,
      af_factual = NA_real_,
      af_counterfactual = NA_real_
    ))
  }

  centred_factual <- sweep(
    basis_factual[valid, , drop = FALSE],
    2L,
    basis_reference,
    "-"
  )
  centred_counterfactual <- sweep(
    basis_counterfactual[valid, , drop = FALSE],
    2L,
    basis_reference,
    "-"
  )

  eta_factual <- drop(centred_factual %*% fit$reduced_coef)
  eta_counterfactual <- drop(
    centred_counterfactual %*% fit$reduced_coef
  )

  daily_af <- data.table(
    year = as.integer(format(dates[valid], "%Y")),
    case_weight = forward_cases[valid],
    attributable_factual =
      forward_cases[valid] * (1 - exp(-eta_factual)),
    attributable_counterfactual =
      forward_cases[valid] * (1 - exp(-eta_counterfactual))
  )

  annual_af <- daily_af[, .(
    case_weight = sum(case_weight, na.rm = TRUE),
    attributable_factual = sum(
      attributable_factual,
      na.rm = TRUE
    ),
    attributable_counterfactual = sum(
      attributable_counterfactual,
      na.rm = TRUE
    )
  ), by = year]

  annual_af[, `:=`(
    af_factual = fifelse(
      case_weight > 0,
      100 * attributable_factual / case_weight,
      NA_real_
    ),
    af_counterfactual = fifelse(
      case_weight > 0,
      100 * attributable_counterfactual / case_weight,
      NA_real_
    )
  )]

  annual_af[, .(
    year,
    case_weight,
    af_factual,
    af_counterfactual
  )]
}

# ==============================================================================
# 4. 计算原始地区编码下的年度AF
# ==============================================================================

af_results <- vector("list", length(gcm_names))

for (gcm_index in seq_along(gcm_names)) {
  current_gcm <- gcm_names[gcm_index]

  message(sprintf(
    "[%02d/%02d] 正在处理：%s",
    gcm_index,
    length(gcm_names),
    current_gcm
  ))

  climate_file <- file.path(
    climate_dir,
    paste0("climate_pair_", current_gcm, ".rds")
  )

  climate_pair <- readRDS(climate_file)
  climate <- acc_complete_climate_pair(
    pair = climate_pair,
    codes = health_codes,
    dates = analysis_dates
  )

  daily_data <- merge(
    climate,
    counts,
    by = c("code", "date"),
    all.x = TRUE
  )
  daily_data[is.na(n), n := 0]
  setorder(daily_data, code, date)

  current_af <- daily_data[, {
    annual_scenario_af_city(
      factual_temperature = tf,
      counterfactual_temperature = tn,
      cases = n,
      dates = date,
      fit = fit
    )
  }, by = code]

  current_af[, gcm := current_gcm]
  af_results[[gcm_index]] <- current_af

  rm(climate_pair, climate, daily_data, current_af)
  gc(verbose = FALSE)
}

af_gcm_raw <- rbindlist(af_results, use.names = TRUE)
af_gcm_raw[, disease := disease_name]
af_gcm_raw[, climate_change_af :=
  af_factual - af_counterfactual]
setorder(af_gcm_raw, gcm, code, year)

# ==============================================================================
# 5. 统一行政区划编码
# ==============================================================================

af_gcm_raw[, code := sprintf("%06d", as.integer(code))]

af_gcm_raw[, code_merge := fcase(
  # 直辖市保留六位区县编码。
  substr(code, 1L, 2L) %chin% municipality_provinces &
    substr(code, 5L, 6L) != "00",
  code,

  # 删除直辖市的上级汇总编码。
  substr(code, 1L, 2L) %chin% municipality_provinces,
  NA_character_,

  # 非直辖市使用前四位 + 00。
  default = paste0(substr(code, 1L, 4L), "00")
)]

af_gcm_raw <- af_gcm_raw[
  !is.na(code_merge) &
    is.finite(case_weight) &
    case_weight > 0 &
    is.finite(af_factual) &
    is.finite(af_counterfactual)
]

# ==============================================================================
# 6. 按code_merge汇总AF
# ==============================================================================

# 多个原始地区合并成一个地级分析单元时，使用case_weight加权AF。
af_gcm <- af_gcm_raw[, .(
  case_weight = sum(case_weight, na.rm = TRUE),
  af_factual = weighted.mean(
    af_factual,
    w = case_weight,
    na.rm = TRUE
  ),
  af_counterfactual = weighted.mean(
    af_counterfactual,
    w = case_weight,
    na.rm = TRUE
  )
), by = .(
  disease,
  code_merge,
  year,
  gcm
)]

af_gcm[, climate_change_af :=
  af_factual - af_counterfactual]
setorder(af_gcm, gcm, code_merge, year)

# 10个GCM的城市年度AF汇总。
af_summary <- af_gcm[, .(
  case_weight = median_or_na(case_weight),
  af_factual_median = median_or_na(af_factual),
  af_factual_q025 = quantile(af_factual, 0.025, na.rm = TRUE),
  af_factual_q975 = quantile(af_factual, 0.975, na.rm = TRUE),
  af_counterfactual_median = median_or_na(af_counterfactual),
  af_counterfactual_q025 = quantile(
    af_counterfactual,
    0.025,
    na.rm = TRUE
  ),
  af_counterfactual_q975 = quantile(
    af_counterfactual,
    0.975,
    na.rm = TRUE
  ),
  climate_change_af_median = median_or_na(climate_change_af),
  climate_change_af_q025 = quantile(
    climate_change_af,
    0.025,
    na.rm = TRUE
  ),
  climate_change_af_q975 = quantile(
    climate_change_af,
    0.975,
    na.rm = TRUE
  ),
  n_gcm = uniqueN(gcm)
), by = .(
  disease,
  code_merge,
  year
)]

setorder(af_summary, code_merge, year)

# ==============================================================================
# 7. 读取2015年人口并统一编码
# ==============================================================================

population_environment <- new.env(parent = emptyenv())
load(population_file, envir = population_environment)
population_raw <- as.data.table(
  population_environment[["pop.total"]]
)

population_raw[, code := sprintf("%06d", as.integer(code))]

population_raw[, code_merge := fcase(
  # 直辖市保留六位区县编码。
  substr(code, 1L, 2L) %chin% municipality_provinces &
    substr(code, 5L, 6L) != "00",
  code,

  # 删除直辖市的上级汇总编码。
  substr(code, 1L, 2L) %chin% municipality_provinces,
  NA_character_,

  # 非直辖市使用前四位 + 00。
  default = paste0(substr(code, 1L, 4L), "00")
)]

population_weight <- population_raw[
  as.integer(year) == population_weight_year &
    !is.na(code_merge) &
    is.finite(as.numeric(pop)) &
    as.numeric(pop) > 0,
  .(
    population = sum(as.numeric(pop), na.rm = TRUE)
  ),
  by = code_merge
]

population_weight[, population_weight_year :=
  population_weight_year]

# ==============================================================================
# 8. 读取GDP并整理成长数据
# ==============================================================================

gdp_wide <- fread(
  gdp_file,
  colClasses = list(character = "code_merge")
)
gdp_wide[, code_merge := sprintf(
  "%06d",
  as.integer(code_merge)
)]

gdp_columns <- paste0("GDP", analysis_years)

gdp_long <- melt(
  gdp_wide,
  id.vars = "code_merge",
  measure.vars = gdp_columns,
  variable.name = "year_column",
  value.name = "GDP"
)

gdp_long[, `:=`(
  year = as.integer(sub("GDP", "", year_column)),
  GDP = as.numeric(GDP)
)]

gdp_long <- gdp_long[, .(
  code_merge,
  year,
  GDP
)]

# ==============================================================================
# 9. 合并AF、GDP和固定人口权重
# ==============================================================================

analysis_data <- merge(
  af_gcm,
  gdp_long,
  by = c("code_merge", "year"),
  all = FALSE
)

analysis_data <- merge(
  analysis_data,
  population_weight,
  by = "code_merge",
  all = FALSE
)

analysis_data <- analysis_data[
  year %in% analysis_years &
    is.finite(GDP) &
    is.finite(population) &
    population > 0 &
    is.finite(af_factual) &
    is.finite(af_counterfactual)
]

setorder(analysis_data, year, gcm, GDP, code_merge)

city_number_check <- analysis_data[, .(
  n_cities = uniqueN(code_merge)
), by = .(
  year,
  gcm
)]

# ==============================================================================
# 10. 根据每年GDP构造人口加权ridit
# ==============================================================================

rank_table <- unique(analysis_data[, .(
  code_merge,
  year,
  GDP,
  population
)])

# GDP完全相同的地区共享同一个ridit。
gdp_groups <- rank_table[, .(
  group_population = sum(population, na.rm = TRUE)
), by = .(
  year,
  GDP
)]

setorder(gdp_groups, year, GDP)

gdp_groups[, cumulative_population_before :=
  data.table::shift(
    cumsum(group_population),
    fill = 0
  ),
by = year]

gdp_groups[, socioeconomic_rank :=
  (
    cumulative_population_before +
      0.5 * group_population
  ) /
  sum(group_population),
by = year]

rank_table <- merge(
  rank_table,
  gdp_groups[, .(
    year,
    GDP,
    socioeconomic_rank
  )],
  by = c("year", "GDP"),
  all.x = TRUE
)

analysis_data <- merge(
  analysis_data,
  rank_table[, .(
    code_merge,
    year,
    socioeconomic_rank
  )],
  by = c("code_merge", "year"),
  all = FALSE
)

# ==============================================================================
# 11. 每个年份、每个GCM计算SII
# ==============================================================================

sii_by_gcm <- analysis_data[, {
  current_data <- .SD[
    is.finite(socioeconomic_rank) &
      is.finite(population) &
      population > 0 &
      is.finite(af_factual) &
      is.finite(af_counterfactual)
  ]

  if (
    nrow(current_data) < minimum_cities ||
      uniqueN(current_data$socioeconomic_rank) < 2L
  ) {
    data.table(
      n_cities = nrow(current_data),
      population_total = sum(current_data$population, na.rm = TRUE),
      sii_factual = NA_real_,
      sii_factual_se = NA_real_,
      sii_factual_lower = NA_real_,
      sii_factual_upper = NA_real_,
      sii_counterfactual = NA_real_,
      sii_counterfactual_se = NA_real_,
      sii_counterfactual_lower = NA_real_,
      sii_counterfactual_upper = NA_real_,
      delta_sii = NA_real_,
      delta_sii_se = NA_real_,
      delta_sii_lower = NA_real_,
      delta_sii_upper = NA_real_,
      delta_abs_sii = NA_real_,
      change_low_gdp_disadvantage = NA_real_
    )
  } else {
    factual_model <- lm(
      af_factual ~ socioeconomic_rank,
      data = current_data,
      weights = population
    )
    counterfactual_model <- lm(
      af_counterfactual ~ socioeconomic_rank,
      data = current_data,
      weights = population
    )
    difference_model <- lm(
      I(af_factual - af_counterfactual) ~ socioeconomic_rank,
      data = current_data,
      weights = population
    )

    factual_sii <- extract_sii(factual_model)
    counterfactual_sii <- extract_sii(counterfactual_model)
    difference_sii <- extract_sii(difference_model)

    data.table(
      n_cities = nrow(current_data),
      population_total = sum(current_data$population),
      sii_factual = factual_sii["estimate"],
      sii_factual_se = factual_sii["standard_error"],
      sii_factual_lower = factual_sii["lower"],
      sii_factual_upper = factual_sii["upper"],
      sii_counterfactual = counterfactual_sii["estimate"],
      sii_counterfactual_se = counterfactual_sii["standard_error"],
      sii_counterfactual_lower = counterfactual_sii["lower"],
      sii_counterfactual_upper = counterfactual_sii["upper"],
      delta_sii = difference_sii["estimate"],
      delta_sii_se = difference_sii["standard_error"],
      delta_sii_lower = difference_sii["lower"],
      delta_sii_upper = difference_sii["upper"],
      delta_abs_sii =
        abs(factual_sii["estimate"]) -
        abs(counterfactual_sii["estimate"]),
      change_low_gdp_disadvantage =
        -difference_sii["estimate"]
    )
  }
}, by = .(
  disease,
  year,
  gcm
)]

sii_by_gcm[, factual_direction := fcase(
  sii_factual < 0,
  "higher_AF_in_low_GDP",
  sii_factual > 0,
  "higher_AF_in_high_GDP",
  default = "no_gradient"
)]

sii_by_gcm[, counterfactual_direction := fcase(
  sii_counterfactual < 0,
  "higher_AF_in_low_GDP",
  sii_counterfactual > 0,
  "higher_AF_in_high_GDP",
  default = "no_gradient"
)]

sii_by_gcm[, absolute_inequality_change := fcase(
  delta_abs_sii > 0,
  "widened",
  delta_abs_sii < 0,
  "narrowed",
  default = "unchanged"
)]

setorder(sii_by_gcm, year, gcm)

# ==============================================================================
# 12. 汇总10个GCM
# ==============================================================================

sii_summary <- sii_by_gcm[, .(
  n_gcm = sum(
    is.finite(sii_factual) &
      is.finite(sii_counterfactual)
  ),
  n_cities = median(n_cities, na.rm = TRUE),

  sii_factual_median = median(sii_factual, na.rm = TRUE),
  sii_factual_q025 = quantile(sii_factual, 0.025, na.rm = TRUE),
  sii_factual_q975 = quantile(sii_factual, 0.975, na.rm = TRUE),

  sii_counterfactual_median = median(
    sii_counterfactual,
    na.rm = TRUE
  ),
  sii_counterfactual_q025 = quantile(
    sii_counterfactual,
    0.025,
    na.rm = TRUE
  ),
  sii_counterfactual_q975 = quantile(
    sii_counterfactual,
    0.975,
    na.rm = TRUE
  ),

  delta_sii_median = median(delta_sii, na.rm = TRUE),
  delta_sii_q025 = quantile(delta_sii, 0.025, na.rm = TRUE),
  delta_sii_q975 = quantile(delta_sii, 0.975, na.rm = TRUE),

  delta_abs_sii_median = median(delta_abs_sii, na.rm = TRUE),
  delta_abs_sii_q025 = quantile(
    delta_abs_sii,
    0.025,
    na.rm = TRUE
  ),
  delta_abs_sii_q975 = quantile(
    delta_abs_sii,
    0.975,
    na.rm = TRUE
  ),

  low_gdp_disadvantage_change_median = median(
    change_low_gdp_disadvantage,
    na.rm = TRUE
  ),
  low_gdp_disadvantage_change_q025 = quantile(
    change_low_gdp_disadvantage,
    0.025,
    na.rm = TRUE
  ),
  low_gdp_disadvantage_change_q975 = quantile(
    change_low_gdp_disadvantage,
    0.975,
    na.rm = TRUE
  ),

  n_gcm_widened = sum(delta_abs_sii > 0, na.rm = TRUE),
  n_gcm_narrowed = sum(delta_abs_sii < 0, na.rm = TRUE)
), by = .(
  disease,
  year
)]

sii_summary[, climate_inequality_conclusion := fcase(
  delta_abs_sii_q025 > 0,
  "widened",
  delta_abs_sii_q975 < 0,
  "narrowed",
  default = "uncertain"
)]

sii_summary[, `:=`(
  ses_indicator = "GDP_total",
  rank_direction = "low_GDP_to_high_GDP",
  outcome = "attributable_fraction_percent",
  population_weight_year = population_weight_year,
  sii_unit = "AF_percentage_points"
)]

setorder(sii_summary, year)

# ==============================================================================
# 13. 查看结果（不保存）
# ==============================================================================

print(city_number_check)
print(sii_summary)
