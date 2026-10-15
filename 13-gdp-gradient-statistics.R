# -*- coding: UTF-8 -*-
# ============================================================================
# Figure 6: quantitative contrasts for GDP-related temperature-attributable AF
# 独立统计脚本，不修改原分析对象或原绘图脚本。
#
# 推荐：在生成原图的 R 会话中运行：
#   source("SII/13-figure6-quantitative-statistics.R", encoding = "UTF-8")
# 必需对象：analysis_data, sii_by_gcm（来自 08-gastrointestinal-gdp-sii.R）。
# 可选对象：panel_F1, panel_F2（原图每个 GCM 的两组估计）。
# 没有内存对象时读取 analysis-inputs.rds；也可通过 SII_STATS_INPUT 指定
# 含上述命名对象的 list RDS。都没有时重跑 08 并缓存 A-D 的输入。
# F1/F2 不从图片或汇总中位数反推；缺失时在状态表中明确标记。
#
# 输出：统计汇总、逐模型结果、年度结果、人口覆盖、输入状态、统计说明。
# 范围为 GCM 估计的 2.5–97.5 百分位（type=7），不是完整统计置信区间。
# 不把年份×GCM当作独立样本，不对十个相关模型做 t / Wilcoxon 显著性检验。
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

sii_stats_root <- function() {
  candidates <- c(Sys.getenv("ATTR_PROJECT_ROOT", ""), getwd(), dirname(getwd()))
  for (p in candidates[nzchar(candidates)]) {
    if (file.exists(file.path(p, "SII", "08-gastrointestinal-gdp-sii.R"))) {
      return(normalizePath(p, winslash = "/", mustWork = TRUE))
    }
  }
  stop("Run from the project root or set ATTR_PROJECT_ROOT to that directory.")
}

sii_stats_require <- function(x, columns, name) {
  missing <- setdiff(columns, names(x))
  if (length(missing)) stop(name, " missing columns: ", paste(missing, collapse = ", "))
}

sii_stats_unique <- function(x, keys, name) {
  if (anyDuplicated(x, by = keys)) stop(name, " contains duplicate keys: ", paste(keys, collapse = ", "))
}

sii_stats_complete_years <- function(x, years, keys, name) {
  coverage <- x[, .(complete = setequal(year, years)), by = keys]
  if (any(!coverage$complete)) stop(name, " has incomplete study years; no partial-period means are allowed.")
}

sii_stats_summarise <- function(x) {
  x[, {
    v <- value[is.finite(value)]
    .(estimate = if (length(v)) median(v) else NA_real_,
      lower = if (length(v)) unname(quantile(v, .025, type = 7)) else NA_real_,
      upper = if (length(v)) unname(quantile(v, .975, type = 7)) else NA_real_,
      n_gcm = length(v), n_positive = sum(v > 0), n_negative = sum(v < 0),
      n_zero = sum(v == 0))
  }, by = .(period, metric, unit)]
}

run_sii_quantitative_statistics <- function(
    analysis_data, sii_by_gcm, panel_F1 = NULL, panel_F2 = NULL,
    disease = "Gastrointestinal and enterovirus", years = 2008:2019,
    output_dir = file.path(sii_stats_root(), "SII", "output", "quantitative-statistics")) {
  # Work on copies so the original figure objects are preserved.
  a <- copy(as.data.table(analysis_data))
  s <- copy(as.data.table(sii_by_gcm))
  sii_stats_require(a, c("disease", "code_merge", "year", "gcm", "GDP",
                        "population", "socioeconomic_rank", "af_factual", "af_counterfactual"), "analysis_data")
  sii_stats_require(s, c("disease", "year", "gcm", "sii_factual", "sii_counterfactual", "delta_sii"), "sii_by_gcm")
  selected_disease <- disease
  a <- a[disease == selected_disease & year %in% years]
  s <- s[disease == selected_disease & year %in% years]
  if (!nrow(a) || !nrow(s)) stop("No matching disease/year observations.")
  a[, gcm := as.character(gcm)]
  s[, gcm := as.character(gcm)]
  sii_stats_unique(a, c("gcm", "year", "code_merge"), "analysis_data")
  sii_stats_unique(s, c("gcm", "year"), "sii_by_gcm")
  if (!setequal(a$gcm, s$gcm)) stop("AF and SII have different GCM sets.")
  if (any(!is.finite(a$population) | a$population <= 0 |
          !is.finite(a$af_factual) | !is.finite(a$af_counterfactual) |
          !is.finite(a$socioeconomic_rank) | a$socioeconomic_rank < 0 |
          a$socioeconomic_rank > 1)) stop("Invalid AF, population or rank values.")
  if (any(!is.finite(s$sii_factual) | !is.finite(s$sii_counterfactual) |
          !is.finite(s$delta_sii))) stop("Non-finite SII values.")
  if (any(abs(s$delta_sii - (s$sii_factual - s$sii_counterfactual)) > 1e-7)) {
    stop("Paired SII identity fails: delta_sii must equal factual minus natural-only.")
  }
  sii_stats_complete_years(s, years, "gcm", "SII")
  model_names <- sort(unique(a$gcm))
  period_label <- paste0(min(years), "-", max(years), " annual mean")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # D: match the original population-rank deciles, then pool cities using their
  # actual population weights, NOT an unweighted average of two decile means.
  a[, gdp_decile := pmin(10L, pmax(1L, floor(socioeconomic_rank * 10) + 1L))]
  a[, gdp_group := fcase(gdp_decile <= 2L, "Lower GDP",
                         gdp_decile >= 9L, "Higher GDP", default = "Middle GDP")]
  a[, delta_af := af_factual - af_counterfactual]
  coverage <- a[, .(n_cities = uniqueN(code_merge), population = sum(population)),
                by = .(gcm, year, gdp_group)]
  coverage[, population_share := population / sum(population), by = .(gcm, year)]
  grouped <- a[gdp_group != "Middle GDP", .(
    delta_af = weighted.mean(delta_af, population),
    factual_af = weighted.mean(af_factual, population),
    natural_af = weighted.mean(af_counterfactual, population)
  ), by = .(gcm, year, gdp_group)]
  d <- dcast(grouped, gcm + year ~ gdp_group, value.var = "delta_af")
  sii_stats_require(d, c("Lower GDP", "Higher GDP"), "GDP group contrasts")
  if (any(!is.finite(d[["Lower GDP"]]) | !is.finite(d[["Higher GDP"]]))) {
    stop("A GCM/year lacks one of the two GDP groups.")
  }
  setnames(d, c("Lower GDP", "Higher GDP"), c("lower_gdp", "higher_gdp"))
  d[, difference := lower_gdp - higher_gdp]
  sii_stats_complete_years(d, years, "gcm", "GDP contrasts")
  dp <- d[, .(lower_gdp = mean(lower_gdp), higher_gdp = mean(higher_gdp),
               difference = mean(difference)), by = gcm]

  # B-C: average paired annual estimates WITHIN each GCM before summarising GCMs.
  sp <- s[, .(sii_factual = mean(sii_factual),
              sii_counterfactual = mean(sii_counterfactual),
              delta_sii = mean(delta_sii)), by = gcm]
  sp[, period := period_label]
  sy <- copy(s)[, period := as.character(year)]
  st <- rbindlist(list(sp, sy), use.names = TRUE, fill = TRUE)
  # Relative change is change in the signed positive GDP gradient, not a
  # percentage change in health outcomes. Suppress a period if any GCM has a
  # non-positive gradient or a denominator numerically indistinguishable from 0.
  st[, relative_eligible := all(sii_counterfactual > 1e-8 & sii_factual > 0), by = period]
  st[, relative_sii_change_pct := fifelse(relative_eligible,
      100 * delta_sii / sii_counterfactual, NA_real_)]
  all_values <- melt(st, id.vars = c("gcm", "period"),
                    measure.vars = c("sii_factual", "sii_counterfactual", "delta_sii", "relative_sii_change_pct"),
                    variable.name = "metric", value.name = "value")
  all_values[, metric := as.character(metric)]
  all_values[, unit := fifelse(metric == "relative_sii_change_pct", "%", "pp")]
  dp[, period := period_label]
  dy <- copy(d)[, period := as.character(year)]
  dt <- rbindlist(list(dp, dy), use.names = TRUE, fill = TRUE)
  dl <- melt(dt, id.vars = c("gcm", "period"),
             measure.vars = c("lower_gdp", "higher_gdp", "difference"),
             variable.name = "metric", value.name = "value")
  dl[, metric := paste0("delta_af_", metric)]
  dl[, unit := "pp"]
  all_values <- rbindlist(list(all_values, dl), use.names = TRUE)

  # F1-F2: consume raw paired GCM estimates from the original figure, preserving
  # its group definitions, study-period aggregation and temperature reference.
  status <- data.table(panel = c("B-C", "D"), status = "computed",
                       detail = c("Paired SII by GCM/year", "Population-weighted GDP groups"))
  paired_f <- function(x, prefix, unit, panel) {
    if (is.null(x)) return(NULL)
    x <- copy(as.data.table(x))
    if ("disease" %in% names(x)) x <- x[disease == selected_disease]
    sii_stats_require(x, c("gcm", "gdp_group", "estimate"), panel)
    x[, gcm := as.character(gcm)]
    x[, gdp_group := as.character(gdp_group)]
    if (!setequal(x$gdp_group, c("Lower GDP", "Higher GDP"))) {
      stop(panel, " requires group labels 'Lower GDP' and 'Higher GDP'.")
    }
    sii_stats_unique(x, c("gcm", "gdp_group"), panel)
    if (!setequal(x$gcm, model_names)) stop(panel, " has a different GCM set from A-D.")
    w <- dcast(x, gcm ~ gdp_group, value.var = "estimate")
    if (any(!is.finite(w[["Lower GDP"]]) | !is.finite(w[["Higher GDP"]]))) {
      stop(panel, " has incomplete or non-finite paired estimates.")
    }
    setnames(w, c("Lower GDP", "Higher GDP"), c("lower_gdp", "higher_gdp"))
    w[, difference := lower_gdp - higher_gdp]
    z <- melt(w, id.vars = "gcm", variable.name = "metric", value.name = "value")
    z[, metric := paste0(prefix, "_", metric)]
    z[, period := period_label]
    z[, unit := unit]
    z[, .(gcm, period, metric, value, unit)]
  }
  f_specs <- list(list(x = panel_F1, prefix = "warming", unit = "deg C", panel = "F1"),
                  list(x = panel_F2, prefix = "response_1C", unit = "pp per +1 deg C", panel = "F2"))
  for (spec in f_specs) {
    z <- paired_f(spec$x, spec$prefix, spec$unit, spec$panel)
    if (is.null(z)) {
      detail <- paste0("Missing panel_", spec$panel, "; run in the original figure session or supply it in SII_STATS_INPUT.")
      warning(detail, call. = FALSE)
      status <- rbind(status, data.table(panel = spec$panel, status = "missing_input", detail = detail))
    } else {
      all_values <- rbindlist(list(all_values, z), use.names = TRUE)
      status <- rbind(status, data.table(panel = spec$panel, status = "computed",
                                        detail = "Paired model-level estimates from original figure object"))
    }
  }
  sii_stats_unique(all_values, c("gcm", "period", "metric"), "Statistics output")
  summary <- sii_stats_summarise(all_values)
  summary[, interval_type := "GCM 2.5th-97.5th percentile range (type 7); not a confidence interval"]
  summary[, label := sprintf("%.2f [%.2f, %.2f] %s", estimate, lower, upper, unit)]
  setorder(summary, period, metric)
  setorder(all_values, period, metric, gcm)
  fwrite(summary, file.path(output_dir, "figure6-statistics-summary.csv"))
  fwrite(all_values, file.path(output_dir, "figure6-statistics-by-gcm.csv"))
  fwrite(d, file.path(output_dir, "figure6-gdp-contrast-by-year-gcm.csv"))
  fwrite(coverage, file.path(output_dir, "figure6-gdp-group-coverage.csv"))
  fwrite(status, file.path(output_dir, "figure6-input-status.csv"))

  notes <- c(
    "# Figure 6 quantitative statistics", "",
    paste0("Disease: ", selected_disease), paste0("Period: ", period_label), "",
    "- Delta AF is factual temperature-attributable AF minus natural-only temperature-attributable AF (percentage points). It is NOT the direct HICC-attributable fraction among all observed cases.",
    "- GDP groups follow the original annual population-weighted ridit deciles (1-2 versus 9-10). Whole cities are retained; actual population shares may differ slightly from 20%. See coverage CSV.",
    "- City AFs are population-weighted within each group/year/GCM. Years are equally weighted within each GCM; GCMs are summarised only after these operations.",
    "- Mean-period relative SII change = 100 * mean annual Delta SII / mean annual natural-only SII, within each GCM. It is not the mean of annual ratios. Relative change is suppressed when a period contains a non-positive gradient.",
    "- F1 and F2 reuse original paired GCM estimates and their original exposure distributions. They are not recalculated from the image or reconstructed from published point estimates.",
    "- Every contrast is paired within GCM. The median of paired differences need not equal the difference of group medians.",
    "- Interval bounds are empirical 2.5th and 97.5th GCM percentiles (R quantile type 7). They exclude exposure-response parameter uncertainty, city sampling uncertainty and other model structural uncertainty.",
    "- Counts of positive/negative GCM estimates describe directional agreement, not independent hypothesis tests. No p-values are manufactured from GCM/year pseudo-replicates.",
    "- A smaller positive SII is attenuation of a GDP gradient. It does not on its own mean health improvement or establish that absolute case burdens are larger in poorer regions.", "",
    "## Input status", capture.output(print(status)), "", "## Period estimates",
    capture.output(print(summary[period == period_label, .(metric, label, n_gcm, n_positive, n_negative)])))
  writeLines(notes, file.path(output_dir, "figure6-statistics-notes.md"), useBytes = TRUE)
  print(summary[period == period_label, .(metric, label, n_gcm, n_positive, n_negative)])
  message("Statistics saved to: ", normalizePath(output_dir, winslash = "/"))
  invisible(list(summary = summary, by_gcm = all_values, annual_gdp = d,
                 coverage = coverage, input_status = status))
}

# Run using the current figure session, or explicit/cached analysis inputs.
figure6_quantitative_statistics <- local({
  root <- sii_stats_root()
  out <- file.path(root, "SII", "output", "quantitative-statistics")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  explicit_input <- Sys.getenv("SII_STATS_INPUT", "")
  current <- parent.env(environment())
  if (nzchar(explicit_input)) {
    inputs <- readRDS(explicit_input)
    if (!is.list(inputs)) stop("SII_STATS_INPUT must contain a named list.")
  } else if (exists("analysis_data", envir = current, inherits = TRUE) &&
             exists("sii_by_gcm", envir = current, inherits = TRUE)) {
    inputs <- list(analysis_data = get("analysis_data", envir = current, inherits = TRUE),
                   sii_by_gcm = get("sii_by_gcm", envir = current, inherits = TRUE))
  } else {
    cache <- file.path(out, "analysis-inputs.rds")
    if (file.exists(cache)) {
      message("Using cached A-D inputs: ", cache)
      inputs <- readRDS(cache)
    } else {
      message("Recalculating A-D inputs with SII/08-gastrointestinal-gdp-sii.R ...")
      old_wd <- getwd()
      on.exit(setwd(old_wd), add = TRUE)
      setwd(root)
      e <- new.env(parent = globalenv())
      sys.source(file.path("SII", "08-gastrointestinal-gdp-sii.R"), envir = e)
      inputs <- list(analysis_data = e$analysis_data, sii_by_gcm = e$sii_by_gcm)
      saveRDS(inputs, cache)
    }
  }
  for (nm in c("panel_F1", "panel_F2")) {
    if (!nzchar(explicit_input) && exists(nm, envir = current, inherits = TRUE)) {
      inputs[[nm]] <- get(nm, envir = current, inherits = TRUE)
    }
  }
  if (is.null(inputs$analysis_data) || is.null(inputs$sii_by_gcm)) {
    stop("Input list must contain analysis_data and sii_by_gcm.")
  }
  selected <- if (!is.null(inputs$plot_disease)) inputs$plot_disease else if (
    !nzchar(explicit_input) && exists("plot_disease", envir = current, inherits = TRUE)) {
    get("plot_disease", envir = current, inherits = TRUE)
  } else "Gastrointestinal and enterovirus"
  run_sii_quantitative_statistics(inputs$analysis_data, inputs$sii_by_gcm,
    panel_F1 = inputs$panel_F1, panel_F2 = inputs$panel_F2,
    disease = selected, output_dir = out)
})
