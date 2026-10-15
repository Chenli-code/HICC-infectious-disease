###############################################################################
## Run the full analysis in order
##
## Usage (from the repository root):
##   Rscript run_all.R
##
## Stages can be selected with ATTR_STAGES, e.g.
##   ATTR_STAGES="models,attribution" Rscript run_all.R
## Available stages: models, attribution, sensitivity, inequality, sanitation.
## The sensitivity stage refits the exposure-response functions under
## alternative lag windows and is slow; it is skipped unless selected
## explicitly or ATTR_RUN_SENSITIVITY=true.
###############################################################################

ROOT <- normalizePath(Sys.getenv("ATTR_PROJECT_ROOT", getwd()), winslash = "/")
setwd(ROOT)
Sys.setenv(ATTR_PROJECT_ROOT = ROOT)
source("00-paths.R")

stages_env <- Sys.getenv("ATTR_STAGES", "")
run_sensitivity <- tolower(Sys.getenv("ATTR_RUN_SENSITIVITY", "false")) %in%
  c("true", "1", "yes")
stages <- if (nzchar(stages_env)) {
  trimws(strsplit(stages_env, ",")[[1]])
} else {
  c("models", "attribution", if (run_sensitivity) "sensitivity",
    "inequality", "sanitation")
}

run_script <- function(file) {
  if (!file.exists(file)) stop("Missing pipeline script: ", file)
  message("\n", strrep("=", 78), "\nRunning: ", file, "\n", strrep("=", 78))
  setwd(ROOT)
  sys.source(file, envir = new.env(parent = globalenv()))
  setwd(ROOT)
}

run_python <- function(file) {
  message("\n", strrep("=", 78), "\nRunning: ", file, "\n", strrep("=", 78))
  python <- Sys.getenv("ATTR_PYTHON", "python")
  status <- system2(python, shQuote(file))
  if (!identical(status, 0L)) stop("Python step failed: ", file)
}

pipeline <- list(
  # Paired climate data and exposure-response functions.
  models = c(
    "01.0-prepare-full-climate-pairs.R",
    "03.1-fit-user-route-pooled-erf.R",
    "01.1-fit-full-major-classes.R",
    "10.1-fit-Figure3-disease-erf.R"
  ),
  # Climate summaries and attribution of disease burden.
  attribution = c(
    "V01-00-climate-summaries-by-gcm.R",
    "V01-01-human-induced-temperature-change.R",
    "V01-03-scenario-temperatures-by-city-region.R",
    "V01-04-exposure-response-summaries.R",
    "V01-07-calculate-direct-counterfactual-AF.R",
    "V01-13-hottest-four-months-AF.R",
    "04.1-run-user-route-annual-acc.R",
    "V01-05-temperature-range-decomposition.R",
    "V01-06-regional-direct-attribution.R",
    "V01-06b-regional-disease-relevant-attribution.R",
    "V01-08-disease-specific-attribution.R",
    "V01-09-supplementary-tables.R"
  ),
  sensitivity = c(
    "V01-10-run-sensitivity-analysis.R",
    "V01-11-sensitivity-summaries.R"
  ),
  # GDP-related inequality (slope index of inequality).
  inequality = c(
    "SII/03-process-gridded-gdp.R",
    "SII/04-process-human-footprint.R",
    "SII/05-process-nighttime-light.R",
    "SII/06-process-gdp-per-capita.R",
    "SII/13-gdp-gradient-statistics.R",
    "SII/09-figure6-plus1c-response.R",
    "SII/10-gdp-gradient-contrasts.R",
    "SII/12-three-disease-and-indicator-sii.R",
    "SII/14-qinling-huaihe-contrast.R"
  ),
  # Exploratory sanitation-system analyses (Supplementary Text).
  sanitation = c(
    "two_stage_analysis/01-calculate-major-class-two-stage-ACC.R",
    "two_stage_analysis/05-fit-disease-specific-two-stage.R",
    "Mechanism/data_external/_metadata/derive_wwtp_capacity.py",
    "Mechanism/02-build-urban-sanitation-stress.R",
    "Mechanism/04-build-sanitation-pressure-capacity.R",
    "Mechanism/09-build-hydro-context.R",
    "Mechanism/05-sanitation-response-and-interaction-models.R",
    "Mechanism/10-sanitation-decomposition-and-robustness.R",
    "Mechanism/06-disease-specific-sanitation-response.R",
    "Mechanism/07-disease-specific-warming-sanitation-interaction.R",
    "Mechanism/14-hu-line-diagnostic.R"
  )
)

unknown <- setdiff(stages, names(pipeline))
if (length(unknown)) stop("Unknown stage(s): ", paste(unknown, collapse = ", "))

for (stage in names(pipeline)[names(pipeline) %in% stages]) {
  message("\n#### Stage: ", stage)
  for (file in pipeline[[stage]]) {
    if (grepl("[.]py$", file)) run_python(file) else run_script(file)
  }
}

message("\nAnalysis pipeline completed.")
