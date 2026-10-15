###############################################################################
## Two-stage DLNM analysis: shared configuration
##
## This analysis is intentionally isolated from the existing pooled-model
## pipeline. Existing figures and results are never overwritten.
###############################################################################

TS_ROOT <- normalizePath(
  file.path(
    Sys.getenv(
      "ATTR_PROJECT_ROOT",
      getwd()
    ),
    "two_stage_analysis"
  ),
  winslash = "/",
  mustWork = FALSE
)
PROJECT_ROOT <- normalizePath(
  file.path(TS_ROOT, ".."),
  winslash = "/",
  mustWork = TRUE
)

source(file.path(PROJECT_ROOT, "03.0-user-route-config.R"))
source(USER_CORE)

TS_OUTPUT <- file.path(TS_ROOT, "output")
TS_MAJOR_OUTPUT <- file.path(TS_OUTPUT, "major_classes")
TS_DISEASE_OUTPUT <- file.path(TS_OUTPUT, "individual_diseases")
TS_FIGURE_ROOT <- file.path(TS_ROOT, "figures")
TS_CACHE <- file.path(TS_OUTPUT, "cache")

invisible(lapply(
  c(
    TS_OUTPUT, TS_MAJOR_OUTPUT, TS_DISEASE_OUTPUT,
    TS_FIGURE_ROOT, TS_CACHE
  ),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))

TS_GCMS <- USER_GCMS
TS_GCM_FILTER <- trimws(Sys.getenv("TS_GCM_FILTER", ""))
if (nzchar(TS_GCM_FILTER)) {
  requested_gcms <- trimws(strsplit(TS_GCM_FILTER, "\\|")[[1L]])
  TS_GCMS <- intersect(TS_GCMS, requested_gcms)
  if (!length(TS_GCMS)) {
    stop("TS_GCM_FILTER did not match a configured GCM.")
  }
}
TS_DATES <- USER_DATES
TS_BASELINE <- Sys.getenv("TS_CASE_BASELINE", "observed_daily")
TS_BOOTSTRAP_REPETITIONS <- 10000L
TS_BOOTSTRAP_SEED <- 20260727L
TS_NONRESP_REFERENCE_PROBABILITY <- 0.80
TS_MIN_CITY_CASES <- as.integer(Sys.getenv("TS_MIN_CITY_CASES", "300"))

TS_TAGS <- c("respiratory", "gastrointestinal", "vectorborne")
TS_LABELS <- c(
  respiratory = "Respiratory infectious diseases",
  gastrointestinal = "Gastrointestinal and enterovirus diseases",
  vectorborne = "Vector-borne infectious diseases"
)
TS_SHORT_LABELS <- c(
  respiratory = "Respiratory",
  gastrointestinal = "Gastrointestinal /\nenterovirus",
  vectorborne = "Vector-borne"
)
TS_COLOURS <- c(
  respiratory = "#5F788D",
  gastrointestinal = "#AD735C",
  vectorborne = "#718473"
)

TS_MAJOR_FIT_FILES <- c(
  respiratory = file.path(
    PROJECT_ROOT, "full_analysis/output/resp_fit_allcities.rds"
  ),
  gastrointestinal = file.path(
    PROJECT_ROOT,
    "full_analysis/output/gastrointestinal_fit_allcities.rds"
  ),
  vectorborne = file.path(
    PROJECT_ROOT, "full_analysis/output/vectorborne_fit_allcities.rds"
  )
)

stopifnot(
  TS_BASELINE %in% c("observed_daily", "day_of_year"),
  all(file.exists(TS_MAJOR_FIT_FILES)),
  all(file.exists(vapply(TS_GCMS, user_pair_file, character(1))))
)

ts_safe_name <- function(x) {
  gsub("^_+|_+$", "", gsub("[^A-Za-z0-9]+", "_", x))
}

ts_major_cache_file <- function(gcm) {
  file.path(
    TS_CACHE,
    paste0("two_stage_major_", user_safe_gcm(gcm), ".rds")
  )
}
