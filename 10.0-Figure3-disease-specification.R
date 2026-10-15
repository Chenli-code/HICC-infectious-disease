###############################################################################
## Figure 3 disease specification
##
## Each disease is fitted separately with the same national pooled
## case-crossover/DLNM strategy used for the three major disease classes.
## Disease-specific lags are copied from the previously published temperature
## analysis. Three very sparse diseases were not fitted in the old scripts;
## they use the corresponding major-class lag as an explicit fallback.
###############################################################################

suppressPackageStartupMessages(library(data.table))
ROOT <- Sys.getenv(
  "ATTR_PROJECT_ROOT",
  getwd()
)
source(file.path(ROOT, "03.0-user-route-config.R"))

FIG3_DIR <- file.path(USER_OUT, "figure3_disease")
FIG3_FIT_DIR <- file.path(FIG3_DIR, "fits")
FIG3_ATTR_DIR <- file.path(FIG3_DIR, "attribution")
dir.create(FIG3_FIT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG3_ATTR_DIR, recursive = TRUE, showWarnings = FALSE)

figure3_safe_name <- function(x) {
  y <- gsub("[^A-Za-z0-9]+", "_", x)
  gsub("^_+|_+$", "", y)
}

FIG3_DISEASE_SPECS <- data.table::data.table(
  classification = c(
    rep("Respiratory diseases", 9),
    rep("Gastrointestinal and enterovirus", 7),
    rep("Vectorborne", 8)
  ),
  class_tag = c(
    rep("respiratory", 9),
    rep("gastrointestinal", 7),
    rep("vectorborne", 8)
  ),
  disease = c(
    "Seasonal influenza", "Mumps", "Tuberculosis", "Scarlet fever",
    "Rubella", "Measles", "Pertussis", "Meningococcal meningitis",
    "Leprosy",
    "Hand, foot, and mouth disease", "Infectious diarrhoea", "Dysentery",
    "Acute haemorrhagic conjunctivitis", "Hepatitis A",
    "Typhoid and paratyphoid", "Cholera",
    "Dengue", "Haemorrhagic fever", "Malaria", "Japanese encephalitis",
    "Typhus", "Schistosomiasis", "Kala-azar", "Plague"
  ),
  lag = c(
    7L, 21L, 21L, 12L, 21L, 14L, 14L, 10L, 21L,
    5L, 3L, 3L, 2L, 30L, 14L, 3L,
    8L, 12L, 14L, 14L, 14L, 21L, 14L, 14L
  ),
  lagnk = c(
    rep(2L, 9),
    2L, 2L, 2L, 1L, 2L, 2L, 2L,
    rep(2L, 8)
  ),
  reference = c(
    rep("mrt", 9),
    rep(USER_NONRESP_REFERENCE_TYPE, 7),
    rep(USER_NONRESP_REFERENCE_TYPE, 8)
  ),
  add_rain = c(
    rep(FALSE, 9),
    rep(FALSE, 7),
    rep(TRUE, 8)
  ),
  lag_source = c(
    rep("published disease-specific temperature model", 9),
    rep("published disease-specific temperature model", 6),
    "major-class fallback; no prior disease-specific model",
    rep("published disease-specific temperature model", 6),
    "major-class fallback; no prior disease-specific model",
    "major-class fallback; no prior disease-specific model"
  )
)

FIG3_DISEASE_SPECS[, disease_tag := figure3_safe_name(disease)]
FIG3_DISEASE_SPECS[, disease_order := seq_len(.N)]

stopifnot(
  nrow(FIG3_DISEASE_SPECS) == 24L,
  uniqueN(FIG3_DISEASE_SPECS$disease) == 24L,
  all(FIG3_DISEASE_SPECS$lag >= 0L)
)

figure3_fit_file <- function(disease_tag) {
  file.path(FIG3_FIT_DIR, paste0("erf_disease_", disease_tag, ".rds"))
}

figure3_attr_file <- function(gcm) {
  file.path(
    FIG3_ATTR_DIR,
    paste0("disease_attribution_", user_safe_gcm(gcm), ".rds")
  )
}
