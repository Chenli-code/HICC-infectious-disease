###############################################################################
## Project configuration
##
## This is the only configuration file used by the final analysis pipeline.
## Final anthropogenic climate-change attribution uses the direct paired
## factual-versus-hist-nat contrast defined in R/attribution-core.R.
###############################################################################

USER_ROOT <- Sys.getenv(
  "ATTR_PROJECT_ROOT",
  getwd()
)
source(file.path(USER_ROOT, "00-paths.R"))
USER_OUT <- file.path(USER_ROOT, "user_route/output")
USER_CACHE <- file.path(USER_ROOT, "full_analysis/cache/climate_pairs")
USER_CORE <- file.path(USER_ROOT, "R", "attribution-core.R")
dir.create(USER_OUT, recursive = TRUE, showWarnings = FALSE)

USER_PATH_HEALTH <- PATH_HEALTH
USER_PATH_OBS <- PATH_OBS_WEATHER
USER_START_DATE <- as.Date("2008-01-01")
USER_END_DATE <- as.Date("2019-12-31")
USER_DATES <- seq(USER_START_DATE, USER_END_DATE, by = "day")

# Gastrointestinal/enterovirus and vector-borne ERFs require a numerical
# centring temperature. P50 is retained for curve estimation and presentation.
# The disease-relevant attribution window is defined separately as each city's
# four hottest consecutive months; it must not be encoded as a percentile
# reference temperature.
USER_NONRESP_REFERENCE_TYPE <- "p50"
USER_NONRESP_REFERENCE_PROBABILITY <- 0.50
USER_WARM_MONTH_COUNT <- 4L

USER_SPECS <- data.frame(
  classification = c(
    "Respiratory diseases",
    "Gastrointestinal and enterovirus",
    "Vectorborne"
  ),
  tag = c("respiratory", "gastrointestinal", "vectorborne"),
  lag = c(7L, 3L, 14L),
  add_rain = c(FALSE, FALSE, TRUE),
  reference = c(
    "mrt",
    USER_NONRESP_REFERENCE_TYPE,
    USER_NONRESP_REFERENCE_TYPE
  ),
  stringsAsFactors = FALSE
)

user_pad_code <- function(x) {
  stringr::str_pad(as.character(x), width = 6L, side = "right", pad = "0")
}

user_safe_gcm <- function(x) gsub("[^A-Za-z0-9_-]", "_", x)

USER_GCMS <- c(
  "ACCESS-CM2", "ACCESS-ESM1-5", "BCC-CSM2-MR", "CanESM5",
  "FGOALS-g3", "GFDL-ESM4", "ipsl-cm6a-lr", "MIROC6",
  "MRI-ESM2-0", "NorESM2-LM"
)

user_pair_file <- function(gcm) {
  file.path(USER_CACHE, paste0("climate_pair_", user_safe_gcm(gcm), ".rds"))
}

USER_DISEASE_LABELS <- c(
  respiratory = "Respiratory infectious diseases",
  gastrointestinal = "Gastrointestinal and enterovirus diseases",
  vectorborne = "Vector-borne diseases"
)

USER_DISEASE_SHORT <- c(
  respiratory = "Respiratory",
  gastrointestinal = "Gastrointestinal / enterovirus",
  vectorborne = "Vector-borne"
)

USER_REGION_BY_PROVINCE <- c(
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

if (!file.exists(USER_CORE)) {
  stop("Missing shared attribution module: ", USER_CORE)
}
