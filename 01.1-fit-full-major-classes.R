###############################################################################
## 01.1 Fit the three published major-class specifications using all cities
###############################################################################

ROOT <- Sys.getenv(
  "ATTR_PROJECT_ROOT",
  getwd()
)
OUT <- file.path(ROOT, "full_analysis/output")
SCRIPT <- file.path(ROOT, "00.2-major-classes-pooled-vs-BLUP.R")
RSCRIPT <- file.path(
  R.home("bin"),
  if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"
)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

classes <- c(
  "Respiratory diseases",
  "Gastrointestinal and enterovirus",
  "Vectorborne"
)

min_cases <- Sys.getenv("FULL_MIN_CASES", "100")
reuse_fits <- tolower(Sys.getenv("FULL_REUSE_FITS", "true"))
status <- data.frame(
  classification = classes,
  started = NA_character_, finished = NA_character_, exit_status = NA_integer_
)

for (i in seq_along(classes)) {
  cls <- classes[i]
  cat("\n", strrep("=", 78), "\n", cls, "\n", sep = "")
  status$started[i] <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  Sys.setenv(
    ATTR_PROJECT_ROOT = ROOT,
    ATTR_OUTPUT_DIR = OUT,
    ATTR_TEST_CLASS = cls,
    RESP_TEST_MAX_CITIES = "0",
    RESP_TEST_MIN_CASES = min_cases,
    RESP_TEST_RUN_ATTRIB = "false",
    RESP_TEST_REUSE_FITS = reuse_fits
  )
  status$exit_status[i] <- system2(RSCRIPT, shQuote(SCRIPT))
  status$finished[i] <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  write.csv(status, file.path(OUT, "01-full-model-fit-status.csv"), row.names = FALSE)
}

print(status)
if (any(status$exit_status != 0L, na.rm = TRUE)) {
  stop("At least one full major-class model failed; inspect logs before attribution.")
}
