###############################################################################
## 01.0 Prepare paired factual / hist-nat temperature files for all 10 GCMs
##
## Reads each large source RDS only once, retains 2008-2019 and all 493 cities,
## and writes one compact paired RDS per GCM for the full attribution analysis.
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
})

ROOT <- Sys.getenv(
  "ATTR_PROJECT_ROOT",
  getwd()
)
CACHE <- file.path(ROOT, "full_analysis/cache/climate_pairs")
OUT <- file.path(ROOT, "full_analysis/output")
dir.create(CACHE, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

source(file.path(ROOT, "00-paths.R"))
PATH_FACT <- PATH_GCM_FACTUAL
PATH_NAT <- PATH_GCM_HISTNAT
REBUILD <- tolower(Sys.getenv("FULL_REBUILD_CLIMATE_CACHE", "false")) == "true"
YEARS <- 2008:2019

pad_code <- function(x) str_pad(as.character(x), 6, side = "right", pad = "0")
safe_name <- function(x) gsub("[^A-Za-z0-9_-]", "_", x)
pair_file <- function(gcm) file.path(CACHE, paste0("climate_pair_", safe_name(gcm), ".rds"))
fact_file <- function(gcm) file.path(CACHE, paste0("_factual_", safe_name(gcm), ".rds"))

existing <- list.files(CACHE, pattern = "^climate_pair_.*\\.rds$", full.names = TRUE)
if (length(existing) == 10L && !REBUILD) {
  message("All 10 paired GCM cache files already exist; set FULL_REBUILD_CLIMATE_CACHE=true to rebuild.")
  quit(save = "no", status = 0L)
}

message("Loading full factual temperature file...")
fact <- as.data.table(readRDS(PATH_FACT))
stopifnot(all(c("CityID", "gcm", "date", "year", "tas_corrected_weighted") %in% names(fact)))
fact <- fact[year %in% YEARS,
             .(code = pad_code(CityID), gcm = as.character(gcm),
               date = as.Date(date), tf = as.numeric(tas_corrected_weighted))]
gcms <- sort(unique(fact$gcm))
if (length(gcms) != 10L) stop("Expected 10 factual GCMs; found ", length(gcms))
if (fact[, anyDuplicated(.SD), .SDcols = c("code", "gcm", "date")] > 0L) {
  stop("Duplicate city-GCM-date keys in factual data.")
}
for (g in gcms) {
  message("Writing factual partition: ", g)
  saveRDS(fact[gcm == g, .(code, date, tf)], fact_file(g), compress = TRUE)
}
rm(fact); invisible(gc())

message("Loading full hist-nat temperature file...")
nat <- as.data.table(readRDS(PATH_NAT))
stopifnot(all(c("CityID", "gcm", "date", "year", "tas_corrected_weighted") %in% names(nat)))
nat <- nat[year %in% YEARS,
           .(code = pad_code(CityID), gcm = as.character(gcm),
             date = as.Date(date), tn = as.numeric(tas_corrected_weighted))]
nat_gcms <- sort(unique(nat$gcm))
if (!identical(gcms, nat_gcms)) stop("Factual and hist-nat GCM names do not match.")
if (nat[, anyDuplicated(.SD), .SDcols = c("code", "gcm", "date")] > 0L) {
  stop("Duplicate city-GCM-date keys in hist-nat data.")
}

qc_rows <- vector("list", length(gcms))
for (i in seq_along(gcms)) {
  g <- gcms[i]
  target <- pair_file(g)
  if (file.exists(target) && !REBUILD) {
    message("Reusing paired cache: ", g)
    pair <- as.data.table(readRDS(target))
  } else {
    message(sprintf("Pairing GCM %d/%d: %s", i, length(gcms), g))
    f <- as.data.table(readRDS(fact_file(g)))
    n <- nat[gcm == g, .(code, date, tn)]
    setkey(f, code, date)
    setkey(n, code, date)
    pair <- merge(f, n, by = c("code", "date"), all = TRUE)
    pair <- pair[is.finite(tf) & is.finite(tn)]
    setorder(pair, code, date)
    saveRDS(pair, target, compress = TRUE)
    rm(f, n); invisible(gc())
  }
  qc_rows[[i]] <- pair[, .(
    gcm = g,
    cities = uniqueN(code),
    paired_rows = .N,
    first_date = min(date),
    last_date = max(date),
    mean_signal = mean(tf - tn),
    median_signal = median(tf - tn),
    min_city_signal = min(tapply(tf - tn, code, mean)),
    max_city_signal = max(tapply(tf - tn, code, mean))
  )]
  rm(pair); invisible(gc())
}
rm(nat); invisible(gc())

temporary <- vapply(gcms, fact_file, character(1))
unlink(temporary[file.exists(temporary)])
qc <- rbindlist(qc_rows)
fwrite(qc, file.path(OUT, "00-climate-pair-cache-qc.csv"))
print(qc)
message("Prepared all GCM pairs in: ", CACHE)
