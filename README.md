# Attributing infectious diseases of children, adolescents and young adults to human-induced climate change in China

Analysis code for the manuscript of the same title. The repository contains
the code that produces every number, table and model estimate reported in the
paper and its Supplementary Information. Code that only draws figures is not
included.

## Requirements

- R 4.4.1 with: data.table (1.18.2.1), dlnm (2.4.10), gnm (1.1.5),
  mixmeta (1.2.2), tsModel (0.6.2), splines, sf (1.0.20), terra (1.9.11),
  exactextractr (0.10.1), lubridate, stringr, dplyr, tidyr, scales, openxlsx,
  archive, parallel.
- Python 3 with pandas, numpy and openpyxl (one data-preparation step).

## Data

The notifiable infectious disease records are confidential and cannot be
shared (see the Data availability statement). All other inputs are public or
derived from public data. Scripts read every input through `00-paths.R`.
By default the inputs are expected under `data/`; set `ATTR_DATA_DIR` to use
another location.

```
data/
  health/DATA.incidence.id.RData          daily case counts by city and disease (restricted)
  population/pop.total.RData              child and adolescent population by city and year (restricted)
  weather/China_cities_total_fire_PM25_O3_weather_2007_2019.rds
                                          observed city-level daily temperature, humidity,
                                          precipitation and PM2.5 (population weighted)
  climate/all_gcms_all_years_factual_pop.rds
  climate/all_gcms_all_years_hist_nat_pop.rds
  climate/all_gcms_all_years_factual_grid.rds     (sensitivity analysis only)
  climate/all_gcms_all_years_hist_nat_grid.rds    (sensitivity analysis only)
                                          bias-adjusted CMIP6 temperatures for 10 GCMs,
                                          aggregated to cities
  boundaries/city_map.shp                 analysis-unit polygons
  boundaries/province.shp                 province polygons
  boundaries/geometryProvince/*.json      city polygons used for the Hu Huanyong line sides
  socioeconomic/gridded_gdp/GDP2008.tif ... GDP2019.tif
                                          annual gridded GDP, https://doi.org/10.5281/zenodo.7898409
  socioeconomic/rast_adm2_gdp_perCapita_1990_2024.tif
                                          gridded GDP per capita, https://doi.org/10.5281/zenodo.10976733
  socioeconomic/human_footprint/hfpYYYY.zip
                                          https://doi.org/10.6084/m9.figshare.16571064
  socioeconomic/nighttime_light/SRUNet_NPP_VIIRS_V2_Like_YYYY.zip
                                          https://doi.org/10.6084/m9.figshare.22262545
  sanitation/water_pollution_raster/<COD|NH3N>/<Urban|Rural|Industry>/chn_1km_pol_<sector>_<year>.tif
                                          https://doi.org/10.5281/zenodo.16930436
                                          (the archive names the folder NH3-N; rename it NH3N)
  sanitation/2006_2019_WWTP_GHG_emissions_in_China.xlsx
                                          https://doi.org/10.6084/m9.figshare.19158206.v1
```

The city-level climate files in `data/climate/` were prepared before this
pipeline: daily CMIP6 historical/SSP2-4.5 (factual) and hist-nat
(natural-only) near-surface temperature from 10 GCMs (r1i1p1f1), bias-adjusted
against ERA5 with an additive monthly correction on a 0.25° grid, then
aggregated to cities with population and boundary-overlap weights (Methods).

## Running

Run everything from the repository root:

```
Rscript run_all.R
```

Stages can be run separately with `ATTR_STAGES`, in this order:
`models`, `attribution`, `sensitivity`, `inequality`, `sanitation`. The
sensitivity stage refits all models under alternative lag windows and is
skipped unless `ATTR_RUN_SENSITIVITY=true` or it is named in `ATTR_STAGES`.
Most scripts cache intermediate results and reuse them when their inputs have
not changed.

Results are written as CSV files to `user_route/output/`, `full_analysis/`,
`Figure2/` to `Figure5/`, `Sensitivity/results/`, `Supplementary/Tables/`,
`SII/output/`, `Mechanism/output/` and `two_stage_analysis/output/`. The
`Figure*` folder names follow the figures whose numbers the scripts compute;
no figures are drawn.

## Where each result is computed

| Result | Script |
|---|---|
| Paired factual and natural-only city temperatures | `01.0-prepare-full-climate-pairs.R` |
| Human-induced temperature change, national and regional (Fig. 1, Supp. Fig. S1) | `V01-00`, `V01-01`, `V01-03` |
| National exposure–response functions, three disease groups | `03.1-fit-user-route-pooled-erf.R` |
| City-specific models, BLUPs and regional pooled RRs | `01.1-fit-full-major-classes.R` (runs `00.2`) |
| Relative risks at temperature percentiles (Supp. Table S2), regional RRs at the regional P95 (Fig. 2f) | `V01-04-exposure-response-summaries.R` |
| Direct HICC-attributable fraction and cases, national and by GCM (Supp. Tables S3, S7) | `V01-07-calculate-direct-counterfactual-AF.R` |
| Disease-relevant warm-season contribution, national (Supp. Table S4) | `V01-13-hottest-four-months-AF.R` |
| Annual estimates and trends (Supp. Table S8) | `04.1-run-user-route-annual-acc.R` |
| Temperature-range decomposition and 1 °C contribution curves (Fig. 3c–f) | `V01-05-temperature-range-decomposition.R` |
| Regional direct attribution and case composition (Supp. Table S5, Fig. 4e) | `V01-06-regional-direct-attribution.R` |
| Regional disease-relevant contribution (Supp. Table S4, Fig. 4) | `V01-06b-regional-disease-relevant-attribution.R` |
| Disease-specific models and attribution (Supp. Table S6) | `10.1-fit-Figure3-disease-erf.R`, `V01-08-disease-specific-attribution.R` |
| Sensitivity analyses (Supp. Fig. S10) | `V01-10-run-sensitivity-analysis.R`, `V01-11-sensitivity-summaries.R` |
| Assembled supplementary tables | `V01-09-supplementary-tables.R` |
| GDP-related SII, ΔSII and decile contrasts (Fig. 6a–c) | `SII/08` (via `SII/13-gdp-gradient-statistics.R`) |
| ΔT and +1 °C response contrasts, Hu Huanyong line, bivariate classes (Fig. 6d–e) | `SII/09`, `SII/10-gdp-gradient-contrasts.R`, `SII/14` |
| SII for three disease groups; night-time light and human footprint (Supp. Figs. S7, S8) | `SII/12-three-disease-and-indicator-sii.R` |
| Sanitation-system stress index | `Mechanism/02`, `Mechanism/04` |
| Sanitation interaction and standardized-response models (Supp. Text, Supp. Fig. S11, Supp. Table S9) | `Mechanism/05`, `Mechanism/10` |
| Hydrological context and interaction (Supp. Fig. S9) | `Mechanism/09`, `Mechanism/10` |
| Disease-specific sanitation rows of Supp. Table S9 | `two_stage_analysis/05`, `Mechanism/06`, `Mechanism/07` |
| Cross-estimator rows of Supp. Table S9 | `Mechanism/14-hu-line-diagnostic.R` |

`V01-09-supplementary-tables.R` numbers its output tables differently from
the paper: its STable3 is Supplementary Table S7, its STable4 is
Supplementary Table S3, and Supplementary Table S4 is the national summary
from `V01-13` combined with the regional summary from `V01-06b`.

## Notes

- The primary estimand is the direct paired contrast
  AF = 1 − RR(natural-only) / RR(factual), using the same exposure–response
  coefficients in both climates (`R/attribution-core.R`). It uses observed
  daily forward case weights; the disease-relevant warm-season measure uses
  day-of-year climatological case weights.
- Uncertainty intervals are percentiles of 10,000 bootstrap resamples of the
  10 GCM point estimates. They do not include exposure–response parameter
  uncertainty.
- `SII/06-process-gdp-per-capita.R` was rewritten for this release because
  the original script was lost. It takes the coverage-weighted mean of the
  per-capita raster within each analysis unit. Its values agree with those
  used in the paper for most units but differ for some units on boundaries
  (by up to about 10%), so rerunning the sanitation stage may change the
  sanitation estimates slightly.
