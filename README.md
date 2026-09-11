# Optimal timing of influenza vaccination

Estimates the MMWR week that maximizes expected direct influenza vaccine protection by state and age group, conditional on receiving one dose during the season. The model does not represent uptake, missed opportunities, transmission, or campaign logistics.

## Reproduce the publication analysis

Requirements are R 4.2 or later with Cairo graphics support, Bash, and network
access to the public CDC and Census sources. On Linux, the `sf`, `tigris`, and
figure-building steps also require system libraries for GDAL, GEOS, PROJ,
SQLite, UDUNITS, Cairo, and fontconfig. Install the R packages from the
repository root:

```sh
Rscript requirements.R
```

`requirements.R` installs missing CRAN packages and enforces the minimum
versions needed by the analysis code. This repository does not include a
package lockfile. R, package, platform, and external-library versions are
recorded with each run but intentionally not locked.

Run the canonical workflow:

```sh
scripts/run_publication_pipeline.sh
```

The runner uses 5,000 simulation draws. It downloads and prepares inputs,
runs the primary and prespecified sensitivity analyses, builds the publication
figures, validates the expected artifacts, and removes unreferenced caches and
exploratory figures. The canonical run freshly downloads the public surveillance,
population, boundary, and archived CDC vaccine-effectiveness tables it uses.
The cited 2010/11 and recent published VE estimates are encoded in the
normalization script, along with an explicit primary-source correction to the
2015/16 overall confidence interval. The pipeline does not rely on a downloaded
data snapshot stored in the repository.

For local development only, an existing set of downloaded inputs can be reused:

```sh
REFRESH_INPUTS=false scripts/run_publication_pipeline.sh
```

This cached mode is not the canonical publication workflow. The same input
validator runs in either mode and stops before model fitting if the inputs do
not match the prespecified analysis.

The sensitivity commands and parameter values are encoded in
`scripts/run_manuscript_sensitivity_suite.sh`. `R/publication_spec.R` lists the
expected analysis manifests used by validation and cache pruning. Each analysis
manifest records its input hashes, configuration, and R version.

The suite includes a no-shift sensitivity, `gam_no_timing_shift_primary`,
which sets `TIMING_SHIFT_SD=0` instead of the primary value of 0.75 weeks.
It uses the same fitted state GAMs, sampled seasons, coefficient draws, and
vaccine-protection draws as the primary analysis. Only the added epidemic
timing shift is removed. Its results are included in the robustness summary
and Figure 4, and publication validation checks its configuration and pairing
with the primary analysis.

Fitted-model cache names fingerprint the modeled surveillance input and the
direct GAM-fitting implementation. Cache reuse is an execution optimization;
a fresh public clone contains no fitted models.

The canonical runner starts a fresh R process for each age group to release
working memory between analyses. `scripts/assemble_age_specific_intervals.R` then combines
the five saved interval tables in their original age order.

## Prespecified analysis scope

The publication pipeline fixes the primary surveillance seasons at 2010/11
through 2018/19 and 2023/24 through 2025/26. It excludes 2019/20 through
2022/23. Influenza-like illness burden is evaluated from MMWR week 36 through
week 22 of the following calendar year. CDC download seasons begin at week 40,
so the download requests include the preceding season to obtain weeks 36--39.
Public state ILINet reporting begins at week 40 in 2010/11; that is the only
permitted later start. The HHS, national ILINet, and national laboratory series
cover the full week 36--22 window in every included season. The analytic
endpoint is fixed at MMWR year 2026, week 22.

The state analysis requires all 50 states and the District of Columbia in every
week. The age-specific input requires all 10 US Department of Health and Human
Services regions and five reported age groups. National ILINet and national
WHO/NREVSS records must align on the same weekly grid. Vaccine-effectiveness
estimates for the all-age analyses cover 2010/11 through 2018/19 and 2023/24
through 2025/26. Each age-specific VE pool includes 2010/11 through 2018/19 and
2023/24. The recent estimates come from the U.S. Flu VE Network: 2023/24
provides the five age strata, whereas the final 2024/25 and preliminary 2025/26
estimates enter only the all-age pool. The 2024/25 input uses the final
Chung et al. report (33%, 95% CI 24–41; DOI 10.1093/cid/ciag437).
The 2023/24 youngest source stratum is
8 months–8 years, an explicit proxy for the historical 6 months–8 years pool.
VE remains sampled independently of epidemic season; burden-restricted
sensitivities retain the same all-age VE pool. Population weights use the
2025 Census vintage, and maps use 2024 Census boundaries. Candidate vaccination
dates are MMWR weeks 36 through 52 and 1 through 12.

`scripts/validate_publication_inputs.R` enforces these conditions immediately
after input preparation. It writes a generated validation record containing the
run time, software and platform versions, and MD5 hashes of the raw and
processed inputs under `outputs/diagnostics/data/`.

## Inputs, outputs, and provenance

- `data/raw/` and `data/processed/` contain downloaded and normalized inputs.
- `models/` contains reusable fitted GAM objects.
- `outputs/primary/` contains primary tables and final figures.
- `outputs/primary/final_tables/` contains the generated supplementary VE and
  empirical-input tables, their companion CSVs, and the VE source bibliography.
- `outputs/sensitivity/` contains prespecified sensitivity results and their
  cross-analysis summary.
- `outputs/diagnostics/` contains data, model, calendar, held-out, and Monte
  Carlo checks.

Every VE input row records its source URL, source age group, source table,
citation key, and final or preliminary status. `references/ve_sources.bib`
contains the source references for VE and the other model inputs.
`scripts/build_ve_input_table.R` generates
`table_s1_ve_inputs.tex`, `table_s1_ve_inputs.csv`, and `ve_sources.bib` in
`outputs/primary/final_tables/`; the canonical runner calls it automatically.
The publication validator checks the exact season–age pairs, checks saved
draws against their source pools, and requires the supplementary CSV to match
the simulation input. The portrait LaTeX table uses `booktabs` and the
provided bibliography. It can be copied into the manuscript package without
re-entering any estimates.
The bibliography is a plain-text source record; running the analysis and
generating these tables does not require LaTeX or BibTeX.

`scripts/build_parameter_input_table.R` generates a concise portrait table of
empirical inputs, published waning scenarios, and key author-specified assumptions
from the completed run. Fitting and computational settings remain in the Methods,
analysis manifests, and code. Its multipage LaTeX output also requires `longtable`
and `array`. `scripts/build_national_curve_figure.R` displays normalized
national seasonal curves for ILI, laboratory positivity, and ILI+.
`scripts/build_waning_curve_figure.R` displays protection distributions from the
saved primary and age-specific draws and checks its calculations against the
decision engine. Both export PDF/PNG figures and CSV source data. These builders
run automatically after the analyses and can also regenerate the supplementary
artifacts from completed outputs without rerunning simulations. No builder
compiles the manuscript or supplement.

Downloaded inputs, fitted models, generated outputs, and local working
materials are intentionally not version-controlled. A fresh clone retains only
the analysis code, source bibliography, dependency declaration, documentation, and empty directory
skeleton. Public CDC and Census records can be revised. A future run therefore
executes the same prespecified code on the records served at that time, and its
numerical results can change if a source changes. The generated validation
record and per-analysis manifests identify the files used for each run.

The ILI+ sensitivity multiplies national ILINet ILI by national NREVSS
positivity. Its output is a standardized outpatient-burden proxy, not a count
of laboratory-confirmed visits.

Code in this repository is available under the [MIT License](LICENSE).
