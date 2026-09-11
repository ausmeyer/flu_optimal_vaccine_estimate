#!/usr/bin/env bash
set -euo pipefail

: "${CLEAN_GENERATED_OUTPUTS:=true}"
: "${REFRESH_INPUTS:=true}"

# Canonical scientific configuration. These assignments are unconditional so
# inherited shell variables cannot alter a publication run. Only the operational
# input-refresh and cleanup switches above are user-configurable here.
N_DRAWS=5000
PRIMARY_SEASON_START_YEARS="2010:2018,2023:2025"
EXCLUDED_SEASON_START_YEARS="2019:2022"
STATE_POPULATION_YEAR="2025"
TIMING_SHIFT_SD="0.75"
VE_REFERENCE_WEEKS="8"
WANING_PRIOR="ray_2019"
IMMUNE_LAG_PRIOR="discrete_about_2w"
GAM_GLOBAL_K="14"
GAM_SEASON_K="12"
GAM_GAMMA="1"
OUTCOME_DATE_PRIOR="none"
SAMPLE_COMPLETE_SEASONS_ONLY="true"
AGE_GROUPS=""
WRITE_BURDEN_DRAWS="false"
WRITE_ANALYSIS_FIGURES="false"
SENSITIVITY_ROOT="outputs/sensitivity"
VALIDATION_OUTPUT_DIR="outputs/diagnostics/held_out_validation"
NREVSS_ANALYSIS_REGION="national"
NREVSS_TARGET="all"
ILI_PLUS_ILI_WEIGHTING="cdc_population_weighted"
ILI_PLUS_LAB_SCOPE="all_available"
CDC_FLUVIEW_API_BASE="https://gis.cdc.gov/flu2"

export N_DRAWS PRIMARY_SEASON_START_YEARS EXCLUDED_SEASON_START_YEARS
export STATE_POPULATION_YEAR TIMING_SHIFT_SD VE_REFERENCE_WEEKS
export WANING_PRIOR IMMUNE_LAG_PRIOR GAM_GLOBAL_K GAM_SEASON_K GAM_GAMMA
export OUTCOME_DATE_PRIOR SAMPLE_COMPLETE_SEASONS_ONLY AGE_GROUPS
export WRITE_BURDEN_DRAWS WRITE_ANALYSIS_FIGURES
export SENSITIVITY_ROOT VALIDATION_OUTPUT_DIR NREVSS_ANALYSIS_REGION NREVSS_TARGET
export ILI_PLUS_ILI_WEIGHTING ILI_PLUS_LAB_SCOPE CDC_FLUVIEW_API_BASE

publication_backup_dir=""
publication_backup_ready=false
publication_outputs_valid=false

is_true() {
  case "$1" in
    1|true|TRUE|yes|YES|y|Y) return 0 ;;
    *) return 1 ;;
  esac
}

finish_publication_run() {
  local status=$?
  trap - EXIT

  if is_true "${publication_backup_ready}" &&
      [[ "${status}" -ne 0 ]] &&
      ! is_true "${publication_outputs_valid}"; then
    echo "Publication pipeline failed; restoring the previous outputs."
    mkdir -p outputs
    find outputs -mindepth 1 -depth -delete
    cp -a "${publication_backup_dir}/outputs/." outputs/
  elif is_true "${publication_backup_ready}" &&
      [[ "${status}" -ne 0 ]]; then
    echo "Publication outputs passed validation and will be retained."
  fi

  if [[ -n "${publication_backup_dir}" && -d "${publication_backup_dir}" ]]; then
    find "${publication_backup_dir}" -depth -delete
  fi
  exit "${status}"
}
trap finish_publication_run EXIT

echo "Running canonical publication pipeline with N_DRAWS=${N_DRAWS}"

if is_true "${CLEAN_GENERATED_OUTPUTS}"; then
  publication_backup_dir="$(mktemp -d)"
  mkdir -p "${publication_backup_dir}/outputs"
  if [[ -d outputs ]]; then
    cp -a outputs/. "${publication_backup_dir}/outputs/"
  fi
  publication_backup_ready=true
  scripts/clean_generated_outputs.sh
fi

if is_true "${REFRESH_INPUTS}"; then
  REDOWNLOAD_ILINET=true Rscript scripts/download_ilinet.R
  REDOWNLOAD_NREVSS=true NREVSS_REGIONS=national \
    Rscript scripts/download_nrevss_clinical_labs.R
  Rscript scripts/download_state_population.R
  Rscript scripts/scrape_cdc_ve_estimates.R
  Rscript scripts/normalize_cdc_ve_estimates.R
  Rscript scripts/update_state_boundaries.R
else
  REDOWNLOAD_ILINET=false Rscript scripts/download_ilinet.R
  REDOWNLOAD_NREVSS=false NREVSS_REGIONS=national \
    Rscript scripts/download_nrevss_clinical_labs.R
  if [[ ! -f data/processed/state_population.csv ]]; then
    Rscript scripts/download_state_population.R
  fi
  if [[ ! -f data/raw/cdc_ve_raw_tables.rds ]]; then
    Rscript scripts/scrape_cdc_ve_estimates.R
  fi
  Rscript scripts/normalize_cdc_ve_estimates.R
  if [[ ! -f data/processed/us_states_2024_cartographic_boundaries.rds ]]; then
    Rscript scripts/update_state_boundaries.R
  fi
fi

Rscript scripts/validate_publication_inputs.R
Rscript scripts/check_model_components.R
Rscript scripts/build_ve_input_table.R

OUTPUT_CLASS=primary \
  OUTPUT_PREFIX=gam_primary \
  REUSE_GAM_MODELS=true \
  Rscript scripts/run_gam_primary_analysis.R

# A fresh R process releases each age group's working memory before the next.
for age_group in "0-4" "5-24" "25-49" "50-64" "65+"; do
  AGE_GROUPS="${age_group}" \
    OUTPUT_CLASS=primary \
    OUTPUT_PREFIX=gam_age \
    REUSE_GAM_MODELS=true \
    Rscript scripts/run_gam_age_specific_analysis.R
done
Rscript scripts/assemble_age_specific_intervals.R

SKIP_EXISTING=false scripts/run_manuscript_sensitivity_suite.sh

Rscript scripts/audit_calendar_guardrails.R
Rscript scripts/build_results_figures.R
Rscript scripts/build_robustness_figure.R
Rscript scripts/build_supplement_figures.R
Rscript scripts/build_national_curve_figure.R
Rscript scripts/build_waning_curve_figure.R
Rscript scripts/build_parameter_input_table.R

VALIDATE_ARTIFACT_HYGIENE=false \
  Rscript scripts/validate_publication_outputs.R
scripts/prune_publication_artifacts.sh
VALIDATE_ARTIFACT_HYGIENE=true \
  Rscript scripts/validate_publication_outputs.R
publication_outputs_valid=true
Rscript scripts/prune_unreferenced_models.R

echo "Publication pipeline complete."
