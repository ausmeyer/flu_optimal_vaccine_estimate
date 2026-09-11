#!/usr/bin/env bash
set -euo pipefail

: "${N_DRAWS:=5000}"
: "${SKIP_EXISTING:=false}"
: "${WRITE_ANALYSIS_FIGURES:=false}"
export N_DRAWS
export WRITE_ANALYSIS_FIGURES

is_true() {
  case "$1" in
    1|true|TRUE|yes|YES|y|Y) return 0 ;;
    *) return 1 ;;
  esac
}

run_analysis() {
  local prefix="$1"
  shift
  local manifest="outputs/sensitivity/${prefix}/tables/${prefix}_analysis_manifest.json"

  if is_true "${SKIP_EXISTING}" && [[ -f "${manifest}" ]]; then
    echo "Skipping ${prefix}; a completed manifest is present."
    return
  fi

  echo "Running ${prefix}"
  env "$@"
}

echo "Running manuscript sensitivities with N_DRAWS=${N_DRAWS}"
echo "Existing completed analyses will be skipped: ${SKIP_EXISTING}"
if is_true "${SKIP_EXISTING}"; then
  echo "Warning: SKIP_EXISTING=true is an unchecked development shortcut; canonical runs set it to false."
fi
echo "Per-analysis exploratory figures will be written: ${WRITE_ANALYSIS_FIGURES}"

run_analysis empirical_primary \
  OUTPUT_CLASS=sensitivity \
  OUTPUT_PREFIX=empirical_primary \
  Rscript scripts/run_primary_analysis.R

run_analysis gam_no_timing_shift_primary \
  TIMING_SHIFT_SD=0 \
  OUTPUT_CLASS=sensitivity \
  OUTPUT_PREFIX=gam_no_timing_shift_primary \
  REUSE_GAM_MODELS=true \
  Rscript scripts/run_gam_primary_analysis.R

run_analysis gam_precovid_primary \
  PRIMARY_SEASON_START_YEARS=2010:2018 \
  OUTPUT_CLASS=sensitivity \
  OUTPUT_PREFIX=gam_precovid_primary \
  REUSE_GAM_MODELS=true \
  Rscript scripts/run_gam_primary_analysis.R

run_analysis gam_fixed_immune_lag_primary \
  IMMUNE_LAG_PRIOR=fixed_2w \
  OUTPUT_CLASS=sensitivity \
  OUTPUT_PREFIX=gam_fixed_immune_lag_primary \
  REUSE_GAM_MODELS=true \
  Rscript scripts/run_gam_primary_analysis.R

run_analysis gam_slow_waning_primary \
  WANING_PRIOR=ray_2019_lower95 \
  OUTPUT_CLASS=sensitivity \
  OUTPUT_PREFIX=gam_slow_waning_primary \
  REUSE_GAM_MODELS=true \
  Rscript scripts/run_gam_primary_analysis.R

run_analysis gam_fast_waning_primary \
  WANING_PRIOR=ray_2019_upper95 \
  OUTPUT_CLASS=sensitivity \
  OUTPUT_PREFIX=gam_fast_waning_primary \
  REUSE_GAM_MODELS=true \
  Rscript scripts/run_gam_primary_analysis.R

run_analysis gam_spencer_ferdinands_fast_waning \
  WANING_PRIOR=spencer_ferdinands_fast \
  OUTPUT_CLASS=sensitivity \
  OUTPUT_PREFIX=gam_spencer_ferdinands_fast_waning \
  REUSE_GAM_MODELS=true \
  Rscript scripts/run_gam_primary_analysis.R

run_analysis gam_spencer_ferdinands_onset_aligned \
  WANING_PRIOR=spencer_ferdinands_fast \
  OUTCOME_DATE_PRIOR=flu_ve_outpatient_onset_to_enrollment \
  OUTPUT_CLASS=sensitivity \
  OUTPUT_PREFIX=gam_spencer_ferdinands_onset_aligned \
  REUSE_GAM_MODELS=true \
  Rscript scripts/run_gam_primary_analysis.R

run_analysis gam_spencer_slow_waning \
  WANING_PRIOR=spencer_slow \
  OUTPUT_CLASS=sensitivity \
  OUTPUT_PREFIX=gam_spencer_slow_waning \
  REUSE_GAM_MODELS=true \
  Rscript scripts/run_gam_primary_analysis.R

run_analysis gam_smooth_k16_k14_primary \
  GAM_GLOBAL_K=16 \
  GAM_SEASON_K=14 \
  OUTPUT_CLASS=sensitivity \
  OUTPUT_PREFIX=gam_smooth_k16_k14_primary \
  REUSE_GAM_MODELS=true \
  Rscript scripts/run_gam_primary_analysis.R

run_analysis gam_ve_reference_4w \
  VE_REFERENCE_WEEKS=4 \
  OUTPUT_CLASS=sensitivity \
  OUTPUT_PREFIX=gam_ve_reference_4w \
  REUSE_GAM_MODELS=true \
  Rscript scripts/run_gam_primary_analysis.R

run_analysis gam_ve_reference_12w \
  VE_REFERENCE_WEEKS=12 \
  OUTPUT_CLASS=sensitivity \
  OUTPUT_PREFIX=gam_ve_reference_12w \
  REUSE_GAM_MODELS=true \
  Rscript scripts/run_gam_primary_analysis.R

nrevss_full_manifest="outputs/sensitivity/gam_nrevss_national/tables/gam_nrevss_national_analysis_manifest.json"
nrevss_a_manifest="outputs/sensitivity/gam_nrevss_post2015_flu_a/tables/gam_nrevss_post2015_flu_a_analysis_manifest.json"
if ! is_true "${SKIP_EXISTING}" || [[ ! -f "${nrevss_full_manifest}" || ! -f "${nrevss_a_manifest}" ]]; then
  echo "Preparing WHO/NREVSS clinical laboratory data"
  PRIMARY_SEASON_START_YEARS=2010:2018,2023:2025 \
    NREVSS_REGIONS=national \
    Rscript scripts/download_nrevss_clinical_labs.R
fi

echo "Building aligned national ILI+ input"
Rscript scripts/build_national_ili_plus.R

run_analysis gam_nrevss_national \
  PRIMARY_SEASON_START_YEARS=2010:2018,2023:2025 \
  NREVSS_ANALYSIS_REGION=national \
  NREVSS_TARGET=all \
  OUTPUT_CLASS=sensitivity \
  OUTPUT_PREFIX=gam_nrevss_national \
  REUSE_GAM_MODELS=true \
  Rscript scripts/run_nrevss_clinical_lab_sensitivity.R

run_analysis gam_nrevss_post2015_flu_a \
  PRIMARY_SEASON_START_YEARS=2015:2018,2023:2025 \
  NREVSS_ANALYSIS_REGION=national \
  NREVSS_TARGET=a \
  OUTPUT_CLASS=sensitivity \
  OUTPUT_PREFIX=gam_nrevss_post2015_flu_a \
  REUSE_GAM_MODELS=true \
  Rscript scripts/run_nrevss_clinical_lab_sensitivity.R

run_analysis empirical_ili_plus_all_seasons_weighted \
  PRIMARY_SEASON_START_YEARS=2010:2018,2023:2025 \
  ILI_PLUS_ILI_WEIGHTING=cdc_population_weighted \
  ILI_PLUS_LAB_SCOPE=all_available \
  OUTPUT_CLASS=sensitivity \
  OUTPUT_PREFIX=empirical_ili_plus_all_seasons_weighted \
  Rscript scripts/run_ili_plus_sensitivity.R

run_analysis empirical_ili_plus_all_seasons_unweighted \
  PRIMARY_SEASON_START_YEARS=2010:2018,2023:2025 \
  ILI_PLUS_ILI_WEIGHTING=national_unweighted \
  ILI_PLUS_LAB_SCOPE=all_available \
  OUTPUT_CLASS=sensitivity \
  OUTPUT_PREFIX=empirical_ili_plus_all_seasons_unweighted \
  Rscript scripts/run_ili_plus_sensitivity.R

run_analysis empirical_ili_plus_clinical_only_weighted \
  PRIMARY_SEASON_START_YEARS=2016:2018,2023:2025 \
  ILI_PLUS_ILI_WEIGHTING=cdc_population_weighted \
  ILI_PLUS_LAB_SCOPE=clinical_only \
  OUTPUT_CLASS=sensitivity \
  OUTPUT_PREFIX=empirical_ili_plus_clinical_only_weighted \
  Rscript scripts/run_ili_plus_sensitivity.R

echo "Running leave-one-season-out stability analyses"
OUTPUT_CLASS=sensitivity \
  SKIP_EXISTING="${SKIP_EXISTING}" \
  Rscript scripts/run_gam_leave_one_season_out_sensitivity.R

echo "Running held-out-season decision validation"
Rscript scripts/run_gam_held_out_validation.R

Rscript scripts/audit_gam_fit_diagnostics.R
Rscript scripts/audit_monte_carlo_precision.R
Rscript scripts/summarise_robustness_results.R

echo "Manuscript sensitivity suite complete."
