#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${repo_root}"

if [[ -d outputs/primary/figures ]]; then
  find outputs/primary/figures -mindepth 1 -type f -delete
fi

if [[ -d outputs/sensitivity ]]; then
  while IFS= read -r figure_dir; do
    find "${figure_dir}" -mindepth 1 -type f -delete
  done < <(find outputs/sensitivity -type d -name figures -print)
fi

for obsolete_figure in \
  outputs/primary/final_figures/figure_3a_state_timing_choropleth.pdf \
  outputs/primary/final_figures/figure_3a_state_timing_choropleth.png
do
  rm -f "${obsolete_figure}"
done

find outputs -type f \( \
  -name ".DS_Store" -o \
  -name "Rplots.pdf" -o \
  -name "*.tmp" -o \
  -name "*.log" \
\) -delete
for output_root in \
  outputs/primary \
  outputs/sensitivity \
  outputs/diagnostics
do
  if [[ -d "${output_root}" ]]; then
    find "${output_root}" -mindepth 1 -type d -empty -delete
  fi
done

if [[ -d data/processed ]]; then
  find data/processed -type f -name ".DS_Store" -delete
fi

for obsolete_artifact in \
  data/processed/ilinet_state_reconstructed_age_primary_seasons.csv \
  data/processed/ilinet_state_age_composition_primary_seasons.csv \
  outputs/diagnostics/data/ilinet_state_age_coverage_by_season.csv
do
  rm -f "${obsolete_artifact}"
done

echo "Pruned per-analysis and legacy figure files, temporary files, logs, obsolete state-age artifacts, and empty subdirectories."
