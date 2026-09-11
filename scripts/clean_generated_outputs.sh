#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${repo_root}"

if [[ -d outputs ]]; then
  find outputs -mindepth 1 -type f \
    ! -path "outputs/.gitkeep" \
    -delete
  find outputs -mindepth 1 -type d -empty -delete
fi

echo "Removed all generated output artifacts while preserving outputs/.gitkeep."
