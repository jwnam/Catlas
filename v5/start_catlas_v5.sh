#!/usr/bin/env bash
set -euo pipefail

umask 077

readonly default_project_dir="/home/minwook/Shiny_CRC_atlas"
readonly project_dir="${CATLAS_PROJECT_DIR:-$default_project_dir}"
readonly app_dir="${CATLAS_APP_DIR:-$project_dir/Catlas/v5}"
readonly data_path="${CATLAS_DATA_PATH:-$project_dir/crc_shiny_app_seurat.rds}"
readonly tmp_root="${CATLAS_TMP_ROOT:-$app_dir/tmp_for_catlas}"
readonly catlas_host="${CATLAS_HOST:-127.0.0.1}"
readonly catlas_port="${CATLAS_PORT:-4511}"

if [[ ! "$catlas_port" =~ ^[0-9]+$ ]] ||
  ((catlas_port < 1 || catlas_port > 65535)); then
  echo "Invalid CATLAS_PORT: $catlas_port" >&2
  exit 64
fi

if [[ ! -d "$app_dir" ]]; then
  echo "CRC Atlas v5 app directory does not exist: $app_dir" >&2
  exit 66
fi

if [[ ! -r "$app_dir/app.R" ]]; then
  echo "CRC Atlas v5 app entrypoint is not readable: $app_dir/app.R" >&2
  exit 66
fi

if [[ ! -r "$data_path" ]]; then
  echo "CRC Atlas input RDS is not readable: $data_path" >&2
  exit 66
fi

mkdir -p -- "$tmp_root"
chmod 0700 -- "$tmp_root"

if [[ ! -d "$tmp_root" || ! -w "$tmp_root" || ! -x "$tmp_root" ]]; then
  echo "CRC Atlas temporary root is not writable/searchable: $tmp_root" >&2
  exit 73
fi

if [[ -n "${CATLAS_R_BIN:-}" ]]; then
  r_bin="$CATLAS_R_BIN"
elif [[ -x /home/minwook/miniconda3/envs/crc_shiny/bin/R ]]; then
  r_bin="/home/minwook/miniconda3/envs/crc_shiny/bin/R"
else
  echo "No CRC Shiny R executable found; set CATLAS_R_BIN." >&2
  exit 69
fi

if [[ ! -x "$r_bin" ]]; then
  echo "CATLAS_R_BIN is not executable: $r_bin" >&2
  exit 69
fi

export CATLAS_PROJECT_DIR="$project_dir"
export CATLAS_APP_DIR="$app_dir"
export CATLAS_DATA_PATH="$data_path"
export CATLAS_TMP_ROOT="$tmp_root"
export CATLAS_HOST="$catlas_host"
export CATLAS_PORT="$catlas_port"
export TMPDIR="$tmp_root"
export TMP="$tmp_root"
export TEMP="$tmp_root"

cd -- "$app_dir"

echo "Starting CRC Atlas v5 with R: $r_bin" >&2
echo "TMPDIR=$TMPDIR; app_dir=$CATLAS_APP_DIR; host=$CATLAS_HOST; port=$CATLAS_PORT" >&2

exec "$r_bin" --no-echo --no-restore --no-save -e \
  'shiny::runApp(appDir = Sys.getenv("CATLAS_APP_DIR"), host = Sys.getenv("CATLAS_HOST"), port = as.integer(Sys.getenv("CATLAS_PORT")), launch.browser = FALSE)'
