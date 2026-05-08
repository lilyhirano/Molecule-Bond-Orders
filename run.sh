#!/bin/sh
set -e

ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"

echo "=== Building molecule_bond_order and calculate_stats ==="
mkdir -p "$ROOT/build"
(
  cd "$ROOT/build"
  cmake .. -DCMAKE_BUILD_TYPE=Debug
  make -j4
)

BIN="$ROOT/build/molecule_bond_order"
STATS_BIN="$ROOT/build/calculate_stats"
for exe in "$BIN" "$STATS_BIN"; do
  if [ ! -x "$exe" ]; then
    echo "Build did not produce executable: $exe" >&2
    exit 1
  fi
done

if [ ! -d "$ROOT/atoms" ]; then
  echo "atoms/ directory not found." >&2
  exit 1
fi

OUT_STD="$ROOT/project_output/final_proj"
OUT_HDF5="$ROOT/project_output/final_proj_hdf5"
RUN_CFG="$ROOT/project_output/run_configs"
STATS_OUT="$ROOT/project_output/stats.stdout"
mkdir -p "$OUT_STD" "$OUT_HDF5" "$RUN_CFG" "$ROOT/project_output"
: >"$STATS_OUT"

ran=0
for xyz in "$ROOT"/atoms/*.xyz; do
  if [ ! -f "$xyz" ]; then
    continue
  fi
  stem="$(basename "$xyz" .xyz)"
  cfg="$ROOT/sample_input/final_proj/${stem}.json"
  if [ ! -f "$cfg" ]; then
    echo "Skipping $stem: no sample_input/final_proj/${stem}.json" >&2
    continue
  fi
  run_cfg="$RUN_CFG/${stem}.json"
  python3 - "$ROOT" "$stem" "$run_cfg" "$cfg" <<'PY'
import json
import sys
from pathlib import Path

root = Path(sys.argv[1])
stem = sys.argv[2]
dst = Path(sys.argv[3])
src = Path(sys.argv[4])
cfg = json.loads(src.read_text(encoding="utf-8"))
cfg["atoms_file_path"] = str((root / "atoms" / f"{stem}.xyz").resolve())
cfg["output_file_path"] = str(
    (root / "project_output/final_proj_hdf5" / f"{stem}.hdf5").resolve()
)
dst.parent.mkdir(parents=True, exist_ok=True)
dst.write_text(json.dumps(cfg, indent=4), encoding="utf-8")
PY
  echo "=== $stem ==="
  "$BIN" "$run_cfg" >"$OUT_STD/${stem}.stdout" 2>&1
  "$STATS_BIN" "$run_cfg" >>"$STATS_OUT" 2>&1
  ran=$((ran + 1))
done

if [ "$ran" -eq 0 ]; then
  echo "No .xyz files processed in atoms/." >&2
  exit 1
fi