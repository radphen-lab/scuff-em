#!/bin/bash
set -euo pipefail

# Support source, out-of-tree build, and installed test directories.
SOURCE_DIR=$(cd -- "${srcdir:-$(dirname -- "$0")}" && pwd)
CODE=${SCUFF_AVESCATTER:-scuff-avescatter}
if ! CODE=$(command -v "$CODE"); then
  echo "Could not find scuff-avescatter. Set SCUFF_AVESCATTER or check PATH." >&2
  exit 1
fi
CODE=$(realpath "$CODE")
PYTHON=${PYTHON:-python3}
if ! command -v "$PYTHON" >/dev/null; then
  echo "SKIP AveScatter requires Python 3." >&2
  exit 77
fi
THREADS=${OMP_NUM_THREADS:-4}
export OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-4}
RESULTS=$(mktemp -d "${TMPDIR:-/tmp}/TestAveScatter.XXXXXX")
echo "AveScatter results in $RESULTS"
cd "$RESULTS"

for GEOM in Oblate BiSphere CoreShell LossySphere; do
  ARGS=(--geometry "$SOURCE_DIR/$GEOM.scuffgeo" --OmegaFile "$SOURCE_DIR/OmegaFile")
  if [ "$GEOM" = Oblate ]; then
    ARGS=(--geometry "$SOURCE_DIR/$GEOM.scuffgeo" --OmegaFile "$SOURCE_DIR/Oblate.OmegaFile")
  fi
  if [ "$GEOM" = BiSphere ]; then
    ARGS+=(--TransFile "$SOURCE_DIR/BiSphere.trans")
  fi
  "$CODE" "${ARGS[@]}" --Threads "$THREADS" --FileBase "$GEOM" </dev/null
  "$PYTHON" "$SOURCE_DIR/CheckAveScatter.py" \
    --data "$GEOM.SIFlux.AVSCAT" \
    --reference "$SOURCE_DIR/reference/$GEOM.SIFlux.AVSCAT" \
    --checklist "$SOURCE_DIR/Checklist.SIFlux"
done
