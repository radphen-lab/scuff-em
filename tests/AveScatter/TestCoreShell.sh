#!/bin/bash
set -euo pipefail

SOURCE_DIR=$(cd -- "${srcdir:-$(dirname -- "$0")}" && pwd)
CODE=${SCUFF_AVESCATTER:-scuff-avescatter}
if ! CODE=$(command -v "$CODE"); then
  echo "Could not find scuff-avescatter. Set SCUFF_AVESCATTER or check PATH." >&2
  exit 1
fi
CODE=$(realpath "$CODE")
PYTHON=${PYTHON:-python3}
if ! "$PYTHON" -c 'import numpy; import scipy' >/dev/null 2>&1; then
  echo "SKIP CoreShell requires Python 3, NumPy, and SciPy." >&2
  exit 77
fi
THREADS=${OMP_NUM_THREADS:-4}
export OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-4}
RESULTS=$(mktemp -d "${TMPDIR:-/tmp}/TestCoreShell.XXXXXX")
echo "CoreShell results in $RESULTS"
cd "$RESULTS"

"$CODE" --geometry "$SOURCE_DIR/CoreShellMie.scuffgeo" \
  --OmegaFile "$SOURCE_DIR/OmegaFile" --Threads "$THREADS" \
  --FileBase CoreShellMie </dev/null
"$PYTHON" "$SOURCE_DIR/ValidateCoreShell.py" CoreShellMie.SIFlux.AVSCAT --json comparison.json
