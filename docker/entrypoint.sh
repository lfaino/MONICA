#!/bin/bash
# Entry point of the monica image: run the pipeline, or an arbitrary command
# when the first argument is `bash`, `sh` or `python`.
set -euo pipefail

case "${1:-}" in
    bash|sh|python|blastn|ktImportTaxonomy) exec "$@" ;;
esac

exec python "${MONICA_SCRIPT:-/opt/monica/monica_v5.py}" "$@"
