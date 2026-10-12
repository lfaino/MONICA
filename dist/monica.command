#!/bin/bash
# MONICA - double-click to start. macOS and Linux.
set -u
IMAGE="${MONICA_IMAGE:-ghcr.io/lfaino/monica:latest}"
PORT="${MONICA_PORT:-8765}"
HOME_DIR="${MONICA_HOME:-$HOME/MONICA}"

cd "$(dirname "$0")"
echo "MONICA - avvio in corso"

if ! command -v docker >/dev/null 2>&1; then
    echo
    echo "Docker non è installato."
    echo "Scaricalo da https://www.docker.com/products/docker-desktop/ , installalo,"
    echo "avvialo e poi riapri questo file."
    read -r -p "Premi invio per chiudere." _
    exit 1
fi
if ! docker info >/dev/null 2>&1; then
    echo "Docker è installato ma non è in esecuzione: avvia Docker Desktop e riprova."
    read -r -p "Premi invio per chiudere." _
    exit 1
fi

mkdir -p "$HOME_DIR/reads" "$HOME_DIR/database" "$HOME_DIR/results"
echo "I tuoi dati stanno in $HOME_DIR"
echo "  reads/     metti qui le cartelle fastq"
echo "  database/  database blast e tassonomia"
echo "  results/   risultati delle analisi"

docker rm -f monica >/dev/null 2>&1
echo "Scarico l'immagine (solo la prima volta, qualche minuto)..."
docker pull "$IMAGE" || { echo "Download fallito."; read -r -p "Invio per chiudere." _; exit 1; }

docker run -d --name monica -p "$PORT:8000" \
    -v "$HOME_DIR/reads:/data:ro" \
    -v "$HOME_DIR/database:/db" \
    -v "$HOME_DIR/results:/jobs" \
    "$IMAGE" >/dev/null || { echo "Avvio fallito."; read -r -p "Invio per chiudere." _; exit 1; }

for _ in $(seq 1 60); do
    if curl -fs "http://localhost:$PORT/health" >/dev/null 2>&1; then break; fi
    sleep 1
done
URL="http://localhost:$PORT/"
echo "Pronto: $URL"
if command -v open >/dev/null 2>&1; then open "$URL"; elif command -v xdg-open >/dev/null 2>&1; then xdg-open "$URL"; fi
echo
echo "MONICA resta attivo in background. Per fermarlo:  docker stop monica"
read -r -p "Premi invio per chiudere questa finestra." _
