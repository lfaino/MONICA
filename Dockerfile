# MONICA - detection and quantification of pathogens from Nanopore amplicons.
# Everything in one image: the pipeline, the command line tools it calls, and
# the web interface. This is the image published to the registry.
#
#   docker build -t monica .
#   docker run -d -p 8765:8000 -v $HOME/MONICA/reads:/data:ro \
#              -v $HOME/MONICA/database:/db -v $HOME/MONICA/results:/jobs monica
#
# Mount points:
#   /data   reads (fastq), read only
#   /db     blast databases and the ete4/Krona taxonomies
#   /jobs   one directory per analysis: parameters, log and results
FROM condaforge/miniforge3:24.11.3-2

LABEL org.opencontainers.image.title="MONICA" \
      org.opencontainers.image.description="Detection and quantification of pathogens from Nanopore amplicon reads" \
      org.opencontainers.image.licenses="SEE LICENSE FILE"

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential ca-certificates procps wget \
    && rm -rf /var/lib/apt/lists/*

# --- external tools and python runtime, in two environments ------------------
# (the pinned interpreters of the bioconda packages clash with the analysis stack)
COPY docker/env_tools_multiarch.yml /tmp/env_tools.yml
RUN mamba env create -f /tmp/env_tools.yml && mamba clean -afy

COPY docker/env_monica.yml /tmp/env_monica.yml
RUN mamba env create -f /tmp/env_monica.yml && mamba clean -afy

RUN pip_bin=/opt/conda/envs/monica/bin/pip \
    && $pip_bin install --no-cache-dir fastapi "uvicorn[standard]" jinja2 python-multipart

# Recent bioconda builds ship only `spades.py`; the pipeline calls `spades`.
RUN ln -sf /opt/conda/envs/tools/bin/spades.py /opt/conda/envs/tools/bin/spades

ENV PATH=/opt/conda/envs/monica/bin:/opt/conda/envs/tools/bin:$PATH \
    MONICA_SCRIPT=/opt/monica/monica_v5.py \
    MONICA_JOBS=/jobs \
    MONICA_SEQ_ROOT=/data \
    BLASTDB=/db \
    XDG_DATA_HOME=/db \
    ETE_DATA_HOME=/db \
    MPLBACKEND=Agg

COPY monica_v5.py /opt/monica/monica_v5.py
COPY docker/web/app.py /opt/web/app.py
COPY docker/web/templates /opt/web/templates
RUN mkdir -p /data /db /jobs && chmod +x /opt/monica/monica_v5.py

EXPOSE 8000
HEALTHCHECK --interval=30s --timeout=5s --start-period=20s \
    CMD python -c "import urllib.request;urllib.request.urlopen('http://localhost:8000/health')"
ENTRYPOINT ["uvicorn", "app:app", "--app-dir", "/opt/web", "--host", "0.0.0.0", "--port", "8000"]
