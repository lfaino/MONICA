# Runtime image for monica_v4.py (Nanopore amplicon pathogen detection).
#
#   docker build -t monica:v4 -f docker/Dockerfile .
#
# The pipeline needs an external blast database and, optionally, the NCBI
# taxonomy dump used by ete4 and Krona: mount them on /mnt/dataQNAP/blastdb.
FROM condaforge/miniforge3:24.11.3-2

LABEL org.opencontainers.image.title="monica" \
      org.opencontainers.image.description="Nanopore amplicon pathogen detection and quantification"

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        ca-certificates \
        procps \
        wget \
    && rm -rf /var/lib/apt/lists/*

# --- conda environments ----------------------------------------------------
COPY docker/env_tools.yml /tmp/env_tools.yml
RUN mamba env create -f /tmp/env_tools.yml && mamba clean -afy

COPY docker/env_monica.yml /tmp/env_monica.yml
RUN mamba env create -f /tmp/env_monica.yml && mamba clean -afy

# The script calls `spades`, while recent bioconda builds only ship `spades.py`.
RUN ln -s /opt/conda/envs/tools/bin/spades.py /opt/conda/envs/tools/bin/spades

# The script carries the shebang of the original deployment: keep it valid.
RUN mkdir -p /opt/conda_envs && ln -s /opt/conda/envs/monica /opt/conda_envs/ampliseq_wooey

# The Python environment comes first: the tools only need to be on the PATH.
ENV PATH=/opt/conda/envs/monica/bin:/opt/conda/envs/tools/bin:$PATH

# --- paths hardcoded in the script ------------------------------------------
# BLASTDB and ETE_DATA_HOME are also set by the script itself; the mount points
# are created here so that a read-only bind mount always succeeds.
# ete4 >= 4 resolves its data directory from XDG_DATA_HOME (the ETE_DATA_HOME
# variable set by the script is read too late to have any effect), so the
# taxonomy ends up in <XDG_DATA_HOME>/ete, i.e. on the same mounted volume.
ENV BLASTDB=/mnt/dataQNAP/blastdb/ \
    ETE_DATA_HOME=/mnt/dataQNAP/blastdb/ \
    XDG_DATA_HOME=/mnt/dataQNAP/blastdb \
    MPLBACKEND=Agg
RUN mkdir -p /mnt/dataQNAP/blastdb \
             /mnt/dataZEUSPV/wooey/tmp \
             /mnt/data3ZEUSPV/SequenceData \
             /etc/slurm \
             /data

# --- pipeline ----------------------------------------------------------------
COPY monica_v4.py /opt/monica/monica_v4.py
COPY monica_v5.py /opt/monica/monica_v5.py
COPY docker/entrypoint.sh /usr/local/bin/entrypoint.sh
RUN chmod +x /opt/monica/monica_v4.py /opt/monica/monica_v5.py /usr/local/bin/entrypoint.sh

WORKDIR /data
ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
CMD ["--help"]
