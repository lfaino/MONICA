# MONICA with its web interface, in one image.
#
#   docker build -t monica:web -f docker/web/Dockerfile .
#
# Three mount points, nothing else to configure:
#   /data   reads (fastq) - read only
#   /db     blast databases + ete4/Krona taxonomies
#   /jobs   one directory per analysis, with its log and its results
ARG BASE=monica:v5
FROM ${BASE}

RUN pip install --no-cache-dir \
        fastapi "uvicorn[standard]" jinja2 python-multipart \
    && rm -rf /root/.cache/pip

COPY monica_v5.py /opt/monica/monica_v5.py
COPY docker/web/app.py /opt/web/app.py
COPY docker/web/templates /opt/web/templates

ENV MONICA_SCRIPT=/opt/monica/monica_v5.py \
    MONICA_JOBS=/jobs \
    MONICA_SEQ_ROOT=/data \
    BLASTDB=/db \
    XDG_DATA_HOME=/db \
    ETE_DATA_HOME=/db \
    MPLBACKEND=Agg
RUN mkdir -p /data /db /jobs

EXPOSE 8000
HEALTHCHECK --interval=30s --timeout=5s --start-period=20s \
    CMD python -c "import urllib.request;urllib.request.urlopen('http://localhost:8000/health')"
ENTRYPOINT ["uvicorn", "app:app", "--app-dir", "/opt/web", "--host", "0.0.0.0", "--port", "8000"]
