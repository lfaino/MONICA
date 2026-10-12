# Web interface

A single FastAPI process: it serves the form, queues the analyses and runs them
one at a time as sub-processes of the container. No broker, no database - each
analysis is a directory under `/jobs` holding `state.json`, `monica.log` and
every file the pipeline wrote.

The advanced half of the form is generated from the pipeline's own argument
parser at request time (`monica_v5.build_parser()`), so an option added to
`monica_v5.py` appears in the page with no registration step and no duplicated
list of parameters.

Environment variables: `MONICA_SCRIPT`, `MONICA_JOBS`, `MONICA_SEQ_ROOT`,
`BLASTDB`, `MONICA_DEFAULT_EMAIL`, `MONICA_DEFAULT_THREADS`, `MONICA_SITE_NAME`,
`MONICA_SCAN_DEPTH`.

Endpoints: `/` form, `/jobs` list, `/jobs/{id}` status, log and results,
`/jobs/{id}/files/{name}` download, `/jobs/{id}/log` full log,
`/jobs/{id}/cancel`, `/jobs/{id}/delete`, `/health`.

There is no authentication: it is meant to be reachable from the machine it runs
on (`localhost`) or from a trusted network. Put a reverse proxy in front of it
before exposing it more widely.
