# MONICA

Detection and quantification of pathogens from Nanopore amplicon reads.

MONICA takes the basecalled reads of a sequencing run, filters and trims them,
blasts them against a nucleotide database, assigns each read to a taxon with the
NCBI taxonomy, and reports the composition of every barcode as a table, a Krona
chart, a stacked histogram, a clustermap and - on request - alpha/beta
diversity.

It runs in three ways, and they all use the same code:

| | for whom |
|---|---|
| **Web interface in Docker** | anyone, on Windows, macOS or Linux - no installation beyond Docker |
| **Command line** | `python monica_v5.py --barcode barcode01-12 ...` |
| **SLURM cluster** | the blast step is submitted as a single job array |

---

## 1. Use it from the browser (recommended)

**Step 1 - install Docker.** Download [Docker
Desktop](https://www.docker.com/products/docker-desktop/), install it, start it.
This is the only installation step, and it is the same on every operating
system.

**Step 2 - start MONICA.** Download this repository (green *Code* button → *Download ZIP*), unzip it, then:

* **Windows**: double-click `dist/monica.cmd`
* **macOS / Linux**: double-click `dist/monica.command`

The first start downloads the image (about 1.5 GB) and may take several minutes.
The browser then opens on <http://localhost:8765/>.

**Step 3 - put your data where MONICA can see it.** The launcher creates a
`MONICA` folder in your home directory:

```
MONICA/
  reads/      <- copy here the folder of your sequencing run (the one with barcode01, barcode02, ...)
  database/   <- blast databases and the taxonomy
  results/    <- one folder per analysis, with the log and all result files
```

Choose the run from the drop-down, tick the barcodes, pick the database, press
*Avvia l'analisi*. The page follows the run and offers every result file for
download when it is done.

To stop it: `docker stop monica`. Your results stay in `MONICA/results`.

### What you need in `database/`

* **the NCBI taxonomy** used to name the taxa - about 800 MB, built once:

  ```bash
  docker exec monica python -c "from ete4 import NCBITaxa; NCBITaxa()"
  docker exec monica ktUpdateTaxonomy.sh /db
  ```

* **at least one blast database**. Either copy yours into `MONICA/database/`
  (it then appears in the drop-down), or fetch a public one, for example the
  16S ribosomal RNA database:

  ```bash
  docker exec -w /db monica update_blastdb.pl --decompress 16S_ribosomal_RNA
  ```

Without a taxonomy the run stops immediately and says so; it does not waste the
blast.

---

## 2. Use it from the command line

```bash
python monica_v5.py \
    --barcode barcode01-12 \
    --folder_fastq /path/to/fastq_pass \
    -d nt.nal \
    -e you@institute.it \
    -t 16 \
    -o my_run
```

`--help` lists every option. The ones worth knowing:

| option | what it does |
|---|---|
| `--barcode` | `barcode01-12`, or a comma-separated list |
| `-d` / `-dbe` | a built-in database name, or the path of your own |
| `--chunk_size` | reads per blast call (default 1000); the single lever on blast speed |
| `--slurm` | submit the blast chunks as one SLURM job array |
| `--taxonomy_db` | path of `taxa.sqlite`, if it is not where ete4 looks |
| `--lca`, `--blast_filter` | the two assignment rules that can be changed; defaults reproduce version 4 |
| `--verbose` | keeps the temporary files and the SLURM logs next to the results |

Results land in the current directory: `<name>.xlsx` (every species),
`<name>.html` (Krona), `<name>.histogram.pdf`, `<name>.heatmap.pdf`,
`<sample>_rank.txt`, `<sample>_blastn_results.txt` and
`<name>.run_manifest.json`, which records the parameters, the database and the
version of every tool used.

Requirements, if you do not use Docker: python ≥ 3.10 with biopython, pandas,
matplotlib, seaborn, scikit-bio, ete4, simple_slurm, tqdm, openpyxl, plus
minimap2, racon, porechop, jellyfish, SPAdes, nanopolish, bwa, samtools,
bcftools, BLAST+ and Krona on the PATH. `docker/env_*.yml` list exact versions.

---

## 3. Run it on a cluster

`--slurm` submits the blast chunks as a single job array instead of one job per
chunk. The temporary directory must be visible from the compute nodes (`--temp`),
and the run stops with an explanation if it is not. The state of every task is
read from `sacct`, and the logs of the failed ones are copied next to the
results.

---

## Repository layout

```
monica_v5.py                 the pipeline
CHANGES.md                   what changed from version 4, and why
Dockerfile                   the published image (pipeline + tools + web interface)
dist/monica.command|.cmd     double-click launchers
docker/Dockerfile            pipeline-only image, for advanced use
docker/web/                  the web interface (FastAPI)
docker/env_*.yml             conda environments
.github/workflows/image.yml  builds and publishes the image for Intel and ARM
```

The web form builds its advanced section by reading the pipeline's own argument
parser, so adding an option to `monica_v5.py` adds the field to the page with no
other change.

## Citing

See `CITATION.cff`.

## Licence

See `LICENSE`.
