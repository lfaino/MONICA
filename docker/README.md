# Immagine Docker per `monica_v4.py`

Immagine che contiene il runtime Python e tutti i tool esterni chiamati dalla
pipeline (minimap2, racon, porechop, jellyfish, SPAdes, nanopolish, bwa,
samtools/bcftools/htslib, BLAST+, Krona, freebayes).

## Build

```bash
# dalla directory che contiene monica_v4.py e la cartella docker/
docker build -t monica:v4 -f docker/Dockerfile .
```

## Percorsi attesi dentro il container

Lo script ha dei percorsi scritti nel codice; nell'immagine esistono già come
punti di mount vuoti:

| percorso                      | contenuto                                             |
|-------------------------------|-------------------------------------------------------|
| `/mnt/dataQNAP/blastdb`       | database BLAST (`$BLASTDB`), tassonomia ete4 e Krona   |
| `/mnt/data3ZEUSPV/SequenceData` | fastq/fast5 di input (default di `--folder_fastq`)   |
| `/mnt/dataZEUSPV/wooey/tmp`   | directory temporanea (default di `--temp`)             |
| `/data`                       | working directory: qui vengono scritti gli output      |

I risultati (`.xlsx`, `.histogram.pdf`, `.heatmap.pdf`, `.diversity.PCA.pdf`,
`*_rank.txt`, `*_blastn_results.txt`, `.html` di Krona) vengono scritti nella
directory di lavoro del container, quindi conviene montarla.

## Preparazione (una tantum) della tassonomia

Il modulo `ete4` costruisce `taxa.sqlite` (~800 MB) al primo avvio scaricandolo
da NCBI. Per farlo una volta sola e riusarlo, tienilo sul volume montato:

```bash
mkdir -p /percorso/blastdb
# tassonomia ete4  ->  /percorso/blastdb/ete/taxa.sqlite
docker run --rm -v /percorso/blastdb:/mnt/dataQNAP/blastdb monica:v4 \
    python -c "from ete4 import NCBITaxa; NCBITaxa()"

# tassonomia Krona ->  /percorso/blastdb/taxonomy.tab  (serve per il report .html)
docker run --rm -v /percorso/blastdb:/mnt/dataQNAP/blastdb --entrypoint bash monica:v4 \
    -lc 'ktUpdateTaxonomy.sh /mnt/dataQNAP/blastdb'
```

Dopo questo passaggio il container gira anche senza rete (`--network none`),
a patto di non usare `--search`/`--species`, che interrogano NCBI e PubMLST.

## Esecuzione

```bash
docker run --rm \
    -v /percorso/blastdb:/mnt/dataQNAP/blastdb \
    -v /percorso/sequencing:/data/reads \
    -v $(pwd)/risultati:/data/out \
    -w /data/out \
    monica:v4 \
        --barcode barcode01-12 \
        --folder_fastq /data/reads \
        -d nt.nal \
        -e nome.cognome@ente.it \
        -t 16 \
        -o run_2026_09
```

Per una shell dentro il container: `docker run --rm -it ... monica:v4 bash`.

## Note

* **SLURM (`--slurm`)**: il client Slurm non è nell'immagine (`sbatch`/`squeue`
  non sono pacchettizzati in conda). Dentro il container usa la modalità locale
  (process pool). Se vuoi comunque sottomettere alla griglia, monta i binari e
  la configurazione dell'host:
  `-v /usr/bin/sbatch:/usr/bin/sbatch:ro -v /usr/bin/squeue:/usr/bin/squeue:ro -v /etc/slurm:/etc/slurm:ro`
  (insieme alle librerie Munge/Slurm), oppure lancia lo script fuori dal
  container e usa l'immagine solo per i tool.
* **`spades`**: le build recenti di bioconda installano solo `spades.py`;
  l'immagine crea il link `spades` richiesto dallo script.
* **`jellyfish`**: installato da bioconda come `kmer-jellyfish` (il pacchetto
  `jellyfish` di conda-forge è un'altra libreria Python).
* **fast5**: servono solo con `-a/-ua/-oa` (nanopolish); monta la cartella e
  passala con `-f5`.
* **Permessi**: il container gira come root e i file di output risultano di
  root. Per evitarlo: `docker run --user $(id -u):$(id -g) ...` (assicurati che
  le directory montate siano scrivibili da quell'utente).
