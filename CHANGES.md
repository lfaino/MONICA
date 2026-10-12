# monica_v5.py - what changed with respect to monica_v4.py

The command line of version 4 still works unchanged. Default settings reproduce
the taxonomic assignment of version 4 exactly: on 4,000 randomly generated blast
tables (including missing taxids, `;`-separated taxids, bitscore ties and taxids
encoded in the sequence id) the two implementations return the same assignment
for every table, and batching the queries does not change any result either.

## Fixed defects

| # | Defect in v4 | Fix |
|---|---|---|
| 1 | `--only_assemble` crashed at startup: with no database resolved the script still ran `blastdbcmd -db "" -info` to compute `nice_slurm`, which was never used anyway | `nice_slurm` removed; the database is only probed when it is needed |
| 2 | Alpha/beta diversity and the PCoA were computed on the **species** instead of the samples (the matrix was transposed before being handed to skbio) | `diversity_analysis` works on a samples x species matrix; Shannon from counts, Bray-Curtis from relative abundances |
| 3 | Species colours were computed reproducibly and then destroyed: the list was fed to a `LinearSegmentedColormap`, which resampled it into a continuous ramp | the colour of each species is passed directly to the bar plot, so it is stable across runs and figures |
| 4 | `species_pd.fillna(0)` was never assigned, so species absent from a barcode stayed NaN in the histogram | the abundance matrix is filled once, when it is built |
| 5 | The `--database_external` fasta branch read `--database`, wrote into `$BLASTDB` (usually read-only) and never ran `makeblastdb` | it reads the right argument, writes next to the outputs and builds the database |
| 6 | Output names were inconsistent (`barcode01_rank.txt` next to `barcode01raw_blastn_results.txt`) because the barcode was renamed halfway through the iteration | the sample name, with or without the `raw` suffix, is decided once and used for every file of that sample |
| 7 | The MLST branch parsed `outputRacon`, which was always empty, and `patternsST`/`orderST` were computed and never used | the MLST search runs on the sequences that are actually being analysed and writes `<sample>.mlst_loci.tsv`; the unused profile parsing is gone |
| 8 | Sub-directories of a barcode folder were lost: file names were rebuilt from `cwd/barcode` instead of the directory being walked | the paths returned by `os.walk` are used |

## Performance

| # | Problem | Fix |
|---|---|---|
| 9 | One fasta file and one `blastn` process **per read**: with 100k reads, 200k temporary files and 100k invocations, each re-reading the database index | reads are batched into chunks (`--chunk_size`, default 1000) and each chunk is one `blastn` call with `-num_threads` (`--blast_threads`) |
| 10 | On SLURM, one `sbatch` per read | a single job array (`--max_job_submit` becomes the `%N` throttle) |
| 11 | Recovering the reads racon could not correct was `any(... startswith ...)` inside a loop, i.e. O(N^2) - about 10^10 string comparisons on 100k reads | the corrected ids go into a set and the fastq is streamed once |
| 12 | A `NamedTemporaryFile` object was kept alive for every read just to use its name, exhausting the file descriptor limit | `temp_path()` returns a path and leaves no open handle |
| 13 | Temporary files named after `record.id` (a read id containing `/` crashed the run; identical ids in different barcodes overwrote each other) | temporary names are generated, read ids are never used as file names |
| 14 | The temporary directory was never removed | it is removed in a `finally` block, unless `--keep_temp` or `--verbose` |

## Robustness

| # | Problem | Fix |
|---|---|---|
| 15 | No external command's exit status was ever checked: a failing minimap2, racon or SPAdes silently produced an empty file and the pipeline carried on | `run_cmd`/`run_pipe` check every stage, log stderr and raise `CommandError` |
| 16 | `assembly()` returned its output path even after an error, and the caller parsed a file that did not exist | `assemble()` returns `None` on failure; a failure of the optional polishing steps degrades to the previous valid assembly |
| 17 | Every command was built by string interpolation and run with `shell=True` (paths with spaces broke the command; values came from the command line) | commands are argument lists, no shell involved |
| 18 | A SLURM job was considered successful as soon as it left `squeue`, which cannot distinguish a completed job from a cancelled or out-of-memory one | the final state of every array task is read from `sacct`, and missing outputs are reported |
| 19 | No input validation | `--interval`, `--min`, `--threads`, `--chunk_size`, the input folders and the email address are checked before any work starts |
| 20 | The logger had no handler: log lines were only visible because `ete4` happens to call `logging.basicConfig` on import | logging is configured explicitly, with `--log_file` to also write to a file |
| 21 | `NCBITaxa()` was instantiated at import time, so even `--help` triggered the download of the ~800 MB taxonomy | the taxonomy is opened on first use |
| 22 | `ETE_DATA_HOME` was set after importing ete4 and had no effect (ete4 >= 4 reads `XDG_DATA_HOME`) | both variables are set before the import |
| 23 | Module-level globals (`dirpathTEMP`, `verbose`, `qos`, `od`, `cwd`) | a `Run` dataclass passed explicitly |
| 24 | Nothing recorded how a run was produced | `<output>.run_manifest.json`: parameters, resolved database, samples, tool versions, timestamp, command line |

## Deliberately kept

* **The workbook holds every species, the figures only those above `--min`.**
  This is intentional, so the legend does not fill up with hundreds of species.
  The workbook now makes it explicit: sheet `counts` is the complete table,
  sheet `plotted` is what the figures show.
* **The bitscore threshold scales with the shorter of query and subject**
  (`legacy_min_bitscore`): a hit is kept only if it scores at least as much as
  an alignment spanning `min(qlen, slen)` at 80% identity. This is the default.
* **The tie-break of version 4**, including its consequence described below.

## Two rules that can now be changed (defaults keep v4 behaviour)

* `--blast_filter identity --min_pident X [--min_qcov Y]` replaces the
  Karlin-Altschul threshold with explicit percent identity and query coverage
  cut-offs (`qcovhsp` is added to the blast output format only when requested).
* `--lca lca [--lca_tolerance 0.01]` assigns the last common ancestor of every
  hit within a relative tolerance of the best bitscore, instead of the
  sequential best-hit rule.

  **Worth knowing before choosing.** The v4 tie-break compares the lineage of
  the current best hit with itself, because the variable holding "the other
  lineage" is only updated when a hit is promoted. The character-wise common
  prefix of a string with itself is the whole string, and the rule then drops
  the last token - so when two hits tie on bitscore and the new one does not win
  on identity, the read is assigned to the **parent** of the current best taxid
  (e.g. genus instead of species). This is reproduced exactly by the default
  `--lca best`; `--lca lca` does not have this behaviour.

## Not exercised by the test run

SLURM job arrays (no scheduler in the container), `--debarcode`, `--search`
(needs NCBI), MLST (needs PubMLST) and nanopolish with real fast5 files.

## Keeping the script usable from Wooey

Wooey builds its form with clinto, which `exec`s the script source with
`__name__` set to `"__main__"` and a patched `parse_args` that raises as soon as
it is called, capturing the parser. Three consequences for future edits:

* `main()` must still reach `parse_args()` as its first action, and the
  `if __name__ == "__main__"` block must still call it. Moving the call behind a
  condition, or building the parser only inside a branch, makes the form empty.
* Module-level code runs during that introspection, so it must stay cheap and
  side-effect free. Version 4 instantiated `NCBITaxa()` at import: every
  registration or re-registration checked, and potentially downloaded, the
  taxonomy. Version 5 opens it lazily.
* `--version` must print a PEP 440 string (`5.0.0`, not `monica 5.0.0`):
  Wooey parses that output to label the script version, and a non-conforming
  string raises `InvalidVersion` during registration.

## Taxonomy resolution (fixed after the first production run)

Version 5.0.0 set `XDG_DATA_HOME` to the blast database directory
unconditionally, to make `ETE_DATA_HOME` effective. On a machine whose
`taxa.sqlite` sits in ete4's default location (`~/.local/share/ete/`) this
HIDES the taxonomy: blast runs, every hit is found, every lineage lookup fails,
and the run ends with `NO READ CLASSIFIED` and no error. Three changes:

* `XDG_DATA_HOME` is only redirected when `<BLASTDB>/ete/taxa.sqlite` actually
  exists, and never when the caller has already set the variable.
* `--taxonomy_db /path/to/taxa.sqlite` selects the file explicitly.
* The taxonomy is checked before the first blast is submitted: a missing or
  unreadable `taxa.sqlite` stops the run immediately, naming the path it looked
  for, instead of wasting the blast and reporting an empty result. If hits are
  found but nothing can be assigned, the error names the taxonomy file rather
  than saying that no read was classified.

## SLURM logs and polling

* The per-task `blast_<jobid>_<task>.out/.err` files and the job script are
  written to the temporary directory and removed with it. They are noise in a
  successful run - one empty pair per chunk - and a web front end would offer
  every one of them as a result file.
* `--verbose` (or `--keep_temp`) keeps them in `<output>/slurm_logs` as before.
* When tasks fail, the logs of the first 20 failed tasks and the job script are
  copied to `<output>/slurm_logs` before the temporary directory is removed, and
  the error message points there: a failure is still fully diagnosable.
* A rejected `sbatch` used to surface as a bare `AssertionError` from
  simple_slurm; it now names the job script that was submitted.
* The queue is polled after 2 s, then 4, 8, 10 s. With the fixed 10 s interval,
  a small array that finished immediately still cost 10 s per barcode.

## Where the intermediate files go

The query chunks (`<sample>.NNNNNN.fasta`) and the per-chunk blast tables
(`<sample>.NNNNNN.blastn`) are written to `<--temp>/monica.XXXXXX/` and removed
when the run ends; the output directory only receives results. Two additions
make that guarantee hold:

* the location is logged at INFO (`INTERMEDIATE FILES IN ...`), so a run can be
  audited from its own log;
* if `--temp` points inside the output directory the chunks would land next to
  the results, so the system temporary directory is used instead, with a
  warning. Under `--slurm` that fallback would be node-local and the compute
  nodes could not read the chunks, so the run stops and asks for a shared
  directory outside the output folder.

The per-sample report `<sample>_blastn_results.txt` is a result, not an
intermediate, and stays in the output directory.
