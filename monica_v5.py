#!/usr/bin/env python3
"""MONICA - detection and quantification of pathogens from Nanopore amplicons.

Per barcode the pipeline:

1. collects and filters the fastq reads (length interval, optional subsampling);
2. trims adapters with porechop (optionally demultiplexing with it);
3. optionally corrects the reads (minimap2 + racon) and assembles them
   (jellyfish + SPAdes, polished with nanopolish and bwa/bcftools);
4. blasts the reads (or the contigs) against a nucleotide database, locally in a
   process pool or on SLURM as a single job array;
5. assigns taxonomy with the NCBI taxonomy dump (ete4) and writes the abundance
   table, the Krona report, the stacked histogram, the clustermap and - on
   request - alpha/beta diversity.

Version 5 keeps the command line and the analysis of version 4 and fixes the
defects listed in CHANGES.md. The two changes that would alter the numbers -
the blast hit filter and the last-common-ancestor rule - are available as
options and are disabled by default.
"""

from __future__ import annotations

import argparse
import gzip
import json
import logging
import math
import os
import random
import re
import shutil
import string
import subprocess
import sys
import tempfile
import time
import uuid
import xml.etree.ElementTree as ET
from collections import Counter
from dataclasses import asdict, dataclass, field
from multiprocessing import Pool
from typing import Iterable, Iterator, Optional, Sequence

__version__ = "5.0.0"

# --------------------------------------------------------------------------- #
# Environment
# --------------------------------------------------------------------------- #
# ete4 >= 4 resolves its data directory from XDG_DATA_HOME (+ '/ete'); setting
# ETE_DATA_HOME after importing ete4 - as version 4 did - has no effect.
#
# XDG_DATA_HOME is only redirected to the blast database directory when the
# taxonomy really lives there. Forcing it unconditionally would HIDE a taxonomy
# that ete4 would otherwise find in its default location (~/.local/share/ete),
# and every lineage lookup would then fail silently, leaving a run that blasts
# correctly and classifies nothing.

# Paths of the original deployment. They are still the defaults on that machine,
# but only if they exist: the same script has to run on a laptop, inside the
# distributed image, where none of them is present. MONICA_TEMP and
# MONICA_FAST5 override them.
LEGACY_TEMP = "/mnt/dataZEUSPV/wooey/tmp/"
LEGACY_SEQUENCE_DATA = "/mnt/data3ZEUSPV/SequenceData/"


def default_temp() -> str:
    """Directory the run's temporary folder is created in."""
    return (os.environ.get("MONICA_TEMP")
            or (LEGACY_TEMP if os.path.isdir(LEGACY_TEMP) else tempfile.gettempdir()))


def default_fast5() -> str:
    """Default fast5 directory; empty when the deployment path is absent."""
    return (os.environ.get("MONICA_FAST5")
            or (LEGACY_SEQUENCE_DATA if os.path.isdir(LEGACY_SEQUENCE_DATA) else ""))


DEFAULT_BLASTDB_DIR = os.environ.get("BLASTDB") or "/mnt/dataQNAP/blastdb/"
os.environ.setdefault("BLASTDB", DEFAULT_BLASTDB_DIR)
if "XDG_DATA_HOME" not in os.environ and os.path.exists(
        os.path.join(DEFAULT_BLASTDB_DIR, "ete", "taxa.sqlite")):
    os.environ["XDG_DATA_HOME"] = DEFAULT_BLASTDB_DIR.rstrip("/")
os.environ.setdefault("ETE_DATA_HOME", DEFAULT_BLASTDB_DIR)
os.environ.setdefault("SLURM_CONF", "/etc/slurm/slurm.conf")

import matplotlib  # noqa: E402

matplotlib.use("Agg")

import Bio  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
import requests  # noqa: E402
import seaborn as sns  # noqa: E402
from Bio import Entrez  # noqa: E402,F401  (imported for its side effect: Bio.Entrez)
from Bio import SeqIO  # noqa: E402
from ete4 import NCBITaxa  # noqa: E402
from simple_slurm import Slurm  # noqa: E402
from skbio.diversity import alpha_diversity, beta_diversity  # noqa: E402
from skbio.stats.ordination import pcoa  # noqa: E402
from tqdm import tqdm  # noqa: E402

matplotlib.style.use("ggplot")

logger = logging.getLogger("monica")

# Karlin-Altschul parameters of the legacy bitscore filter (see keep_hit).
BLAST_LAMBDA = 1.37
BLAST_K = 0.711
BLAST_MIN_IDENTITY_FRACTION = 0.8
BLAST_MATCH_SCORE = 1
BLAST_MISMATCH_SCORE = -2

BLAST_FIELDS = ["qseqid", "sseqid", "bitscore", "sscinames", "pident",
                "evalue", "staxids", "qlen", "slen"]
BLAST_FIELDS_QCOV = BLAST_FIELDS + ["qcovhsp"]

DESIRED_RANKS = ["species", "genus", "family", "order", "class", "phylum",
                 "clade", "superkingdom", "subspecies", "species subgroup",
                 "species group"]

MLST_DATABASE = "https://pubmlst.org/data/dbases.xml"

DATABASE_CHOICES = [
    "Plant_Bacteria_Funghi.nal", "nt.nal", "ITS_16S_18S_28S_LSU_SSU.nal",
    "Metazoa_Organism_OR_metazoa_All_Fields_.nal", "Xanthomonas_genomes.nal",
    "curatedXylellaDatabase.nal", "customITSreduced.nal",
    "customITSreducedReduced.nal", "BacteriaFungiVitis.nal",
    "50_4000RefSeq_biomol_genomic_ITS_16S_18S_28S_LSU_SSU.nal",
    "bacterialGenome2024_01_24.fna.nal", "emu", "16S_ITS.nal", "uniteMonica",
]

_NCBI: Optional[NCBITaxa] = None
TAXONOMY_DB: Optional[str] = None   # set by main() from --taxonomy_db


def taxonomy_dbfile() -> str:
    """Path of the taxa.sqlite ete4 will use, without opening (or building) it."""
    if TAXONOMY_DB:
        return TAXONOMY_DB
    import ete4

    home = getattr(ete4, "ETE_DATA_HOME", None) or os.path.join(
        os.environ.get("XDG_DATA_HOME", os.path.expanduser("~/.local/share")), "ete")
    return os.path.join(home, "taxa.sqlite")


def ncbi() -> NCBITaxa:
    """Return the NCBI taxonomy handle, opening it on first use.

    Version 4 instantiated NCBITaxa at import time, so even `--help` triggered
    the download of the ~800 MB taxonomy dump.
    """
    global _NCBI
    if _NCBI is None:
        _NCBI = NCBITaxa(dbfile=TAXONOMY_DB) if TAXONOMY_DB else NCBITaxa()
        logger.debug("TAXONOMY DATABASE %s", _NCBI.dbfile)
    return _NCBI


def temp_root(requested: str, out_dir: str, shared: bool = False) -> str:
    """Directory the run's temporary folder is created in.

    The query chunks and the per-chunk blast tables are written there, so it
    must not sit inside the output directory: those intermediates would end up
    next to the results (and, under a web front end, be offered as result files
    of the job).
    """
    root = os.path.abspath(requested) if requested else tempfile.gettempdir()
    out_dir = os.path.abspath(out_dir)
    try:
        inside = os.path.commonpath([root, out_dir]) == out_dir
    except ValueError:      # different drives: cannot be inside
        inside = False
    if inside:
        if shared:
            # The fallback would be node-local, and the compute nodes of the array
            # would not see the chunks at all.
            raise CommandError(
                f"--temp {root} is inside the output directory {out_dir}: the blast "
                "chunks would be written next to the results. With --slurm they must "
                "also be on a filesystem the compute nodes can read, so pass --temp "
                "with a shared directory outside the output folder.")
        fallback = tempfile.gettempdir()
        logger.warning(
            "--temp %s is inside the output directory %s: the blast chunks would be "
            "written next to the results; using %s instead", root, out_dir, fallback)
        return fallback
    return root


def check_taxonomy() -> str:
    """Verify the taxonomy before any blast runs; return the file in use.

    A missing or unreadable taxa.sqlite used to surface only as "NO READ
    CLASSIFIED" at the end of a barcode, after the blast had already been paid
    for. It is now a startup error that names the file it looked for.
    """
    dbfile = taxonomy_dbfile()
    if not os.path.exists(dbfile):
        raise CommandError(
            f"the ete4 taxonomy was not found at {dbfile}.\n"
            "Either pass --taxonomy_db /path/to/taxa.sqlite, or set XDG_DATA_HOME to the "
            "directory that contains ete/taxa.sqlite, or build it once with --update_ete.\n"
            "Hint: `locate taxa.sqlite` or look in ~/.local/share/ete/")
    if get_lineage("9606") is None:
        raise CommandError(
            f"{dbfile} exists but no lineage can be read from it "
            f"({os.path.getsize(dbfile) / 1e6:.1f} MB); rebuild it with --update_ete")
    return dbfile


class CommandError(RuntimeError):
    """An external command exited with a non-zero status."""


# --------------------------------------------------------------------------- #
# Run configuration
# --------------------------------------------------------------------------- #


@dataclass
class Run:
    """Everything the worker functions need; replaces the module globals of v4."""

    temp_dir: str
    out_dir: str
    fastq_dir: str
    blastdb_dir: str
    threads: int = 10
    blast_threads: int = 1
    chunk_size: int = 1000
    verbose: bool = False
    keep_temp: bool = False
    cleanup: bool = False
    qos: str = "phyto"
    max_array_jobs: int = 5000
    blast_filter: str = "legacy"
    min_pident: float = 0.0
    min_qcov: float = 0.0
    lca_mode: str = "best"
    lca_tolerance: float = 0.0

    @property
    def outfmt_fields(self) -> Sequence[str]:
        """Query coverage is only requested when a filter actually needs it."""
        if self.blast_filter == "identity" and self.min_qcov > 0:
            return BLAST_FIELDS_QCOV
        return BLAST_FIELDS


@dataclass
class BlastJob:
    """One blastn invocation over a chunk of query sequences."""

    query: str
    out: str
    db: str
    threads: int = 1
    outfmt: str = " ".join(BLAST_FIELDS)
    resume: bool = False      # reuse a non-empty output from a previous run
    remove_query: bool = True  # drop the query file once blast is done


@dataclass
class Hit:
    """One blastn tabular row."""

    qseqid: str
    sseqid: str
    bitscore: float
    sscinames: str
    pident: float
    evalue: float
    staxids: str
    qlen: int
    slen: int
    qcovhsp: Optional[float] = None


@dataclass
class Assignment:
    """Taxonomic assignment of one query sequence."""

    taxid: str
    qseqid: str
    bitscore: float
    pident: float
    sseqid: str
    lineage: list = field(default_factory=list)


# --------------------------------------------------------------------------- #
# Subprocess helpers
# --------------------------------------------------------------------------- #


def run_cmd(args: Sequence[str], cwd: Optional[str] = None, stdout=None,
            check: bool = True, log: bool = True) -> subprocess.CompletedProcess:
    """Run an external command given as an argument list (no shell).

    Version 4 built every command by string interpolation and ran it through
    `shell=True`, never looking at the exit status: a failing minimap2 or SPAdes
    silently produced an empty file and the pipeline carried on. Here stderr is
    captured and a non-zero exit raises CommandError.
    """
    if log:
        logger.debug("RUN %s", " ".join(args))
    proc = subprocess.run(
        list(args), cwd=cwd,
        stdout=stdout if stdout is not None else subprocess.PIPE,
        stderr=subprocess.PIPE, check=False,
    )
    stderr = (proc.stderr or b"").decode(errors="replace").strip()
    if proc.returncode != 0:
        if check:
            raise CommandError(
                f"{args[0]} exited with {proc.returncode}\n"
                f"command: {' '.join(args)}\n{stderr[-2000:]}"
            )
        logger.warning("%s exited with %s: %s", args[0], proc.returncode, stderr[-500:])
    elif stderr:
        logger.debug("%s stderr: %s", args[0], stderr[-500:])
    return proc


def run_pipe(cmds: Sequence[Sequence[str]], stdout=None, cwd: Optional[str] = None,
             check: bool = True) -> None:
    """Run `cmd1 | cmd2 | ...`, checking the exit status of every stage."""
    logger.debug("RUN %s", " | ".join(" ".join(c) for c in cmds))
    err_files = [tempfile.TemporaryFile() for _ in cmds]
    procs: list[subprocess.Popen] = []
    prev_stdout = None
    for index, cmd in enumerate(cmds):
        last = index == len(cmds) - 1
        proc = subprocess.Popen(
            list(cmd), cwd=cwd, stdin=prev_stdout,
            stdout=(stdout if last else subprocess.PIPE),
            stderr=err_files[index],
        )
        if prev_stdout is not None:
            prev_stdout.close()
        prev_stdout = proc.stdout
        procs.append(proc)
    for proc in procs:
        proc.wait()

    failures = []
    for cmd, proc, err_file in zip(cmds, procs, err_files):
        err_file.seek(0)
        message = err_file.read().decode(errors="replace").strip()
        err_file.close()
        if proc.returncode != 0:
            failures.append(f"{cmd[0]} exited with {proc.returncode}: {message[-1000:]}")
        elif message:
            logger.debug("%s stderr: %s", cmd[0], message[-500:])
    if failures and check:
        raise CommandError("; ".join(failures))
    for failure in failures:
        logger.warning(failure)


def tool_version(cmd: Sequence[str]) -> str:
    """Best-effort one-line version string of an external tool."""
    try:
        proc = subprocess.run(list(cmd), stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT, timeout=20, check=False)
        return proc.stdout.decode(errors="replace").strip().splitlines()[0]
    except Exception:
        return "unknown"


def collect_tool_versions() -> dict:
    """Versions of the external tools, for the run manifest."""
    return {
        "minimap2": tool_version(["minimap2", "--version"]),
        "racon": tool_version(["racon", "--version"]),
        "porechop": tool_version(["porechop", "--version"]),
        "jellyfish": tool_version(["jellyfish", "--version"]),
        "spades": tool_version(["spades", "--version"]),
        "nanopolish": tool_version(["nanopolish", "--version"]),
        "bwa": tool_version(["bwa"]),
        "samtools": tool_version(["samtools", "--version"]),
        "bcftools": tool_version(["bcftools", "--version"]),
        "blastn": tool_version(["blastn", "-version"]),
    }


def temp_path(run: Run, suffix: str = "", prefix: str = "monica") -> str:
    """Return a unique path in the run temp directory, leaving no open handle.

    Version 4 kept a NamedTemporaryFile object alive for every read just to use
    its name, so long runs ran out of file descriptors.
    """
    fd, path = tempfile.mkstemp(dir=run.temp_dir, suffix=suffix, prefix=prefix)
    os.close(fd)
    return path


# --------------------------------------------------------------------------- #
# Blast: filtering and taxonomic assignment
# --------------------------------------------------------------------------- #


def legacy_min_bitscore(qlen: int, slen: int) -> float:
    """Bitscore floor used by version 4.

    The threshold is set by the SHORTER of query and subject: a hit is kept only
    if it scores at least as much as an alignment spanning min(qlen, slen) at 80%
    identity (+1/-2, lambda=1.37, K=0.711). Coupling the threshold to the shorter
    sequence is what stops a short subject from passing on a long query, and vice
    versa. This is the default; `--blast_filter identity` expresses the same
    intent with explicit percent identity and query coverage thresholds.
    """
    length = min(slen, qlen)
    matches = round(length * BLAST_MIN_IDENTITY_FRACTION)
    mismatches = length - matches
    raw = matches * BLAST_MATCH_SCORE + mismatches * BLAST_MISMATCH_SCORE
    return (BLAST_LAMBDA * raw - math.log(BLAST_K)) / math.log(2)


def keep_hit(hit: Hit, run: Run) -> bool:
    """Decide whether a blast hit is good enough to carry a taxonomic label."""
    if run.blast_filter == "identity":
        if hit.pident < run.min_pident:
            return False
        if run.min_qcov > 0 and (hit.qcovhsp is None or hit.qcovhsp < run.min_qcov):
            return False
        return True
    return hit.bitscore > legacy_min_bitscore(hit.qlen, hit.slen)


def hit_taxid(hit: Hit) -> tuple:
    """Return (taxid, from_sequence_id) for a hit, or (None, False).

    The flag matters for the tie-break: when a taxid taken from the sequence id
    is missing from the taxonomy dump, version 4 cleared the last known lineage
    (disabling the tie-break for the next hit), while a missing staxids left it
    in place. That asymmetry is preserved.
    """
    if hit.staxids != "0" and ";" not in hit.staxids:
        return hit.staxids, False
    # Custom databases encode the taxid in the sequence identifier.
    if hit.staxids == "0" and hit.sscinames == "N/A" and ":" in hit.sseqid:
        return hit.sseqid.split(":")[0], True
    return None, False


def get_lineage(taxid: str) -> Optional[list]:
    try:
        return ncbi().get_lineage(taxid)
    except Exception:
        logger.debug("TAXID %s NOT FOUND", taxid)
        return None


def common_prefix(*sequences) -> list:
    """Longest common prefix of the given sequences."""
    if not sequences:
        return []
    shortest = min(len(s) for s in sequences)
    prefix = []
    for index in range(shortest):
        value = sequences[0][index]
        if all(s[index] == value for s in sequences[1:]):
            prefix.append(value)
        else:
            break
    return prefix


def legacy_common_taxid(lineage_a: Sequence, lineage_b: Sequence) -> Optional[str]:
    """Tie-break taxid of version 4, reproduced exactly.

    The two lineages are joined into strings, their common prefix is taken and
    the LAST token is dropped, because a prefix taken on characters can cut a
    taxid in half. Note the consequence, kept here for comparability: when the
    two lineages are identical - which is what happens when the rule compares
    the current best with itself - the assignment moves one rank up, to the
    PARENT of the best taxid. Use `--lca lca` for a rule that instead takes the
    last common ancestor of all the hits near the top score.
    """
    joined = {"_".join(str(taxid) for taxid in lineage_a),
              "_".join(str(taxid) for taxid in lineage_b)}
    shortest = min(len(item) for item in joined)
    common = ""
    for index in range(shortest):
        character = next(iter(joined))[index]
        if all(item[index] == character for item in joined):
            common += character
        else:
            break
    if "_" in common and len(common.split("_")) >= 2:
        return common.split("_")[-2]
    return None


def assign_best(hits: Sequence[Hit], run: Run) -> Optional[Assignment]:
    """Version 4 rule: keep the best hit, break ties on identity, else the LCA.

    `current` is the hit every candidate is compared against and `result` is what
    the function returns: version 4 kept them distinct (a tie-break updated the
    reported assignment but not the comparison state) and so does this.
    """
    current: Optional[Assignment] = None   # comparison state
    result: Optional[Assignment] = None    # reported assignment
    last_lineage: Optional[list] = None

    for hit in hits:
        if not keep_hit(hit, run):
            continue
        taxid, from_sequence_id = hit_taxid(hit)
        if taxid is None:
            continue

        def promote() -> None:
            nonlocal current, result, last_lineage
            lineage = get_lineage(taxid)
            if lineage is None:
                if from_sequence_id:
                    last_lineage = []
                return
            last_lineage = lineage
            current = Assignment(taxid, hit.qseqid, hit.bitscore, hit.pident,
                                 hit.sseqid, lineage)
            result = current

        if current is None:
            promote()
        elif current.bitscore < hit.bitscore:
            promote()
        elif (current.bitscore == hit.bitscore and current.taxid != taxid
                and last_lineage is not None and len(last_lineage) > 1):
            if current.pident < hit.pident:
                promote()
            else:
                ancestor = legacy_common_taxid(current.lineage, last_lineage)
                if ancestor is None:
                    result = None
                else:
                    result = Assignment(ancestor, hit.qseqid, hit.bitscore,
                                        hit.pident, hit.sseqid,
                                        get_lineage(ancestor) or [])
    return result


def assign_lca(hits: Sequence[Hit], run: Run) -> Optional[Assignment]:
    """Assign the last common ancestor of every hit near the top bitscore.

    Unlike `assign_best` this does not depend on the order of the rows: all hits
    within `--lca_tolerance` of the best bitscore contribute.
    """
    usable = [(h, t) for h in hits if keep_hit(h, run)
              for t in [hit_taxid(h)[0]] if t is not None]
    if not usable:
        return None

    top = max(h.bitscore for h, _ in usable)
    threshold = top * (1.0 - run.lca_tolerance)
    selected = [(h, t) for h, t in usable if h.bitscore >= threshold]

    lineages, seen = [], set()
    for _, taxid in selected:
        if taxid in seen:
            continue
        seen.add(taxid)
        lineage = get_lineage(taxid)
        if lineage:
            lineages.append(lineage)
    if not lineages:
        return None

    shared = common_prefix(*lineages) if len(lineages) > 1 else lineages[0]
    if not shared:
        return None

    best_hit = max(selected, key=lambda item: (item[0].bitscore, item[0].pident))[0]
    ancestor = str(shared[-1])
    return Assignment(ancestor, best_hit.qseqid, best_hit.bitscore,
                      best_hit.pident, best_hit.sseqid, list(shared))


def parse_hit_line(line: str, fields: Sequence[str]) -> Optional[Hit]:
    values = line.rstrip("\n").split("\t")
    if len(values) != len(fields):
        return None
    record = dict(zip(fields, values))
    try:
        return Hit(
            qseqid=record["qseqid"], sseqid=record["sseqid"],
            bitscore=float(record["bitscore"]), sscinames=record["sscinames"],
            pident=float(record["pident"]), evalue=float(record["evalue"]),
            staxids=record["staxids"], qlen=int(record["qlen"]),
            slen=int(record["slen"]),
            qcovhsp=float(record["qcovhsp"]) if "qcovhsp" in record else None,
        )
    except ValueError:
        return None


def parse_blast_file(job: tuple) -> list:
    """Assign every query of one blast table; returns a list of Assignment.

    Version 4 wrote one fasta and one blastn call per read, so this function
    handled a single query. Queries are now batched, and the hits are grouped by
    query id before applying the very same rule.
    """
    path, run = job
    grouped: dict[str, list] = {}
    fields = run.outfmt_fields
    with open(path) as handle:
        for line in handle:
            hit = parse_hit_line(line, fields)
            if hit is not None:
                grouped.setdefault(hit.qseqid, []).append(hit)

    assign = assign_lca if run.lca_mode == "lca" else assign_best
    results = [assign(hits, run) for hits in grouped.values()]

    if not run.keep_temp:
        query_file = path[: -len(".blastn")] + ".fasta" if path.endswith(".blastn") else None
        for candidate in (query_file, path):
            if candidate and os.path.exists(candidate):
                os.remove(candidate)
    return [r for r in results if r is not None]


# --------------------------------------------------------------------------- #
# Blast: execution
# --------------------------------------------------------------------------- #


def blastn_args(job: BlastJob) -> list:
    return ["blastn", "-db", job.db, "-query", job.query, "-max_target_seqs", "10",
            "-evalue", "0.001", "-num_threads", str(job.threads),
            "-out", job.out, "-outfmt", f"6 {job.outfmt}"]


def run_blast_job(job: BlastJob) -> str:
    """Run one blastn chunk locally; returns the output path."""
    if job.resume and os.path.isfile(job.out) and os.path.getsize(job.out) > 0:
        return job.out
    try:
        run_cmd(blastn_args(job), log=False)
    except CommandError as error:
        logger.error("%s", error)
        open(job.out, "a").close()
    if job.remove_query and os.path.exists(job.query):
        os.remove(job.query)
    return job.out


def write_query_chunks(records: Iterable, run: Run, db: str, prefix: str) -> tuple:
    """Split the sequences into chunks of `run.chunk_size` and build the jobs.

    One blastn per read - the version 4 behaviour - re-opens the database index
    for every single read, which dominates the runtime on a large database, and
    named the temporary files after the read id (a crash waiting for an id with
    a slash in it).
    """
    jobs: list[BlastJob] = []
    written = 0
    total = 0
    handle = None
    outfmt = " ".join(run.outfmt_fields)
    try:
        for record in records:
            if len(record.seq) == 0:
                continue
            if handle is None:
                base = os.path.join(run.temp_dir, f"{prefix}.{len(jobs):06d}")
                handle = open(base + ".fasta", "w")
                jobs.append(BlastJob(
                    query=base + ".fasta", out=base + ".blastn", db=db,
                    threads=run.blast_threads, outfmt=outfmt,
                    resume=run.verbose and not run.cleanup,
                    remove_query=not run.keep_temp,
                ))
                written = 0
            SeqIO.write(record, handle, "fasta")
            written += 1
            total += 1
            if written >= run.chunk_size:
                handle.close()
                handle = None
    finally:
        if handle is not None:
            handle.close()
    return jobs, total


def blast_locally(jobs: Sequence[BlastJob], run: Run) -> list:
    """Run the blast chunks in a process pool."""
    workers = max(1, run.threads // max(1, run.blast_threads))
    logger.info("RUNNING BLAST ON %s CHUNKS USING %s WORKERS", len(jobs), workers)
    outputs = []
    with Pool(processes=workers) as pool:
        for out in tqdm(pool.imap(run_blast_job, jobs), total=len(jobs),
                        disable=not sys.stderr.isatty()):
            outputs.append(out)
    return outputs


def slurm_job_states(job_id: str) -> dict:
    """Final state of every task of a job, from sacct (empty if unavailable)."""
    if shutil.which("sacct") is None:
        return {}
    proc = subprocess.run(["sacct", "-j", str(job_id), "-n", "-P", "-X",
                           "--format=JobID,State"],
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    if proc.returncode != 0:
        return {}
    states = {}
    for line in proc.stdout.decode(errors="replace").splitlines():
        if "|" in line:
            task, state = line.split("|", 1)
            states[task.strip()] = state.strip()
    return states


def wait_for_slurm(job_id: str, poll_seconds: int = 10) -> None:
    """Block until no task of the array is left in the queue.

    `squeue` fails with "Invalid job id" once the job has left the queue, which
    is the normal exit condition; any other failure is transient (a busy
    controller) and is retried rather than taken as completion.
    """
    consecutive_errors = 0
    interval = 2  # short at first: a small array can be over in seconds
    while True:
        proc = subprocess.run(["squeue", "-h", "-j", str(job_id), "-o", "%A"],
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        if proc.returncode == 0:
            consecutive_errors = 0
            if not proc.stdout.decode().strip():
                return
        else:
            message = proc.stderr.decode(errors="replace").strip()
            if "invalid job id" in message.lower():
                return
            consecutive_errors += 1
            logger.warning("squeue failed (%s/5): %s", consecutive_errors, message)
            if consecutive_errors >= 5:
                logger.error("GIVING UP ON squeue; ASSUMING THE ARRAY IS DONE")
                return
        time.sleep(interval)
        interval = min(poll_seconds, interval * 2)


def slurm_array_script(jobs: Sequence[BlastJob], manifest: str, run: Run) -> str:
    """Body of the array job: one blastn per task, chunk taken from the manifest.

    Deliberately POSIX-only and free of bash-isms: simple_slurm writes the
    script with a `#!/bin/sh` shebang unless told otherwise, and under dash
    `set -o pipefail` or `${var%%<tab>*}` would kill the task on its first line.
    """
    outfmt = " ".join(run.outfmt_fields)
    return "\n".join([
        "set -eu",
        f'export BLASTDB="{run.blastdb_dir}"',
        f'entry=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "{manifest}")',
        'query=$(printf "%s" "$entry" | cut -f 1)',
        'out=$(printf "%s" "$entry" | cut -f 2)',
        'if [ -z "$query" ] || [ -z "$out" ]; then',
        '    echo "no manifest entry for task $SLURM_ARRAY_TASK_ID" >&2',
        "    exit 1",
        "fi",
        'if ! command -v blastn >/dev/null 2>&1; then',
        '    echo "blastn is not on the PATH of the compute node" >&2',
        "    exit 127",
        "fi",
        f'exec blastn -db "{jobs[0].db}" -query "$query" -max_target_seqs 10 '
        f'-evalue 0.001 -num_threads {run.blast_threads} -out "$out" '
        f'-outfmt "6 {outfmt}"',
        "",
    ])


def preserve_slurm_logs(log_dir: str, run: Run, job_id: str, failed: dict,
                        job_script: str, limit: int = 20) -> str:
    """Copy the logs of the failed tasks into the output directory.

    In a successful run the per-task logs are noise - they are empty, there is
    one pair per chunk, and a web front end would offer them all as results -
    so they are written to the temporary directory and go away with it. When
    something fails they are the only evidence, so the ones that matter are
    copied next to the results.
    """
    kept = os.path.join(run.out_dir, "slurm_logs")
    if os.path.abspath(log_dir) == os.path.abspath(kept):
        return log_dir
    os.makedirs(kept, exist_ok=True)
    if os.path.exists(job_script):
        shutil.copy2(job_script, kept)
    for task in sorted(failed)[:limit]:
        index = task.split("_")[-1]
        for suffix in (".err", ".out"):
            source = os.path.join(log_dir, f"blast_{job_id}_{index}{suffix}")
            if os.path.exists(source):
                shutil.copy2(source, kept)
    logger.error("SLURM LOGS OF THE FAILED TASKS COPIED TO %s", kept)
    return kept


def report_slurm_failures(failed: dict, job_id: str, log_dir: str,
                          job_script: str) -> str:
    """Log the stderr of the first failed task; return its content."""
    logger.warning("%s SLURM TASKS DID NOT COMPLETE: %s",
                   len(failed), sorted(failed.items())[:10])
    logger.error("JOB SCRIPT KEPT AT %s", job_script)
    for task in sorted(failed):
        index = task.split("_")[-1]
        for candidate in (os.path.join(log_dir, f"blast_{job_id}_{index}.err"),
                          os.path.join(log_dir, f"blast_{job_id}_{index}.out")):
            if os.path.exists(candidate) and os.path.getsize(candidate) > 0:
                with open(candidate, errors="replace") as handle:
                    tail = handle.read()[-2000:]
                logger.error("TASK %s WROTE TO %s:\n%s", task, candidate, tail)
                return tail
    logger.error("NO SLURM LOG FOUND IN %s", log_dir)
    return ""


def blast_on_slurm(jobs: Sequence[BlastJob], run: Run) -> list:
    """Submit the blast chunks as a single SLURM job array.

    Version 4 submitted one sbatch per read and inferred success from the job
    disappearing from squeue, which cannot distinguish a completed job from one
    that was cancelled or hit the memory limit; sacct is queried instead, and
    the stderr of a failed task is reported.
    """
    for tool in ("sbatch", "squeue"):
        if shutil.which(tool) is None:
            raise CommandError(f"{tool} is not available: run without --slurm")

    manifest = os.path.join(run.temp_dir, "blast_chunks.tsv")
    with open(manifest, "w") as handle:
        for job in jobs:
            handle.write(f"{job.query}\t{job.out}\n")

    # Per-task logs and the job script live in the temporary directory, which is
    # removed at the end of the run; --verbose (or --keep_temp) puts them next to
    # the results instead, and a failure copies the relevant ones there anyway.
    log_root = run.out_dir if (run.verbose or run.keep_temp) else run.temp_dir
    log_dir = os.path.join(log_root, "slurm_logs")
    os.makedirs(log_dir, exist_ok=True)
    script = slurm_array_script(jobs, manifest, run)
    logger.debug("SLURM ARRAY SCRIPT:\n%s", script)

    slurm = Slurm(
        ntasks=1,
        cpus_per_task=run.blast_threads,
        array=f"0-{len(jobs) - 1}%{run.max_array_jobs}",
        qos=run.qos,
        job_name="monica_blast",
        # ALL propagates the submitting environment, so the compute node sees
        # the same PATH (and therefore the same blastn) as this process.
        export="ALL",
        output=os.path.join(log_dir, "blast_%A_%a.out"),
        error=os.path.join(log_dir, "blast_%A_%a.err"),
    )
    logger.info("SUBMITTING %s BLAST CHUNKS AS ONE SLURM ARRAY", len(jobs))
    # job_file + convert=False writes the script to disk and submits it with
    # `sbatch <file>`: no here-document, so nothing re-quotes the script, and
    # the exact file that ran stays available for inspection.
    job_script = os.path.join(log_dir, "blast_array.sh")
    try:
        job_id = slurm.sbatch(script, convert=False, job_file=job_script,
                              shell="/bin/bash", verbose=run.verbose)
    except Exception as failure:
        # simple_slurm asserts on the sbatch output, so a rejected submission
        # surfaces as a bare AssertionError with no context.
        raise CommandError(
            f"sbatch refused the job ({type(failure).__name__}: {failure}); "
            f"the script it was given is {job_script}") from failure
    logger.info("SLURM ARRAY JOB %s SUBMITTED", job_id)

    wait_for_slurm(str(job_id))

    states = slurm_job_states(str(job_id))
    failed = {task: state for task, state in states.items()
              if state.split()[0] != "COMPLETED"}
    if failed:
        kept_dir = preserve_slurm_logs(log_dir, run, str(job_id), failed, job_script)
        tail = report_slurm_failures(failed, str(job_id), kept_dir,
                                     os.path.join(kept_dir, os.path.basename(job_script)))
        if states and len(failed) == len(states):
            raise CommandError(
                f"every task of SLURM array {job_id} failed; last error:\n{tail}")

    missing = [job.out for job in jobs if not os.path.exists(job.out)]
    if missing:
        logger.warning("%s BLAST OUTPUTS ARE MISSING (FIRST: %s)",
                       len(missing), missing[0])
    return [job.out for job in jobs if os.path.exists(job.out)]


def concatenate_blast_tables(paths: Sequence[str], destination: str) -> None:
    """Concatenate the per-chunk blast tables into the per-sample report."""
    with open(destination, "w") as out:
        for path in paths:
            if os.path.exists(path):
                with open(path) as handle:
                    shutil.copyfileobj(handle, out)


# --------------------------------------------------------------------------- #
# Taxonomy tables
# --------------------------------------------------------------------------- #


def taxid_to_name(taxid: str) -> str:
    try:
        return list(ncbi().get_taxid_translator([taxid]).values())[0]
    except Exception:
        return str(taxid)


def ranks_of_taxid(taxid: str, abundance: float) -> Optional[dict]:
    """Full rank record of one taxid, or None when it is not in the dump."""
    try:
        lineage = ncbi().get_lineage(taxid)
        lineage2ranks = ncbi().get_rank(lineage)
        ranks2lineage = {rank: tid for tid, rank in lineage2ranks.items()}
        record = {"tax_id": taxid, "abundance": abundance}
        for rank in DESIRED_RANKS:
            tid = ranks2lineage.get(rank)
            record[rank] = (list(ncbi().get_taxid_translator([tid]).values())[0]
                            if tid is not None else "<not present>")
        return record
    except Exception:
        logger.debug("TAXID %s NOT FOUND", taxid)
        return None


def write_rank_table(counts: dict, path: str, min_fraction: float) -> None:
    """Write the per-sample rank table, filtered at `min_fraction`."""
    total = sum(counts.values())
    if total == 0:
        return
    rows = []
    for taxid, count in counts.items():
        fraction = count / total
        if fraction < min_fraction:
            continue
        record = ranks_of_taxid(taxid, fraction)
        if record is not None:
            rows.append(record)
    columns = ["tax_id", "abundance"] + DESIRED_RANKS
    pd.DataFrame(rows, columns=columns).to_csv(path, sep="\t", index=False)


# --------------------------------------------------------------------------- #
# Read handling
# --------------------------------------------------------------------------- #


def open_maybe_gzip(path: str):
    return gzip.open(path, "rt") if path.endswith("gz") else open(path)


def filtered_records(handle, min_len: int, max_len: int, subset: int) -> Iterator:
    """Yield the reads inside the length interval, at most `subset` (0 = all)."""
    kept = 0
    for record in SeqIO.parse(handle, "fastq"):
        if not min_len < len(record.seq) < max_len:
            continue
        if subset > 0 and kept >= subset:
            break
        kept += 1
        yield record


def collect_barcode_reads(source: str, run: Run, min_len: int, max_len: int,
                          subset: int) -> Optional[str]:
    """Return a single filtered fastq for the barcode, or None if there is none."""
    destination = temp_path(run, suffix=".fastq", prefix="reads")
    total = 0
    with open(destination, "w") as out:
        if os.path.isfile(source):
            if not (source.endswith("gz") or source.endswith("fastq")):
                logger.error("UNKNOWN INPUT TYPE: %s (expected .fastq or .fastq.gz)", source)
                return None
            with open_maybe_gzip(source) as handle:
                total += SeqIO.write(filtered_records(handle, min_len, max_len, subset),
                                     out, "fastq")
        elif os.path.isdir(source):
            # Version 4 rebuilt the file names from cwd/barcode instead of the
            # directory it was walking, so sub-directories were silently lost.
            files = sorted(
                os.path.join(root, name)
                for root, _dirs, names in os.walk(source)
                for name in names
                if name.endswith("fastq") or name.endswith("fastq.gz") or name.endswith(".gz")
            )
            remaining = subset
            for path in files:
                if subset > 0 and remaining <= 0:
                    break
                with open_maybe_gzip(path) as handle:
                    written = SeqIO.write(
                        filtered_records(handle, min_len, max_len, remaining), out, "fastq")
                total += written
                if subset > 0:
                    remaining -= written
        else:
            logger.error("UNKNOWN INPUT TYPE: %s (expected a file or a folder)", source)
            return None

    if total == 0:
        logger.warning("NO READS PASSED THE LENGTH FILTER FOR %s", source)
        return None
    logger.info("COLLECTED %s READS FROM %s", total, source)
    return destination


def porechop(input_path: str, threads: int, output_path: str, demultiplex: bool,
             run: Run) -> None:
    """Trim adapters, or demultiplex into `output_path` when `demultiplex`."""
    flag = "-b" if demultiplex else "-o"
    run_cmd(["porechop", "-i", input_path, "-t", str(threads), flag, output_path],
            cwd=run.fastq_dir)


# --------------------------------------------------------------------------- #
# Correction and assembly
# --------------------------------------------------------------------------- #


def best_align_match(paf: str, max_overlaps: int) -> str:
    """Keep at most `max_overlaps` extra alignments per query in a PAF file."""
    if max_overlaps <= 0:
        return paf
    filtered = paf + ".reduced.paf"
    seen: set = set()
    count = 0
    with open(paf) as source, open(filtered, "w") as destination:
        for line in source:
            query = line.split("\t", 1)[0]
            if query in seen:
                if count < max_overlaps:
                    destination.write(line)
                    count += 1
            else:
                seen.add(query)
                destination.write(line)
                count = 0
    return filtered


def racon_round(target: str, reads: str, run: Run, max_overlaps: int) -> str:
    """One minimap2 + racon polishing round; returns the polished fasta."""
    paf = temp_path(run, suffix=".paf")
    logger.info("RUNNING MINIMAP")
    with open(paf, "wb") as handle:
        run_cmd(["minimap2", "-x", "ava-ont", "-t", str(run.threads),
                 "-N", str(max_overlaps), target, reads], stdout=handle,
                cwd=run.fastq_dir)
    filtered = best_align_match(paf, max_overlaps)

    logger.info("RUNNING RACON")
    polished = temp_path(run, suffix=".fasta")
    with open(polished, "wb") as handle:
        run_cmd(["racon", "-t", str(run.threads), "-u", "-f", reads, filtered, target],
                stdout=handle, cwd=run.fastq_dir)

    if not run.keep_temp:
        for path in {paf, filtered}:
            if os.path.exists(path):
                os.remove(path)
    return polished


def racon_correct(reads: str, run: Run, num_reads_cor: int) -> str:
    """Two rounds of minimap2 + racon correction.

    racon is called with -u, so the returned fasta holds every input read: the
    ones it could polish in corrected form, the others unchanged. Nothing has to
    be merged back afterwards.
    """
    max_overlaps = num_reads_cor if num_reads_cor > 0 else 100000
    first = racon_round(reads, reads, run, max_overlaps)
    corrected = racon_round(first, reads, run, max_overlaps)
    n_out = sum(1 for line in open(corrected) if line.startswith(">"))
    logger.info("CORRECTION RETURNED %s READS", n_out)
    return corrected


# merge_corrected_reads() used to re-add the reads racon had not corrected, so
# that no read was lost before blast. It was both redundant and harmful:
#
#   * racon runs with -u (--include-unpolished), so its output already contains
#     every input read - corrected where it could be, unchanged elsewhere.
#     Measured on a 22,215-read sample: racon returned 22,215 sequences, none
#     missing, of which 14,480 were actually modified.
#   * it compared read ids verbatim, normalising only the "<id>:0-1200" form of
#     racon suffix. Each racon round appends a bare "r" to the name instead, so
#     after the two rounds of racon_correct() no id matched and *every* read was
#     re-added: the blast query doubled (44,430 sequences for 22,215 reads),
#     doubling the blast cost and giving each read two independent assignments.
#
# The reads now reach blast exactly as racon returns them.


def assemble(corrected_reads: str, raw_reads: str, run: Run, sample: str,
             fast5: str, kmer: int) -> Optional[str]:
    """Assemble and polish the reads; returns the polished fasta or None.

    Version 4 wrapped the whole function in one try/except and returned the
    output path even after a failure, so the caller parsed a file that was never
    written. Each step is checked here and a failure of the optional polishing
    steps degrades to the previous, still valid, assembly.
    """
    logger.info("STARTED TO RUN ASSEMBLY")
    try:
        counts = temp_path(run, suffix=".jfc.fasta")
        logger.info("RUNNING JELLYFISH")
        with open(counts, "wb") as handle:
            run_pipe([
                ["jellyfish", "count", "-m", str(kmer), "-s", "100M", "-t", "10",
                 "-o", "/dev/stdout", "-C", corrected_reads],
                ["jellyfish", "dump", "-L", "4", "/dev/stdin"],
            ], stdout=handle, cwd=run.fastq_dir)

        kmer_file = temp_path(run, suffix=".kmers.fasta")
        n_kmers = 0
        with open(kmer_file, "w") as handle:
            for record in SeqIO.parse(counts, "fasta"):
                if len(record.seq) == kmer:
                    n_kmers += SeqIO.write(record, handle, "fasta")
        if n_kmers == 0:
            logger.error("NO %s-MERS LEFT AFTER JELLYFISH: SKIPPING ASSEMBLY", kmer)
            return None

        assembly_dir = tempfile.mkdtemp(dir=run.temp_dir, prefix="spades.")
        logger.info("RUNNING SPADES")
        run_cmd(["spades", "-s", kmer_file, "-o", assembly_dir, "--only-assembler"],
                cwd=run.fastq_dir)
        contigs = os.path.join(assembly_dir, "contigs.fasta")
        if not os.path.exists(contigs) or os.path.getsize(contigs) == 0:
            logger.error("SPADES PRODUCED NO CONTIGS")
            return None
        spades_contigs = SeqIO.to_dict(SeqIO.parse(contigs, "fasta"))
        logger.info("SPADES ASSEMBLED %s CONTIGS", len(spades_contigs))
    except CommandError as error:
        logger.error("ASSEMBLY FAILED: %s", error)
        return None

    polished = contigs

    # --- optional signal-level polishing (needs the fast5 files) -------------
    if fast5 and os.path.isdir(fast5):
        try:
            bam = temp_path(run, suffix=".bam")
            logger.info("RUNNING MINIMAP")
            with open(bam, "wb") as handle:
                run_pipe([
                    ["minimap2", "-a", "-x", "map-ont", "--secondary=no",
                     "-t", str(run.threads), contigs, raw_reads],
                    ["samtools", "sort", "-O", "BAM"],
                ], stdout=handle, cwd=run.fastq_dir)
            run_cmd(["samtools", "index", bam], cwd=run.fastq_dir)

            logger.info("RUNNING NANOPOLISH")
            run_cmd(["nanopolish", "index", "-d", fast5, raw_reads], cwd=run.fastq_dir)
            regions = []
            for record in SeqIO.parse(contigs, "fasta"):
                vcf = temp_path(run, suffix=".vcf", prefix="polished.")
                run_cmd(["nanopolish", "variants", "--consensus", "-o", vcf,
                         "-w", f"{record.id}:0-{len(record.seq)}", "-r", raw_reads,
                         "-b", bam, "-g", contigs, "-t", str(run.threads),
                         "--min-candidate-frequency", "0.1", "-p", "1",
                         "--fix-homopolymers"], cwd=run.fastq_dir)
                if os.path.isfile(vcf) and os.path.getsize(vcf) > 0:
                    regions.append(vcf)

            if regions:
                nanopolished = os.path.join(assembly_dir, "contigs.nanopolish.fasta")
                with open(nanopolished, "wb") as handle:
                    run_cmd(["nanopolish", "vcf2fasta", "-g", contigs, *regions],
                            stdout=handle, cwd=run.fastq_dir)
                corrected = SeqIO.to_dict(SeqIO.parse(nanopolished, "fasta"))
                n_corrected = len(corrected)
                for name, record in spades_contigs.items():
                    corrected.setdefault(name, record)
                merged = os.path.join(assembly_dir, "contigs.nanopolish.all.fasta")
                SeqIO.write(corrected.values(), merged, "fasta")
                logger.info("%s OUT OF %s CONTIGS CORRECTED BY NANOPOLISH",
                            n_corrected, len(corrected))
                polished = merged
        except CommandError as error:
            logger.warning("NANOPORE CORRECTION FAILED, USING THE SPADES CONTIGS: %s", error)
    else:
        logger.info("NO FAST5 FOLDER: SKIPPING NANOPOLISH")

    # --- consensus from the read alignment ------------------------------------
    try:
        bam = temp_path(run, suffix=".bam")
        run_cmd(["bwa", "index", polished], cwd=run.temp_dir)
        with open(bam, "wb") as handle:
            run_pipe([
                ["bwa", "mem", "-t", str(run.threads), polished, corrected_reads],
                ["samtools", "sort", "-O", "BAM"],
            ], stdout=handle, cwd=run.temp_dir)
        run_cmd(["samtools", "index", bam], cwd=run.temp_dir)

        vcf = temp_path(run, suffix=".vcf")
        with open(vcf, "wb") as handle:
            run_pipe([
                ["bcftools", "mpileup", "-Ou", "-f", polished, bam],
                ["bcftools", "call", "-mv", "--ploidy", "1"],
            ], stdout=handle, cwd=run.temp_dir)
        run_cmd(["bgzip", "-f", vcf], cwd=run.temp_dir)
        run_cmd(["tabix", "-f", "-p", "vcf", vcf + ".gz"], cwd=run.temp_dir)

        final = os.path.join(run.out_dir, sample + ".MLST.done.fasta")
        run_cmd(["bcftools", "consensus", "-f", polished, "-o", final, vcf + ".gz"],
                cwd=run.temp_dir)
        logger.info("FINISHED TO RUN ASSEMBLY")
        return final
    except CommandError as error:
        logger.warning("CONSENSUS STEP FAILED, RETURNING THE UNPOLISHED ASSEMBLY: %s", error)
        fallback = os.path.join(run.out_dir, sample + ".assembly.fasta")
        shutil.copyfile(polished, fallback)
        return fallback


# --------------------------------------------------------------------------- #
# MLST
# --------------------------------------------------------------------------- #


def fetch_mlst_alleles(species: str, run: Run) -> Optional[str]:
    """Download the PubMLST alleles of a species and build a blast database."""
    response = requests.get(MLST_DATABASE, timeout=120)
    response.raise_for_status()
    root = ET.fromstring(response.text)

    alleles = []
    for child in root:
        if child.text and species in child.text:
            for locus in child.findall("mlst/database/loci/locus/url"):
                page = requests.get(locus.text, timeout=120)
                page.raise_for_status()
                alleles.extend(SeqIO.parse(io_string(page.text), "fasta"))
    if not alleles:
        logger.error("NO MLST SCHEME FOUND FOR SPECIES %r", species)
        return None

    path = os.path.join(run.temp_dir, uuid.uuid4().hex + ".mlst.fasta")
    SeqIO.write(alleles, path, "fasta")
    run_cmd(["makeblastdb", "-dbtype", "nucl", "-in", path])
    logger.info("MLST DATABASE BUILT WITH %s ALLELES", len(alleles))
    return path


def io_string(text: str):
    from io import StringIO

    return StringIO(text)


def mlst_locus_counts(reads: str, reads_format: str, mlst_db: str, run: Run) -> dict:
    """Count the best matching MLST locus for every read."""
    jobs, _ = write_query_chunks(SeqIO.parse(reads, reads_format), run, mlst_db, "mlst")
    if not jobs:
        return {}
    counts: Counter = Counter()
    for out in blast_locally(jobs, run):
        best: dict = {}
        with open(out) as handle:
            for line in handle:
                hit = parse_hit_line(line, run.outfmt_fields)
                if hit is None:
                    continue
                if hit.qseqid not in best or hit.bitscore > best[hit.qseqid][0]:
                    best[hit.qseqid] = (hit.bitscore, hit.sseqid)
        counts.update(locus for _, locus in best.values())
        if not run.keep_temp and os.path.exists(out):
            os.remove(out)
    return dict(counts)


# --------------------------------------------------------------------------- #
# Plots and tables
# --------------------------------------------------------------------------- #


def species_colors(names: Sequence[str]) -> dict:
    """Reproducible colour per species name.

    The colour only depends on the name, so a species keeps its colour across
    runs and across figures. Version 4 computed the same colours and then fed
    them to a LinearSegmentedColormap, which resampled them into a continuous
    ramp: the colours actually drawn were interpolations that changed with the
    number of species.
    """
    values = {letter: index + 1 for index, letter in enumerate(string.ascii_letters)}
    colors = {}
    for name in names:
        if " " in name:
            parts = name.split(" ")
            genus, species = parts[0], "".join(parts[1:])
            full = name
        else:
            genus, species, full = name, "unknown", name + "unknown"
        channels = []
        for component in (genus, full, species):
            letters = re.sub("[^a-z]+", "", component.lower())
            seed = sum(values[letter] for letter in letters) if letters else 0
            channels.append(random.Random(seed).randint(0, 255) / 255)
        colors[name] = tuple(channels)
    return colors


def plot_krona(krona_files: Sequence[str], taxonomy_dir: str, output: str,
               out_dir: str) -> None:
    """Build the Krona html report; a missing Krona taxonomy is not fatal."""
    try:
        run_cmd(["ktImportTaxonomy", "-tax", taxonomy_dir, "-o", output + ".html",
                 *krona_files], cwd=out_dir)
    except (CommandError, FileNotFoundError) as error:
        logger.warning("KRONA REPORT NOT PRODUCED: %s", error)


def plot_stacked_histogram(relative: pd.DataFrame, colors: dict, path: str) -> None:
    """Stacked bar plot of the relative abundance per sample."""
    figure, axes = plt.subplots(figsize=(12, 8))
    relative.plot(kind="bar", stacked=True, ax=axes,
                  color=[colors[name] for name in relative.columns])
    axes.legend(loc="center left", bbox_to_anchor=(1.0, 0.5), fontsize=8)
    axes.set_title("Abundance species")
    axes.set_xlabel("barcode")
    axes.set_ylabel("Relative abundance (%)")
    axes.set_ylim(0, 100)
    figure.tight_layout()
    figure.savefig(path, bbox_inches="tight", dpi=300)
    plt.close(figure)


def plot_heatmap(relative: pd.DataFrame, path: str, top: int = 30) -> None:
    """Clustermap of the `top` most abundant species (species x samples)."""
    by_species = relative.T
    selected = by_species.loc[
        by_species.mean(axis=1).sort_values(ascending=False).head(top).index
    ]
    if selected.shape[0] < 2 or selected.shape[1] < 2:
        logger.warning("NOT ENOUGH DATA FOR THE CLUSTERMAP")
        return
    cluster_map = sns.clustermap(selected, cmap="viridis", annot=True, fmt=".1f",
                                 metric="euclidean", method="average", figsize=(12, 10))
    cluster_map.fig.suptitle(
        f"CLUSTERMAP OF THE HIGHER {top} SPECIES FOR EACH SAMPLE", y=1.02)
    cluster_map.savefig(path)
    plt.close(cluster_map.fig)


def diversity_analysis(counts: pd.DataFrame, figure_path: str, shannon_path: str,
                       braycurtis_path: str) -> None:
    """Shannon and Bray-Curtis over the SAMPLES, with the PCoA of the latter.

    Version 4 transposed the matrix before calling skbio, so both indices - and
    the PCoA - described the species instead of the samples.
    """
    if counts.shape[0] < 2:
        logger.warning("ALPHA/BETA DIVERSITY NEEDS AT LEAST 2 SAMPLES")
        return

    shannon = alpha_diversity(metric="shannon", counts=counts.values,
                              ids=list(counts.index))
    shannon.to_csv(shannon_path, sep="\t", header=["shannon"])

    relative = counts.div(counts.sum(axis=1), axis=0).fillna(0)
    braycurtis = beta_diversity(metric="braycurtis", counts=relative.values,
                                ids=list(relative.index))
    braycurtis.to_data_frame().to_csv(braycurtis_path, sep="\t")

    if counts.shape[0] < 3:
        logger.warning("PCoA PLOT NEEDS AT LEAST 3 SAMPLES; CSV FILES WRITTEN")
        return

    coordinates = pcoa(braycurtis).samples
    figure, axes = plt.subplots(figsize=(7, 7))
    axes.scatter(coordinates["PC1"], coordinates["PC2"])
    for sample in coordinates.index:
        axes.text(coordinates.loc[sample, "PC1"], coordinates.loc[sample, "PC2"], sample)
    axes.set_xlabel("PC1")
    axes.set_ylabel("PC2")
    axes.set_title("PCoA - Bray Curtis")
    figure.tight_layout()
    figure.savefig(figure_path, bbox_inches="tight", dpi=300)
    plt.close(figure)


def apply_min_abundance(counts: pd.DataFrame, min_percent: float) -> pd.DataFrame:
    """Zero the species below `min_percent` of a sample and drop empty columns."""
    if min_percent <= 0:
        return counts
    totals = counts.sum(axis=1)
    thresholds = totals * (min_percent / 100.0)
    filtered = counts.where(counts.ge(thresholds, axis=0), 0.0)
    return filtered.loc[:, filtered.sum(axis=0) > 0]


# --------------------------------------------------------------------------- #
# Database resolution
# --------------------------------------------------------------------------- #


def blastdb_exists(name: str) -> bool:
    if not name:
        return False
    proc = subprocess.run(["blastdbcmd", "-db", name, "-info"],
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    return proc.returncode == 0


def strip_extension(path: str) -> str:
    return os.path.splitext(path)[0]


def download_gi_list(args, blastdb_dir: str) -> str:
    """Download the GIs matching an NCBI query into a temporary file."""
    Bio.Entrez.email = args.email
    if args.ncbi_api_key:
        Bio.Entrez.api_key = args.ncbi_api_key
    # Without a key NCBI allows 3 requests/s, with one 10/s.
    pause = 0.34 if args.ncbi_api_key else 1.0

    handle = Bio.Entrez.esearch(db="nucleotide", term=args.search)
    total = int(Bio.Entrez.read(handle)["Count"])
    logger.info("FOUND %s SEQUENCES", total)
    pages = max(1, math.ceil(total / args.retmax))

    identifiers = []
    for page in tqdm(range(pages), disable=not sys.stderr.isatty()):
        for attempt in range(5):
            try:
                handle = Bio.Entrez.esearch(db="nucleotide", term=args.search,
                                            retmax=args.retmax,
                                            retstart=page * args.retmax)
                identifiers.extend(Bio.Entrez.read(handle)["IdList"])
                break
            except Exception as error:  # transient NCBI failures are common
                wait = 2 ** attempt
                logger.warning("NCBI ESEARCH FAILED (%s), RETRY IN %ss", error, wait)
                time.sleep(wait)
        else:
            raise RuntimeError("NCBI esearch kept failing")
        time.sleep(pause)

    fd, path = tempfile.mkstemp(dir=blastdb_dir, prefix="gilist.")
    with os.fdopen(fd, "w") as handle:
        handle.write("\n".join(identifiers))
    return path


def build_fasta_database(fasta: str, run: Run) -> Optional[str]:
    """Clean a fasta and turn it into a blast database.

    Version 4 read `--database` instead of `--database_external` here, wrote the
    cleaned file into $BLASTDB (usually read-only) and never called makeblastdb,
    so the resulting "database" could not be searched.
    """
    kept = []
    for record in SeqIO.parse(fasta, "fasta"):
        if not 200 < len(record.seq) < 10000:
            continue
        parts = record.description.split(" ")
        if len(parts) >= 3 and "sp." not in parts[2]:
            record.description = " ".join([parts[1], parts[2]])
            kept.append(record)
    if not kept:
        logger.error("NO USABLE SEQUENCE IN %s", fasta)
        return None

    destination = os.path.join(run.out_dir,
                               os.path.basename(fasta) + ".clean.fasta")
    SeqIO.write(kept, destination, "fasta")
    run_cmd(["makeblastdb", "-dbtype", "nucl", "-in", destination, "-parse_seqids"])
    logger.info("BUILT BLAST DATABASE %s FROM %s SEQUENCES", destination, len(kept))
    return destination


def resolve_database(args, run: Run) -> Optional[str]:
    """Return the blast database prefix to query, or None when none is usable."""
    external = args.database_external

    if blastdb_exists(external):
        logger.info("USING DATABASE NAME %s", external)
        return external

    if args.search:
        name = os.path.join(run.blastdb_dir, re.sub("[^0-9a-zA-Z]+", "_", args.search))
        if not os.path.exists(name + ".nal") or args.update:
            gi_list = download_gi_list(args, run.blastdb_dir)
            run_cmd(["blastdb_aliastool", "-db", "nt", "-dbtype", "nucl",
                     "-gilist", gi_list, "-out", name])
        logger.info("USING DATABASE NAME %s", name)
        return name

    if args.database:
        candidate = os.path.join(run.blastdb_dir, args.database)
        if os.path.isfile(candidate) and candidate.lower().endswith(".nal"):
            name = strip_extension(candidate)
            logger.info("USING DATABASE NAME %s", name)
            return name
        logger.error("DATABASE %s NOT FOUND IN %s", args.database, run.blastdb_dir)
        return None

    if external:
        if external.endswith("fasta") or external.endswith(".fa"):
            return build_fasta_database(external, run)
        if external.lower().endswith(".nal"):
            for candidate in (os.path.abspath(external),
                              os.path.join(run.blastdb_dir, external)):
                if os.path.isfile(candidate):
                    name = strip_extension(candidate)
                    logger.info("USING DATABASE NAME %s", name)
                    return name
            name = strip_extension(os.path.basename(external))
            if blastdb_exists(name):
                logger.info("USING DATABASE NAME %s", name)
                return name
        candidate = os.path.join(run.blastdb_dir, external)
        if blastdb_exists(candidate):
            logger.info("USING DATABASE NAME %s", candidate)
            return candidate
        logger.error("DATABASE %s NOT FOUND", external)
        return None

    name = os.path.join(run.blastdb_dir, "nt")
    logger.info("USING DATABASE NAME %s", name)
    return name


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="quantify and detect pathogens",
                                     usage="%(prog)s [options]")
    parser.add_argument("-o", "--output", nargs="?", default="output", help="output name")
    parser.add_argument("-r", "--retmax", nargs="?", type=int, default=100000,
                        help="number of id to retrieve per iteration")
    parser.add_argument("--barcode", help="barcode to analyse", required=True)
    parser.add_argument("-f5", "--folder_fast5", default=default_fast5(),
                        help="folder containing fast5 files")
    parser.add_argument("-cwd", "--working_directory", default="", help="working directory")
    parser.add_argument("-t", "--threads", nargs="?", type=int, default=10,
                        help="number of threads")
    parser.add_argument("-d", "--database", nargs="?", choices=DATABASE_CHOICES, default="",
                        help="database name; a subset of the ncbi database")
    parser.add_argument("-q", "--qos", nargs="?", choices=["phyto", "blastn"],
                        default="phyto", help="QOS to use in SLURM")
    parser.add_argument("-dbe", "--database_external", nargs="?", default="",
                        help="blast database prefix, .nal alias or fasta file")
    parser.add_argument("-tmp", "--temp", nargs="?", default=default_temp(),
                        help="directory for the intermediate files (query chunks and "
                             "per-chunk blast tables); with --slurm it must be readable "
                             "from the compute nodes")
    parser.add_argument("--folder_fastq", required=True,
                        default=default_fast5(),
                        help="folder where fastq file are located")
    parser.add_argument("-e", "--email", nargs="?", default="", required=True,
                        help="email for blast search")
    parser.add_argument("-s", "--search", nargs="?", default="",
                        help="ncbi search to retrieve GIs")
    parser.add_argument("-sp", "--species", nargs="?", default="",
                        help="PUBMLST species to use")
    parser.add_argument("-m", "--min", default=1.0, type=float,
                        help="minimum abundance of a species to be reported, "
                             "as percentage of the classified reads [0-100]")
    parser.add_argument("-mr", "--min_reads", type=int, default=100,
                        help="minimum number of reads to consider a sample")
    parser.add_argument("-km", "--kmer", type=int, default=100, help="kmer to use in JELLYFISH")
    parser.add_argument("-int", "--interval", default="100-20000",
                        help="minimum and maximum read length; default [100-20000]")
    parser.add_argument("-su", "--subset", default=0, type=int,
                        help="number of reads to use for each barcode [0 = all]")
    parser.add_argument("-mrc", "--num_reads_cor", type=int, default=0,
                        help="maximum overlaps per read used in the correction step "
                             "[0 = no filtering]")
    parser.add_argument("-a", "--assemble", action="store_true", help="assemble reads")
    parser.add_argument("-ab", "--alpha_beta", action="store_true",
                        help="compute alpha and beta diversity")
    parser.add_argument("-oa", "--only_assemble", action="store_true",
                        help="assemble reads only - no mapping")
    parser.add_argument("-deb", "--debarcode", action="store_true",
                        help="debarcode reads using porechop")
    parser.add_argument("-sl", "--slurm", action="store_true", help="activate slurm grid")
    parser.add_argument("-ua", "--use_assembled", action="store_true",
                        help="use assembled reads for discovery")
    parser.add_argument("-dh", "--database_history", action="store_true",
                        help=argparse.SUPPRESS)
    parser.add_argument("-c", "--correct", action="store_true",
                        help="correct reads via racon")
    parser.add_argument("-u", "--update", action="store_true", help="update database")
    parser.add_argument(
        "--taxonomy_db", nargs="?", default="",
        help="path of the ete4 taxa.sqlite to use; by default the one ete4 resolves "
             "from XDG_DATA_HOME")
    parser.add_argument("-u_e3", "--update_ete", action="store_true",
                        help="update the ete taxonomy database")
    parser.add_argument("-v", "--verbose", action="store_true", help="be verbose")
    parser.add_argument("-cu", "--cleanup", action="store_true",
                        help="ignore the blast output of a previous run")
    parser.add_argument("-mjs", "--max_job_submit", nargs="?", type=int, default=5000,
                        help="maximum number of SLURM array tasks running at once")

    # --- new in version 5 -----------------------------------------------------
    parser.add_argument("--chunk_size", type=int, default=1000,
                        help="number of sequences per blastn call [1000]")
    parser.add_argument("--blast_threads", type=int, default=1,
                        help="threads given to each blastn call [1]")
    parser.add_argument("--keep_temp", action="store_true",
                        help="keep the temporary directory (implied by --verbose)")
    parser.add_argument("--log_file", default="", help="also write the log to this file")
    parser.add_argument("--ncbi_api_key", default=os.environ.get("NCBI_API_KEY", ""),
                        help="NCBI API key for --search (or NCBI_API_KEY)")
    parser.add_argument("--blast_filter", choices=["legacy", "identity"], default="legacy",
                        help="hit filter: 'legacy' reproduces version 4, 'identity' uses "
                             "--min_pident/--min_qcov [legacy]")
    parser.add_argument("--min_pident", type=float, default=90.0,
                        help="minimum percent identity with --blast_filter identity [90]")
    parser.add_argument("--min_qcov", type=float, default=0.0,
                        help="minimum query coverage with --blast_filter identity "
                             "[0 = not used]")
    parser.add_argument("--lca", choices=["best", "lca"], default="best",
                        help="assignment rule: 'best' reproduces version 4, 'lca' takes the "
                             "last common ancestor of the hits near the top score [best]")
    parser.add_argument("--lca_tolerance", type=float, default=0.01,
                        help="relative bitscore tolerance of --lca lca [0.01]")
    parser.add_argument("--version", action="version", version=__version__)
    return parser


def configure_logging(verbose: bool, log_file: str) -> None:
    """Attach handlers explicitly.

    Version 4 created a Formatter and never used it; the log was only visible
    because ete4 happens to call logging.basicConfig on import.
    """
    logger.setLevel(logging.DEBUG if verbose else logging.INFO)
    logger.propagate = False
    formatter = logging.Formatter("[{asctime}] {levelname} - {message}",
                                  style="{", datefmt="%Y-%m-%d %H:%M:%S")
    stream = logging.StreamHandler()
    stream.setFormatter(formatter)
    logger.addHandler(stream)
    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)


def validate_args(args) -> tuple:
    """Check the arguments up front and return the read length interval."""
    errors = []
    try:
        min_len, max_len = (int(value) for value in args.interval.split("-"))
        if min_len >= max_len or min_len < 0:
            errors.append(f"--interval {args.interval}: expected MIN-MAX with MIN < MAX")
    except ValueError:
        min_len, max_len = 0, 0
        errors.append(f"--interval {args.interval}: expected two integers, e.g. 100-20000")

    if not 0 <= args.min <= 100:
        errors.append(f"--min {args.min}: expected a percentage between 0 and 100")
    if args.threads < 1:
        errors.append(f"--threads {args.threads}: expected at least 1")
    if args.chunk_size < 1:
        errors.append(f"--chunk_size {args.chunk_size}: expected at least 1")
    if not os.path.isdir(args.folder_fastq):
        errors.append(f"--folder_fastq {args.folder_fastq}: not a directory")
    if not os.path.isdir(args.temp):
        errors.append(f"--temp {args.temp}: not a directory")
    if "@" not in args.email:
        errors.append(f"--email {args.email}: not an email address")
    if args.use_assembled and not args.folder_fast5:
        errors.append("--use_assembled requires --folder_fast5")

    if errors:
        for error in errors:
            logger.error("%s", error)
        sys.exit(1)
    return min_len, max_len


def expand_barcodes(spec: str) -> list:
    """Expand "barcode01,barcode02", "barcode01-12" or a single name."""
    if "," in spec:
        return [name for name in spec.replace(" ", "").split(",") if name]
    if "-" in spec and spec.startswith("barcode"):
        first, last = spec.replace("barcode", "").split("-")
        return [f"barcode{number:02d}" for number in range(int(first), int(last) + 1)]
    return [spec]


def write_manifest(path: str, args, run: Run, database: str, samples: Sequence[str]) -> None:
    """Record parameters, versions and inputs so a run can be reconstructed."""
    manifest = {
        "monica_version": __version__,
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "command_line": " ".join(sys.argv),
        "arguments": vars(args),
        "run": asdict(run),
        "database": database,
        "samples": list(samples),
        "tool_versions": collect_tool_versions(),
    }
    with open(path, "w") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)


# --------------------------------------------------------------------------- #
# Per-sample processing
# --------------------------------------------------------------------------- #


def process_barcode(barcode: str, args, run: Run, database: str,
                    min_len: int, max_len: int, mlst_db: Optional[str]) -> Optional[tuple]:
    """Process one barcode; returns (sample name, {taxid: reads}, krona file)."""
    source = (os.path.join(run.temp_dir, barcode + ".fastq") if args.debarcode
              else os.path.join(run.fastq_dir, barcode))
    logger.info("EXECUTING BARCODE: %s", barcode)
    if not os.path.exists(source):
        logger.warning("PATH FOR BARCODE %s DOES NOT EXIST: %s", barcode, source)
        return None

    reads = collect_barcode_reads(source, run, min_len, max_len, args.subset)
    if reads is None:
        return None

    # --- adapter trimming ------------------------------------------------------
    if args.debarcode:
        trimmed = reads
    else:
        logger.info("RUNNING PORECHOP")
        trimmed = temp_path(run, suffix=".fastq", prefix="trimmed")
        try:
            porechop(reads, args.threads, trimmed, False, run)
        except CommandError as error:
            logger.error("PORECHOP FAILED FOR %s: %s", barcode, error)
            return None

    # --- correction and assembly -----------------------------------------------
    corrected: Optional[str] = None
    assembled: Optional[str] = None
    if args.correct or args.assemble or args.use_assembled or args.only_assemble:
        try:
            corrected = racon_correct(trimmed, run, args.num_reads_cor)
        except CommandError as error:
            logger.error("READ CORRECTION FAILED FOR %s: %s", barcode, error)
            corrected = None

        if corrected and (args.assemble or args.use_assembled or args.only_assemble):
            assembled = assemble(corrected, trimmed, run, barcode,
                                 args.folder_fast5, args.kmer)
            if assembled:
                logger.info("ASSEMBLY OF %s IN %s", barcode, assembled)

    if args.only_assemble:
        return None

    # --- choose what goes into blast --------------------------------------------
    if args.use_assembled:
        if not assembled:
            logger.error("NO ASSEMBLY AVAILABLE FOR %s: SKIPPING", barcode)
            return None
        sample = barcode
        query_file, query_format = assembled, "fasta"
    elif corrected:
        sample = barcode
        query_file, query_format = corrected, "fasta"
    else:
        # Uncorrected reads: the sample keeps the historical "raw" suffix, and
        # every output of this barcode uses it (version 4 renamed the barcode
        # halfway through, producing barcode01_rank.txt next to
        # barcode01raw_blastn_results.txt).
        sample = barcode if barcode.endswith("raw") else barcode.rstrip("/").split("/")[-1] + "raw"
        query_file, query_format = trimmed, "fastq"

    # --- optional MLST -----------------------------------------------------------
    if mlst_db and query_file:
        counts = mlst_locus_counts(query_file, query_format, mlst_db, run)
        if counts:
            path = os.path.join(run.out_dir, sample + ".mlst_loci.tsv")
            pd.Series(counts, name="reads").rename_axis("locus").to_csv(path, sep="\t")
            logger.info("MLST LOCUS COUNTS WRITTEN TO %s", path)

    # --- blast --------------------------------------------------------------------
    jobs, n_queries = write_query_chunks(
        SeqIO.parse(query_file, query_format), run, database, sample)
    if not jobs:
        logger.warning("NO SEQUENCE TO BLAST FOR %s", barcode)
        return None

    outputs = blast_on_slurm(jobs, run) if args.slurm else blast_locally(jobs, run)

    report = os.path.join(run.out_dir, sample + "_blastn_results.txt")
    concatenate_blast_tables(outputs, report)

    logger.info("PARSING BLAST RESULTS")
    assignments: list = []
    with Pool(processes=run.threads) as pool:
        for chunk in tqdm(pool.imap(parse_blast_file, [(out, run) for out in outputs]),
                          total=len(outputs), disable=not sys.stderr.isatty()):
            assignments.extend(chunk)

    if not assignments:
        if os.path.exists(report) and os.path.getsize(report) > 0:
            logger.error(
                "BLAST FOUND HITS FOR %s BUT NONE COULD BE ASSIGNED: the taxonomy in "
                "%s resolves no lineage for those taxids", barcode, taxonomy_dbfile())
        else:
            logger.warning("NO READ CLASSIFIED FOR %s", barcode)
        return None

    logger.info("NUMBER OF READS CLASSIFIED %s REPRESENTING %.2f%% OF THE %s READS",
                len(assignments), 100.0 * len(assignments) / max(1, n_queries), n_queries)

    counts: Counter = Counter(assignment.taxid for assignment in assignments)

    krona_file = os.path.join(run.temp_dir, sample + ".krona.txt")
    with open(krona_file, "w") as handle:
        for assignment in assignments:
            handle.write(f"{assignment.qseqid}\t{assignment.taxid}\n")

    if run.verbose:
        path = os.path.join(run.out_dir, sample + ".full_ranks.txt")
        with open(path, "w") as handle:
            for assignment in assignments:
                handle.write(f"{taxid_to_name(assignment.taxid)}\t{assignment.sseqid}\n")

    write_rank_table(counts, os.path.join(run.out_dir, sample + "_rank.txt"),
                     args.min / 100.0)
    return sample, dict(counts), krona_file


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    configure_logging(args.verbose, args.log_file)
    logger.info("MONICA %s STARTED", __version__)

    min_len, max_len = validate_args(args)

    out_dir = os.getcwd()
    blastdb_dir = os.environ.get("BLASTDB", DEFAULT_BLASTDB_DIR)
    temp_dir = tempfile.mkdtemp(
        dir=temp_root(args.temp, out_dir, shared=args.slurm), prefix="monica.")
    run = Run(
        temp_dir=temp_dir, out_dir=out_dir, fastq_dir=args.folder_fastq,
        blastdb_dir=blastdb_dir, threads=args.threads,
        blast_threads=max(1, args.blast_threads), chunk_size=args.chunk_size,
        verbose=args.verbose, keep_temp=args.keep_temp or args.verbose,
        cleanup=args.cleanup, qos=args.qos, max_array_jobs=args.max_job_submit,
        blast_filter=args.blast_filter, min_pident=args.min_pident,
        min_qcov=args.min_qcov, lca_mode=args.lca, lca_tolerance=args.lca_tolerance,
    )
    logger.info("INTERMEDIATE FILES IN %s (REMOVED AT THE END%s)",
                temp_dir, "" if not (args.keep_temp or args.verbose) else "; KEPT: --keep_temp/--verbose")

    try:
        if args.taxonomy_db:
            global TAXONOMY_DB
            TAXONOMY_DB = args.taxonomy_db
        if args.update_ete:
            ncbi().update_taxonomy_database()
        if not args.only_assemble:
            logger.info("USING TAXONOMY DATABASE %s", check_taxonomy())

        mlst_db = fetch_mlst_alleles(args.species, run) if args.species else None

        database = ""
        if args.only_assemble:
            logger.info("PERFORMING ONLY ASSEMBLY")
        else:
            database = resolve_database(args, run)
            if not database:
                logger.error("DATABASE NOT FOUND")
                return 2

        barcodes = expand_barcodes(args.barcode)

        if args.debarcode:
            logger.info("RUNNING PORECHOP (DEMULTIPLEXING)")
            porechop(run.fastq_dir, args.threads, temp_dir, True, run)
            for name in os.listdir(temp_dir):
                if name.lower().startswith("bc"):
                    renamed = "barcode" + name[2:]
                    os.rename(os.path.join(temp_dir, name),
                              os.path.join(temp_dir, renamed))
                    logger.info("RENAMED FILE: %s -> %s", name, renamed)

        table: dict = {}
        krona_files: list = []
        for barcode in barcodes:
            logger.debug("STARTING TO PROCESS BARCODE: %s", barcode)
            result = process_barcode(barcode, args, run, database, min_len, max_len, mlst_db)
            if result is None:
                continue
            sample, counts, krona_file = result
            if sum(counts.values()) < args.min_reads and not args.use_assembled:
                logger.warning("NOT ANALYSING %s: LESS THAN %s CLASSIFIED READS",
                               sample, args.min_reads)
                continue
            table[sample] = {taxid_to_name(taxid): value for taxid, value in counts.items()}
            krona_files.append(krona_file)

        if not table:
            logger.warning("NOTHING TO PLOT")
            return 0

        # --- outputs --------------------------------------------------------------
        counts = pd.DataFrame.from_dict(table, orient="index").fillna(0.0)  # samples x species
        filtered = apply_min_abundance(counts, args.min)

        # The workbook intentionally carries the COMPLETE table: the --min
        # threshold only applies to the figures, whose legend would otherwise
        # list hundreds of species. Sheet 'counts' = every species, sheet
        # 'plotted' = what the figures show.
        excel_path = os.path.join(out_dir, args.output + ".xlsx")
        with pd.ExcelWriter(excel_path) as writer:
            counts.T.to_excel(writer, sheet_name="counts")
            filtered.T.to_excel(writer, sheet_name="plotted")

        relative = filtered.div(filtered.sum(axis=1), axis=0).fillna(0.0) * 100
        colors = species_colors(list(relative.columns))

        logger.info("PREPARING PLOTS")
        plot_krona(krona_files, blastdb_dir, args.output, out_dir)
        plot_stacked_histogram(relative, colors,
                               os.path.join(out_dir, args.output + ".histogram.pdf"))
        if len(table) > 1:
            plot_heatmap(relative, os.path.join(out_dir, args.output + ".heatmap.pdf"))
            if args.alpha_beta:
                diversity_analysis(
                    counts,
                    os.path.join(out_dir, args.output + ".diversity.PCA.pdf"),
                    os.path.join(out_dir, args.output + ".shannon.csv"),
                    os.path.join(out_dir, args.output + ".beta_diversity.csv"),
                )

        write_manifest(os.path.join(out_dir, args.output + ".run_manifest.json"),
                       args, run, database, sorted(table))
        logger.info("DONE")
        return 0
    finally:
        if run.keep_temp:
            logger.info("TEMPORARY FILES KEPT IN %s", temp_dir)
        else:
            shutil.rmtree(temp_dir, ignore_errors=True)


if __name__ == "__main__":
    try:
        sys.exit(main())
    except CommandError as failure:
        # Configuration and external-command failures are reported as a message,
        # not as a traceback: the cause is in the message, not in the stack.
        logger.error("%s", failure)
        sys.exit(2)
