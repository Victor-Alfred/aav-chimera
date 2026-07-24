#!/usr/bin/env python3
"""
AAV-Chimera: Nanopore AAV Vector QC & Chimeric Read Detection
===============================
A standalone tool for quality control of recombinant AAV preparations
sequenced on Oxford Nanopore platforms.

Features:
- Coordinate-based genome structure classification (ssAAV/scAAV/snapback/GDM/ICG/backbone)
- Corrected ITR truncation logic using explicit plasmid annotations
- Backbone contamination detection (reads extending beyond ITR-to-ITR region)
- Single-pass BAM analysis for performance
- Memory-efficient unmapped read handling
- Piped alignment (no intermediate SAM)
- Per-base transgene coverage
- Strand bias analysis
- Concatemer/over-packaging detection
- HTML report generation
- Combined reference caching (build once, reuse thereafter)
- Flexible transgene naming (use YOUR construct names)
- Pre- and post-filtering NanoPlot QC
- Dry-run mode and checkpointing
- Input validation and dependency checking
- Parallel sample processing
"""

import argparse
import csv
import gzip
import json
import logging
import re
import shutil
import subprocess
import sys
import textwrap
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path
from typing import Optional

try:
    import pysam
except ImportError:
    sys.exit("ERROR: pysam is required. Install via: pip install pysam")


# ============================================================
# Configuration
# ============================================================

@dataclass
class PipelineConfig:
    """Pipeline configuration with validation."""

    raw_fastq_dir: str
    work_dir: str
    refs_dir: str

    # Reference handling
    reuse_combined_ref: bool = False
    rebuild_combined_ref: bool = False

    # Processing toggles
    run_porechop: bool = True
    run_nanofilt: bool = True
    run_nanoplot: bool = True
    run_mapping: bool = True

    # Filtering parameters
    min_quality: int = 10
    min_length: int = 500
    max_length: int = 6000

    # Alignment parameters
    minimap2_preset: str = "map-ont"
    minimap2_extra: str = "-Y --secondary=no"

    # Chimeric detection parameters
    min_chimeric_segment_length: int = 30
    min_chimeric_mapq: int = 10

    # Transgene/plasmid identification — use YOUR specific construct name
    transgene_name: str = ""

    # Coordinate annotations (0-based, within the transgene/plasmid reference)
    itr_5_start: int = 0
    itr_5_end: int = 145
    transgene_start: int = 145
    transgene_end: int = 4600
    itr_3_start: int = 4600
    itr_3_end: int = 4745

    # ITR masking (positions relative to each ITR start)
    mask_itr_variable_region: bool = True
    itr_variable_start: int = 42
    itr_variable_end: int = 85

    # ITR truncation threshold
    itr_full_length_threshold: int = 100

    # Host genome patterns
    host_ref_patterns: tuple = (
        r"^chr([1-9]|1[0-9]|2[0-2]|X|Y|M|MT)$",
        r"^(?:[1-9]|1[0-9]|2[0-2]|X|Y|M|MT)$",
        r"^hg38.*$",
        r"^GRCh38.*$",
    )

    # Helper/Rep-Cap patterns
    helper_ref_patterns: tuple = (
        r".*rep.*cap.*",
        r".*helper.*",
        r".*pHelper.*",
    )

    # Execution parameters
    threads: int = 4
    parallel_samples: int = 1
    resume: bool = False
    dry_run: bool = False
    generate_html_report: bool = True

    @property
    def itr_length(self) -> int:
        """ITR length derived from coordinates."""
        return self.itr_5_end - self.itr_5_start

    @property
    def expected_genome_start(self) -> int:
        """Start of expected AAV genome (ITR-to-ITR)."""
        return self.itr_5_start

    @property
    def expected_genome_end(self) -> int:
        """End of expected AAV genome (ITR-to-ITR)."""
        return self.itr_3_end

    @property
    def expected_genome_length(self) -> int:
        """Expected full-length AAV genome size."""
        return self.itr_3_end - self.itr_5_start

    def backbone_regions(self, ref_length: int) -> list:
        """
        Return backbone coordinate ranges for a given reference length.
        Backbone = everything outside the ITR-to-ITR region.
        """
        regions = []
        if self.itr_5_start > 0:
            regions.append((0, self.itr_5_start))
        if self.itr_3_end < ref_length:
            regions.append((self.itr_3_end, ref_length))
        return regions

    def validate(self) -> list:
        """Return list of validation errors (empty if valid)."""
        errors = []
        if not Path(self.raw_fastq_dir).is_dir():
            errors.append(f"raw_fastq_dir does not exist: {self.raw_fastq_dir}")
        if not Path(self.refs_dir).is_dir():
            errors.append(f"refs_dir does not exist: {self.refs_dir}")
        if self.min_quality < 0:
            errors.append(f"min_quality must be >= 0, got {self.min_quality}")
        if self.min_length < 0:
            errors.append(f"min_length must be >= 0, got {self.min_length}")
        if self.max_length <= self.min_length:
            errors.append(
                f"max_length ({self.max_length}) must be > min_length ({self.min_length})"
            )
        if self.threads < 1:
            errors.append(f"threads must be >= 1, got {self.threads}")
        if self.itr_5_end <= self.itr_5_start:
            errors.append(
                f"itr_5_end ({self.itr_5_end}) must be > itr_5_start ({self.itr_5_start})"
            )
        if self.itr_3_end <= self.itr_3_start:
            errors.append(
                f"itr_3_end ({self.itr_3_end}) must be > itr_3_start ({self.itr_3_start})"
            )
        if self.itr_3_start < self.itr_5_end:
            errors.append(
                f"itr_3_start ({self.itr_3_start}) must be >= itr_5_end ({self.itr_5_end})"
            )
        if self.transgene_start < self.itr_5_end:
            errors.append(
                f"transgene_start ({self.transgene_start}) must be >= itr_5_end ({self.itr_5_end})"
            )
        if self.transgene_end > self.itr_3_start:
            errors.append(
                f"transgene_end ({self.transgene_end}) must be <= itr_3_start ({self.itr_3_start})"
            )
        return errors


# ============================================================
# Utility functions
# ============================================================

def setup_logging(log_file: Path) -> None:
    """Configure logging to file and stdout."""
    log_file.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout),
        ],
    )


def check_tool_available(tool_name: str) -> bool:
    """Check if a command-line tool is available on PATH."""
    return shutil.which(tool_name) is not None


def validate_dependencies(cfg: PipelineConfig) -> list:
    """Check that all required external tools are installed."""
    missing = []
    required = ["minimap2", "samtools"]
    if cfg.run_porechop:
        required.append("porechop")
    if cfg.run_nanofilt:
        required.append("NanoFilt")
    if cfg.run_nanoplot:
        required.append("NanoPlot")

    for tool in required:
        if not check_tool_available(tool):
            missing.append(tool)
    return missing


def run_command(cmd, log_prefix=None, shell=False) -> str:
    """Execute a shell command with logging. Raises RuntimeError on failure."""
    if isinstance(cmd, list):
        cmd_display = " ".join(cmd)
    else:
        cmd_display = cmd

    if log_prefix:
        logging.info("%s: %s", log_prefix, cmd_display)
    else:
        logging.info("Running: %s", cmd_display)

    process = subprocess.Popen(
        cmd,
        shell=shell,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        bufsize=1,
    )

    stdout_lines = []
    stderr_lines = []

    for line in process.stdout:
        line = line.rstrip("\n")
        stdout_lines.append(line)
        logging.debug(line)

    for line in process.stderr:
        line = line.rstrip("\n")
        stderr_lines.append(line)
        logging.debug(line)

    process.wait()

    if process.returncode != 0:
        error_output = "\n".join(stderr_lines[-20:])
        raise RuntimeError(
            f"Command failed (exit {process.returncode}): {cmd_display}\n"
            f"Last stderr:\n{error_output}"
        )

    return "\n".join(stdout_lines + stderr_lines)


def count_fastq_reads(fastq_path: Path) -> int:
    """Count reads in a FASTQ file efficiently."""
    count = 0
    with pysam.FastxFile(str(fastq_path)) as fh:
        for _ in fh:
            count += 1
    return count


def write_json(data, path: Path) -> None:
    """Write data to a JSON file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as handle:
        json.dump(data, handle, indent=2, default=str)


def is_host_reference(ref_name: str, patterns) -> bool:
    """Check if a reference name matches host genome patterns."""
    if ref_name is None:
        return False
    for pattern in patterns:
        if re.match(pattern, ref_name, flags=re.IGNORECASE):
            return True
    return False


def is_helper_reference(ref_name: str, patterns) -> bool:
    """Check if a reference name matches helper/Rep-Cap plasmid patterns."""
    if ref_name is None:
        return False
    for pattern in patterns:
        if re.match(pattern, ref_name, flags=re.IGNORECASE):
            return True
    return False


def parse_cigar_ops(cigarstring: str) -> list:
    """Parse a CIGAR string into (length, op) tuples."""
    if not cigarstring:
        return []
    return [
        (int(length), op)
        for length, op in re.findall(r"(\d+)([MIDNSHP=X])", cigarstring)
    ]


def query_consuming_length(cigarstring: str) -> int:
    """Total query-consuming bases from a CIGAR string."""
    total = 0
    for length, op in parse_cigar_ops(cigarstring):
        if op in {"M", "I", "S", "=", "X"}:
            total += length
    return total


def aligned_reference_length(cigarstring: str) -> int:
    """Total reference-consuming bases from a CIGAR string."""
    total = 0
    for length, op in parse_cigar_ops(cigarstring):
        if op in {"M", "D", "N", "=", "X"}:
            total += length
    return total


def softclip_left(cigarstring: str) -> int:
    """Get left soft-clip length."""
    ops = parse_cigar_ops(cigarstring)
    if ops and ops[0][1] == "S":
        return ops[0][0]
    return 0


def softclip_right(cigarstring: str) -> int:
    """Get right soft-clip length."""
    ops = parse_cigar_ops(cigarstring)
    if ops and ops[-1][1] == "S":
        return ops[-1][0]
    return 0


def get_query_span_from_cigar(cigarstring: str) -> tuple:
    """
    Get query start and end (0-based, half-open) from CIGAR.
    Accounts for soft-clipping at both ends.
    """
    total_query = query_consuming_length(cigarstring)
    left = softclip_left(cigarstring)
    right = softclip_right(cigarstring)
    qstart = left
    qend = max(left, total_query - right)
    return qstart, qend


def gc_percent(seq: str) -> float:
    """Calculate GC content percentage."""
    if not seq:
        return 0.0
    gc = seq.count("G") + seq.count("C") + seq.count("g") + seq.count("c")
    return (gc / len(seq)) * 100.0


def open_maybe_gzip(path: Path, mode="rt"):
    """Open a file, transparently handling gzip compression."""
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


# ============================================================
# Reference handling
# ============================================================

def list_available_references(refs_dir: Path) -> dict:
    """
    Scan all FASTA files in the reference directory and list
    all sequence IDs with their lengths. Helps the user identify
    the correct --transgene-name value.
    """
    ref_files = sorted([
        p for p in refs_dir.iterdir()
        if p.suffix.lower() in {".fa", ".fasta"}
    ])

    all_sequences = {}

    for ref_file in ref_files:
        current_name = None
        current_len = 0
        with open(ref_file) as handle:
            for line in handle:
                if line.startswith(">"):
                    if current_name is not None:
                        all_sequences[current_name] = {
                            "length": current_len,
                            "source_file": ref_file.name,
                        }
                    current_name = line[1:].strip().split()[0]
                    current_len = 0
                else:
                    current_len += len(line.strip())
        if current_name is not None:
            all_sequences[current_name] = {
                "length": current_len,
                "source_file": ref_file.name,
            }

    return all_sequences


def fasta_lengths(fasta_path: Path) -> dict:
    """Parse a FASTA file and return {name: length} dict."""
    lengths = {}
    current_name = None
    current_len = 0
    with open(fasta_path) as handle:
        for line in handle:
            if line.startswith(">"):
                if current_name is not None:
                    lengths[current_name] = current_len
                current_name = line[1:].strip().split()[0]
                current_len = 0
            else:
                current_len += len(line.strip())
    if current_name is not None:
        lengths[current_name] = current_len
    return lengths


def mask_itr_regions(
    sequence: str,
    cfg: PipelineConfig,
) -> str:
    """
    Mask the variable (flip/flop) region of ITRs with N characters.
    Uses absolute coordinates from cfg to mask both 5' and 3' ITR variable regions.
    """
    seq_list = list(sequence)
    seq_len = len(seq_list)

    # Mask 5' ITR variable region (absolute coordinates)
    mask_5_start = cfg.itr_5_start + cfg.itr_variable_start
    mask_5_end = cfg.itr_5_start + cfg.itr_variable_end
    for i in range(mask_5_start, min(mask_5_end, seq_len)):
        seq_list[i] = "N"

    # Mask 3' ITR variable region (absolute coordinates)
    mask_3_start = cfg.itr_3_start + cfg.itr_variable_start
    mask_3_end = cfg.itr_3_start + cfg.itr_variable_end
    for i in range(mask_3_start, min(mask_3_end, seq_len)):
        seq_list[i] = "N"

    return "".join(seq_list)


def combine_references(refs_dir: Path, output_fasta: Path, cfg: PipelineConfig) -> dict:
    """
    Combine .fa/.fasta files into one FASTA.
    - Human genome files: only conventional chromosomes are kept.
    - Transgene plasmid: optionally mask variable ITR regions.
    Returns {name: length} dict.
    """
    conventional = re.compile(
        r"^chr([1-9]|1[0-9]|2[0-2]|X|Y|M|MT)$", re.IGNORECASE
    )
    ref_files = sorted([
        p for p in refs_dir.iterdir()
        if p.suffix.lower() in {".fa", ".fasta"}
        and p.name != "combined_reference.fasta"
    ])
    if not ref_files:
        raise FileNotFoundError(f"No FASTA files found in {refs_dir}")

    logging.info("Combining %d reference FASTA files into %s", len(ref_files), output_fasta)
    output_fasta.parent.mkdir(parents=True, exist_ok=True)

    with open(output_fasta, "w") as out:
        for ref_file in ref_files:
            is_human = (
                "hg38" in ref_file.name.lower()
                or "grch38" in ref_file.name.lower()
            )
            keep_current = False
            current_name = None
            current_seq_lines = []

            with open(ref_file) as inp:
                for line in inp:
                    if line.startswith(">"):
                        # Write previous sequence if kept
                        if keep_current and current_name:
                            seq = "".join(current_seq_lines)
                            # Apply ITR masking if applicable
                            if (cfg.mask_itr_variable_region
                                    and cfg.transgene_name
                                    and current_name == cfg.transgene_name):
                                seq = mask_itr_regions(seq, cfg)
                            out.write(f">{current_name}\n")
                            for i in range(0, len(seq), 80):
                                out.write(seq[i:i + 80] + "\n")

                        current_name = line[1:].strip().split()[0]
                        current_seq_lines = []
                        if is_human:
                            keep_current = bool(conventional.match(current_name))
                        else:
                            keep_current = True
                    else:
                        if keep_current:
                            current_seq_lines.append(line.strip())

            # Write last sequence in file
            if keep_current and current_name:
                seq = "".join(current_seq_lines)
                if (cfg.mask_itr_variable_region
                        and cfg.transgene_name
                        and current_name == cfg.transgene_name):
                    seq = mask_itr_regions(seq, cfg)
                out.write(f">{current_name}\n")
                for i in range(0, len(seq), 80):
                    out.write(seq[i:i + 80] + "\n")

    return fasta_lengths(output_fasta)


# ============================================================
# FASTQ processing
# ============================================================

def run_porechop(input_fastq: Path, output_fastq: Path) -> dict:
    """Run Porechop for adapter trimming with middle-adapter discard."""
    output_fastq.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        "porechop",
        "-i", str(input_fastq),
        "-o", str(output_fastq),
        "--discard_middle",
    ]
    output = run_command(cmd, log_prefix=f"Porechop {input_fastq.name}")

    stats = {
        "input_reads": None,
        "discarded_middle_adapter_reads": None,
        "retained_reads": None,
    }

    m = re.search(
        r"(\d+)\s*/\s*(\d+)\s+reads were discarded based on middle adapters",
        output,
    )
    if m:
        discarded = int(m.group(1))
        total = int(m.group(2))
        stats["input_reads"] = total
        stats["discarded_middle_adapter_reads"] = discarded
        stats["retained_reads"] = total - discarded
    else:
        try:
            stats["input_reads"] = count_fastq_reads(input_fastq)
            stats["retained_reads"] = count_fastq_reads(output_fastq)
            stats["discarded_middle_adapter_reads"] = (
                stats["input_reads"] - stats["retained_reads"]
            )
        except Exception:
            pass

    return stats


def run_nanofilt(
    input_fastq: Path,
    output_fastq: Path,
    min_q: int,
    min_len: int,
    max_len: int,
) -> dict:
    """Run NanoFilt for quality and length filtering."""
    output_fastq.parent.mkdir(parents=True, exist_ok=True)

    cmd = (
        f"NanoFilt -q {min_q} -l {min_len} --maxlength {max_len} "
        f"< '{input_fastq}' > '{output_fastq}'"
    )
    run_command(cmd, log_prefix=f"NanoFilt {input_fastq.name}", shell=True)

    input_reads = count_fastq_reads(input_fastq)
    output_reads = count_fastq_reads(output_fastq)

    return {
        "input_reads": input_reads,
        "retained_reads": output_reads,
        "discarded_reads": input_reads - output_reads,
        "min_quality": min_q,
        "min_length": min_len,
        "max_length": max_len,
    }


def run_nanoplot(input_fastq: Path, output_dir: Path, title: str) -> None:
    """Run NanoPlot for read quality visualisation."""
    output_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        "NanoPlot",
        "--plots", "dot",
        "--fastq", str(input_fastq),
        "--outdir", str(output_dir),
        "--title", title,
    ]
    run_command(cmd, log_prefix=f"NanoPlot {input_fastq.name}")


# ============================================================
# Alignment (optimised — piped, no intermediate SAM)
# ============================================================

def align_fastq_to_reference(
    fastq_file: Path,
    combined_ref: Path,
    output_dir: Path,
    cfg: PipelineConfig,
) -> Path:
    """
    Align FASTQ to combined reference using minimap2 piped into samtools sort.
    No intermediate SAM/unsorted BAM is written to disc.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    prefix = fastq_file.stem
    sorted_bam = output_dir / f"{prefix}.sorted.bam"

    pipe_cmd = (
        f"minimap2 -t {cfg.threads} -ax {cfg.minimap2_preset} {cfg.minimap2_extra} "
        f"'{combined_ref}' '{fastq_file}' "
        f"| samtools sort -@ {cfg.threads} -o '{sorted_bam}'"
    )
    run_command(pipe_cmd, log_prefix=f"Align {fastq_file.name}", shell=True)

    pysam.index(str(sorted_bam))

    return sorted_bam


# ============================================================
# SA tag and chimeric read utilities
# ============================================================

def parse_sa_tag(sa_tag: str) -> list:
    """
    Parse SA auxiliary tag.
    Format: ref,pos,strand,CIGAR,mapQ,NM;...
    """
    entries = []
    for raw in sa_tag.strip(";").split(";"):
        if not raw:
            continue
        parts = raw.split(",")
        if len(parts) < 6:
            continue
        ref, pos, strand, cigar, mapq, nm = parts[:6]
        entries.append({
            "ref": ref,
            "pos1": int(pos),
            "strand": strand,
            "cigar": cigar,
            "mapq": int(mapq),
            "nm": int(nm),
        })
    return entries


def analyse_junction(left_seq: str, right_seq: str, max_window: int = 20) -> dict:
    """
    Analyse breakpoint junction for microhomology.
    Reports longest suffix/prefix match.
    """
    mh = ""
    max_len = min(max_window, len(left_seq), len(right_seq))
    for i in range(1, max_len + 1):
        if left_seq[-i:] == right_seq[:i]:
            mh = left_seq[-i:]

    return {
        "microhomology_length": len(mh),
        "microhomology_seq": mh,
        "left_flank": left_seq,
        "right_flank": right_seq,
    }


def merge_same_ref_neighbours(segments: list) -> list:
    """Merge adjacent segments on the same reference and strand."""
    if not segments:
        return []

    merged = [segments[0].copy()]
    for seg in segments[1:]:
        prev = merged[-1]
        same_ref = prev["ref"] == seg["ref"]
        same_strand = prev["strand"] == seg["strand"]

        if same_ref and same_strand and seg["qstart"] <= prev["qend"] + 20:
            prev["qstart"] = min(prev["qstart"], seg["qstart"])
            prev["qend"] = max(prev["qend"], seg["qend"])
            prev["ref_start"] = min(prev["ref_start"], seg["ref_start"])
            prev["ref_end"] = max(prev["ref_end"], seg["ref_end"])
            prev["mapq"] = max(prev["mapq"], seg["mapq"])
            prev["aligned_length"] += seg["aligned_length"]
        else:
            merged.append(seg.copy())

    return merged


# ============================================================
# Genome structure classification
# ============================================================

class GenomeCategory:
    """Constants for genome structure categories."""
    FULL_SS_AAV = "full_ssAAV"
    ICG = "incomplete_genome"
    GDM = "genome_deletion_mutant"
    FULL_SC_AAV = "full_scAAV"
    SNAPBACK_SYMMETRIC = "snapback_symmetric"
    SNAPBACK_ASYMMETRIC = "snapback_asymmetric"
    BACKBONE = "backbone_contamination"
    ITR_ONLY = "itr_only"
    HOST_CONTAMINATION = "host_contamination"
    HELPER_CONTAMINATION = "helper_contamination"
    CHIMERIC = "chimeric"
    UNCLASSIFIED = "unclassified"


def classify_genome_structure(
    segments: list,
    read_length: int,
    ref_lengths: dict,
    cfg: PipelineConfig,
) -> dict:
    """
    Classify a read into an AAV genome structure category based on
    alignment coordinates relative to the annotated plasmid regions.

    Regions:
      [itr_5_start, itr_5_end)              = 5' ITR
      [transgene_start, transgene_end)       = Transgene cassette
      [itr_3_start, itr_3_end)              = 3' ITR
      [0, itr_5_start) + [itr_3_end, ref_len) = Backbone
    """
    if not segments:
        return {
            "category": GenomeCategory.UNCLASSIFIED,
            "subcategory": "no_alignments",
            "confidence": 0.0,
        }

    transgene_ref = cfg.transgene_name

    # Separate segments by reference type
    transgene_segments = [s for s in segments if s["ref"] == transgene_ref]
    host_segments = [
        s for s in segments
        if is_host_reference(s["ref"], cfg.host_ref_patterns)
    ]
    helper_segments = [
        s for s in segments
        if is_helper_reference(s["ref"], cfg.helper_ref_patterns)
    ]

    refs_involved = set(s["ref"] for s in segments)

    # --- Pure host contamination ---
    if host_segments and not transgene_segments:
        return {
            "category": GenomeCategory.HOST_CONTAMINATION,
            "subcategory": "host_only",
            "confidence": 0.9,
        }

    # --- Pure helper contamination ---
    if helper_segments and not transgene_segments:
        return {
            "category": GenomeCategory.HELPER_CONTAMINATION,
            "subcategory": "helper_only",
            "confidence": 0.9,
        }

    # --- Chimeric (multi-reference) ---
    if len(refs_involved) > 1 and transgene_segments:
        if host_segments:
            return {
                "category": GenomeCategory.CHIMERIC,
                "subcategory": "host_vector_chimera",
                "confidence": 0.85,
            }
        return {
            "category": GenomeCategory.CHIMERIC,
            "subcategory": "multi_reference",
            "confidence": 0.7,
        }

    # --- All segments on transgene/plasmid reference ---
    if transgene_segments and len(refs_involved) == 1:
        ref_len = ref_lengths.get(transgene_ref, 0)
        if ref_len == 0:
            return {
                "category": GenomeCategory.UNCLASSIFIED,
                "subcategory": "missing_ref_length",
                "confidence": 0.0,
            }

        # Determine the total aligned span on the reference
        all_ref_starts = [s["ref_start"] for s in transgene_segments]
        all_ref_ends = [s["ref_end"] for s in transgene_segments]
        leftmost = min(all_ref_starts)
        rightmost = max(all_ref_ends)

        # --- Backbone involvement ---
        backbone_regions = cfg.backbone_regions(ref_len)
        backbone_bp = 0
        for bb_start, bb_end in backbone_regions:
            for seg in transgene_segments:
                overlap_start = max(seg["ref_start"], bb_start)
                overlap_end = min(seg["ref_end"], bb_end)
                if overlap_end > overlap_start:
                    backbone_bp += (overlap_end - overlap_start)

        if backbone_bp > 50:
            return {
                "category": GenomeCategory.BACKBONE,
                "subcategory": f"backbone_{backbone_bp}bp",
                "confidence": 0.9,
            }

        # --- Classify within the expected genome region ---
        genome_length = cfg.expected_genome_length

        # ITR coverage calculations
        if leftmost >= cfg.itr_5_end:
            itr_5_missing = cfg.itr_length
        else:
            itr_5_missing = max(0, leftmost - cfg.itr_5_start)

        if rightmost <= cfg.itr_3_start:
            itr_3_missing = cfg.itr_length
        else:
            itr_3_missing = max(0, cfg.itr_3_end - rightmost)

        covers_5_itr = itr_5_missing < cfg.itr_full_length_threshold
        covers_3_itr = itr_3_missing < cfg.itr_full_length_threshold

        # Aligned fraction of expected genome
        genome_start = cfg.expected_genome_start
        genome_end = cfg.expected_genome_end
        aligned_within_genome = 0
        for seg in transgene_segments:
            seg_start = max(seg["ref_start"], genome_start)
            seg_end = min(seg["ref_end"], genome_end)
            if seg_end > seg_start:
                aligned_within_genome += (seg_end - seg_start)

        genome_coverage_frac = (
            aligned_within_genome / genome_length if genome_length > 0 else 0
        )

        # Check for self-complementary / snapback (inverted segments)
        if len(transgene_segments) >= 2:
            strands = [s["strand"] for s in transgene_segments]
            has_inversion = len(set(strands)) > 1

            if has_inversion:
                total_aligned = sum(
                    s["aligned_length"] for s in transgene_segments
                )
                coverage_ratio = (
                    total_aligned / genome_length if genome_length > 0 else 0
                )

                if coverage_ratio >= 1.5:
                    return {
                        "category": GenomeCategory.FULL_SC_AAV,
                        "subcategory": "self_complementary",
                        "confidence": 0.8,
                    }

                plus_len = sum(
                    s["aligned_length"] for s in transgene_segments
                    if s["strand"] == "+"
                )
                minus_len = sum(
                    s["aligned_length"] for s in transgene_segments
                    if s["strand"] == "-"
                )
                denom = max(plus_len, minus_len)
                ratio = min(plus_len, minus_len) / denom if denom > 0 else 0

                if ratio >= 0.8:
                    return {
                        "category": GenomeCategory.SNAPBACK_SYMMETRIC,
                        "subcategory": f"strand_ratio_{ratio:.2f}",
                        "confidence": 0.75,
                    }
                else:
                    return {
                        "category": GenomeCategory.SNAPBACK_ASYMMETRIC,
                        "subcategory": f"strand_ratio_{ratio:.2f}",
                        "confidence": 0.75,
                    }

        # ITR-only fragment
        total_aligned = sum(s["aligned_length"] for s in transgene_segments)
        if total_aligned < cfg.itr_length * 3:
            if covers_5_itr or covers_3_itr:
                return {
                    "category": GenomeCategory.ITR_ONLY,
                    "subcategory": "itr_fragment",
                    "confidence": 0.7,
                }

        # Full ssAAV
        if genome_coverage_frac >= 0.85 and covers_5_itr and covers_3_itr:
            return {
                "category": GenomeCategory.FULL_SS_AAV,
                "subcategory": f"coverage_{genome_coverage_frac:.2f}",
                "confidence": 0.9,
            }

        # Genome Deletion Mutant: both ITRs present but large internal gap
        if covers_5_itr and covers_3_itr and genome_coverage_frac < 0.85:
            return {
                "category": GenomeCategory.GDM,
                "subcategory": f"coverage_{genome_coverage_frac:.2f}",
                "confidence": 0.75,
            }

        # Incomplete Genome: missing one or both ITRs
        if not covers_5_itr or not covers_3_itr:
            missing = []
            if not covers_5_itr:
                missing.append("5prime")
            if not covers_3_itr:
                missing.append("3prime")
            return {
                "category": GenomeCategory.ICG,
                "subcategory": f"missing_{'_and_'.join(missing)}",
                "confidence": 0.8,
            }

    return {
        "category": GenomeCategory.UNCLASSIFIED,
        "subcategory": "unresolved",
        "confidence": 0.0,
    }


# ============================================================
# ITR truncation detection (coordinate-based)
# ============================================================

def detect_itr_truncation(
    ref_start: int,
    ref_end: int,
    ref_name: str,
    read_name: str,
    read_length: int,
    mapq: int,
    ref_lengths: dict,
    cfg: PipelineConfig,
) -> Optional[dict]:
    """
    Detect ITR truncation using explicit coordinate annotations.

    Uses cfg.itr_5_start, cfg.itr_5_end, cfg.itr_3_start, cfg.itr_3_end
    to determine how much of each ITR is covered by the alignment.
    Also detects backbone read-through.
    """
    if ref_name != cfg.transgene_name:
        return None

    ref_len = ref_lengths.get(ref_name)
    if not ref_len:
        return None

    # 5' ITR coverage
    if ref_start >= cfg.itr_5_end:
        itr_5_missing = cfg.itr_length
        itr_5_coverage = 0
    elif ref_start <= cfg.itr_5_start:
        itr_5_missing = 0
        itr_5_coverage = min(cfg.itr_5_end, ref_end) - cfg.itr_5_start
    else:
        itr_5_missing = ref_start - cfg.itr_5_start
        itr_5_coverage = min(cfg.itr_5_end, ref_end) - ref_start

    # 3' ITR coverage
    if ref_end <= cfg.itr_3_start:
        itr_3_missing = cfg.itr_length
        itr_3_coverage = 0
    elif ref_end >= cfg.itr_3_end:
        itr_3_missing = 0
        itr_3_coverage = cfg.itr_3_end - max(cfg.itr_3_start, ref_start)
    else:
        itr_3_missing = cfg.itr_3_end - ref_end
        itr_3_coverage = ref_end - max(cfg.itr_3_start, ref_start)

    # Backbone involvement
    backbone_regions = cfg.backbone_regions(ref_len)
    backbone_bp = 0
    for bb_start, bb_end in backbone_regions:
        overlap_start = max(ref_start, bb_start)
        overlap_end = min(ref_end, bb_end)
        if overlap_end > overlap_start:
            backbone_bp += (overlap_end - overlap_start)

    # Categorise
    threshold = cfg.itr_full_length_threshold
    is_5_truncated = itr_5_missing > threshold
    is_3_truncated = itr_3_missing > threshold

    if is_5_truncated and is_3_truncated:
        category = "both_ends_truncated"
    elif is_5_truncated:
        category = "5_prime_truncated"
    elif is_3_truncated:
        category = "3_prime_truncated"
    else:
        category = "full_length"

    # Add backbone flag
    if backbone_bp > 0:
        category += "_with_backbone"

    return {
        "read_name": read_name,
        "reference": ref_name,
        "read_length": read_length,
        "ref_start": ref_start,
        "ref_end": ref_end,
        "ref_length": ref_len,
        "itr_5_missing": itr_5_missing,
        "itr_3_missing": itr_3_missing,
        "itr_5_coverage": max(0, itr_5_coverage),
        "itr_3_coverage": max(0, itr_3_coverage),
        "backbone_bp": backbone_bp,
        "category": category,
        "full_length": int(category == "full_length"),
        "mapq": mapq,
    }


# ============================================================
# Single-pass BAM analysis
# ============================================================

def single_pass_bam_analysis(
    bam_path: Path,
    ref_lengths: dict,
    cfg: PipelineConfig,
) -> dict:
    """
    Perform all BAM-based analyses in a single pass for efficiency.

    Returns a comprehensive results dict containing:
    - Mapping statistics
    - Chimeric reads and breakpoints
    - ITR truncation data
    - Genome structure classification
    - Per-base coverage (for transgene)
    - Strand bias
    """
    mapped_reads = set()
    unmapped_read_names = set()
    ref_counts = Counter()
    mapq_buckets = defaultdict(lambda: Counter())
    coverage_rows = []
    truncation_rows = []
    chimeric_reads = {}
    host_chimeric_reads = {}
    breakpoint_rows = []
    genome_classifications = []
    strand_counts = Counter()

    # Per-base coverage for transgene
    transgene_ref_len = ref_lengths.get(cfg.transgene_name, 0)
    per_base_coverage = [0] * transgene_ref_len if transgene_ref_len > 0 else []

    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam:
            if read.is_secondary:
                continue

            if read.is_unmapped:
                unmapped_read_names.add(read.query_name)
                continue

            if read.is_supplementary:
                continue

            if not read.query_sequence:
                continue

            read_name = read.query_name
            mapped_reads.add(read_name)
            ref = read.reference_name
            ref_counts[ref] += 1

            # MAPQ bucketing
            mq = read.mapping_quality
            if mq >= 30:
                mapq_buckets[ref]["high"] += 1
            elif mq >= 10:
                mapq_buckets[ref]["medium"] += 1
            else:
                mapq_buckets[ref]["low"] += 1

            # Coverage calculation
            qlen = read.query_length or 0
            aln_len = read.query_alignment_length or 0
            qcov = (aln_len / qlen) if qlen > 0 else 0
            rlen = ref_lengths.get(ref, 0)
            rcov = (aln_len / rlen) if rlen > 0 else 0

            coverage_rows.append({
                "read_name": read_name,
                "reference": ref,
                "query_coverage": round(qcov, 4),
                "reference_coverage": round(rcov, 4),
                "mapq": mq,
            })

            # Per-base coverage for transgene reference
            if ref == cfg.transgene_name and per_base_coverage:
                for pos in read.get_reference_positions():
                    if 0 <= pos < transgene_ref_len:
                        per_base_coverage[pos] += 1

            # Strand bias (for transgene)
            if ref == cfg.transgene_name:
                strand_counts["-" if read.is_reverse else "+"] += 1

            # ITR truncation detection
            trunc = detect_itr_truncation(
                ref_start=read.reference_start,
                ref_end=read.reference_end,
                ref_name=ref,
                read_name=read_name,
                read_length=qlen,
                mapq=mq,
                ref_lengths=ref_lengths,
                cfg=cfg,
            )
            if trunc is not None:
                truncation_rows.append(trunc)

            # --- Chimeric read detection ---
            segments = []
            if read.cigarstring:
                qstart, qend = get_query_span_from_cigar(read.cigarstring)
                seg_aligned_len = aln_len if aln_len > 0 else (qend - qstart)
                if (seg_aligned_len >= cfg.min_chimeric_segment_length
                        and mq >= cfg.min_chimeric_mapq):
                    segments.append({
                        "ref": ref,
                        "ref_start": read.reference_start,
                        "ref_end": read.reference_end,
                        "qstart": qstart,
                        "qend": qend,
                        "strand": "-" if read.is_reverse else "+",
                        "mapq": mq,
                        "cigar": read.cigarstring,
                        "aligned_length": seg_aligned_len,
                    })

            # Parse SA tag for supplementary alignments
            if read.has_tag("SA"):
                for sa in parse_sa_tag(read.get_tag("SA")):
                    qstart, qend = get_query_span_from_cigar(sa["cigar"])
                    sa_aligned_len = sum(
                        length for length, op in parse_cigar_ops(sa["cigar"])
                        if op in {"M", "=", "X"}
                    )
                    if (sa_aligned_len >= cfg.min_chimeric_segment_length
                            and sa["mapq"] >= cfg.min_chimeric_mapq):
                        ref_start0 = sa["pos1"] - 1
                        ref_end0 = ref_start0 + aligned_reference_length(sa["cigar"])
                        segments.append({
                            "ref": sa["ref"],
                            "ref_start": ref_start0,
                            "ref_end": ref_end0,
                            "qstart": qstart,
                            "qend": qend,
                            "strand": sa["strand"],
                            "mapq": sa["mapq"],
                            "cigar": sa["cigar"],
                            "aligned_length": sa_aligned_len,
                        })

            # Genome structure classification
            classification = classify_genome_structure(
                segments=segments,
                read_length=qlen,
                ref_lengths=ref_lengths,
                cfg=cfg,
            )
            classification["read_name"] = read_name
            classification["read_length"] = qlen
            genome_classifications.append(classification)

            # Continue chimeric analysis only if multi-segment
            if len(segments) < 2:
                continue

            segments = sorted(segments, key=lambda x: (x["qstart"], x["qend"]))
            segments = merge_same_ref_neighbours(segments)

            if len(segments) < 2:
                continue

            refs = [s["ref"] for s in segments]
            if len(set(refs)) < 2:
                continue

            query_seq = read.query_sequence
            junctions = []

            for i in range(len(segments) - 1):
                left = segments[i]
                right = segments[i + 1]

                left_break = left["qend"]
                right_break = right["qstart"]

                if right_break > left_break:
                    inserted_seq = query_seq[left_break:right_break]
                else:
                    inserted_seq = ""

                left_flank = query_seq[max(0, left_break - 20):left_break]
                right_flank = query_seq[
                    right_break:min(len(query_seq), right_break + 20)
                ]
                junction = analyse_junction(left_flank, right_flank, max_window=20)

                junction_row = {
                    "read_name": read_name,
                    "junction_index": i + 1,
                    "left_ref": left["ref"],
                    "left_ref_start": left["ref_start"],
                    "left_ref_end": left["ref_end"],
                    "left_query_end": left_break,
                    "left_strand": left["strand"],
                    "right_ref": right["ref"],
                    "right_ref_start": right["ref_start"],
                    "right_ref_end": right["ref_end"],
                    "right_query_start": right_break,
                    "right_strand": right["strand"],
                    "query_gap_bp": max(0, right_break - left_break),
                    "inserted_seq": inserted_seq,
                    "microhomology_length": junction["microhomology_length"],
                    "microhomology_seq": junction["microhomology_seq"],
                    "left_flank_20bp": junction["left_flank"],
                    "right_flank_20bp": junction["right_flank"],
                }

                junctions.append(junction_row)
                breakpoint_rows.append(junction_row)

            summary = {
                "read_name": read_name,
                "read_length": len(query_seq),
                "num_segments": len(segments),
                "refs_order": [s["ref"] for s in segments],
                "segment_lengths": [s["aligned_length"] for s in segments],
                "ref_starts": [s["ref_start"] for s in segments],
                "ref_ends": [s["ref_end"] for s in segments],
                "query_starts": [s["qstart"] for s in segments],
                "query_ends": [s["qend"] for s in segments],
                "strands": [s["strand"] for s in segments],
                "mapqs": [s["mapq"] for s in segments],
                "host_involved": any(
                    is_host_reference(s["ref"], cfg.host_ref_patterns)
                    for s in segments
                ),
                "segments": segments,
                "junctions": junctions,
                "sequence": query_seq,
            }

            chimeric_reads[read_name] = summary
            if summary["host_involved"]:
                host_chimeric_reads[read_name] = summary

    # Aggregate genome classification counts
    classification_counts = Counter(c["category"] for c in genome_classifications)

    # Strand bias summary
    plus = strand_counts.get("+", 0)
    minus = strand_counts.get("-", 0)
    total_stranded = plus + minus
    strand_bias = {
        "plus_strand": plus,
        "minus_strand": minus,
        "total": total_stranded,
        "plus_fraction": round(plus / total_stranded, 4) if total_stranded > 0 else 0,
        "minus_fraction": round(minus / total_stranded, 4) if total_stranded > 0 else 0,
        "balanced": (
            abs(plus - minus) / total_stranded < 0.1 if total_stranded > 0 else None
        ),
    }

    return {
        "mapped_reads": mapped_reads,
        "unmapped_read_names": unmapped_read_names,
        "ref_counts": ref_counts,
        "mapq_buckets": mapq_buckets,
        "coverage_rows": coverage_rows,
        "truncation_rows": truncation_rows,
        "chimeric_reads": chimeric_reads,
        "host_chimeric_reads": host_chimeric_reads,
        "breakpoint_rows": breakpoint_rows,
        "genome_classifications": genome_classifications,
        "classification_counts": classification_counts,
        "per_base_coverage": per_base_coverage,
        "strand_bias": strand_bias,
    }


# ============================================================
# Unmapped read analysis (memory-efficient)
# ============================================================

def analyse_sequence_flags(seq: str) -> list:
    """Flag potential quality issues in a sequence."""
    flags = []
    if not seq:
        return ["no_sequence"]

    gc = gc_percent(seq)
    if gc > 65:
        flags.append("high_gc")
    elif gc < 35:
        flags.append("low_gc")

    if len(seq) < 100:
        flags.append("very_short")
    if len(seq) > 10000:
        flags.append("very_long")

    upper_seq = seq.upper()
    for base in "ACGT":
        if base * 10 in upper_seq:
            flags.append("homopolymer_10bp")
            break

    return flags if flags else ["no_issues"]


def analyse_unmapped_reads(
    fastq_path: Path,
    bam_path: Path,
    out_dir: Path,
) -> dict:
    """
    Memory-efficient unmapped read analysis.
    First pass: collect unmapped read names from BAM.
    Second pass: stream FASTQ and extract only those reads.
    """
    out_dir.mkdir(parents=True, exist_ok=True)

    mapped_names = set()
    unmapped_names = set()

    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam:
            if read.is_secondary or read.is_supplementary:
                continue
            if read.is_unmapped:
                unmapped_names.add(read.query_name)
            else:
                mapped_names.add(read.query_name)

    total_fastq = 0
    issue_counts = Counter()
    per_read_rows = []
    unmapped_fasta_path = out_dir / f"{fastq_path.stem}_unmapped.fasta"

    with pysam.FastxFile(str(fastq_path)) as fh, \
         open(unmapped_fasta_path, "w") as fasta_out:
        for rec in fh:
            total_fastq += 1
            if rec.name in unmapped_names:
                seq = rec.sequence or ""
                flags = analyse_sequence_flags(seq)
                for flag in flags:
                    issue_counts[flag] += 1
                per_read_rows.append({
                    "read_name": rec.name,
                    "length": len(seq),
                    "gc_percent": round(gc_percent(seq), 2),
                    "flags": "|".join(flags),
                })
                fasta_out.write(f">{rec.name}\n{seq}\n")

    all_bam_names = mapped_names | unmapped_names
    missing_count = max(0, total_fastq - len(all_bam_names))

    unmapped_csv = out_dir / f"{fastq_path.stem}_unmapped_reads.csv"
    with open(unmapped_csv, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["read_name", "length", "gc_percent", "flags"]
        )
        writer.writeheader()
        writer.writerows(per_read_rows)

    summary_data = {
        "total_fastq_reads": total_fastq,
        "mapped_reads": len(mapped_names),
        "unmapped_reads": len(unmapped_names),
        "missing_reads": missing_count,
        "issue_counts": dict(issue_counts),
    }
    summary_json = out_dir / f"{fastq_path.stem}_unmapped_summary.json"
    write_json(summary_data, summary_json)

    return {
        "unmapped_count": len(unmapped_names),
        "missing_count": missing_count,
        "issue_counts": dict(issue_counts),
        "unmapped_csv": str(unmapped_csv),
        "unmapped_fasta": str(unmapped_fasta_path),
        "summary_json": str(summary_json),
    }


# ============================================================
# Output writers
# ============================================================

def write_chimeric_outputs(
    sample_name: str,
    chimeric_reads: dict,
    host_chimeric_reads: dict,
    breakpoint_rows: list,
    output_dir: Path,
) -> dict:
    """Write all chimeric read analysis outputs."""
    output_dir.mkdir(parents=True, exist_ok=True)

    chimeric_csv = output_dir / f"{sample_name}_chimeric_reads.csv"
    breakpoints_csv = output_dir / f"{sample_name}_chimeric_breakpoints.csv"
    host_csv = output_dir / f"{sample_name}_host_chimeric_reads.csv"
    chimeric_fasta = output_dir / f"{sample_name}_chimeric_reads.fasta"

    fieldnames_chimeric = [
        "read_name", "read_length", "num_segments", "refs_order",
        "segment_lengths", "ref_starts", "ref_ends", "query_starts",
        "query_ends", "strands", "mapqs", "host_involved",
    ]
    with open(chimeric_csv, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames_chimeric)
        writer.writeheader()
        for row in chimeric_reads.values():
            writer.writerow({
                "read_name": row["read_name"],
                "read_length": row["read_length"],
                "num_segments": row["num_segments"],
                "refs_order": "|".join(map(str, row["refs_order"])),
                "segment_lengths": "|".join(map(str, row["segment_lengths"])),
                "ref_starts": "|".join(map(str, row["ref_starts"])),
                "ref_ends": "|".join(map(str, row["ref_ends"])),
                "query_starts": "|".join(map(str, row["query_starts"])),
                "query_ends": "|".join(map(str, row["query_ends"])),
                "strands": "|".join(map(str, row["strands"])),
                "mapqs": "|".join(map(str, row["mapqs"])),
                "host_involved": int(row["host_involved"]),
            })

    fieldnames_host = [
        "read_name", "read_length", "num_segments", "refs_order", "host_involved",
    ]
    with open(host_csv, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames_host)
        writer.writeheader()
        for row in host_chimeric_reads.values():
            writer.writerow({
                "read_name": row["read_name"],
                "read_length": row["read_length"],
                "num_segments": row["num_segments"],
                "refs_order": "|".join(map(str, row["refs_order"])),
                "host_involved": 1,
            })

    bp_fieldnames = [
        "read_name", "junction_index", "left_ref", "left_ref_start",
        "left_ref_end", "left_query_end", "left_strand", "right_ref",
        "right_ref_start", "right_ref_end", "right_query_start",
        "right_strand", "query_gap_bp", "inserted_seq",
        "microhomology_length", "microhomology_seq",
        "left_flank_20bp", "right_flank_20bp",
    ]
    with open(breakpoints_csv, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=bp_fieldnames)
        writer.writeheader()
        if breakpoint_rows:
            writer.writerows(breakpoint_rows)

    with open(chimeric_fasta, "w") as handle:
        for row in chimeric_reads.values():
            handle.write(f">{row['read_name']}\n{row['sequence']}\n")

    return {
        "chimeric_csv": str(chimeric_csv),
        "host_chimeric_csv": str(host_csv),
        "breakpoints_csv": str(breakpoints_csv),
        "chimeric_fasta": str(chimeric_fasta),
    }


def write_chimeric_bam(
    sample_name: str,
    input_bam: Path,
    chimeric_read_names: set,
    output_dir: Path,
) -> Path:
    """Extract chimeric reads into a separate sorted BAM."""
    output_dir.mkdir(parents=True, exist_ok=True)
    out_bam = output_dir / f"{sample_name}_chimeric.unsorted.bam"
    sorted_bam = output_dir / f"{sample_name}_chimeric.sorted.bam"

    with pysam.AlignmentFile(str(input_bam), "rb") as inp, \
         pysam.AlignmentFile(str(out_bam), "wb", template=inp) as out:
        for read in inp:
            if read.query_name in chimeric_read_names:
                out.write(read)

    pysam.sort("-o", str(sorted_bam), str(out_bam))
    pysam.index(str(sorted_bam))

    if out_bam.exists():
        out_bam.unlink()

    return sorted_bam


def write_per_base_coverage(
    sample_name: str,
    per_base_coverage: list,
    output_dir: Path,
) -> Path:
    """Write per-base coverage as a BEDGraph-like TSV."""
    output_dir.mkdir(parents=True, exist_ok=True)
    out_path = output_dir / f"{sample_name}_transgene_coverage.tsv"

    with open(out_path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["position", "coverage"])
        for i, cov in enumerate(per_base_coverage):
            writer.writerow([i, cov])

    return out_path


def write_genome_classifications(
    sample_name: str,
    classifications: list,
    output_dir: Path,
) -> Path:
    """Write per-read genome structure classification."""
    output_dir.mkdir(parents=True, exist_ok=True)
    out_path = output_dir / f"{sample_name}_genome_classifications.csv"

    fieldnames = ["read_name", "read_length", "category", "subcategory", "confidence"]
    with open(out_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in classifications:
            writer.writerow({
                "read_name": row.get("read_name", ""),
                "read_length": row.get("read_length", ""),
                "category": row["category"],
                "subcategory": row["subcategory"],
                "confidence": row["confidence"],
            })

    return out_path


# ============================================================
# HTML report generation
# ============================================================

def generate_html_report(
    sample_name: str,
    summary: dict,
    output_dir: Path,
) -> Path:
    """Generate a self-contained HTML report for a sample."""
    output_dir.mkdir(parents=True, exist_ok=True)
    html_path = output_dir / f"{sample_name}_report.html"

    total = summary.get("total_reads", 0)
    mapped = summary.get("mapped_reads", 0)
    unmapped = summary.get("unmapped_reads", 0)
    chimeric = summary.get("chimeric_reads", 0)
    host_chimeric = summary.get("host_chimeric_reads", 0)
    classification_counts = summary.get("classification_counts", {})
    strand_bias = summary.get("strand_bias", {})
    ref_counts = summary.get("references", {})
    truncation_summary = summary.get("truncation_summary", {})

    class_rows = ""
    for cat, count in sorted(classification_counts.items(), key=lambda x: -x[1]):
        pct = (count / total * 100) if total > 0 else 0
        class_rows += f"<tr><td>{cat}</td><td>{count:,}</td><td>{pct:.1f}%</td></tr>\n"

    ref_rows = ""
    for ref, count in sorted(ref_counts.items(), key=lambda x: -x[1]):
        pct = (count / total * 100) if total > 0 else 0
        ref_rows += f"<tr><td>{ref}</td><td>{count:,}</td><td>{pct:.1f}%</td></tr>\n"

    trunc_rows = ""
    for cat, count in sorted(truncation_summary.items(), key=lambda x: -x[1]):
        trunc_rows += f"<tr><td>{cat}</td><td>{count:,}</td></tr>\n"

    plus_frac = strand_bias.get("plus_fraction", 0)
    minus_frac = strand_bias.get("minus_fraction", 0)
    strand_status = "Balanced" if strand_bias.get("balanced") else "Biased"

    html_content = textwrap.dedent(f"""\
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>AAV-Chimera Report: {sample_name}</title>
        <style>
            * {{ margin: 0; padding: 0; box-sizing: border-box; }}
            body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
                   background: #f5f7fa; color: #2d3748; padding: 2rem; line-height: 1.6; }}
            .container {{ max-width: 1200px; margin: 0 auto; }}
            h1 {{ color: #1a365d; border-bottom: 3px solid #3182ce; padding-bottom: 0.5rem; margin-bottom: 1.5rem; }}
            h2 {{ color: #2c5282; margin: 1.5rem 0 0.75rem 0; border-left: 4px solid #3182ce; padding-left: 0.75rem; }}
            .stats-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 1rem; margin: 1rem 0; }}
            .stat-card {{ background: white; border-radius: 8px; padding: 1.25rem; box-shadow: 0 1px 3px rgba(0,0,0,0.1); text-align: center; }}
            .stat-card .value {{ font-size: 2rem; font-weight: 700; color: #2b6cb0; }}
            .stat-card .label {{ font-size: 0.875rem; color: #718096; text-transform: uppercase; }}
            table {{ width: 100%; border-collapse: collapse; background: white; border-radius: 8px; overflow: hidden; box-shadow: 0 1px 3px rgba(0,0,0,0.1); margin: 0.75rem 0; }}
            th, td {{ padding: 0.75rem 1rem; text-align: left; border-bottom: 1px solid #e2e8f0; }}
            th {{ background: #edf2f7; font-weight: 600; color: #4a5568; }}
            tr:hover {{ background: #f7fafc; }}
            .bar-container {{ background: #e2e8f0; border-radius: 4px; height: 24px; margin: 0.5rem 0; overflow: hidden; display: flex; }}
            .bar-segment {{ height: 100%; display: flex; align-items: center; justify-content: center; font-size: 0.75rem; color: white; font-weight: 600; }}
            .bar-mapped {{ background: #48bb78; }}
            .bar-unmapped {{ background: #fc8181; }}
            .note {{ background: #ebf8ff; border-left: 4px solid #3182ce; padding: 1rem; margin: 1rem 0; border-radius: 0 4px 4px 0; }}
            .timestamp {{ color: #a0aec0; font-size: 0.875rem; margin-top: 2rem; text-align: center; }}
        </style>
    </head>
    <body>
        <div class="container">
            <h1>AAV-Chimera Report: {sample_name}</h1>
            <h2>Summary Statistics</h2>
            <div class="stats-grid">
                <div class="stat-card"><div class="value">{total:,}</div><div class="label">Total Reads</div></div>
                <div class="stat-card"><div class="value">{mapped:,}</div><div class="label">Mapped Reads</div></div>
                <div class="stat-card"><div class="value">{unmapped:,}</div><div class="label">Unmapped Reads</div></div>
                <div class="stat-card"><div class="value">{chimeric:,}</div><div class="label">Chimeric Reads</div></div>
                <div class="stat-card"><div class="value">{host_chimeric:,}</div><div class="label">Host Chimeric</div></div>
            </div>
            <h2>Mapping Overview</h2>
            <div class="bar-container">
                <div class="bar-segment bar-mapped" style="width: {mapped/total*100 if total else 0:.1f}%">Mapped ({mapped/total*100 if total else 0:.1f}%)</div>
                <div class="bar-segment bar-unmapped" style="width: {unmapped/total*100 if total else 0:.1f}%">Unmapped ({unmapped/total*100 if total else 0:.1f}%)</div>
            </div>
            <h2>Genome Structure Classification</h2>
            <table><thead><tr><th>Category</th><th>Count</th><th>Percentage</th></tr></thead><tbody>{class_rows}</tbody></table>
            <h2>Reference Mapping Distribution</h2>
            <table><thead><tr><th>Reference</th><th>Count</th><th>Percentage</th></tr></thead><tbody>{ref_rows}</tbody></table>
            <h2>ITR Truncation Summary</h2>
            <table><thead><tr><th>Category</th><th>Count</th></tr></thead><tbody>{trunc_rows}</tbody></table>
            <h2>Strand Bias (Transgene)</h2>
            <div class="note">
                <strong>Plus strand:</strong> {strand_bias.get('plus_strand', 0):,} ({plus_frac*100:.1f}%) |
                <strong>Minus strand:</strong> {strand_bias.get('minus_strand', 0):,} ({minus_frac*100:.1f}%) |
                <strong>Status:</strong> {strand_status}
            </div>
            <div class="timestamp">Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</div>
        </div>
    </body>
    </html>
    """)

    with open(html_path, "w") as handle:
        handle.write(html_content)

    return html_path


# ============================================================
# Per-sample orchestration
# ============================================================

def calculate_sample_statistics(
    sample_name: str,
    fastq_path: Path,
    bam_path: Path,
    ref_lengths: dict,
    cfg: PipelineConfig,
    sample_out_dir: Path,
) -> dict:
    """Run all analyses for a single sample and write outputs."""
    sample_out_dir.mkdir(parents=True, exist_ok=True)

    total_reads = count_fastq_reads(fastq_path)

    logging.info("Running single-pass BAM analysis for %s", sample_name)
    results = single_pass_bam_analysis(bam_path, ref_lengths, cfg)

    mapped_reads = results["mapped_reads"]
    chimeric_reads = results["chimeric_reads"]
    host_chimeric_reads = results["host_chimeric_reads"]
    breakpoint_rows = results["breakpoint_rows"]
    truncation_rows = results["truncation_rows"]
    coverage_rows = results["coverage_rows"]
    genome_classifications = results["genome_classifications"]
    per_base_coverage = results["per_base_coverage"]

    # Write chimeric outputs
    chimeric_dir = sample_out_dir / "chimeric_reads"
    chim_files = write_chimeric_outputs(
        sample_name=sample_name,
        chimeric_reads=chimeric_reads,
        host_chimeric_reads=host_chimeric_reads,
        breakpoint_rows=breakpoint_rows,
        output_dir=chimeric_dir,
    )

    chimeric_bam = write_chimeric_bam(
        sample_name=sample_name,
        input_bam=bam_path,
        chimeric_read_names=set(chimeric_reads.keys()),
        output_dir=chimeric_dir,
    )

    # Unmapped read analysis
    unmapped_stats = analyse_unmapped_reads(
        fastq_path=fastq_path,
        bam_path=bam_path,
        out_dir=sample_out_dir / "unmapped_reads",
    )

    # Write coverage CSV
    coverage_csv = sample_out_dir / f"{sample_name}_coverage.csv"
    with open(coverage_csv, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["read_name", "reference", "query_coverage", "reference_coverage", "mapq"],
        )
        writer.writeheader()
        writer.writerows(coverage_rows)

    # Write per-base transgene coverage
    coverage_tsv = write_per_base_coverage(
        sample_name=sample_name,
        per_base_coverage=per_base_coverage,
        output_dir=sample_out_dir,
    )

    # Write genome classifications
    class_csv = write_genome_classifications(
        sample_name=sample_name,
        classifications=genome_classifications,
        output_dir=sample_out_dir,
    )

    # Write truncation CSV (updated fieldnames to match new dict keys)
    truncation_csv = sample_out_dir / f"{sample_name}_itr_truncation.csv"
    trunc_fieldnames = [
        "read_name", "reference", "read_length", "ref_start", "ref_end",
        "ref_length", "itr_5_missing", "itr_3_missing", "itr_5_coverage",
        "itr_3_coverage", "backbone_bp", "category", "full_length", "mapq",
    ]
    with open(truncation_csv, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=trunc_fieldnames)
        writer.writeheader()
        if truncation_rows:
            writer.writerows(truncation_rows)

    # Truncation summary
    truncation_summary = Counter(r["category"] for r in truncation_rows)

    # Concatemer/over-packaging detection
    concatemer_reads = [
        row for row in coverage_rows
        if row["reference"] == cfg.transgene_name
        and row["reference_coverage"] > 1.3
    ]

    # Build summary
    summary = {
        "sample": sample_name,
        "transgene": cfg.transgene_name,
        "fastq": str(fastq_path),
        "bam": str(bam_path),
        "total_reads": total_reads,
        "mapped_reads": len(mapped_reads),
        "unmapped_reads": total_reads - len(mapped_reads),
        "mapping_rate": round(len(mapped_reads) / total_reads, 4) if total_reads > 0 else 0,
        "references": dict(results["ref_counts"]),
        "mapq_buckets": {k: dict(v) for k, v in results["mapq_buckets"].items()},
        "chimeric_reads": len(chimeric_reads),
        "host_chimeric_reads": len(host_chimeric_reads),
        "breakpoints": len(breakpoint_rows),
        "classification_counts": dict(results["classification_counts"]),
        "truncation_summary": dict(truncation_summary),
        "strand_bias": results["strand_bias"],
        "concatemer_candidates": len(concatemer_reads),
        "unmapped_analysis": unmapped_stats,
        "outputs": {
            "coverage_csv": str(coverage_csv),
            "per_base_coverage_tsv": str(coverage_tsv),
            "truncation_csv": str(truncation_csv),
            "genome_classifications_csv": str(class_csv),
            "chimeric_bam": str(chimeric_bam),
            **chim_files,
        },
    }

    detailed_json = sample_out_dir / f"{sample_name}_summary.json"
    write_json(summary, detailed_json)

    if cfg.generate_html_report:
        html_path = generate_html_report(
            sample_name=sample_name,
            summary=summary,
            output_dir=sample_out_dir,
        )
        summary["outputs"]["html_report"] = str(html_path)
        logging.info("HTML report: %s", html_path)

    return summary


# ============================================================
# Checkpoint management
# ============================================================

def get_checkpoint_path(work_dir: Path, sample_base: str) -> Path:
    """Get the checkpoint file path for a sample."""
    return work_dir / "checkpoints" / f"{sample_base}.done"


def is_sample_complete(work_dir: Path, sample_base: str) -> bool:
    """Check if a sample has been successfully processed."""
    return get_checkpoint_path(work_dir, sample_base).exists()


def mark_sample_complete(work_dir: Path, sample_base: str) -> None:
    """Mark a sample as successfully processed."""
    cp = get_checkpoint_path(work_dir, sample_base)
    cp.parent.mkdir(parents=True, exist_ok=True)
    cp.write_text(datetime.now().isoformat())


# ============================================================
# Pipeline orchestration
# ============================================================

def find_fastq_files(directory: Path) -> list:
    """Find all FASTQ files in a directory."""
    return sorted([
        p for p in directory.iterdir()
        if p.name.endswith(".fastq") or p.name.endswith(".fastq.gz")
    ])


def get_sample_base(fastq_path: Path) -> str:
    """Extract sample base name from a FASTQ filename."""
    name = fastq_path.name
    for suffix in [".fastq.gz", ".fastq"]:
        if name.endswith(suffix):
            return name[:-len(suffix)]
    return name


def process_one_sample(
    raw_fastq: Path,
    cfg: PipelineConfig,
    combined_ref: Path,
    ref_lengths: dict,
    dirs: dict,
) -> dict:
    """Process a single sample through the full pipeline."""
    sample_base = get_sample_base(raw_fastq)
    work_dir = Path(cfg.work_dir)

    # Check checkpoint (resume support)
    if cfg.resume and is_sample_complete(work_dir, sample_base):
        logging.info("Skipping %s (already complete)", sample_base)
        summary_path = (
            dirs["samples"] / sample_base / "analysis" / f"{sample_base}_summary.json"
        )
        if summary_path.exists():
            with open(summary_path) as f:
                return json.load(f)
        return {"sample": sample_base, "status": "skipped_checkpoint"}

    sample_dir = dirs["samples"] / sample_base
    sample_dir.mkdir(parents=True, exist_ok=True)

    current_fastq = raw_fastq

    # --- Pre-filtering NanoPlot (raw reads baseline) ---
    if cfg.run_nanoplot:
        nanoplot_raw_dir = sample_dir / "nanoplot_raw"
        run_nanoplot(
            current_fastq,
            nanoplot_raw_dir,
            f"Raw Reads (Pre-filter): {sample_base}",
        )
        logging.info("Pre-filtering NanoPlot complete: %s", nanoplot_raw_dir)

    # --- Adapter trimming ---
    if cfg.run_porechop:
        trimmed_fastq = sample_dir / f"{sample_base}.trimmed.fastq"
        porechop_stats = run_porechop(current_fastq, trimmed_fastq)
        write_json(porechop_stats, sample_dir / f"{sample_base}_porechop_stats.json")
        current_fastq = trimmed_fastq

    # --- Quality/length filtering ---
    if cfg.run_nanofilt:
        filtered_fastq = sample_dir / f"{sample_base}.filtered.fastq"
        filt_stats = run_nanofilt(
            current_fastq, filtered_fastq,
            cfg.min_quality, cfg.min_length, cfg.max_length,
        )
        write_json(filt_stats, sample_dir / f"{sample_base}_nanofilt_stats.json")
        current_fastq = filtered_fastq

    # --- Post-filtering NanoPlot (filtered reads) ---
    if cfg.run_nanoplot:
        nanoplot_filtered_dir = sample_dir / "nanoplot_filtered"
        run_nanoplot(
            current_fastq,
            nanoplot_filtered_dir,
            f"Filtered Reads (Post-filter): {sample_base}",
        )
        logging.info("Post-filtering NanoPlot complete: %s", nanoplot_filtered_dir)

    # --- Alignment ---
    if cfg.run_mapping:
        bam_dir = sample_dir / "mapping"
        sorted_bam = align_fastq_to_reference(current_fastq, combined_ref, bam_dir, cfg)
    else:
        raise ValueError("Mapping must be enabled for downstream analyses.")

    # --- Full analysis ---
    summary = calculate_sample_statistics(
        sample_name=sample_base,
        fastq_path=current_fastq,
        bam_path=sorted_bam,
        ref_lengths=ref_lengths,
        cfg=cfg,
        sample_out_dir=sample_dir / "analysis",
    )

    mark_sample_complete(work_dir, sample_base)

    return summary


def write_run_summary(summaries: list, out_csv: Path) -> None:
    """Write aggregated run-level summary CSV."""
    out_csv.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "sample", "total_reads", "mapped_reads", "unmapped_reads",
        "mapping_rate", "chimeric_reads", "host_chimeric_reads",
        "breakpoints", "concatemer_candidates",
    ]

    with open(out_csv, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for s in summaries:
            writer.writerow({
                "sample": s.get("sample", ""),
                "total_reads": s.get("total_reads", ""),
                "mapped_reads": s.get("mapped_reads", ""),
                "unmapped_reads": s.get("unmapped_reads", ""),
                "mapping_rate": s.get("mapping_rate", ""),
                "chimeric_reads": s.get("chimeric_reads", ""),
                "host_chimeric_reads": s.get("host_chimeric_reads", ""),
                "breakpoints": s.get("breakpoints", ""),
                "concatemer_candidates": s.get("concatemer_candidates", ""),
            })


def build_dirs(work_dir: Path) -> dict:
    """Create the directory structure."""
    return {
        "root": work_dir,
        "logs": work_dir / "logs",
        "refs": work_dir / "refs",
        "samples": work_dir / "samples",
        "summaries": work_dir / "summaries",
        "checkpoints": work_dir / "checkpoints",
    }


# ============================================================
# Main entry point
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description="AAV-Chimera: Nanopore AAV Vector QC & Chimeric Read Detection",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent("""\
            Example usage:
              python aav_chimera.py \\
                --raw-fastq-dir /data/fastq \\
                --work-dir /data/output \\
                --refs-dir /data/references \\
                --transgene-name "pAAV-CMV-eGFP" \\
                --itr-5-start 0 --itr-5-end 145 \\
                --itr-3-start 4600 --itr-3-end 4745 \\
                --threads 8

            To see available reference names:
              python aav_chimera.py --list-references --refs-dir /data/references
        """),
    )

    # Required arguments
    parser.add_argument("--raw-fastq-dir", required="--list-references" not in sys.argv,
                        help="Directory containing input FASTQ files")
    parser.add_argument("--work-dir", required="--list-references" not in sys.argv,
                        help="Output working directory")
    parser.add_argument("--refs-dir", required=True,
                        help="Directory containing reference FASTA files")

    # Reference handling
    parser.add_argument("--reuse-combined-ref", action="store_true",
                        help="Reuse existing combined reference without prompting")
    parser.add_argument("--rebuild-combined-ref", action="store_true",
                        help="Force rebuild of combined reference even if it exists")
    parser.add_argument("--list-references", action="store_true",
                        help="List all sequence names in reference directory and exit")

    # Processing toggles
    parser.add_argument("--skip-porechop", action="store_true",
                        help="Skip adapter trimming")
    parser.add_argument("--skip-nanofilt", action="store_true",
                        help="Skip quality/length filtering")
    parser.add_argument("--skip-nanoplot", action="store_true",
                        help="Skip NanoPlot visualisation")
    parser.add_argument("--skip-html-report", action="store_true",
                        help="Skip HTML report generation")

    # Filtering parameters
    parser.add_argument("--min-quality", type=int, default=10,
                        help="Minimum read quality score (default: 10)")
    parser.add_argument("--min-length", type=int, default=500,
                        help="Minimum read length (default: 500)")
    parser.add_argument("--max-length", type=int, default=6000,
                        help="Maximum read length (default: 6000)")

    # Alignment parameters
    parser.add_argument("--threads", type=int, default=4,
                        help="Number of threads (default: 4)")
    parser.add_argument("--parallel-samples", type=int, default=1,
                        help="Number of samples to process in parallel (default: 1)")

    # Transgene identification
    parser.add_argument("--transgene-name", default="",
                        help=(
                            "Sequence ID of your vector plasmid/construct as it appears "
                            "in the FASTA header (e.g., 'pAAV-CMV-eGFP'). "
                            "Use --list-references to see available names."
                        ))

    # Coordinate annotations
    parser.add_argument("--itr-5-start", type=int, default=0,
                        help="Start coordinate of 5' ITR (0-based, default: 0)")
    parser.add_argument("--itr-5-end", type=int, default=145,
                        help="End coordinate of 5' ITR (default: 145)")
    parser.add_argument("--transgene-start", type=int, default=145,
                        help="Start of transgene cassette (default: 145)")
    parser.add_argument("--transgene-end", type=int, default=4600,
                        help="End of transgene cassette (default: 4600)")
    parser.add_argument("--itr-3-start", type=int, default=4600,
                        help="Start coordinate of 3' ITR (default: 4600)")
    parser.add_argument("--itr-3-end", type=int, default=4745,
                        help="End coordinate of 3' ITR (default: 4745)")

    # ITR parameters
    parser.add_argument("--itr-full-length-threshold", type=int, default=100,
                        help="Max missing ITR bp to be called full-length (default: 100)")
    parser.add_argument("--no-itr-masking", action="store_true",
                        help="Disable ITR variable region masking")
    parser.add_argument("--itr-variable-start", type=int, default=42,
                        help="Start of ITR variable region relative to ITR start (default: 42)")
    parser.add_argument("--itr-variable-end", type=int, default=85,
                        help="End of ITR variable region relative to ITR start (default: 85)")

    # Chimeric detection parameters
    parser.add_argument("--min-chimeric-segment-length", type=int, default=30,
                        help="Minimum aligned segment length for chimeric detection (default: 30)")
    parser.add_argument("--min-chimeric-mapq", type=int, default=10,
                        help="Minimum MAPQ for chimeric segments (default: 10)")

    # Workflow control
    parser.add_argument("--resume", action="store_true",
                        help="Resume from checkpoints (skip completed samples)")
    parser.add_argument("--dry-run", action="store_true",
                        help="Validate inputs and configuration without running")

    args = parser.parse_args()

    # --- Handle --list-references mode ---
    if args.list_references:
        refs_dir = Path(args.refs_dir)
        if not refs_dir.is_dir():
            print(f"ERROR: refs_dir does not exist: {refs_dir}")
            sys.exit(1)

        all_refs = list_available_references(refs_dir)

        print(f"\n{'='*70}")
        print(f"Available sequences in: {refs_dir}")
        print(f"{'='*70}")
        print(f"{'Sequence ID':<45} {'Length (bp)':<12} {'Source File'}")
        print(f"{'-'*45} {'-'*12} {'-'*20}")

        for seq_id, info in sorted(all_refs.items(), key=lambda x: x[1]["source_file"]):
            print(f"{seq_id:<45} {info['length']:<12,} {info['source_file']}")

        print(f"\n{'='*70}")
        print(
            f"Total: {len(all_refs)} sequences across "
            f"{len(set(v['source_file'] for v in all_refs.values()))} files"
        )
        if all_refs:
            non_host = [
                k for k, v in all_refs.items()
                if not re.match(r"^chr([1-9]|1[0-9]|2[0-2]|X|Y|M|MT)$", k)
            ]
            if non_host:
                print(f"\nLikely transgene candidates: {', '.join(non_host[:5])}")
                print(f"Example: --transgene-name \"{non_host[0]}\"")
        print(f"{'='*70}\n")
        sys.exit(0)

    # Build configuration
    cfg = PipelineConfig(
        raw_fastq_dir=args.raw_fastq_dir,
        work_dir=args.work_dir,
        refs_dir=args.refs_dir,
        reuse_combined_ref=args.reuse_combined_ref,
        rebuild_combined_ref=args.rebuild_combined_ref,
        run_porechop=not args.skip_porechop,
        run_nanofilt=not args.skip_nanofilt,
        run_nanoplot=not args.skip_nanoplot,
        min_quality=args.min_quality,
        min_length=args.min_length,
        max_length=args.max_length,
        threads=args.threads,
        parallel_samples=args.parallel_samples,
        transgene_name=args.transgene_name,
        itr_5_start=args.itr_5_start,
        itr_5_end=args.itr_5_end,
        transgene_start=args.transgene_start,
        transgene_end=args.transgene_end,
        itr_3_start=args.itr_3_start,
        itr_3_end=args.itr_3_end,
        itr_full_length_threshold=args.itr_full_length_threshold,
        mask_itr_variable_region=not args.no_itr_masking,
        itr_variable_start=args.itr_variable_start,
        itr_variable_end=args.itr_variable_end,
        min_chimeric_segment_length=args.min_chimeric_segment_length,
        min_chimeric_mapq=args.min_chimeric_mapq,
        resume=args.resume,
        dry_run=args.dry_run,
        generate_html_report=not args.skip_html_report,
    )

    # Create directory structure
    work_dir = Path(cfg.work_dir)
    dirs = build_dirs(work_dir)
    for d in dirs.values():
        d.mkdir(parents=True, exist_ok=True)

    # Set up logging
    log_file = dirs["logs"] / f"pipeline_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
    setup_logging(log_file)

    logging.info("=" * 60)
    logging.info("AAV-Chimera: Nanopore AAV Vector QC & Chimeric Read Detection")
    logging.info("=" * 60)

    # Log config (asdict won't include properties, so add them manually)
    config_dict = asdict(cfg)
    config_dict["itr_length_derived"] = cfg.itr_length
    config_dict["expected_genome_length"] = cfg.expected_genome_length
    logging.info("Config: %s", json.dumps(config_dict, indent=2))

    # --- Validate configuration ---
    config_errors = cfg.validate()
    if config_errors:
        for err in config_errors:
            logging.error("Config error: %s", err)
        sys.exit(1)

    # --- Validate dependencies ---
    missing_tools = validate_dependencies(cfg)
    if missing_tools:
        for tool in missing_tools:
            logging.error("Missing required tool: %s", tool)
        sys.exit(1)

    logging.info("All dependencies validated successfully")

    # --- Find input files ---
    raw_fastq_dir = Path(cfg.raw_fastq_dir)
    refs_dir = Path(cfg.refs_dir)

    fastqs = find_fastq_files(raw_fastq_dir)
    if not fastqs:
        logging.error("No FASTQ files found in %s", raw_fastq_dir)
        sys.exit(1)

    logging.info("Found %d FASTQ file(s) to process", len(fastqs))
    for fq in fastqs:
        logging.info("  - %s", fq.name)

    # --- Combine references (build once, reuse on subsequent runs) ---
    combined_ref = refs_dir / "combined_reference.fasta"
    combined_ref_index = Path(str(combined_ref) + ".fai")
    ref_lengths_cache = refs_dir / "combined_reference.lengths.json"

    needs_rebuild = False

    if not combined_ref.exists():
        logging.info("No combined reference found at %s — will build one.", combined_ref)
        needs_rebuild = True
    elif cfg.rebuild_combined_ref:
        logging.info("--rebuild-combined-ref flag set — forcing rebuild.")
        needs_rebuild = True
    elif cfg.reuse_combined_ref:
        logging.info("Reusing existing combined reference: %s", combined_ref)
        needs_rebuild = False
    else:
        # Interactive prompt
        file_size_mb = combined_ref.stat().st_size / (1024 * 1024)
        print(f"\n{'='*60}")
        print("Combined reference already exists:")
        print(f"  Path: {combined_ref}")
        print(f"  Size: {file_size_mb:.1f} MB")
        print(f"{'='*60}")
        response = input("Reuse existing combined reference? [Y/n]: ").strip().lower()
        if response in ("n", "no"):
            needs_rebuild = True
            logging.info("User chose to rebuild combined reference.")
        else:
            needs_rebuild = False
            logging.info("User chose to reuse existing combined reference.")

    if needs_rebuild:
        if combined_ref.exists():
            combined_ref.unlink()
        if combined_ref_index.exists():
            combined_ref_index.unlink()
        if ref_lengths_cache.exists():
            ref_lengths_cache.unlink()

        ref_lengths = combine_references(refs_dir, combined_ref, cfg)
        write_json(ref_lengths, ref_lengths_cache)
        logging.info(
            "Combined reference built: %d sequences, total %d bp",
            len(ref_lengths), sum(ref_lengths.values()),
        )
    else:
        if ref_lengths_cache.exists():
            with open(ref_lengths_cache) as fh:
                ref_lengths = json.load(fh)
            logging.info(
                "Loaded cached reference lengths: %d sequences, total %d bp",
                len(ref_lengths), sum(ref_lengths.values()),
            )
        else:
            logging.info("Parsing reference lengths from existing FASTA...")
            ref_lengths = fasta_lengths(combined_ref)
            write_json(ref_lengths, ref_lengths_cache)
            logging.info(
                "Parsed reference lengths: %d sequences, total %d bp",
                len(ref_lengths), sum(ref_lengths.values()),
            )

    # --- Resolve transgene name ---
    if cfg.transgene_name:
        # User explicitly specified
        if cfg.transgene_name not in ref_lengths:
            logging.error(
                "Transgene '%s' not found in combined reference!",
                cfg.transgene_name,
            )
            logging.error(
                "Available sequences: %s",
                ", ".join(sorted(ref_lengths.keys())),
            )
            logging.error("Use --list-references to see all available names.")
            sys.exit(1)
        logging.info(
            "Transgene: '%s' (%d bp)",
            cfg.transgene_name, ref_lengths[cfg.transgene_name],
        )
    else:
        # Auto-detect: find sequences that are NOT host and NOT helper
        candidate_transgenes = [
            name for name in ref_lengths.keys()
            if not is_host_reference(name, cfg.host_ref_patterns)
            and not is_helper_reference(name, cfg.helper_ref_patterns)
        ]

        if len(candidate_transgenes) == 1:
            cfg.transgene_name = candidate_transgenes[0]
            logging.info(
                "Auto-detected transgene: '%s' (%d bp)",
                cfg.transgene_name, ref_lengths[cfg.transgene_name],
            )
        elif len(candidate_transgenes) > 1:
            logging.warning(
                "Multiple potential transgene references detected: %s",
                ", ".join(f"'{n}'" for n in candidate_transgenes),
            )
            cfg.transgene_name = candidate_transgenes[0]
            logging.warning(
                "Using '%s' as default. Specify with --transgene-name to override.",
                cfg.transgene_name,
            )
        else:
            logging.error(
                "No transgene reference could be identified! "
                "All sequences appear to be host or helper. "
                "Use --transgene-name to specify explicitly."
            )
            sys.exit(1)

    for name, length in ref_lengths.items():
        logging.info("  %s: %d bp", name, length)

    # Log coordinate annotations
    logging.info("Plasmid annotations for '%s':", cfg.transgene_name)
    logging.info("  5' ITR: [%d, %d) (%d bp)", cfg.itr_5_start, cfg.itr_5_end, cfg.itr_length)
    logging.info("  Transgene: [%d, %d)", cfg.transgene_start, cfg.transgene_end)
    logging.info("  3' ITR: [%d, %d) (%d bp)", cfg.itr_3_start, cfg.itr_3_end, cfg.itr_length)
    logging.info("  Expected genome (ITR-to-ITR): %d bp", cfg.expected_genome_length)
    backbone = cfg.backbone_regions(ref_lengths.get(cfg.transgene_name, 0))
    if backbone:
        bb_total = sum(end - start for start, end in backbone)
        logging.info("  Backbone regions: %s (%d bp total)", backbone, bb_total)
    else:
        logging.info("  Backbone regions: none (plasmid = ITR-to-ITR only)")

    # --- Dry run: stop here ---
    if cfg.dry_run:
        logging.info("=" * 60)
        logging.info("DRY RUN complete. All inputs and dependencies validated.")
        logging.info("Would process %d sample(s).", len(fastqs))
        logging.info("=" * 60)
        return

    # --- Process samples ---
    summaries = []

    if cfg.parallel_samples > 1 and len(fastqs) > 1:
        logging.info(
            "Processing %d samples with %d parallel workers",
            len(fastqs), cfg.parallel_samples,
        )
        with ProcessPoolExecutor(max_workers=cfg.parallel_samples) as executor:
            futures = {
                executor.submit(
                    process_one_sample, fq, cfg, combined_ref, ref_lengths, dirs
                ): fq
                for fq in fastqs
            }
            for future in as_completed(futures):
                fq = futures[future]
                try:
                    summary = future.result()
                    summaries.append(summary)
                    logging.info("Completed: %s", fq.name)
                except Exception as e:
                    logging.exception("Failed on sample %s: %s", fq.name, str(e))
    else:
        for fastq in fastqs:
            logging.info("Processing sample: %s", fastq.name)
            try:
                summary = process_one_sample(fastq, cfg, combined_ref, ref_lengths, dirs)
                summaries.append(summary)
            except Exception as e:
                logging.exception("Failed on sample %s: %s", fastq.name, str(e))

    # --- Write run-level summaries ---
    if summaries:
        write_run_summary(summaries, dirs["summaries"] / "run_summary.csv")
        write_json(summaries, dirs["summaries"] / "run_summary.json")
        logging.info("Run summary written to %s", dirs["summaries"])
    else:
        logging.warning("No samples were successfully processed")

    logging.info("=" * 60)
    logging.info(
        "Workflow complete. Processed %d/%d samples successfully.",
        len(summaries), len(fastqs),
    )
    logging.info("=" * 60)


if __name__ == "__main__":
    main()
