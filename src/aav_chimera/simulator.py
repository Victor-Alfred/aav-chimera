#!/usr/bin/env python3
"""
AAV Chimeric Read Simulator & Benchmarking Tool
================================================
Generates simulated nanopore reads with known ground truth for
validating chimeric read detection and genome classification in
the AAV-Chimera pipeline.

Features:
- Full plasmid reference with coordinate annotations
- Backbone contamination simulation (read-through, pure backbone)
- Realistic nanopore error model (context-dependent)
- Chimeric read generation with microhomology/insertion junctions
- Ground truth output in same format as pipeline output
- Dedicated benchmarking module (precision/recall/F1)
- Configurable via command-line arguments
"""

import argparse
import csv
import json
import logging
import random
import re
import sys
import textwrap
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from statistics import mean
from typing import Optional

try:
    from Bio import SeqIO
except ImportError:
    sys.exit("ERROR: Biopython is required. Install via: pip install biopython")


# ============================================================
# Configuration
# ============================================================

@dataclass
class SimulationConfig:
    """Configuration for read simulation."""

    ref_dir: str
    output_dir: str

    # The full plasmid reference (including backbone) — must match FASTA header
    transgene_name: str = "pAAV-CMV-eGFP"

    # Coordinate annotations (0-based, same as pipeline)
    itr_5_start: int = 0
    itr_5_end: int = 145
    transgene_start: int = 145
    transgene_end: int = 4331
    itr_3_start: int = 4331
    itr_3_end: int = 4472

    # Reference files (basenames within ref_dir)
    ref_files: list = field(default_factory=lambda: [
        "pAAV-CMV-eGFP.fasta",
        "ref_rep_cap.fasta",
        "ref_helper.fasta",
        "hg38.fa",
    ])

    # Host genome patterns (for identifying host references)
    host_ref_patterns: tuple = (
        r"^chr([1-9]|1[0-9]|2[0-2]|X|Y|M|MT)$",
        r"^(?:[1-9]|1[0-9]|2[0-2]|X|Y|M|MT)$",
    )

    # Helper patterns
    helper_ref_patterns: tuple = (
        r".*rep.*cap.*",
        r".*helper.*",
        r".*pHelper.*",
    )

    # Simulation parameters
    num_fastq_files: int = 2
    reads_per_fastq_min: int = 5000
    reads_per_fastq_max: int = 10000
    min_read_length: int = 500
    max_read_length: int = 6000
    random_seed: int = 42

    # Read category proportions
    chimeric_proportion: float = 0.10
    backbone_proportion: float = 0.05
    host_dna_proportion: float = 0.05
    unmapped_proportion: float = 0.08

    # Error model parameters
    base_error_rate: float = 0.05
    homopolymer_error_rate: float = 0.15
    cluster_probability: float = 0.08
    max_cluster_length: int = 5

    # Chimeric read parameters
    min_segments: int = 2
    max_segments: int = 4
    min_segment_length: int = 300
    microhomology_probability: float = 0.4
    insertion_probability: float = 0.3
    max_microhomology_length: int = 15
    max_insertion_length: int = 20

    # Truncation parameters
    truncation_5prime_proportion: float = 0.10
    truncation_3prime_proportion: float = 0.10
    truncation_both_proportion: float = 0.05

    @property
    def itr_length(self) -> int:
        """ITR length derived from coordinates."""
        return self.itr_5_end - self.itr_5_start

    @property
    def expected_genome_length(self) -> int:
        """Expected full-length AAV genome size (ITR-to-ITR)."""
        return self.itr_3_end - self.itr_5_start

    def backbone_regions(self, ref_length: int) -> list:
        """Return backbone coordinate ranges."""
        regions = []
        if self.itr_5_start > 0:
            regions.append((0, self.itr_5_start))
        if self.itr_3_end < ref_length:
            regions.append((self.itr_3_end, ref_length))
        return regions


# ============================================================
# Logging
# ============================================================

def setup_logging(output_dir: Path) -> None:
    """Configure logging."""
    output_dir.mkdir(parents=True, exist_ok=True)
    log_file = output_dir / "simulation.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout),
        ],
    )


# ============================================================
# Sequence utilities
# ============================================================

def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    complement = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
    return "".join(complement.get(base, "N") for base in reversed(seq.upper()))


def gc_content(seq: str) -> float:
    """Calculate GC content as a percentage."""
    if not seq:
        return 0.0
    gc = seq.upper().count("G") + seq.upper().count("C")
    return (gc / len(seq)) * 100.0


def generate_quality_scores(length: int, mean_q: float = 15.0, std_q: float = 5.0) -> str:
    """
    Generate realistic nanopore quality scores.
    Nanopore Q-scores typically range from 5–25 with a mode around 12–15.
    """
    scores = []
    for _ in range(length):
        q = max(2, min(40, int(random.gauss(mean_q, std_q))))
        scores.append(chr(q + 33))
    return "".join(scores)


def is_host_reference(ref_name: str, patterns) -> bool:
    """Check if a reference name matches host genome patterns."""
    if ref_name is None:
        return False
    for pattern in patterns:
        if re.match(pattern, ref_name, flags=re.IGNORECASE):
            return True
    return False


def is_helper_reference(ref_name: str, patterns) -> bool:
    """Check if a reference name matches helper patterns."""
    if ref_name is None:
        return False
    for pattern in patterns:
        if re.match(pattern, ref_name, flags=re.IGNORECASE):
            return True
    return False


# ============================================================
# Error model
# ============================================================

def simulate_nanopore_errors(
    sequence: str,
    cfg: SimulationConfig,
) -> str:
    """
    Apply realistic nanopore sequencing errors to a sequence.

    Error types:
    - Substitutions (~2–3% of errors)
    - Deletions (~3–4% of errors)
    - Insertions (~1–2% of errors)
    - Homopolymer length errors (context-dependent)
    - Clustered errors (burst error model)
    """
    if not sequence:
        return sequence

    bases = "ACGT"
    error_read = []
    i = 0
    seq_len = len(sequence)
    in_cluster = False
    cluster_remaining = 0

    while i < seq_len:
        base = sequence[i].upper()
        if base not in "ACGT":
            error_read.append(base)
            i += 1
            continue

        # Determine if we're entering an error cluster
        if not in_cluster and random.random() < cfg.cluster_probability:
            in_cluster = True
            cluster_remaining = random.randint(2, cfg.max_cluster_length)

        # Determine effective error rate
        if in_cluster:
            effective_rate = min(0.5, cfg.base_error_rate * 3)
            cluster_remaining -= 1
            if cluster_remaining <= 0:
                in_cluster = False
        else:
            effective_rate = cfg.base_error_rate

        # Check for homopolymer context
        is_homopolymer = False
        homopolymer_length = 1
        if i + 1 < seq_len and sequence[i + 1].upper() == base:
            is_homopolymer = True
            j = i + 1
            while j < seq_len and sequence[j].upper() == base:
                j += 1
            homopolymer_length = j - i

        # Apply homopolymer-specific errors
        if is_homopolymer and homopolymer_length >= 3:
            if random.random() < cfg.homopolymer_error_rate:
                length_change = random.choice([-2, -1, 1, 1, 2])
                new_length = max(1, homopolymer_length + length_change)
                error_read.extend([base] * new_length)
                i += homopolymer_length
                continue

        # Apply standard errors
        if random.random() < effective_rate:
            error_type = random.random()
            if error_type < 0.35:
                # Substitution
                error_read.append(random.choice([b for b in bases if b != base]))
                i += 1
            elif error_type < 0.75:
                # Deletion (skip base)
                i += 1
            else:
                # Insertion (add random base then keep original)
                error_read.append(random.choice(bases))
                error_read.append(base)
                i += 1
        else:
            error_read.append(base)
            i += 1

    return "".join(error_read)


# ============================================================
# Data classes for ground truth
# ============================================================

@dataclass
class JunctionInfo:
    """Information about a chimeric junction."""
    junction_type: str
    microhomology_length: int = 0
    microhomology_seq: str = ""
    insertion_length: int = 0
    insertion_seq: str = ""
    left_ref: str = ""
    right_ref: str = ""
    left_ref_end: int = 0
    right_ref_start: int = 0
    left_strand: str = "+"
    right_strand: str = "+"


@dataclass
class SegmentInfo:
    """Information about a chimeric segment."""
    ref_name: str
    ref_seq_id: str
    ref_start: int
    ref_end: int
    strand: str
    aligned_length: int
    query_start: int = 0
    query_end: int = 0


@dataclass
class SimulatedRead:
    """Complete information about a simulated read."""
    read_name: str
    sequence: str
    quality: str
    read_length: int
    is_chimeric: bool
    is_host_involved: bool
    has_backbone: bool
    source_ref: str
    orientation: str

    # Chimeric-specific fields
    num_segments: int = 1
    segments: list = field(default_factory=list)
    junctions: list = field(default_factory=list)
    refs_order: list = field(default_factory=list)

    # Classification info
    read_category: str = "normal"  # normal, chimeric, backbone, host, unmapped
    truncation_category: str = "not_applicable"
    backbone_bp: int = 0

    # Quality metadata
    gc_content: float = 0.0
    mean_quality: float = 0.0


# ============================================================
# Read generation: Normal AAV reads
# ============================================================

def generate_normal_read(
    ref_sequences: dict,
    read_name: str,
    target_length: int,
    cfg: SimulationConfig,
) -> Optional[SimulatedRead]:
    """
    Generate a non-chimeric read from the packaged genome region
    (within [itr_5_start, itr_3_end] of the transgene plasmid).
    """
    # Get the transgene/plasmid sequence
    plasmid_seq = None
    for seq_id, seq in ref_sequences.get("transgene", []):
        if seq_id == cfg.transgene_name:
            plasmid_seq = seq
            break

    if not plasmid_seq:
        return None

    # Only sample from the packaged genome region
    genome_start = cfg.itr_5_start
    genome_end = cfg.itr_3_end
    genome_length = genome_end - genome_start

    if genome_length < cfg.min_segment_length:
        return None

    # Determine read length and position within the genome
    read_len = min(target_length, genome_length)
    max_start = max(0, genome_length - read_len)
    relative_start = random.randint(0, max_start)

    # Convert to absolute reference coordinates
    abs_start = genome_start + relative_start
    abs_end = min(abs_start + read_len, genome_end)

    fragment = plasmid_seq[abs_start:abs_end]

    # Determine orientation
    orientation = "forward" if random.random() < 0.5 else "reverse"
    strand = "+"
    if orientation == "reverse":
        fragment = reverse_complement(fragment)
        strand = "-"

    # Apply errors
    error_fragment = simulate_nanopore_errors(fragment, cfg)

    # Generate quality
    quality = generate_quality_scores(len(error_fragment))

    # Determine truncation category
    itr_5_missing = max(0, abs_start - cfg.itr_5_start)
    itr_3_missing = max(0, cfg.itr_3_end - abs_end)

    threshold = 100  # itr_full_length_threshold
    is_5_truncated = itr_5_missing > threshold
    is_3_truncated = itr_3_missing > threshold

    if is_5_truncated and is_3_truncated:
        truncation_category = "both_ends_truncated"
    elif is_5_truncated:
        truncation_category = "5_prime_truncated"
    elif is_3_truncated:
        truncation_category = "3_prime_truncated"
    else:
        truncation_category = "full_length"

    segment = SegmentInfo(
        ref_name=cfg.transgene_name,
        ref_seq_id=cfg.transgene_name,
        ref_start=abs_start,
        ref_end=abs_end,
        strand=strand,
        aligned_length=len(error_fragment),
        query_start=0,
        query_end=len(error_fragment),
    )

    return SimulatedRead(
        read_name=read_name,
        sequence=error_fragment,
        quality=quality,
        read_length=len(error_fragment),
        is_chimeric=False,
        is_host_involved=False,
        has_backbone=False,
        source_ref=cfg.transgene_name,
        orientation=orientation,
        num_segments=1,
        segments=[segment],
        junctions=[],
        refs_order=[cfg.transgene_name],
        read_category="normal",
        truncation_category=truncation_category,
        backbone_bp=0,
        gc_content=gc_content(error_fragment),
        mean_quality=sum(ord(c) - 33 for c in quality) / len(quality) if quality else 0,
    )


# ============================================================
# Read generation: Backbone contamination
# ============================================================

def generate_backbone_read(
    ref_sequences: dict,
    read_name: str,
    target_length: int,
    cfg: SimulationConfig,
) -> Optional[SimulatedRead]:
    """
    Generate a read that contains backbone sequence.

    Three modes:
    - pure_backbone: entirely within backbone region
    - readthrough_3prime: transgene extending past 3' ITR into backbone
    - readthrough_5prime: backbone extending into 5' ITR region
    """
    plasmid_seq = None
    for seq_id, seq in ref_sequences.get("transgene", []):
        if seq_id == cfg.transgene_name:
            plasmid_seq = seq
            break

    if not plasmid_seq:
        return None

    ref_len = len(plasmid_seq)
    backbone_regions = cfg.backbone_regions(ref_len)

    if not backbone_regions:
        return None

    # Calculate total backbone available
    total_backbone = sum(end - start for start, end in backbone_regions)
    if total_backbone < 100:
        return None

    # Choose mode based on what's available
    available_modes = []
    if cfg.itr_3_end < ref_len and (ref_len - cfg.itr_3_end) >= 50:
        available_modes.append("readthrough_3prime")
        available_modes.append("pure_backbone")
    if cfg.itr_5_start > 0 and cfg.itr_5_start >= 50:
        available_modes.append("readthrough_5prime")
        available_modes.append("pure_backbone")

    if not available_modes:
        return None

    mode = random.choice(available_modes)

    if mode == "pure_backbone":
        # Read entirely within backbone
        # Choose a backbone region
        bb_start, bb_end = random.choice(backbone_regions)
        bb_length = bb_end - bb_start
        read_len = min(target_length, bb_length)

        if read_len < 200:
            return None

        start = random.randint(bb_start, max(bb_start, bb_end - read_len))
        end = min(start + read_len, bb_end)
        fragment = plasmid_seq[start:end]
        backbone_bp = len(fragment)
        abs_start = start
        abs_end = end

    elif mode == "readthrough_3prime":
        # Starts in transgene/3'ITR, extends into backbone
        max_backbone_ext = min(500, ref_len - cfg.itr_3_end)
        backbone_extension = random.randint(50, max(51, max_backbone_ext))
        transgene_portion = target_length - backbone_extension

        abs_start = max(cfg.itr_5_start, cfg.itr_3_end - transgene_portion)
        abs_end = min(ref_len, cfg.itr_3_end + backbone_extension)
        fragment = plasmid_seq[abs_start:abs_end]
        backbone_bp = abs_end - cfg.itr_3_end

    elif mode == "readthrough_5prime":
        # Starts in backbone before 5' ITR, extends into transgene
        max_backbone_ext = min(500, cfg.itr_5_start)
        backbone_extension = random.randint(50, max(51, max_backbone_ext))

        abs_start = cfg.itr_5_start - backbone_extension
        abs_end = min(abs_start + target_length, cfg.itr_3_end)
        fragment = plasmid_seq[abs_start:abs_end]
        backbone_bp = cfg.itr_5_start - abs_start

    else:
        return None

    if len(fragment) < 200:
        return None

    # Apply errors and orientation
    error_fragment = simulate_nanopore_errors(fragment, cfg)
    orientation = random.choice(["forward", "reverse"])
    strand = "+"
    if orientation == "reverse":
        error_fragment = reverse_complement(error_fragment)
        strand = "-"

    quality = generate_quality_scores(len(error_fragment))

    segment = SegmentInfo(
        ref_name=cfg.transgene_name,
        ref_seq_id=cfg.transgene_name,
        ref_start=abs_start,
        ref_end=abs_end,
        strand=strand,
        aligned_length=len(error_fragment),
        query_start=0,
        query_end=len(error_fragment),
    )

    # Backbone reads carry the contamination mode as their truncation label.
    truncation_category = f"backbone_{mode}"

    return SimulatedRead(
        read_name=read_name,
        sequence=error_fragment,
        quality=quality,
        read_length=len(error_fragment),
        is_chimeric=False,
        is_host_involved=False,
        has_backbone=True,
        source_ref=cfg.transgene_name,
        orientation=orientation,
        num_segments=1,
        segments=[segment],
        junctions=[],
        refs_order=[cfg.transgene_name],
        read_category="backbone",
        truncation_category=truncation_category,
        backbone_bp=backbone_bp,
        gc_content=gc_content(error_fragment),
        mean_quality=sum(ord(c) - 33 for c in quality) / len(quality) if quality else 0,
    )


# ============================================================
# Read generation: Host DNA contamination
# ============================================================

def generate_host_read(
    ref_sequences: dict,
    read_name: str,
    target_length: int,
    cfg: SimulationConfig,
) -> Optional[SimulatedRead]:
    """Generate a pure host DNA contamination read."""
    host_seqs = ref_sequences.get("host", [])
    if not host_seqs:
        return None

    seq_id, seq = random.choice(host_seqs)

    read_len = min(target_length, len(seq) - 1)
    if read_len < cfg.min_segment_length:
        return None

    max_start = max(0, len(seq) - read_len)
    start_pos = random.randint(0, max_start)
    end_pos = start_pos + read_len

    fragment = seq[start_pos:end_pos]

    orientation = random.choice(["forward", "reverse"])
    strand = "+"
    if orientation == "reverse":
        fragment = reverse_complement(fragment)
        strand = "-"

    error_fragment = simulate_nanopore_errors(fragment, cfg)
    quality = generate_quality_scores(len(error_fragment))

    segment = SegmentInfo(
        ref_name="host",
        ref_seq_id=seq_id,
        ref_start=start_pos,
        ref_end=end_pos,
        strand=strand,
        aligned_length=len(error_fragment),
        query_start=0,
        query_end=len(error_fragment),
    )

    return SimulatedRead(
        read_name=read_name,
        sequence=error_fragment,
        quality=quality,
        read_length=len(error_fragment),
        is_chimeric=False,
        is_host_involved=True,
        has_backbone=False,
        source_ref=seq_id,
        orientation=orientation,
        num_segments=1,
        segments=[segment],
        junctions=[],
        refs_order=[seq_id],
        read_category="host",
        truncation_category="not_applicable",
        backbone_bp=0,
        gc_content=gc_content(error_fragment),
        mean_quality=sum(ord(c) - 33 for c in quality) / len(quality) if quality else 0,
    )


# ============================================================
# Read generation: Chimeric reads
# ============================================================

def construct_junction(
    left_seq: str,
    right_seq: str,
    cfg: SimulationConfig,
) -> tuple:
    """
    Construct a realistic junction between two segments.
    Returns (junction_sequence, JunctionInfo).
    """
    junction_seq = ""
    r = random.random()

    if r < cfg.microhomology_probability:
        mh_len = random.randint(
            2, min(cfg.max_microhomology_length, len(left_seq), len(right_seq))
        )
        mh_seq = left_seq[-mh_len:]
        junction_info = JunctionInfo(
            junction_type="microhomology",
            microhomology_length=mh_len,
            microhomology_seq=mh_seq,
        )
        junction_seq = ""

    elif r < cfg.microhomology_probability + cfg.insertion_probability:
        ins_len = random.randint(1, cfg.max_insertion_length)
        ins_seq = "".join(random.choice("ACGT") for _ in range(ins_len))
        junction_seq = ins_seq
        junction_info = JunctionInfo(
            junction_type="insertion",
            insertion_length=ins_len,
            insertion_seq=ins_seq,
        )

    else:
        junction_seq = ""
        junction_info = JunctionInfo(junction_type="blunt")

    return junction_seq, junction_info


def generate_chimeric_read(
    ref_sequences: dict,
    read_name: str,
    target_length: int,
    cfg: SimulationConfig,
    force_host: bool = False,
) -> Optional[SimulatedRead]:
    """
    Generate a chimeric read composed of segments from different references.
    """
    num_segments = random.randint(cfg.min_segments, cfg.max_segments)

    # Build pool of available sequences
    all_seq_entries = []  # (category, seq_id, sequence)

    # Transgene plasmid — sample from the genome region only
    for seq_id, seq in ref_sequences.get("transgene", []):
        if seq_id == cfg.transgene_name:
            genome_seq = seq[cfg.itr_5_start:cfg.itr_3_end]
            all_seq_entries.append(("transgene", seq_id, genome_seq))

    # Host sequences
    for seq_id, seq in ref_sequences.get("host", []):
        all_seq_entries.append(("host", seq_id, seq))

    # Helper sequences
    for seq_id, seq in ref_sequences.get("helper", []):
        all_seq_entries.append(("helper", seq_id, seq))

    if len(all_seq_entries) < 2:
        return None

    # Choose segment sources (ensure different categories for chimeric)
    segment_sources = []
    for seg_idx in range(num_segments):
        if force_host and seg_idx == 0:
            host_entries = [e for e in all_seq_entries if e[0] == "host"]
            if host_entries:
                chosen = random.choice(host_entries)
            else:
                chosen = random.choice(all_seq_entries)
        elif seg_idx == 0:
            chosen = random.choice(all_seq_entries)
        else:
            # Try to pick from a different category
            prev_category = segment_sources[-1][0]
            different = [e for e in all_seq_entries if e[0] != prev_category]
            if different:
                chosen = random.choice(different)
            else:
                chosen = random.choice(all_seq_entries)
        segment_sources.append(chosen)

    # Calculate segment lengths
    base_size = max(cfg.min_segment_length, target_length // num_segments)
    segment_lengths = []
    for _ in range(num_segments):
        variation = random.randint(-base_size // 4, base_size // 4)
        seg_len = max(cfg.min_segment_length, base_size + variation)
        segment_lengths.append(seg_len)

    # Build the chimeric sequence
    chimeric_parts = []
    segments = []
    junctions = []
    query_pos = 0

    for seg_idx in range(num_segments):
        category, seq_id, ref_seq = segment_sources[seg_idx]

        seg_len = min(segment_lengths[seg_idx], len(ref_seq) - 1)
        if seg_len < cfg.min_segment_length:
            seg_len = min(cfg.min_segment_length, len(ref_seq) - 1)

        max_start = max(0, len(ref_seq) - seg_len)
        start_pos = random.randint(0, max_start)
        end_pos = start_pos + seg_len

        fragment = ref_seq[start_pos:end_pos]

        strand = random.choice(["+", "-"])
        if strand == "-":
            fragment = reverse_complement(fragment)

        # Handle junction with previous segment
        if seg_idx > 0 and chimeric_parts:
            prev_fragment = chimeric_parts[-1]
            junction_seq, junction_info = construct_junction(
                prev_fragment, fragment, cfg
            )

            junction_info.left_ref = segment_sources[seg_idx - 1][1]
            junction_info.right_ref = seq_id
            junction_info.left_ref_end = segments[-1].ref_end
            junction_info.right_ref_start = start_pos
            junction_info.left_strand = segments[-1].strand
            junction_info.right_strand = strand

            # Handle microhomology (trim right fragment start)
            if junction_info.junction_type == "microhomology":
                mh_len = junction_info.microhomology_length
                fragment = fragment[mh_len:]
                if strand == "+":
                    start_pos += mh_len
                seg_len -= mh_len

            # Add junction sequence if any
            if junction_seq:
                chimeric_parts.append(junction_seq)
                query_pos += len(junction_seq)

            junctions.append(junction_info)

        # Record segment info
        # For transgene segments, convert back to absolute plasmid coordinates
        if category == "transgene":
            abs_start = cfg.itr_5_start + start_pos
            abs_end = cfg.itr_5_start + end_pos
        else:
            abs_start = start_pos
            abs_end = end_pos

        seg_info = SegmentInfo(
            ref_name=category,
            ref_seq_id=seq_id,
            ref_start=abs_start,
            ref_end=abs_end,
            strand=strand,
            aligned_length=len(fragment),
            query_start=query_pos,
            query_end=query_pos + len(fragment),
        )
        segments.append(seg_info)

        chimeric_parts.append(fragment)
        query_pos += len(fragment)

    # Combine all parts
    full_sequence = "".join(chimeric_parts)

    # Apply nanopore errors
    error_sequence = simulate_nanopore_errors(full_sequence, cfg)

    # Trim to target length if needed
    if len(error_sequence) > target_length:
        error_sequence = error_sequence[:target_length]

    # Apply read-level orientation
    orientation = "forward" if random.random() < 0.5 else "reverse"
    final_sequence = error_sequence
    if orientation == "reverse":
        final_sequence = reverse_complement(error_sequence)
        segments = segments[::-1]
        junctions = junctions[::-1]
        for seg in segments:
            seg.strand = "-" if seg.strand == "+" else "+"

    quality = generate_quality_scores(len(final_sequence))

    # Check if host is involved
    is_host = any(seg.ref_name == "host" for seg in segments)
    refs_order = [seg.ref_seq_id for seg in segments]

    return SimulatedRead(
        read_name=read_name,
        sequence=final_sequence,
        quality=quality,
        read_length=len(final_sequence),
        is_chimeric=True,
        is_host_involved=is_host,
        has_backbone=False,
        source_ref=refs_order[0] if refs_order else "",
        orientation=orientation,
        num_segments=len(segments),
        segments=segments,
        junctions=junctions,
        refs_order=refs_order,
        read_category="chimeric",
        truncation_category="not_applicable",
        backbone_bp=0,
        gc_content=gc_content(final_sequence),
        mean_quality=sum(ord(c) - 33 for c in quality) / len(quality) if quality else 0,
    )


# ============================================================
# Main simulation
# ============================================================

def load_references(cfg: SimulationConfig) -> dict:
    """
    Load all reference sequences, categorised by type.
    Returns: {"transgene": [(id, seq)], "host": [(id, seq)], "helper": [(id, seq)]}
    """
    ref_sequences = {"transgene": [], "host": [], "helper": []}
    ref_dir = Path(cfg.ref_dir)

    # Only keep conventional chromosomes from hg38
    conventional_chr = re.compile(
        r"^chr([1-9]|1[0-9]|2[0-2]|X|Y|M|MT)$", re.IGNORECASE
    )

    for ref_file in cfg.ref_files:
        ref_path = ref_dir / ref_file
        if not ref_path.exists():
            logging.warning("Reference file not found: %s", ref_path)
            continue

        # Skip the combined reference if it's in the folder
        if ref_file == "combined_reference.fasta":
            logging.info("Skipping combined reference: %s", ref_file)
            continue

        is_human = (
            "hg38" in ref_file.lower()
            or "grch38" in ref_file.lower()
        )

        logging.info("Loading reference: %s", ref_file)

        for record in SeqIO.parse(str(ref_path), "fasta"):
            seq_id = record.id
            seq = str(record.seq).upper()

            # Categorise
            if seq_id == cfg.transgene_name:
                ref_sequences["transgene"].append((seq_id, seq))
                logging.info("  Transgene: %s (%d bp)", seq_id, len(seq))

            elif is_human:
                # Only keep conventional chromosomes
                if conventional_chr.match(seq_id):
                    ref_sequences["host"].append((seq_id, seq))
                    logging.info("  Host: %s (%d bp)", seq_id, len(seq))
                else:
                    logging.debug("  Skipping non-conventional: %s", seq_id)

            elif is_host_reference(seq_id, cfg.host_ref_patterns):
                ref_sequences["host"].append((seq_id, seq))
                logging.info("  Host: %s (%d bp)", seq_id, len(seq))

            elif is_helper_reference(seq_id, cfg.helper_ref_patterns):
                ref_sequences["helper"].append((seq_id, seq))
                logging.info("  Helper: %s (%d bp)", seq_id, len(seq))

            else:
                logging.debug("  Skipping unknown: %s (%d bp)", seq_id, len(seq))

    logging.info(
        "Loaded: %d transgene, %d host, %d helper sequences",
        len(ref_sequences["transgene"]),
        len(ref_sequences["host"]),
        len(ref_sequences["helper"]),
    )

    return ref_sequences


def run_simulation(cfg: SimulationConfig) -> dict:
    """
    Run the full simulation, generating FASTQ files and ground truth.
    Returns a summary dict with statistics.
    """
    random.seed(cfg.random_seed)
    output_dir = Path(cfg.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    setup_logging(output_dir)
    logging.info("Starting simulation with config:")
    logging.info("  Transgene: %s", cfg.transgene_name)
    logging.info("  5' ITR: [%d, %d)", cfg.itr_5_start, cfg.itr_5_end)
    logging.info("  Transgene region: [%d, %d)", cfg.transgene_start, cfg.transgene_end)
    logging.info("  3' ITR: [%d, %d)", cfg.itr_3_start, cfg.itr_3_end)
    logging.info("  Expected genome: %d bp", cfg.expected_genome_length)

    # Load references
    ref_sequences = load_references(cfg)
    if not ref_sequences["transgene"]:
        logging.error("No transgene reference loaded! Check transgene_name matches FASTA header.")
        return {}

    all_summaries = []

    for fastq_idx in range(cfg.num_fastq_files):
        fastq_name = f"simulated_reads_{fastq_idx + 1}"
        fastq_path = output_dir / f"{fastq_name}.fastq"
        sample_dir = output_dir / fastq_name
        sample_dir.mkdir(parents=True, exist_ok=True)

        num_reads = random.randint(cfg.reads_per_fastq_min, cfg.reads_per_fastq_max)
        logging.info("Generating %s: %d reads", fastq_name, num_reads)

        # Determine read category counts
        num_chimeric = int(num_reads * cfg.chimeric_proportion)
        num_backbone = int(num_reads * cfg.backbone_proportion)
        num_host = int(num_reads * cfg.host_dna_proportion)
        num_host_chimeric = int(num_chimeric * 0.3)
        num_normal = num_reads - num_chimeric - num_backbone - num_host

        # Generate all reads
        all_reads = []
        ground_truth_chimeric = []
        ground_truth_backbone = []
        read_counter = 0

        # Generate chimeric reads
        logging.info("  Generating %d chimeric reads...", num_chimeric)
        for i in range(num_chimeric):
            read_name = f"read_{read_counter:06d}"
            read_counter += 1
            target_length = random.randint(cfg.min_read_length, cfg.max_read_length)
            force_host = i < num_host_chimeric

            sim_read = generate_chimeric_read(
                ref_sequences=ref_sequences,
                read_name=read_name,
                target_length=target_length,
                cfg=cfg,
                force_host=force_host,
            )
            if sim_read is not None:
                all_reads.append(sim_read)
                ground_truth_chimeric.append(sim_read)

        # Generate backbone reads
        logging.info("  Generating %d backbone reads...", num_backbone)
        for _ in range(num_backbone):
            read_name = f"read_{read_counter:06d}"
            read_counter += 1
            target_length = random.randint(cfg.min_read_length, cfg.max_read_length)

            sim_read = generate_backbone_read(
                ref_sequences=ref_sequences,
                read_name=read_name,
                target_length=target_length,
                cfg=cfg,
            )
            if sim_read is not None:
                all_reads.append(sim_read)
                ground_truth_backbone.append(sim_read)

        # Generate host reads
        logging.info("  Generating %d host reads...", num_host)
        for _ in range(num_host):
            read_name = f"read_{read_counter:06d}"
            read_counter += 1
            target_length = random.randint(cfg.min_read_length, cfg.max_read_length)

            sim_read = generate_host_read(
                ref_sequences=ref_sequences,
                read_name=read_name,
                target_length=target_length,
                cfg=cfg,
            )
            if sim_read is not None:
                all_reads.append(sim_read)

        # Generate normal reads
        logging.info("  Generating %d normal reads...", num_normal)
        for _ in range(num_normal):
            read_name = f"read_{read_counter:06d}"
            read_counter += 1
            target_length = random.randint(cfg.min_read_length, cfg.max_read_length)

            sim_read = generate_normal_read(
                ref_sequences=ref_sequences,
                read_name=read_name,
                target_length=target_length,
                cfg=cfg,
            )
            if sim_read is not None:
                all_reads.append(sim_read)

        # Shuffle reads
        random.shuffle(all_reads)

        # Write FASTQ
        logging.info("Writing FASTQ: %s", fastq_path)
        with open(fastq_path, "w") as fh:
            for read in all_reads:
                fh.write(f"@{read.read_name}\n")
                fh.write(f"{read.sequence}\n")
                fh.write("+\n")
                fh.write(f"{read.quality}\n")

        # Write ground truth files
        write_ground_truth(
            all_reads=all_reads,
            chimeric_reads=ground_truth_chimeric,
            backbone_reads=ground_truth_backbone,
            output_dir=sample_dir,
            sample_name=fastq_name,
            cfg=cfg,
        )

        # Compute statistics
        summary = compute_simulation_stats(
            all_reads=all_reads,
            chimeric_reads=ground_truth_chimeric,
            backbone_reads=ground_truth_backbone,
            sample_name=fastq_name,
        )
        all_summaries.append(summary)

        logging.info(
            "  Generated: %d total, %d chimeric, %d backbone, %d host-chimeric",
            len(all_reads),
            len(ground_truth_chimeric),
            len(ground_truth_backbone),
            sum(1 for r in ground_truth_chimeric if r.is_host_involved),
        )

    # Write overall summary
    summary_path = output_dir / "simulation_summary.json"
    with open(summary_path, "w") as fh:
        json.dump(all_summaries, fh, indent=2, default=str)

    logging.info("Simulation complete. Output: %s", output_dir)
    return {"summaries": all_summaries, "output_dir": str(output_dir)}


# ============================================================
# Ground truth output (matches pipeline output format)
# ============================================================

def write_ground_truth(
    all_reads: list,
    chimeric_reads: list,
    backbone_reads: list,
    output_dir: Path,
    sample_name: str,
    cfg: SimulationConfig,
) -> None:
    """Write ground truth files in the same format as the pipeline output."""
    output_dir.mkdir(parents=True, exist_ok=True)

    # 1. Chimeric reads CSV
    chimeric_csv = output_dir / f"{sample_name}_ground_truth_chimeric_reads.csv"
    fieldnames_chimeric = [
        "read_name", "read_length", "num_segments", "refs_order",
        "segment_lengths", "ref_starts", "ref_ends", "query_starts",
        "query_ends", "strands", "mapqs", "host_involved",
    ]
    with open(chimeric_csv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames_chimeric)
        writer.writeheader()
        for read in chimeric_reads:
            writer.writerow({
                "read_name": read.read_name,
                "read_length": read.read_length,
                "num_segments": read.num_segments,
                "refs_order": "|".join(read.refs_order),
                "segment_lengths": "|".join(
                    str(s.aligned_length) for s in read.segments
                ),
                "ref_starts": "|".join(
                    str(s.ref_start) for s in read.segments
                ),
                "ref_ends": "|".join(
                    str(s.ref_end) for s in read.segments
                ),
                "query_starts": "|".join(
                    str(s.query_start) for s in read.segments
                ),
                "query_ends": "|".join(
                    str(s.query_end) for s in read.segments
                ),
                "strands": "|".join(s.strand for s in read.segments),
                "mapqs": "|".join(["60"] * read.num_segments),
                "host_involved": int(read.is_host_involved),
            })

    # 2. Chimeric breakpoints CSV
    breakpoints_csv = output_dir / f"{sample_name}_ground_truth_breakpoints.csv"
    bp_fieldnames = [
        "read_name", "junction_index", "junction_type",
        "left_ref", "left_ref_end", "left_strand",
        "right_ref", "right_ref_start", "right_strand",
        "microhomology_length", "microhomology_seq",
        "insertion_length", "insertion_seq",
    ]
    with open(breakpoints_csv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=bp_fieldnames)
        writer.writeheader()
        for read in chimeric_reads:
            for j_idx, junction in enumerate(read.junctions):
                writer.writerow({
                    "read_name": read.read_name,
                    "junction_index": j_idx + 1,
                    "junction_type": junction.junction_type,
                    "left_ref": junction.left_ref,
                    "left_ref_end": junction.left_ref_end,
                    "left_strand": junction.left_strand,
                    "right_ref": junction.right_ref,
                    "right_ref_start": junction.right_ref_start,
                    "right_strand": junction.right_strand,
                    "microhomology_length": junction.microhomology_length,
                    "microhomology_seq": junction.microhomology_seq,
                    "insertion_length": junction.insertion_length,
                    "insertion_seq": junction.insertion_seq,
                })

    # 3. All reads classification (ground truth labels)
    labels_csv = output_dir / f"{sample_name}_ground_truth_labels.csv"
    label_fieldnames = [
        "read_name", "is_chimeric", "is_host_involved", "has_backbone",
        "num_segments", "refs_order", "source_ref", "read_category",
        "truncation_category", "backbone_bp", "read_length", "gc_content",
    ]
    with open(labels_csv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=label_fieldnames)
        writer.writeheader()
        for read in all_reads:
            writer.writerow({
                "read_name": read.read_name,
                "is_chimeric": int(read.is_chimeric),
                "is_host_involved": int(read.is_host_involved),
                "has_backbone": int(read.has_backbone),
                "num_segments": read.num_segments,
                "refs_order": "|".join(read.refs_order),
                "source_ref": read.source_ref,
                "read_category": read.read_category,
                "truncation_category": read.truncation_category,
                "backbone_bp": read.backbone_bp,
                "read_length": read.read_length,
                "gc_content": round(read.gc_content, 2),
            })

    # 4. Backbone reads CSV
    backbone_csv = output_dir / f"{sample_name}_ground_truth_backbone_reads.csv"
    bb_fieldnames = [
        "read_name", "read_length", "backbone_bp", "truncation_category",
        "ref_start", "ref_end",
    ]
    with open(backbone_csv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=bb_fieldnames)
        writer.writeheader()
        for read in backbone_reads:
            seg = read.segments[0] if read.segments else None
            writer.writerow({
                "read_name": read.read_name,
                "read_length": read.read_length,
                "backbone_bp": read.backbone_bp,
                "truncation_category": read.truncation_category,
                "ref_start": seg.ref_start if seg else 0,
                "ref_end": seg.ref_end if seg else 0,
            })

    # 5. Summary statistics JSON
    summary = {
        "sample_name": sample_name,
        "transgene_name": cfg.transgene_name,
        "coordinates": {
            "itr_5": [cfg.itr_5_start, cfg.itr_5_end],
            "transgene": [cfg.transgene_start, cfg.transgene_end],
            "itr_3": [cfg.itr_3_start, cfg.itr_3_end],
        },
        "total_reads": len(all_reads),
        "chimeric_reads": len(chimeric_reads),
        "backbone_reads": len(backbone_reads),
        "host_chimeric_reads": sum(1 for r in chimeric_reads if r.is_host_involved),
        "non_chimeric_reads": len(all_reads) - len(chimeric_reads),
        "junction_types": dict(Counter(
            j.junction_type for r in chimeric_reads for j in r.junctions
        )),
        "segment_count_distribution": dict(Counter(
            r.num_segments for r in chimeric_reads
        )),
        "category_distribution": dict(Counter(
            r.read_category for r in all_reads
        )),
        "mean_read_length": mean(r.read_length for r in all_reads) if all_reads else 0,
        "mean_gc_content": mean(r.gc_content for r in all_reads) if all_reads else 0,
    }

    summary_json = output_dir / f"{sample_name}_ground_truth_summary.json"
    with open(summary_json, "w") as fh:
        json.dump(summary, fh, indent=2, default=str)

    logging.info("Ground truth written to: %s", output_dir)


def compute_simulation_stats(
    all_reads: list,
    chimeric_reads: list,
    backbone_reads: list,
    sample_name: str,
) -> dict:
    """Compute summary statistics for a simulated dataset."""
    category_counts = Counter(r.read_category for r in all_reads)
    truncation_counts = Counter(
        r.truncation_category for r in all_reads
        if r.truncation_category != "not_applicable"
    )

    return {
        "sample_name": sample_name,
        "total_reads": len(all_reads),
        "chimeric_reads": len(chimeric_reads),
        "chimeric_proportion": len(chimeric_reads) / len(all_reads) if all_reads else 0,
        "backbone_reads": len(backbone_reads),
        "backbone_proportion": len(backbone_reads) / len(all_reads) if all_reads else 0,
        "host_chimeric_reads": sum(1 for r in chimeric_reads if r.is_host_involved),
        "category_distribution": dict(category_counts),
        "truncation_distribution": dict(truncation_counts),
        "mean_read_length": round(
            mean(r.read_length for r in all_reads), 1
        ) if all_reads else 0,
        "mean_segments_per_chimeric": round(
            mean(r.num_segments for r in chimeric_reads), 2
        ) if chimeric_reads else 0,
    }


# ============================================================
# Benchmarking module
# ============================================================

@dataclass
class BenchmarkResults:
    """Results from comparing pipeline output to ground truth."""

    # Chimeric detection metrics
    chimeric_tp: int = 0
    chimeric_fp: int = 0
    chimeric_fn: int = 0
    chimeric_tn: int = 0
    chimeric_sensitivity: float = 0.0
    chimeric_specificity: float = 0.0
    chimeric_precision: float = 0.0
    chimeric_recall: float = 0.0
    chimeric_f1: float = 0.0

    # Host chimeric metrics
    host_tp: int = 0
    host_fp: int = 0
    host_fn: int = 0

    # Backbone detection metrics
    backbone_tp: int = 0
    backbone_fp: int = 0
    backbone_fn: int = 0
    backbone_sensitivity: float = 0.0
    backbone_precision: float = 0.0

    # Reference pair accuracy
    ref_pair_correct: int = 0
    ref_pair_total: int = 0
    ref_pair_accuracy: float = 0.0

    # Overall accuracy
    accuracy: float = 0.0

    # Error details
    fp_read_names: list = field(default_factory=list)
    fn_read_names: list = field(default_factory=list)


def load_ground_truth_labels(ground_truth_dir: Path) -> dict:
    """Load ground truth labels CSV into a dict keyed by read_name."""
    labels = {}
    labels_files = list(ground_truth_dir.glob("*_ground_truth_labels.csv"))

    for labels_file in labels_files:
        with open(labels_file) as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                labels[row["read_name"]] = {
                    "is_chimeric": bool(int(row["is_chimeric"])),
                    "is_host_involved": bool(int(row["is_host_involved"])),
                    "has_backbone": bool(int(row["has_backbone"])),
                    "num_segments": int(row["num_segments"]),
                    "refs_order": row["refs_order"].split("|") if row["refs_order"] else [],
                    "read_category": row["read_category"],
                    "backbone_bp": int(row["backbone_bp"]),
                }

    return labels


def load_pipeline_chimeric_output(pipeline_output_dir: Path) -> dict:
    """Load pipeline chimeric reads CSV into a dict keyed by read_name."""
    detected = {}
    # Search recursively — chimeric CSV is in a subdirectory
    chimeric_files = list(pipeline_output_dir.rglob("*_chimeric_reads.csv"))

    for chimeric_file in chimeric_files:
        logging.info("Loading pipeline chimeric output: %s", chimeric_file)
        with open(chimeric_file) as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                detected[row["read_name"]] = {
                    "num_segments": int(row["num_segments"]),
                    "refs_order": row["refs_order"].split("|") if row["refs_order"] else [],
                    "host_involved": bool(int(row["host_involved"])),
                }

    return detected


def load_pipeline_classifications(pipeline_output_dir: Path) -> dict:
    """Load pipeline genome classification CSV."""
    classifications = {}
    # Also search recursively for safety
    class_files = list(pipeline_output_dir.rglob("*_genome_classifications.csv"))

    for class_file in class_files:
        logging.info("Loading pipeline classifications: %s", class_file)
        with open(class_file) as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                classifications[row["read_name"]] = {
                    "category": row["category"],
                    "subcategory": row["subcategory"],
                }

    return classifications


def benchmark_pipeline(
    ground_truth_dir: Path,
    pipeline_output_dir: Path,
    output_dir: Path,
) -> BenchmarkResults:
    """
    Compare pipeline detection results against ground truth.

    Evaluates:
    - Chimeric read detection (TP/FP/FN/TN)
    - Host chimeric detection
    - Backbone contamination detection
    - Reference pair accuracy
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load data
    ground_truth = load_ground_truth_labels(ground_truth_dir)
    pipeline_chimeric = load_pipeline_chimeric_output(pipeline_output_dir)
    pipeline_classifications = load_pipeline_classifications(pipeline_output_dir)

    if not ground_truth:
        logging.error("No ground truth labels found in %s", ground_truth_dir)
        return BenchmarkResults()

    logging.info(
        "Benchmarking: %d ground truth reads, %d pipeline chimeric detections, "
        "%d pipeline classifications",
        len(ground_truth), len(pipeline_chimeric), len(pipeline_classifications),
    )

    results = BenchmarkResults()

    all_read_names = set(ground_truth.keys())
    detected_chimeric_names = set(pipeline_chimeric.keys())

    # --- Chimeric detection metrics ---
    for read_name in all_read_names:
        gt = ground_truth[read_name]
        is_detected_chimeric = read_name in detected_chimeric_names

        if gt["is_chimeric"] and is_detected_chimeric:
            results.chimeric_tp += 1
            # Check reference pair accuracy
            gt_refs = set(gt["refs_order"])
            det_refs = set(pipeline_chimeric[read_name]["refs_order"])
            if gt_refs == det_refs:
                results.ref_pair_correct += 1
            results.ref_pair_total += 1

        elif gt["is_chimeric"] and not is_detected_chimeric:
            results.chimeric_fn += 1
            results.fn_read_names.append(read_name)

        elif not gt["is_chimeric"] and is_detected_chimeric:
            results.chimeric_fp += 1
            results.fp_read_names.append(read_name)

        else:
            results.chimeric_tn += 1

        # Host-specific metrics
        if gt["is_host_involved"] and gt["is_chimeric"]:
            if is_detected_chimeric and pipeline_chimeric.get(read_name, {}).get("host_involved", False):
                results.host_tp += 1
            elif not is_detected_chimeric:
                results.host_fn += 1

        if is_detected_chimeric and pipeline_chimeric.get(read_name, {}).get("host_involved", False):
            if not gt.get("is_host_involved", False):
                results.host_fp += 1

    # --- Backbone detection metrics ---
    for read_name in all_read_names:
        gt = ground_truth[read_name]
        classification = pipeline_classifications.get(read_name, {})
        detected_as_backbone = classification.get("category", "") == "backbone_contamination"

        if gt["has_backbone"] and detected_as_backbone:
            results.backbone_tp += 1
        elif gt["has_backbone"] and not detected_as_backbone:
            results.backbone_fn += 1
        elif not gt["has_backbone"] and detected_as_backbone:
            results.backbone_fp += 1

    # Extra detections not in ground truth
    extra_detections = detected_chimeric_names - all_read_names
    results.chimeric_fp += len(extra_detections)
    results.fp_read_names.extend(list(extra_detections))

    # --- Compute derived metrics ---
    tp = results.chimeric_tp
    fp = results.chimeric_fp
    fn = results.chimeric_fn
    tn = results.chimeric_tn

    results.chimeric_sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    results.chimeric_recall = results.chimeric_sensitivity
    results.chimeric_specificity = tn / (tn + fp) if (tn + fp) > 0 else 0.0
    results.chimeric_precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    results.accuracy = (tp + tn) / (tp + tn + fp + fn) if (tp + tn + fp + fn) > 0 else 0.0

    if results.chimeric_precision + results.chimeric_recall > 0:
        results.chimeric_f1 = (
            2 * results.chimeric_precision * results.chimeric_recall
            / (results.chimeric_precision + results.chimeric_recall)
        )

    results.ref_pair_accuracy = (
        results.ref_pair_correct / results.ref_pair_total
        if results.ref_pair_total > 0 else 0.0
    )

    # Backbone metrics
    bb_tp = results.backbone_tp
    bb_fp = results.backbone_fp
    bb_fn = results.backbone_fn
    results.backbone_sensitivity = bb_tp / (bb_tp + bb_fn) if (bb_tp + bb_fn) > 0 else 0.0
    results.backbone_precision = bb_tp / (bb_tp + bb_fp) if (bb_tp + bb_fp) > 0 else 0.0

    # --- Write results ---
    results_dict = {
        "chimeric_detection": {
            "true_positives": tp,
            "false_positives": fp,
            "false_negatives": fn,
            "true_negatives": tn,
            "sensitivity": round(results.chimeric_sensitivity, 4),
            "specificity": round(results.chimeric_specificity, 4),
            "precision": round(results.chimeric_precision, 4),
            "recall": round(results.chimeric_recall, 4),
            "f1_score": round(results.chimeric_f1, 4),
            "accuracy": round(results.accuracy, 4),
        },
        "host_detection": {
            "true_positives": results.host_tp,
            "false_positives": results.host_fp,
            "false_negatives": results.host_fn,
        },
        "backbone_detection": {
            "true_positives": results.backbone_tp,
            "false_positives": results.backbone_fp,
            "false_negatives": results.backbone_fn,
            "sensitivity": round(results.backbone_sensitivity, 4),
            "precision": round(results.backbone_precision, 4),
        },
        "reference_pair_metrics": {
            "correct": results.ref_pair_correct,
            "total": results.ref_pair_total,
            "accuracy": round(results.ref_pair_accuracy, 4),
        },
        "error_analysis": {
            "num_false_positives": fp,
            "num_false_negatives": fn,
            "fp_examples": results.fp_read_names[:20],
            "fn_examples": results.fn_read_names[:20],
        },
    }

    results_json = output_dir / "benchmark_results.json"
    with open(results_json, "w") as fh:
        json.dump(results_dict, fh, indent=2)

    # Write detailed per-read CSV
    results_csv = output_dir / "benchmark_per_read.csv"
    with open(results_csv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=[
            "read_name", "ground_truth_chimeric", "ground_truth_backbone",
            "pipeline_detected_chimeric", "pipeline_classified_backbone",
            "classification", "ground_truth_refs", "detected_refs",
            "ground_truth_category",
        ])
        writer.writeheader()
        for read_name in sorted(all_read_names):
            gt = ground_truth[read_name]
            is_detected = read_name in detected_chimeric_names
            classification = pipeline_classifications.get(read_name, {})
            is_backbone = classification.get("category", "") == "backbone_contamination"

            if gt["is_chimeric"] and is_detected:
                cls = "TP"
            elif gt["is_chimeric"] and not is_detected:
                cls = "FN"
            elif not gt["is_chimeric"] and is_detected:
                cls = "FP"
            else:
                cls = "TN"

            writer.writerow({
                "read_name": read_name,
                "ground_truth_chimeric": int(gt["is_chimeric"]),
                "ground_truth_backbone": int(gt["has_backbone"]),
                "pipeline_detected_chimeric": int(is_detected),
                "pipeline_classified_backbone": int(is_backbone),
                "classification": cls,
                "ground_truth_refs": "|".join(gt["refs_order"]),
                "detected_refs": "|".join(
                    pipeline_chimeric.get(read_name, {}).get("refs_order", [])
                ),
                "ground_truth_category": gt["read_category"],
            })

    # Print summary
    logging.info("=" * 60)
    logging.info("BENCHMARK RESULTS")
    logging.info("=" * 60)
    logging.info("CHIMERIC DETECTION:")
    logging.info("  True Positives:  %d", tp)
    logging.info("  False Positives: %d", fp)
    logging.info("  False Negatives: %d", fn)
    logging.info("  True Negatives:  %d", tn)
    logging.info("  Sensitivity:     %.4f", results.chimeric_sensitivity)
    logging.info("  Specificity:     %.4f", results.chimeric_specificity)
    logging.info("  Precision:       %.4f", results.chimeric_precision)
    logging.info("  F1 Score:        %.4f", results.chimeric_f1)
    logging.info("  Accuracy:        %.4f", results.accuracy)
    logging.info("  Ref Pair Acc:    %.4f", results.ref_pair_accuracy)
    logging.info("")
    logging.info("BACKBONE DETECTION:")
    logging.info("  True Positives:  %d", results.backbone_tp)
    logging.info("  False Positives: %d", results.backbone_fp)
    logging.info("  False Negatives: %d", results.backbone_fn)
    logging.info("  Sensitivity:     %.4f", results.backbone_sensitivity)
    logging.info("  Precision:       %.4f", results.backbone_precision)
    logging.info("")
    logging.info("HOST CHIMERIC DETECTION:")
    logging.info("  True Positives:  %d", results.host_tp)
    logging.info("  False Positives: %d", results.host_fp)
    logging.info("  False Negatives: %d", results.host_fn)
    logging.info("=" * 60)

    return results


# ============================================================
# CLI entry point
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description="AAV Chimeric Read Simulator & Benchmarking Tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent("""\
Examples:
  # Simulate chimeric reads (10%% chimeric, 5%% backbone)
  python aav_simulator.py simulate \\
    --ref-dir /path/to/refs \\
    --output-dir /path/to/output/chimeric \\
    --transgene-name "pAAV-CMV-eGFP" \\
    --itr-5-start 0 --itr-5-end 145 \\
    --itr-3-start 4331 --itr-3-end 4472 \\
    --chimeric-proportion 0.10 \\
    --backbone-proportion 0.05

  # Simulate non-chimeric reads (negative control)
  python aav_simulator.py simulate \\
    --ref-dir /path/to/refs \\
    --output-dir /path/to/output/nonchimeric \\
    --transgene-name "pAAV-CMV-eGFP" \\
    --itr-5-start 0 --itr-5-end 145 \\
    --itr-3-start 4331 --itr-3-end 4472 \\
    --chimeric-proportion 0.0 \\
    --backbone-proportion 0.0

  # Benchmark pipeline results against ground truth
  python aav_simulator.py benchmark \\
    --ground-truth-dir /path/to/output/chimeric/simulated_reads_1 \\
    --pipeline-output-dir /path/to/pipeline_output/samples/simulated_reads_1/analysis \\
    --output-dir /path/to/benchmark_results
        """),
    )

    subparsers = parser.add_subparsers(dest="command", help="Command to run")

    # --- Simulate subcommand ---
    sim_parser = subparsers.add_parser("simulate", help="Simulate nanopore reads")
    sim_parser.add_argument("--ref-dir", required=True,
                            help="Reference FASTA directory")
    sim_parser.add_argument("--output-dir", required=True,
                            help="Output directory")
    sim_parser.add_argument("--transgene-name", required=True,
                            help="Sequence ID of transgene plasmid in FASTA header")

    # Coordinate annotations
    sim_parser.add_argument("--itr-5-start", type=int, default=0,
                            help="Start of 5' ITR (0-based, default: 0)")
    sim_parser.add_argument("--itr-5-end", type=int, default=145,
                            help="End of 5' ITR (default: 145)")
    sim_parser.add_argument("--transgene-start", type=int, default=145,
                            help="Start of transgene cassette (default: 145)")
    sim_parser.add_argument("--transgene-end", type=int, default=4331,
                            help="End of transgene cassette (default: 4331)")
    sim_parser.add_argument("--itr-3-start", type=int, default=4331,
                            help="Start of 3' ITR (default: 4331)")
    sim_parser.add_argument("--itr-3-end", type=int, default=4472,
                            help="End of 3' ITR (default: 4472)")

    # Simulation parameters
    sim_parser.add_argument("--num-files", type=int, default=2,
                            help="Number of FASTQ files (default: 2)")
    sim_parser.add_argument("--reads-min", type=int, default=5000,
                            help="Min reads per file (default: 5000)")
    sim_parser.add_argument("--reads-max", type=int, default=10000,
                            help="Max reads per file (default: 10000)")
    sim_parser.add_argument("--chimeric-proportion", type=float, default=0.10,
                            help="Proportion of chimeric reads (default: 0.10)")
    sim_parser.add_argument("--backbone-proportion", type=float, default=0.05,
                            help="Proportion of backbone reads (default: 0.05)")
    sim_parser.add_argument("--host-proportion", type=float, default=0.05,
                            help="Proportion of host DNA reads (default: 0.05)")
    sim_parser.add_argument("--seed", type=int, default=42,
                            help="Random seed (default: 42)")
    sim_parser.add_argument("--error-rate", type=float, default=0.05,
                            help="Base error rate (default: 0.05)")
    sim_parser.add_argument("--min-length", type=int, default=500,
                            help="Min read length (default: 500)")
    sim_parser.add_argument("--max-length", type=int, default=6000,
                            help="Max read length (default: 6000)")

    # Reference file list
    sim_parser.add_argument("--ref-files", nargs="+", default=None,
                            help="Reference FASTA filenames (within --ref-dir)")

    # --- Benchmark subcommand ---
    bench_parser = subparsers.add_parser("benchmark",
                                         help="Benchmark pipeline vs ground truth")
    bench_parser.add_argument("--ground-truth-dir", required=True,
                              help="Directory with ground truth CSVs")
    bench_parser.add_argument("--pipeline-output-dir", required=True,
                              help="Directory with pipeline output (analysis folder)")
    bench_parser.add_argument("--output-dir", required=True,
                              help="Directory for benchmark results")

    args = parser.parse_args()

    if args.command == "simulate":
        # Build ref_files list
        if args.ref_files:
            ref_files = args.ref_files
        else:
            # Auto-detect FASTA files in ref-dir
            ref_dir = Path(args.ref_dir)
            ref_files = [
                p.name for p in sorted(ref_dir.iterdir())
                if p.suffix.lower() in {".fa", ".fasta"}
            ]

        cfg = SimulationConfig(
            ref_dir=args.ref_dir,
            output_dir=args.output_dir,
            transgene_name=args.transgene_name,
            itr_5_start=args.itr_5_start,
            itr_5_end=args.itr_5_end,
            transgene_start=args.transgene_start,
            transgene_end=args.transgene_end,
            itr_3_start=args.itr_3_start,
            itr_3_end=args.itr_3_end,
            ref_files=ref_files,
            num_fastq_files=args.num_files,
            reads_per_fastq_min=args.reads_min,
            reads_per_fastq_max=args.reads_max,
            chimeric_proportion=args.chimeric_proportion,
            backbone_proportion=args.backbone_proportion,
            host_dna_proportion=args.host_proportion,
            random_seed=args.seed,
            base_error_rate=args.error_rate,
            min_read_length=args.min_length,
            max_read_length=args.max_length,
        )
        run_simulation(cfg)

    elif args.command == "benchmark":
        # Set up basic logging for benchmark mode
        output_dir = Path(args.output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        logging.basicConfig(
            level=logging.INFO,
            format="%(asctime)s | %(levelname)s | %(message)s",
            handlers=[
                logging.FileHandler(output_dir / "benchmark.log"),
                logging.StreamHandler(sys.stdout),
            ],
        )

        benchmark_pipeline(
            ground_truth_dir=Path(args.ground_truth_dir),
            pipeline_output_dir=Path(args.pipeline_output_dir),
            output_dir=output_dir,
        )

    else:
        parser.print_help()


if __name__ == "__main__":
    main()
