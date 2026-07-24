"""Unit tests for pure functions in aav_chimera.simulator."""

import random

import pytest

from aav_chimera.simulator import (
    SimulationConfig,
    gc_content,
    generate_quality_scores,
    is_helper_reference,
    is_host_reference,
    reverse_complement,
    simulate_nanopore_errors,
)


def test_reverse_complement_basic():
    assert reverse_complement("ATGC") == "GCAT"
    assert reverse_complement("AATTCCGG") == "CCGGAATT"


def test_reverse_complement_is_involution():
    seq = "ACGTACGTTTGGCCA"
    assert reverse_complement(reverse_complement(seq)) == seq


def test_reverse_complement_handles_n_and_case():
    assert reverse_complement("nnAC") == "GTNN"


def test_gc_content():
    assert gc_content("GGCC") == 100.0
    assert gc_content("ATAT") == 0.0
    assert gc_content("ATGC") == 50.0
    assert gc_content("") == 0.0


@pytest.mark.parametrize(
    "name,expected",
    [
        ("chr1", True),
        ("chr22", True),
        ("chrX", True),
        ("chrM", True),
        ("22", True),
        ("pAAV-CMV-eGFP", False),
        ("ref_rep_cap", False),
        (None, False),
    ],
)
def test_is_host_reference(name, expected):
    patterns = SimulationConfig(ref_dir=".", output_dir=".").host_ref_patterns
    assert is_host_reference(name, patterns) is expected


def test_is_helper_reference():
    patterns = SimulationConfig(ref_dir=".", output_dir=".").helper_ref_patterns
    assert is_helper_reference("ref_rep_cap", patterns) is True
    assert is_helper_reference("some_helper_plasmid", patterns) is True
    assert is_helper_reference("chr1", patterns) is False


def test_quality_scores_length_and_range():
    q = generate_quality_scores(200)
    assert len(q) == 200
    # Phred+33 encoding must stay in a sane ASCII band.
    assert all(35 <= ord(c) <= 73 for c in q)


def test_error_model_preserves_nonempty_and_alphabet():
    random.seed(0)
    cfg = SimulationConfig(ref_dir=".", output_dir=".")
    seq = "ACGT" * 250
    out = simulate_nanopore_errors(seq, cfg)
    assert len(out) > 0
    assert set(out) <= set("ACGTN")


def test_error_model_empty_input():
    cfg = SimulationConfig(ref_dir=".", output_dir=".")
    assert simulate_nanopore_errors("", cfg) == ""


def test_config_derived_properties():
    cfg = SimulationConfig(ref_dir=".", output_dir=".")
    assert cfg.itr_length == cfg.itr_5_end - cfg.itr_5_start
    assert cfg.expected_genome_length == cfg.itr_3_end - cfg.itr_5_start


def test_backbone_regions():
    cfg = SimulationConfig(ref_dir=".", output_dir=".")
    regions = cfg.backbone_regions(ref_length=6000)
    # 5' region [0, itr_5_start) is empty by default (itr_5_start == 0);
    # 3' region [itr_3_end, ref_length) should be present.
    assert (cfg.itr_3_end, 6000) in regions
