"""Unit tests for pure functions in aav_chimera.qc (CIGAR / SA-tag / junction logic)."""

from aav_chimera.qc import (
    aligned_reference_length,
    analyse_junction,
    gc_percent,
    get_query_span_from_cigar,
    merge_same_ref_neighbours,
    parse_cigar_ops,
    parse_sa_tag,
    query_consuming_length,
    softclip_left,
    softclip_right,
)


def test_parse_cigar_ops():
    assert parse_cigar_ops("10S50M2D48M") == [
        (10, "S"), (50, "M"), (2, "D"), (48, "M"),
    ]
    assert parse_cigar_ops("") == []


def test_query_consuming_length():
    # M + I + S + =/X consume query; D/N/H do not.
    assert query_consuming_length("50M") == 50
    assert query_consuming_length("10S50M10S") == 70
    assert query_consuming_length("50M2D48M") == 98


def test_aligned_reference_length():
    # M + D + N + =/X consume reference; I/S do not.
    assert aligned_reference_length("50M") == 50
    assert aligned_reference_length("50M2D48M") == 100
    assert aligned_reference_length("10S50M10S") == 50


def test_softclip_left_right():
    assert softclip_left("10S50M") == 10
    assert softclip_left("50M10S") == 0
    assert softclip_right("50M10S") == 10
    assert softclip_right("10S50M") == 0


def test_query_span_from_cigar():
    # 10S then 50M then 5S over a 65bp read -> aligned query span [10, 60).
    qstart, qend = get_query_span_from_cigar("10S50M5S")
    assert qstart == 10
    assert qend == 60


def test_gc_percent():
    assert gc_percent("GGCC") == 100.0
    assert gc_percent("atat") == 0.0
    assert gc_percent("ATGC") == 50.0
    assert gc_percent("") == 0.0


def test_parse_sa_tag():
    sa = "chr20,1000,+,50M100S,60,2;ref_rep_cap,55,-,100S50M,55,1;"
    entries = parse_sa_tag(sa)
    assert len(entries) == 2
    assert entries[0]["ref"] == "chr20"
    assert entries[0]["pos1"] == 1000
    assert entries[0]["strand"] == "+"
    assert entries[0]["mapq"] == 60
    assert entries[1]["ref"] == "ref_rep_cap"
    assert entries[1]["nm"] == 1


def test_parse_sa_tag_skips_malformed():
    assert parse_sa_tag(";;bad,entry;") == []


def test_analyse_junction_microhomology():
    # left ends with "AGCT", right starts with "AGCT" -> 4bp microhomology.
    res = analyse_junction("TTTTAGCT", "AGCTgggg".upper())
    assert res["microhomology_length"] == 4
    assert res["microhomology_seq"] == "AGCT"


def test_analyse_junction_no_microhomology():
    res = analyse_junction("AAAA", "TTTT")
    assert res["microhomology_length"] == 0


def test_merge_same_ref_neighbours():
    segs = [
        {"ref": "chr1", "strand": "+", "qstart": 0, "qend": 100,
         "ref_start": 500, "ref_end": 600, "mapq": 60, "aligned_length": 100},
        {"ref": "chr1", "strand": "+", "qstart": 105, "qend": 200,
         "ref_start": 600, "ref_end": 695, "mapq": 55, "aligned_length": 95},
        {"ref": "ref_rep_cap", "strand": "-", "qstart": 210, "qend": 300,
         "ref_start": 10, "ref_end": 100, "mapq": 50, "aligned_length": 90},
    ]
    merged = merge_same_ref_neighbours(segs)
    # First two (same ref/strand, adjacent) merge; third stays separate.
    assert len(merged) == 2
    assert merged[0]["qend"] == 200
    assert merged[0]["aligned_length"] == 195
    assert merged[1]["ref"] == "ref_rep_cap"


def test_merge_same_ref_neighbours_empty():
    assert merge_same_ref_neighbours([]) == []
