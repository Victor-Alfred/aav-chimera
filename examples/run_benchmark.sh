#!/usr/bin/env bash
#
# End-to-end reproduction of the benchmark:  simulate -> qc -> benchmark
#
# Requires real reference FASTAs in $REF_DIR:
#   - your transgene/plasmid (e.g. pAAV-CMV-eGFP.fasta)  [ITR coords below must match]
#   - rep/cap and helper plasmids
#   - a host genome (e.g. hg38.fa)
# plus external tools on PATH: minimap2, samtools, porechop, NanoFilt, NanoPlot.
#
# Usage:
#   REF_DIR=/path/to/refs ./examples/run_benchmark.sh
#
set -euo pipefail

REF_DIR="${REF_DIR:?Set REF_DIR to the directory containing reference FASTAs}"
OUT="${OUT:-./benchmark_run}"
TRANSGENE="${TRANSGENE:-pAAV-CMV-eGFP}"

SIM_DIR="$OUT/simulated"
PIPE_DIR="$OUT/pipeline_output"
BENCH_DIR="$OUT/benchmark_results"

echo ">> 1/3  Simulating reads with known ground truth"
aav-chimera-sim simulate \
  --ref-dir "$REF_DIR" \
  --output-dir "$SIM_DIR" \
  --transgene-name "$TRANSGENE" \
  --itr-5-start 0 --itr-5-end 145 \
  --itr-3-start 4331 --itr-3-end 4472 \
  --chimeric-proportion 0.10 \
  --backbone-proportion 0.05

echo ">> 2/3  Running the AAV-Chimera pipeline on the simulated FASTQs"
aav-chimera \
  --raw-fastq-dir "$SIM_DIR" \
  --work-dir "$PIPE_DIR" \
  --refs-dir "$REF_DIR" \
  --transgene-name "$TRANSGENE" \
  --itr-5-start 0 --itr-5-end 145 \
  --itr-3-start 4331 --itr-3-end 4472

echo ">> 3/3  Scoring the pipeline against ground truth"
aav-chimera-sim benchmark \
  --ground-truth-dir "$SIM_DIR/simulated_reads_1" \
  --pipeline-output-dir "$PIPE_DIR/samples/simulated_reads_1/analysis" \
  --output-dir "$BENCH_DIR"

echo ">> Regenerating figures"
cp "$BENCH_DIR"/benchmark_results.json results/ 2>/dev/null || true
cp "$BENCH_DIR"/benchmark_per_read.csv  results/ 2>/dev/null || true
python scripts/make_figures.py

echo "Done. See $BENCH_DIR and results/figures/"
