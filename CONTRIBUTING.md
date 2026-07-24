# Contributing

Thanks for your interest in improving AAV-Chimera.

## Development setup

```bash
git clone https://github.com/Victor-Alfred/aav-chimera
cd aav-chimera
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
```

## Before opening a PR

```bash
ruff check src tests scripts   # lint
pytest                          # unit tests
```

Both run in CI on Python 3.9, 3.11 and 3.12. Please add or update tests for any behaviour
change, and keep new functions small and independently testable — the CIGAR/SA-tag/junction
helpers in `aav_chimera.qc` are good models.

## Scope of good first issues

- Improving reference-pair accuracy on multi-segment chimeras (currently ~83%).
- Additional error-model presets for other Nanopore chemistries.
- More genome-structure edge cases in the classifier test coverage.
