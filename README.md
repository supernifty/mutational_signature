# Mutational Signature

Command-line tools for counting mutational contexts, decomposing samples into reference signatures, and generating publication-oriented plots.

The repository now supports a modern `uv` workflow for development while remaining installable with standard `pip`.

## What This Repository Provides

- `count`: derive SBS, DBS, or indel context counts from VCF or MAF input
- `decompose`: fit a sample's counts against a reference signature matrix
- `plot-counts`: render SBS or indel context plots
- `plot-components`: render signature exposure panels across samples
- utility commands for bootstrap analysis, similarity checks, context annotation, and signature conversion

Reference signature files are bundled under [`mutational_signature/data`](mutational_signature/data).

## Quickstart

### Recommended: `uv`

```bash
uv sync --extra dev
uv run count --genome example/mini_reference.fa --vcf example/mini_sample.vcf > mini.count
uv run decompose --signatures example/mini_signatures.tsv --counts mini.count > mini.exposures
uv run plot-counts --target mini.png < mini.count
```

### Standard `pip`

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install .
count --genome example/mini_reference.fa --vcf example/mini_sample.vcf > mini.count
decompose --signatures example/mini_signatures.tsv --counts mini.count > mini.exposures
plot-counts --target mini.png < mini.count
```

If you need the optional bootstrap boxplots backed by `plotme`, install the extra:

```bash
uv sync --extra plotme
```

or

```bash
python -m pip install .[plotme]
```

## End-to-End Example

The [`example`](example) directory contains a small self-contained workflow:

- [`mini_reference.fa`](example/mini_reference.fa)
- [`mini_sample.vcf`](example/mini_sample.vcf)
- [`mini_signatures.tsv`](example/mini_signatures.tsv)
- [`README.md`](example/README.md)

That example is also exercised in automated tests.

## Common Workflows

### Count contexts from a VCF

```bash
uv run count --genome genome.fa --vcf sample.vcf.gz > sample.count
```

### Decompose counts against bundled COSMIC signatures

```bash
uv run decompose \
  --signatures mutational_signature/data/signatures_cosmic_v3_sbs.txt \
  --counts sample.count > sample.exposures
```

### Plot SBS context counts

```bash
uv run plot-counts --target sample.sbs.png < sample.count
```

### Plot signature exposures across multiple samples

```bash
uv run plot-components --target cohort.exposures.png < combined.exposures.tsv
```

### Count context opportunity in target regions

```bash
uv run count-contexts --genome genome.fa --bed regions.bed > opportunity.tsv
```

### Adjust counts for capture opportunity

```bash
uv run adjust-counts \
  --adjust_from opportunity.wgs.tsv \
  --adjust_to opportunity.panel.tsv < sample.count
```

### Convert downloaded signatures into repository format

```bash
uv run convert-signatures --conversion sbs32 < downloaded_signatures.csv > converted.tsv
```

## Development

```bash
uv sync --extra dev
uv run pytest
uv build
```

The compatibility path for environments that do not use `uv` is:

```bash
python -m pip install -e .[dev]
pytest
python -m build
```

## Project Layout

- [`mutational_signature`](mutational_signature): package source and CLI modules
- [`mutational_signature/data`](mutational_signature/data): bundled signature resources
- [`example`](example): runnable sample data and walkthrough
- [`tests`](tests): smoke tests and small fixtures
- [`docs`](docs): contributor, CLI, and release notes

## Documentation

- [`docs/cli.md`](docs/cli.md)
- [`docs/development.md`](docs/development.md)
- [`docs/releasing.md`](docs/releasing.md)

## Current Scope and Caveats

- The codebase is CLI-first and still contains research-oriented utilities with uneven interface polish.
- Some advanced commands assume domain familiarity or specific upstream file conventions.
- Large bundled reference files should be updated carefully and with provenance notes.

## License

See [`LICENSE.md`](LICENSE.md).
