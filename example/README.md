# Example Workflow

This directory contains a minimal self-contained workflow that does not require a large reference genome.

## Files

- `mini_reference.fa`: tiny FASTA with a single chromosome
- `mini_sample.vcf`: one SBS variant on that chromosome
- `mini_signatures.tsv`: a small signature definition containing the expected context

## Run

```bash
uv run count --genome example/mini_reference.fa --vcf example/mini_sample.vcf > mini.count
uv run decompose --signatures example/mini_signatures.tsv --counts mini.count > mini.exposures
uv run plot-counts --target mini.png < mini.count
```

## Expected Outcome

- `mini.count` contains one `ACA>A` event
- `mini.exposures` reports a dominant `SBS_Test_A`
- `mini.png` is a non-empty SBS plot image
