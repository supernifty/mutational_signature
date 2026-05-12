# CLI Reference

## Core Commands

- `count`: count mutation contexts from one or more VCF or MAF inputs
- `decompose`: fit a count table to a reference signature definition
- `plot-counts`: plot SBS or indel count distributions
- `plot-components`: plot exposure profiles across samples

## Analysis Utilities

- `bootstrap`: bootstrap decomposition uncertainty
- `bootstrap-simple`: downsample count files and compare stability
- `compare-signatures`: compare signature definitions by similarity
- `compare-exposures`: compare exposure outputs across samples
- `compare-exposures-to-signatures`: compare exposure tables back to signature definitions
- `generate-signatures`: derive de novo signatures from multiple count files
- `refine-signatures`: backward selection based on cosine error
- `reduce-similarity`: merge or simplify highly similar signatures
- `linear-dependence`: detect linearly dependent signature combinations
- `stability`: estimate sensitivity to loss of individual contexts
- `positive-by-chance`: estimate how often a signature combination appears by chance

## Annotation and Opportunity Utilities

- `annotate-context`: annotate VCF records with local sequence context
- `annotate-context-summary`: summarize context annotations across files
- `annotate-maf`: annotate MAF records with contexts
- `count-contexts`: count context opportunity across genomic regions
- `context-best-sig`: identify the best-matching signature per context
- `context-likelihood`: estimate likelihood by context
- `extended-context`: search for longer or custom sequence contexts
- `primary-context`: extract the primary context representation

## Format Conversion and Aggregation

- `adjust-counts`: adjust counts based on opportunity tables
- `combine-counts`: combine multiple count files into one table
- `combine-signatures`: combine exposure files into one table
- `convert-signatures`: convert external signature files into repository format
- `count-indels`: count annotated indels
- `count-maf-indels`: count indels from MAF input
- `decompose-bulk`: run decomposition over multiple count files
- `decompose-likelihood`: alternative decomposition workflow
- `exclude-contexts`: remove contexts from an analysis
- `histograms`: plot signature distributions
- `plot-sigs`: render reference signatures
- `simulate-signatures`: generate simulated signature mixtures

## Help

Each command exposes built-in help:

```bash
uv run count --help
uv run decompose --help
uv run plot-counts --help
```
