# CLI Reference

This package installs the commands listed in `pyproject.toml`. Each command is a thin wrapper around a module in `mutational_signature/`, so the same behavior is available with either the installed command or `python -m mutational_signature.<module>`.

Most tools read and write tab-separated text. Commands that say "stdin" expect redirected input, for example:

```bash
uv run plot-counts --target sample.png < sample.count
```

## Count and Annotate Variants

| Command | Module | Purpose |
| --- | --- | --- |
| `count` | `count.py` | Count SBS contexts from one or more VCF files against `--genome`. With `--maf_sample`, treats `--vcf` inputs as MAF files and filters to that sample. Optional modes include `--doublets`, `--indels`, `--just_indels`, `--transcripts`, wider SBS context with `--mer`, and VCF INFO weighting with `--weight_field`. Writes count TSVs with `Variation`, `Count`, `Probability`, and strand-bias columns. |
| `annotate-context` | `annotate_context.py` | Add `sig_context` and optional surrounding sequence to each VCF or MAF row using `--genome` and `--vcf`. With `--signatures`, also annotates the signature with the highest probability for each context. Writes annotated VCF or MAF TSV to stdout. |
| `annotate-context-summary` | `annotate_context_summary.py` | Summarize `surrounding` annotations already present in one or more annotated VCFs. `--sequence` is the flank length originally annotated; `--lengths` selects shorter central windows to count. Writes contexts as rows and input VCFs as columns. |
| `annotate-maf` | `annotate_maf.py` | MAF-specific context annotation/count helper. Reads `--maf`, resolves contexts against `--genome`, appends a `context` column, and supports `--doublets`, `--indels`, and `--just_indels`. |
| `count-contexts` | `count_contexts.py` | Count reference-context opportunity in regions from `--bed` against `--genome`. Supports wider context with `--context_length`, optional custom context names, context filtering, and `--write_contexts` for per-position context records. Writes `Variation`, `Count`, `Probability`. The `--indels` and `--doublets` flags are parsed but marked unsupported in the code. |
| `count-indels` | `count_indels.py` | Count indel lengths by repeat annotation. Reads VCF from stdin and an annotation BED from `--annotation` whose fourth field contains `repeat=...`. Writes `Type`, `Count`, `Proportion`. |
| `count-maf-indels` | `count_maf_indels.py` | MAF equivalent of `count-indels`. Reads gzipped `--maf`, filters `--sample`, overlaps `--annotation`, and writes indel repeat/length counts. |
| `extended-context` | `extended_context.py` | Test custom sequence-context rules against a VCF or MAF. Rules look like `T>*,-3=A` or `A>*,+3=T`; context clauses support equals, membership with `~`, and exclusion with `!`. Writes per-rule pass counts and optional per-variant annotations with `--output`. |

## Decompose and Model Signatures

| Command | Module | Purpose |
| --- | --- | --- |
| `decompose` | `decompose.py` | Fit a count table from `--counts` to reference signatures from `--signatures`. Supports `--metric` (`cosine`, `euclidean`, `kl`, `l1`), `--solver` (`basin`, `grid`), `--max_sigs`, context cutoffs, strand signatures, and `--count_column`. Writes one signature/value row per exposure plus error and mutation total rows. |
| `decompose-bulk` | `decompose_bulk.py` | Run `decompose` over multiple `--counts` files and combine the results into one TSV with one row per file. Uses the same decomposition options as `decompose`. |
| `decompose-likelihood` | `decompose_likelihood.py` | Score each signature independently by log likelihood for a count table, then convert scores to posterior probabilities under a uniform prior. Useful for low-count comparisons where mixture fitting is not desired. |
| `bootstrap` | `bootstrap.py` | Read a count table from stdin, repeatedly resample contexts, decompose each bootstrap sample, and report confidence intervals for error, count, and signature proportions. Optional `--plot` requires the `plotme` extra. |
| `bootstrap-simple` | `bootstrap_simple.py` | Downsample one or more count files by percentage or fixed mutation count, decompose each replicate, and report similarity to the full-data estimate. Supports median aggregation with `--aggregate median`. |
| `generate-signatures` | `generate.py` | Derive de novo signatures from multiple count files using NMF for each requested `--components` value. Can write component stats, compare generated signatures to a reference file, perform bootstraps, and save generated signatures. |
| `refine-signatures` | `refine_signatures.py` | Read counts from stdin, decompose against `--signatures`, and iteratively remove signatures whose removal does not increase reconstruction error beyond `--minimum_error`. Writes the final decomposition result. |
| `reduce-similarity` | `reduce_similarity.py` | Read a signature matrix from stdin and repeatedly merge the most similar pair until all pairwise cosine similarities are at or below `--max_similarity`. `--merge_strategy first` keeps the first signature; the default averages merged signatures. |
| `linear-dependence` | `linear_dependence.py` | Approximate each signature in `--signatures` as a non-negative weighted mixture of other signatures and report the cosine-error residual plus the largest component weights. With `--signatures2`, only signatures in that second file are used as candidate components. This is an approximation/similarity diagnostic, not a formal proof of exact linear dependence. |
| `stability` | `stability.py` | Read counts from stdin, compute a point decomposition, then remove one mutation from each observed context in turn and report how much each signature estimate changes. |
| `positive-by-chance` | `positive_by_chance.py` | Read counts from stdin, sample random subsets of `--variant_count` variants according to the observed context distribution, decompose each sample, and report how often a chosen `--signature_sum` exceeds `--threshold`. |
| `simulate-signatures` | `simulate.py` | Generate synthetic count files from signature definitions. Background exposures are `signature,mean,sd`; injected signals are `signature,proportion`. `NUM` in `--output_template` is replaced by the simulation index. |

## Compare and Summarize Results

| Command | Module | Purpose |
| --- | --- | --- |
| `compare-signatures` | `compare_signatures.py` | Compute pairwise similarity between signatures in one matrix, or between two matrices with `--signatures2`. Supports `--measure cosine`, `pearson`, or `pvalue`, and converts four-base COSMIC-style context headers unless `--no_convert` is set. |
| `compare-exposures` | `compare_exposures.py` | Compare two exposure files by cosine or Pearson correlation, ignoring `Error`, `Mutations`, and `Variation`. Writes a one-line comparison and can draw an annotated scatter plot with `--plot`. |
| `compare-exposures-to-signatures` | `compare_exposures_to_signatures.py` | Compare a count/exposure table with `Variation` and `Count` columns against every signature in a signature definition file. Writes `Sig` and `cosine_similarity`. |
| `combine-counts` | `combine_counts.py` | Combine multiple count files into a matrix with one row per input file and one column per observed context. Uses `Count` by default or `Probability` with `--use_probability`. |
| `combine-signatures` | `combine_signatures.py` | Combine multiple `decompose` output files into a sample-by-signature matrix using a reference signature file to define columns. Adds `Error` and `Mutations` columns when present. |
| `context-best-sig` | `context_best_sig.py` | Read a signature matrix from stdin and report the top `--ranks` signatures for each context. Optional `--pie` writes per-context pie charts by replacing `CTX` in the filename template. |
| `context-likelihood` | `context_likelihood.py` | Read a signature matrix from stdin and write, for each context, the percentage contribution of each signature plus a short `best` summary of the top signatures. |
| `primary-context` | `primary_context.py` | Read a signature matrix from stdin and report the highest-probability context for each signature. Four-base context columns are converted to an explicit substitution representation. |
| `exclude-contexts` | `exclude_contexts.py` | Read a signature matrix from stdin, remove contexts named as positional command-line arguments, renormalize each signature row, and write the reduced matrix. Effects are reported on stderr. |
| `adjust-counts` | `adjust_counts.py` | Read a count table from stdin and adjust each context by the ratio between two opportunity/count tables supplied as `--adjust_from` and `--adjust_to`. By default it matches on the pre-substitution genomic context; `--entire_variation` matches the full `Variation`. |
| `convert-signatures` | `convert.py` | Convert external COSMIC/signature tables from stdin into this repository's `Sig`-by-context TSV format. Supported conversions include `sbs`, `sbs32`, `db`, `db32`, `id`, `id32`, `sbstx`, and `sbs_signal`. |

## Plotting

| Command | Module | Purpose |
| --- | --- | --- |
| `plot-counts` | `plot_counts.py` | Read a count TSV from stdin and draw an SBS or indel context bar plot. Uses the probability column when present; `--normalize` renormalizes displayed values. `--type` supports `sbs` and `id`. |
| `plot-sigs` | `plot_sigs.py` | Read a signature matrix from stdin and generate one context plot per signature. `--category sbs` plots SBS signatures; `--category ids` plots indel signatures. Output filenames are `--prefix` + signature name + `--suffix`. |
| `plot-components` | `plot_components.py` | Read a combined exposure matrix from stdin and draw stacked signature-exposure bars. Supports ordering, custom colors, thresholds, labels, vertical layout, indicator strips, separator lines, and renormalization with an `Other` component. |
| `histograms` | `histograms.py` | Plot histograms of signature proportions across exposure files. Reads separate decompose outputs by default or a combined sample-by-signature matrix with `--combined`. |

## Variant-Level Signature Assignment

| Command | Module | Purpose |
| --- | --- | --- |
| `assign-signatures` | `assign_signatures.py` | Annotate an already context-annotated VCF or MAF with posterior signature likelihoods per variant. `--definitions` supplies one or more signature matrices; `--signatures` supplies prior exposure estimates, otherwise a uniform prior is used. Optional `--artefacts` can add artefact likelihoods. |

## Additional Package Scripts

`spikiness.py` is present in `mutational_signature/` but is not installed as a console command. Run it directly with:

```bash
python mutational_signature/spikiness.py < signatures.tsv > spikiness.tsv
```

It reads a signature matrix from stdin and writes per-signature standard deviation and Gini coefficient.

`cli.py` only contains wrappers used by the installed console commands.

## Bundled Signature Data

Reference signatures and helper tables live under `mutational_signature/data/`. The most commonly used files are:

- `signatures_cosmic_v3_sbs.txt`
- `signatures_cosmic_v3.1_sbs.txt`
- `signatures_cosmic_v3.4_sbs.txt`
- `signatures_cosmic_v3_id.txt`
- `signatures_cosmic_v3.1_id.txt`
- `signatures_cosmic_v3_dbs.txt`

Use the exact reference file that matches the count context type being analyzed: SBS counts with SBS signatures, DBS counts with DBS signatures, and indel counts with ID signatures.

## Help

Every installed command exposes argparse help:

```bash
uv run count --help
uv run decompose --help
uv run linear-dependence --help
```
