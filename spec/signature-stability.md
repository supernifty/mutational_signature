# Signature Stability Framework

## Summary

Implement `signature-stability` as the rigorous replacement for the repository's experimental bootstrap and stability approaches.

The command perturbs one observed mutation catalogue, reruns complete fitting from the full supplied signature matrix for every replicate, and reports bootstrap stability intervals for exposure proportions and mutation counts.

These intervals are conditional on the observed mutation catalogue, supplied reference matrix, fitting algorithm, and perturbation model. They must not be described as intervals that capture all biological or experimental uncertainty.

## Scope

- Add a new CLI command named `signature-stability`.
- Accept exactly one count table per invocation.
- Read counts from `--counts PATH`, or from stdin when `--counts` is omitted.
- Require `--sample` so output can be combined safely across samples.
- Keep the existing `decompose` CLI behavior compatible.
- Do not keep the old experimental bootstrap/stability ideas as independent maintained implementations. Compatibility wrappers are acceptable, but the rigorous implementation should be the recommended path.

## Perturbation Models

Default perturbation: `dirichlet-multinomial`.

Supported methods:

1. `dirichlet-multinomial`

   Draw:

   ```text
   p_boot ~ Dirichlet(x + alpha)
   x_boot ~ Multinomial(N, p_boot)
   ```

   Default `--alpha 0.1`.

2. `multinomial`

   Draw:

   ```text
   x_boot ~ Multinomial(N, x / N)
   ```

3. `downsample`

   Sample observed mutations without replacement using exactly one of:

   - `--fraction`
   - `--mutation-count`

Do not inject random sequencing errors or pseudocount mutations into observed catalogues.

Use deterministic random-number generation through `numpy.random.Generator`. Avoid the global `random` module in the new implementation.

## Fitting

Expose decomposition arguments including:

- `--signatures`
- `--metric`
- `--solver`
- `--max-sigs`
- `--context-cutoff`
- `--seed`

Support fitting modes:

1. `direct`: fit against the full supplied signature matrix.
2. `refine`: use redesigned backward elimination.

Refine mode starts from the full supplied matrix, fits all active signatures, tries removing each active signature, removes the signature whose removal produces the smallest acceptable error increase, and repeats until no removal satisfies:

```text
candidate_error <= current_error + --refine-minimum-error
```

Default `--refine-minimum-error 0.01`.

The point estimate and every replicate must use exactly the same fitting mode and fitting parameters. Replicates must not be restricted to signatures selected in the point estimate.

## Output

Support:

- `--summary-output PATH`
- `--replicates-output PATH`
- stdout for summary output when no summary path is supplied

Use TSV output.

Summary output is one row per sample/signature and includes:

- `Sample`
- `Signature`
- `point_exposure_mutations`
- `point_exposure_proportion`
- `bootstrap_mean_exposure_mutations`
- `bootstrap_mean_exposure_proportion`
- `bootstrap_median_exposure_mutations`
- `bootstrap_median_exposure_proportion`
- `bootstrap_sd_exposure_mutations`
- `bootstrap_sd_exposure_proportion`
- explicit bootstrap percentile bounds for mutations
- explicit bootstrap percentile bounds for proportions
- detection frequency
- detected replicate count
- total replicate count
- conditional detected medians and percentile bounds
- coefficient of variation where defined
- probability exposure exceeds a configurable mutation threshold
- point reconstruction error
- bootstrap reconstruction error summaries
- total mutations
- perturbation method
- random seed
- fitting mode
- metric
- solver

Replicate output is one row per sample/replicate/signature and includes:

- `Sample`
- `Replicate`
- `Signature`
- `exposure_mutations`
- `exposure_proportion`
- `detected`
- `total_mutations`
- `reconstruction_error`
- `fitting_status`

Every reference signature must have an entry for every replicate. If a signature is not selected or receives zero exposure, record zero rather than omitting it.

## Validation

Validate:

- counts are non-negative integers
- signature and count contexts match
- total count is greater than zero
- alpha is positive
- interval lies strictly between 0 and 1
- replicate count is positive
- warn below 200 replicates
- downsample arguments are mutually exclusive and valid
- requested downsampling does not exceed the observed count
- detection thresholds are non-negative

## Tests

Add fast deterministic tests covering:

1. Reproducibility with a fixed seed.
2. Correct multinomial sample size.
3. Correct Dirichlet-multinomial sample size.
4. Correct downsample size.
5. Missing signatures recorded as zero.
6. Detection frequency calculated correctly.
7. Interval calculations on a known synthetic result.
8. A pure synthetic signature is recovered stably at high mutation count.
9. A weak injected signature becomes less stable at low mutation count.
10. Two similar signatures show substitution or instability.
11. Direct and refine modes both run.
12. Existing CLI tests continue to pass.

## Documentation

Update `README.md` and `docs/cli.md` to document:

- what the intervals mean
- what they do not mean
- why every replicate must rerun fitting from the full signature set
- differences between Dirichlet-multinomial, multinomial, and downsampling
- interpretation of point exposure, bootstrap interval, and detection frequency
- a complete CLI example
- that earlier bootstrap/stability commands are deprecated compatibility entry points
- that a narrow interval does not protect against an incomplete or inappropriate signature reference set
