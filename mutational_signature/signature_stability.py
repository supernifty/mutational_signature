#!/usr/bin/env python
'''
  estimate stability intervals for mutational signature exposures
'''

import argparse
import csv
import logging
import sys

import numpy as np

import mutational_signature.decompose


def percentile_label(value):
  text = ('{:.6g}'.format(value)).replace('.', '_').replace('-', 'minus_')
  return text


def validate_counts(counts):
  if np.any(counts < 0):
    raise ValueError('counts must be non-negative')
  if not np.all(np.equal(counts, np.floor(counts))):
    raise ValueError('counts must be integers')
  total = int(np.sum(counts))
  if total <= 0:
    raise ValueError('total count must be greater than zero')
  return total


def perturb_counts(counts, method, rng, alpha=0.5, fraction=None, mutation_count=None):
  counts = np.array(counts, dtype=int)
  total = validate_counts(counts)

  if method == 'multinomial':
    return rng.multinomial(total, counts / total)

  if method == 'dirichlet-multinomial':
    if alpha <= 0:
      raise ValueError('alpha must be positive')
    p_boot = rng.dirichlet(counts + alpha)
    return rng.multinomial(total, p_boot)

  if method == 'downsample':
    if (fraction is None) == (mutation_count is None):
      raise ValueError('downsample requires exactly one of fraction or mutation_count')
    if fraction is not None:
      if fraction <= 0 or fraction > 1:
        raise ValueError('fraction must be greater than zero and no more than one')
      target = int(total * fraction)
    else:
      target = int(mutation_count)
    if target <= 0:
      raise ValueError('downsample mutation count must be greater than zero')
    if target > total:
      raise ValueError('downsample mutation count cannot exceed observed count')
    labels = np.repeat(np.arange(len(counts)), counts)
    selected = rng.choice(labels, size=target, replace=False)
    return np.bincount(selected, minlength=len(counts))

  raise ValueError('unknown perturbation method {}'.format(method))


def fit_direct(signature_matrix, counts, metric, solver, max_sigs):
  return mutational_signature.decompose.fit_signature_matrix(
    signature_matrix['matrix'],
    counts,
    signature_matrix['signature_names'],
    metric=metric,
    solver=solver,
    max_sigs=max_sigs,
  )


def expand_result(result, all_names, active_indices, total_mutations):
  proportions = np.zeros(len(all_names), dtype=float)
  for result_index, active_index in enumerate(active_indices):
    proportions[active_index] = result['signature_proportions'][result_index]
  return {
    'signature_names': np.array(all_names),
    'signature_proportions': proportions,
    'signature_mutations': proportions * total_mutations,
    'error': result['error'],
    'status': result.get('status', 'ok'),
  }


def context_cutoff_active_indices(signature_matrix, counts, context_cutoff):
  if context_cutoff is None:
    return list(range(len(signature_matrix['signature_names'])))
  active = []
  for signature_index in range(len(signature_matrix['signature_names'])):
    exclude = False
    for context_index, count in enumerate(counts):
      if count == 0 and signature_matrix['matrix'][context_index, signature_index] > context_cutoff:
        exclude = True
        break
    if not exclude:
      active.append(signature_index)
  if len(active) == 0:
    raise ValueError('context cutoff excluded all signatures')
  return active


def fit_refine(signature_matrix, counts, metric, solver, max_sigs, refine_minimum_error):
  all_names = signature_matrix['signature_names']
  active = list(range(len(all_names)))

  current = mutational_signature.decompose.fit_signature_matrix(
    signature_matrix['matrix'][:, active],
    counts,
    all_names[active],
    metric=metric,
    solver=solver,
    max_sigs=max_sigs,
  )
  current_error = current['error'][0]

  refining = True
  while refining and len(active) > 1:
    refining = False
    best_candidate = None
    for remove_index in list(active):
      candidate_active = [i for i in active if i != remove_index]
      candidate = mutational_signature.decompose.fit_signature_matrix(
        signature_matrix['matrix'][:, candidate_active],
        counts,
        all_names[candidate_active],
        metric=metric,
        solver=solver,
        max_sigs=max_sigs,
      )
      candidate_error = candidate['error'][0]
      if candidate_error <= current_error + refine_minimum_error:
        if best_candidate is None or candidate_error < best_candidate[0]:
          best_candidate = (candidate_error, candidate_active, candidate)
    if best_candidate is not None:
      current_error, active, current = best_candidate
      refining = True

  return expand_result(current, all_names, active, validate_counts(counts))


def subset_signature_matrix(signature_matrix, indices):
  return {
    'signature_names': signature_matrix['signature_names'][indices],
    'contexts': signature_matrix['contexts'],
    'matrix': signature_matrix['matrix'][:, indices],
  }


def fit_catalogue(signature_matrix, counts, mode, metric, solver, max_sigs, refine_minimum_error, context_cutoff=1e6):
  total_mutations = validate_counts(counts)
  active = context_cutoff_active_indices(signature_matrix, counts, context_cutoff)
  fit_matrix = subset_signature_matrix(signature_matrix, active)
  if mode == 'direct':
    result = fit_direct(fit_matrix, counts, metric, solver, max_sigs)
    return expand_result(result, signature_matrix['signature_names'], active, total_mutations)
  if mode == 'refine':
    refined = fit_refine(fit_matrix, counts, metric, solver, max_sigs, refine_minimum_error)
    proportions = np.zeros(len(signature_matrix['signature_names']), dtype=float)
    for name, prop in zip(refined['signature_names'], refined['signature_proportions']):
      original_index = active[list(fit_matrix['signature_names']).index(name)]
      proportions[original_index] = prop
    return {
      'signature_names': signature_matrix['signature_names'],
      'signature_proportions': proportions,
      'signature_mutations': proportions * total_mutations,
      'error': refined['error'],
      'status': refined.get('status', 'ok'),
    }
  raise ValueError('unknown fitting mode {}'.format(mode))


def make_replicate_rows(sample, replicate, signature_names, fit_result, total_mutations, detection_threshold_mutations, detection_threshold_proportion):
  rows = []
  for name, mutations, proportion in zip(signature_names, fit_result['signature_mutations'], fit_result['signature_proportions']):
    detected = mutations > detection_threshold_mutations and proportion > detection_threshold_proportion
    rows.append({
      'Sample': sample,
      'Replicate': replicate,
      'Signature': name,
      'exposure_mutations': mutations,
      'exposure_proportion': proportion,
      'detected': 'true' if detected else 'false',
      'total_mutations': total_mutations,
      'reconstruction_error': fit_result['error'][0],
      'fitting_status': fit_result.get('status', 'ok'),
    })
  return rows


def summarise(sample, signature_names, point_result, replicate_rows, interval, perturbation_method, seed, mode, metric, solver, exceeds_threshold_mutations):
  lower_pct = (1.0 - interval) / 2.0 * 100.0
  upper_pct = (1.0 + interval) / 2.0 * 100.0
  lower_label = percentile_label(lower_pct)
  upper_label = percentile_label(upper_pct)
  errors = np.array([float(r['reconstruction_error']) for r in replicate_rows if r['Signature'] == signature_names[0]])
  total_replicates = len(errors)
  rows = []

  for index, signature in enumerate(signature_names):
    sig_rows = [r for r in replicate_rows if r['Signature'] == signature]
    muts = np.array([float(r['exposure_mutations']) for r in sig_rows])
    props = np.array([float(r['exposure_proportion']) for r in sig_rows])
    detected = np.array([r['detected'] == 'true' for r in sig_rows])
    detected_muts = muts[detected]
    detected_props = props[detected]
    detected_count = int(np.sum(detected))
    mean_mut = float(np.mean(muts))
    sd_mut = float(np.std(muts))
    cv = '' if mean_mut == 0 else sd_mut / mean_mut

    row = {
      'Sample': sample,
      'Signature': signature,
      'point_exposure_mutations': point_result['signature_mutations'][index],
      'point_exposure_proportion': point_result['signature_proportions'][index],
      'bootstrap_mean_exposure_mutations': mean_mut,
      'bootstrap_mean_exposure_proportion': float(np.mean(props)),
      'bootstrap_median_exposure_mutations': float(np.median(muts)),
      'bootstrap_median_exposure_proportion': float(np.median(props)),
      'bootstrap_sd_exposure_mutations': sd_mut,
      'bootstrap_sd_exposure_proportion': float(np.std(props)),
      'bootstrap_percentile_{}_mutations'.format(lower_label): float(np.percentile(muts, lower_pct)),
      'bootstrap_percentile_{}_mutations'.format(upper_label): float(np.percentile(muts, upper_pct)),
      'bootstrap_percentile_{}_proportion'.format(lower_label): float(np.percentile(props, lower_pct)),
      'bootstrap_percentile_{}_proportion'.format(upper_label): float(np.percentile(props, upper_pct)),
      'detection_frequency': detected_count / total_replicates,
      'detected_replicates': detected_count,
      'total_replicates': total_replicates,
      'conditional_median_detected_mutations': '' if detected_count == 0 else float(np.median(detected_muts)),
      'conditional_median_detected_proportion': '' if detected_count == 0 else float(np.median(detected_props)),
      'conditional_percentile_{}_detected_mutations'.format(lower_label): '' if detected_count == 0 else float(np.percentile(detected_muts, lower_pct)),
      'conditional_percentile_{}_detected_mutations'.format(upper_label): '' if detected_count == 0 else float(np.percentile(detected_muts, upper_pct)),
      'conditional_percentile_{}_detected_proportion'.format(lower_label): '' if detected_count == 0 else float(np.percentile(detected_props, lower_pct)),
      'conditional_percentile_{}_detected_proportion'.format(upper_label): '' if detected_count == 0 else float(np.percentile(detected_props, upper_pct)),
      'coefficient_of_variation': cv,
      'probability_exceeds_threshold': float(np.mean(muts > exceeds_threshold_mutations)),
      'point_reconstruction_error': point_result['error'][0],
      'bootstrap_median_reconstruction_error': float(np.median(errors)),
      'bootstrap_percentile_{}_reconstruction_error'.format(lower_label): float(np.percentile(errors, lower_pct)),
      'bootstrap_percentile_{}_reconstruction_error'.format(upper_label): float(np.percentile(errors, upper_pct)),
      'total_mutations': int(np.sum(point_result['signature_mutations'])),
      'perturbation_method': perturbation_method,
      'seed': '' if seed is None else seed,
      'fitting_mode': mode,
      'metric': metric,
      'solver': solver,
    }
    rows.append(row)

  return rows


def write_tsv(path, rows, out=None):
  if len(rows) == 0:
    return
  fieldnames = list(rows[0].keys())
  if path is None:
    fh = out or sys.stdout
    writer = csv.DictWriter(fh, delimiter='\t', fieldnames=fieldnames, lineterminator='\n')
    writer.writeheader()
    writer.writerows(format_row(r) for r in rows)
  else:
    with open(path, 'w', newline='') as fh:
      writer = csv.DictWriter(fh, delimiter='\t', fieldnames=fieldnames, lineterminator='\n')
      writer.writeheader()
      writer.writerows(format_row(r) for r in rows)


def format_row(row):
  formatted = {}
  for key, value in row.items():
    if isinstance(value, float):
      formatted[key] = '{:.6g}'.format(value)
    else:
      formatted[key] = value
  return formatted


def read_inputs(signatures_path, counts_path, counts_column):
  with open(signatures_path, 'r') as signatures_fh:
    signature_matrix = mutational_signature.decompose.read_signature_matrix(signatures_fh)
  if counts_path is None:
    counts_data = mutational_signature.decompose.read_counts_vector(sys.stdin, signature_matrix['contexts'], counts_column=counts_column, strict=False)
  else:
    with open(counts_path, 'r') as counts_fh:
      counts_data = mutational_signature.decompose.read_counts_vector(counts_fh, signature_matrix['contexts'], counts_column=counts_column, strict=False)
  counts = counts_data['counts']
  validate_counts(counts)
  return signature_matrix, counts


def run(args, out=sys.stdout):
  if args.interval <= 0 or args.interval >= 1:
    raise ValueError('interval must lie strictly between 0 and 1')
  if args.replicates <= 0:
    raise ValueError('replicates must be greater than zero')
  if args.replicates < 200:
    logging.warning('using fewer than 200 replicates; intervals may be unstable')
  if args.alpha <= 0:
    raise ValueError('alpha must be positive')
  if args.detection_threshold_mutations < 0 or args.detection_threshold_proportion < 0:
    raise ValueError('detection thresholds must be non-negative')

  signature_matrix, counts = read_inputs(args.signatures, args.counts, args.count_column)
  if args.perturbation == 'downsample':
    perturb_counts(counts, args.perturbation, np.random.default_rng(args.seed), args.alpha, args.fraction, args.mutation_count)

  rng = np.random.default_rng(args.seed)
  point_result = fit_catalogue(signature_matrix, counts, args.mode, args.metric, args.solver, args.max_sigs, args.refine_minimum_error, args.context_cutoff)
  replicate_rows = []

  for replicate in range(1, args.replicates + 1):
    boot_counts = perturb_counts(counts, args.perturbation, rng, args.alpha, args.fraction, args.mutation_count)
    fit_result = fit_catalogue(signature_matrix, boot_counts, args.mode, args.metric, args.solver, args.max_sigs, args.refine_minimum_error, args.context_cutoff)
    replicate_rows.extend(make_replicate_rows(
      args.sample,
      replicate,
      signature_matrix['signature_names'],
      fit_result,
      int(np.sum(boot_counts)),
      args.detection_threshold_mutations,
      args.detection_threshold_proportion,
    ))

  summary_rows = summarise(
    args.sample,
    signature_matrix['signature_names'],
    point_result,
    replicate_rows,
    args.interval,
    args.perturbation,
    args.seed,
    args.mode,
    args.metric,
    args.solver,
    args.exceeds_threshold_mutations,
  )

  write_tsv(args.summary_output, summary_rows, out=out)
  if args.replicates_output is not None:
    write_tsv(args.replicates_output, replicate_rows)

  return summary_rows, replicate_rows


def build_parser():
  parser = argparse.ArgumentParser(description='Estimate mutational signature exposure stability intervals')
  parser.add_argument('--signatures', required=True, help='mutational signatures e.g. cosmic')
  parser.add_argument('--counts', help='counts file; stdin is used when omitted')
  parser.add_argument('--sample', required=True, help='sample identifier to include in output')
  parser.add_argument('--replicates', type=int, default=1000, help='number of bootstrap replicates')
  parser.add_argument('--perturbation', choices=('dirichlet-multinomial', 'multinomial', 'downsample'), default='dirichlet-multinomial', help='catalogue perturbation method')
  parser.add_argument('--alpha', type=float, default=0.5, help='Dirichlet concentration added to observed counts')
  parser.add_argument('--fraction', type=float, help='fraction of mutations to sample without replacement in downsample mode')
  parser.add_argument('--mutation-count', type=int, help='number of mutations to sample without replacement in downsample mode')
  parser.add_argument('--seed', type=int, help='random number seed for reproducibility')
  parser.add_argument('--metric', default='cosine', choices=('cosine', 'euclidean', 'kl', 'l1'), help='decomposition error metric')
  parser.add_argument('--solver', default='basin', choices=('basin', 'grid'), help='decomposition solver')
  parser.add_argument('--max-sigs', type=int, help='maximum number of signatures for compatible solvers')
  parser.add_argument('--context-cutoff', type=float, default=1e6, help='reserved for compatibility with decompose')
  parser.add_argument('--mode', choices=('direct', 'refine'), default='direct', help='fitting mode')
  parser.add_argument('--refine-minimum-error', type=float, default=0.01, help='maximum acceptable error increase when removing a signature')
  parser.add_argument('--interval', type=float, default=0.95, help='bootstrap percentile interval width')
  parser.add_argument('--detection-threshold-mutations', type=float, default=0, help='mutation threshold for calling a signature detected')
  parser.add_argument('--detection-threshold-proportion', type=float, default=0, help='proportion threshold for calling a signature detected')
  parser.add_argument('--exceeds-threshold-mutations', type=float, default=0, help='mutation threshold for probability_exceeds_threshold')
  parser.add_argument('--summary-output', help='write summary TSV to this path instead of stdout')
  parser.add_argument('--replicates-output', help='write replicate-level TSV to this path')
  parser.add_argument('--count-column', default='Count', help='specify column containing counts')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  return parser


if __name__ == '__main__':
  parser = build_parser()
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)
  try:
    run(args)
  except ValueError as exc:
    raise SystemExit(str(exc)) from exc
