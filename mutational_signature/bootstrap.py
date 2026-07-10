#!/usr/bin/env python
'''
  deprecated compatibility wrapper for signature-stability
'''

import argparse
import logging
import sys
from types import SimpleNamespace

import mutational_signature.signature_stability


def main():
  parser = argparse.ArgumentParser(description='Deprecated compatibility wrapper; prefer signature-stability')
  parser.add_argument('--signatures', required=True, help='signatures to decompose to')
  parser.add_argument('--count', required=False, default=1000, type=int, help='number of stability replicates')
  parser.add_argument('--sample', default='legacy', help='sample identifier for the new stability output')
  parser.add_argument('--seed', type=int, help='random number seed')
  parser.add_argument('--alpha', type=float, default=0.1, help='Dirichlet concentration for signature-stability')
  parser.add_argument('--confidence', required=False, nargs='+', default=None, help='legacy option ignored; use --interval with signature-stability')
  parser.add_argument('--subsample', required=False, type=float, default=1.0, help='legacy option mapped to downsample fraction when below 1')
  parser.add_argument('--subsample_count', required=False, type=int, help='legacy option mapped to downsample mutation count')
  parser.add_argument('--plot', required=False, help='legacy plotting is not supported; write --replicates-output with signature-stability')
  parser.add_argument('--plot_title', required=False, help='legacy plotting is not supported')
  parser.add_argument('--all_contexts_possible', action='store_true', help='legacy option ignored')
  parser.add_argument('--error_probability', required=False, default=0.01, type=float, help='legacy option ignored')
  parser.add_argument('--signature_sum', required=False, nargs='+', help='legacy option ignored')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()

  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  logging.warning('bootstrap is deprecated; use signature-stability for rigorous exposure stability intervals')
  if args.plot is not None:
    logging.warning('legacy bootstrap plotting is not supported by this compatibility wrapper')

  perturbation = 'dirichlet-multinomial'
  fraction = None
  mutation_count = None
  if args.subsample_count is not None:
    perturbation = 'downsample'
    mutation_count = args.subsample_count
  elif args.subsample != 1.0:
    perturbation = 'downsample'
    fraction = args.subsample

  new_args = SimpleNamespace(
    signatures=args.signatures,
    counts=None,
    sample=args.sample,
    replicates=args.count,
    perturbation=perturbation,
    alpha=args.alpha,
    fraction=fraction,
    mutation_count=mutation_count,
    seed=args.seed,
    metric='cosine',
    solver='basin',
    max_sigs=None,
    context_cutoff=1e6,
    mode='direct',
    refine_minimum_error=0.01,
    interval=0.95,
    detection_threshold_mutations=0,
    detection_threshold_proportion=0,
    exceeds_threshold_mutations=0,
    summary_output=None,
    replicates_output=None,
    count_column='Count',
  )
  mutational_signature.signature_stability.run(new_args, out=sys.stdout)


if __name__ == '__main__':
  main()
