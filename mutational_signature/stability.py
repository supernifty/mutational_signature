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
  parser.add_argument('--sample', default='legacy', help='sample identifier for the new stability output')
  parser.add_argument('--replicates', type=int, default=1000, help='number of stability replicates')
  parser.add_argument('--seed', type=int, help='random number seed')
  parser.add_argument('--alpha', type=float, default=0.1, help='Dirichlet concentration')
  parser.add_argument('--signature_sum', required=False, nargs='+', help='legacy option ignored')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()

  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  logging.warning('stability is deprecated; use signature-stability for rigorous exposure stability intervals')
  new_args = SimpleNamespace(
    signatures=args.signatures,
    counts=None,
    sample=args.sample,
    replicates=args.replicates,
    perturbation='dirichlet-multinomial',
    alpha=args.alpha,
    fraction=None,
    mutation_count=None,
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
