#!/usr/bin/env python
'''
  deprecated compatibility wrapper for signature-stability
'''

import argparse
import io
import logging
import sys
from types import SimpleNamespace

import mutational_signature.signature_stability


def main():
  parser = argparse.ArgumentParser(description='Deprecated compatibility wrapper; prefer signature-stability --perturbation downsample')
  parser.add_argument('--definition', required=True, help='signatures to decompose to')
  parser.add_argument('--replicates', required=False, default=1000, type=int, help='number of stability replicates')
  parser.add_argument('--count_files', required=True, nargs='+', help='count files to process')
  parser.add_argument('--downsample_style', required=False, default='percent', help='percent or count')
  parser.add_argument('--downsample_points', nargs='+', type=float, default=[50], help='first downsample point is used')
  parser.add_argument('--aggregate', help='legacy option ignored')
  parser.add_argument('--start', type=int, default=0, help='first count file index to process')
  parser.add_argument('--count', type=int, default=None, help='number of count files to process')
  parser.add_argument('--seed', type=int, help='random number seed')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()

  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  logging.warning('bootstrap-simple is deprecated; use signature-stability --perturbation downsample')

  files = args.count_files[args.start:]
  if args.count is not None:
    files = files[:args.count]
  point = args.downsample_points[0]
  if args.downsample_style == 'percent':
    fraction = point / 100.0
    mutation_count = None
  else:
    fraction = None
    mutation_count = int(point)

  all_rows = []
  for count_file in files:
    new_args = SimpleNamespace(
      signatures=args.definition,
      counts=count_file,
      sample=count_file,
      replicates=args.replicates,
      perturbation='downsample',
      alpha=0.1,
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
    summary_rows, _ = mutational_signature.signature_stability.run(new_args, out=io.StringIO())
    all_rows.extend(summary_rows)

  mutational_signature.signature_stability.write_tsv(None, all_rows, out=sys.stdout)


if __name__ == '__main__':
  main()
