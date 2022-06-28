#!/usr/bin/env python

import argparse
import collections
import csv
import io
import itertools
import logging
import math
import random
import sys

import numpy as np
import scipy.optimize

import mutational_signature.decompose

def decompose_bulk(out, signatures_fh, counts, metric, seed, evaluate, solver, max_sigs, context_cutoff, error_contribution=False, strand=False, counts_column='Count'):
  results = []
  seen = set()
  for i, c in enumerate(counts):
    result = io.StringIO()
    row = {'Filename': c}
    logging.info('decomposing %i of %i: %s...', i+1, len(counts), c)
    mutational_signature.decompose.decompose(open(args.signatures, 'r'), open(c, 'r'), result, args.metric, args.seed, args.evaluate, args.solver, args.max_sigs, args.context_cutoff, args.error_contribution, args.strand, args.count_column)
    for line in result.getvalue().split('\n'):
      #logging.info('processing %s', line)
      if '\t' in line:
        k, v = line.strip().split('\t')
        row[k] = v
        seen.add(k)
      else:
        pass #logging.debug('skipping %s', line)
    results.append(row)

  # now write results
  logging.info('writing %i results...', len(results))
  o = csv.DictWriter(out, delimiter='\t', fieldnames=['Filename'] + sorted(list(seen)))
  o.writeheader()
  for r in results:
    o.writerow(r)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='mutational signature finder')
  parser.add_argument('--signatures', required=True, help='mutational signatures e.g. cosmic')
  parser.add_argument('--counts', required=True, nargs='+', help='counts files')
  parser.add_argument('--metric', required=False, default='cosine', help='metric. cosine, euclidean, or l1') # todo: add hellinger
  parser.add_argument('--solver', required=False, default='basin', help='solver. basin, grid')
  parser.add_argument('--max_sigs', required=False, type=int, help='maximum number of sigs to use')
  parser.add_argument('--context_cutoff', required=False, type=float, default=1e6, help='exclude signatures with contexts above this percent that are not represented in the sample') # deconstructSigs = 0.2
  parser.add_argument('--error_contribution', action='store_true', help='show contribution of each context to error')
  parser.add_argument('--seed', required=False, type=int, help='random number seed for reproducibility')
  parser.add_argument('--evaluate', required=False, help='evaluate a mutational profile instead of calculating')
  parser.add_argument('--strand', action='store_true', help='strand based signatures')
  parser.add_argument('--count_column', required=False, default='Count', help='specify column containing counts')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  decompose_bulk(sys.stdout, open(args.signatures, 'r'), args.counts, args.metric, args.seed, args.evaluate, args.solver, args.max_sigs, args.context_cutoff, args.error_contribution, args.strand, args.count_column)
