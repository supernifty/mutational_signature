#!/usr/bin/env python
'''
  signature decomposition based on log likelihood (suitable for low mutation count)
'''

import argparse
import collections
import itertools
import logging
import math
import random
import sys

import numpy as np
import scipy.optimize
import scipy.special

def log_likelihood(counts, probs, eps=1e-12):
    """
    counts: array of observed counts (length K)
    probs: array of probabilities for a candidate process (length K)
    """
    counts = np.asarray(counts)
    probs = np.asarray(probs)

    # Numerical safety
    probs = np.clip(probs, eps, 1.0)

    return np.sum(counts * np.log(probs))

def solve(A, b):
  '''
    calculate log likelihood for each signature
    A = ctx x sig
    b = counts for each context
  '''
  rs = []
  for probs in A:
    r = log_likelihood(b, probs) # this is a single value for this signature
    rs.append(r)
  return rs

def decompose(sigs_fh, counts_fh, ofh, seed, counts_column):
  logging.debug('starting...')

  if seed is not None:
    logging.info('setting random seed to %i', seed)
    random.seed(seed)

  # read sig defs
  first = True
  names = []
  for line in sigs_fh:
    fields = line.strip('\n\r').split('\t')
    if first:
      first = False
      if '>' in fields[1] or len(fields[1]) != 4:
        signature_classes = fields[1:]
      else: # convert from abcd -> abd>c
        signature_classes = ['{}{}{}>{}'.format(f[0], f[1], f[3], f[2]) for f in fields[1:]]
      logging.debug('%i signature classes: %s...', len(signature_classes), signature_classes[:3])
      A = np.empty((0, len(signature_classes)), float)
      continue

    row = [float(x) for x in fields[1:]]
    names.append(fields[0])
    A = np.vstack([A, row]) # 30x96 -> each row is a set of context values

  if first: # no records
    logging.fatal('no signature definitions found')
  
  #A = np.transpose(A) # -> each col is a signature
  b = [0] * len(signature_classes) # i.e. 96
  logging.info('%i sigs and %i ctxs', len(A), len(A[0]))

  # read counts
  first = True
  total_count = 0
  for idx, line in enumerate(counts_fh):
    if first:
      first = False
      header = line.strip('\n').split('\t')
      continue
    if line.startswith('#'):
      continue
    fields = line.strip('\n').split('\t')
    if fields[0] not in signature_classes:
      logging.debug('context %s not in signature definitions', fields[0])
      continue
    signature_index = signature_classes.index(fields[0])
    b[signature_index] = float(fields[header.index(counts_column)]) # counts
    #b[signature_index] = float(fields[2]) # percentage
    total_count += float(fields[header.index(counts_column)])
  b = np.array(b)
  logL = solve(A, b)
  # convert back to prob
  prior = np.array([1 / len(A)] * len(A))
  logging.debug('prior: %s', prior)
  logPrior = np.log(prior)
  logPosterior = logL + logPrior
  logPosterior -= logPosterior.max()
  posterior = np.exp(logPosterior)
  posterior /= posterior.sum()

  if ofh is not None:
    # sort by name
    for i in sorted(range(len(posterior)), key=lambda k: names[k]):
      ofh.write('{}\t{:.6f}\n'.format(names[i], posterior[i]))

    ofh.write('Mutations\t{}\n'.format(total_count))

  return posterior

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='mutational signature finder')
  parser.add_argument('--signatures', required=True, help='mutational signatures e.g. cosmic')
  parser.add_argument('--counts', required=True, help='counts file')
  parser.add_argument('--seed', required=False, type=int, help='random number seed for reproducibility')
  parser.add_argument('--count_column', required=False, default='Count', help='specify column containing counts')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  decompose(open(args.signatures, 'r'), open(args.counts, 'r'), sys.stdout, args.seed, args.count_column)
