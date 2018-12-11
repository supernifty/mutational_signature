#!/usr/bin/env python

import argparse
import collections
import itertools
import logging
import math
import random
import sys

import numpy as np
import scipy.optimize

def make_distance(A, b, metric):
  '''
    assess the current candidate - how close is Ax to b?
    A = signature definitions: mutation_types x signatures
    b = observed mutation frequencies: mutation_types x 1
    x = predicted exposures: signatures x 1
  '''

  def distance(x):
    if metric == 'euclidean':
      return euclidean(x)
    if metric == 'cosine':
      return cosine(x)
    if metric == 'l1':
      return l1(x)
    raise ValueError('Unknown metric {}'.format(metric))

  def euclidean(x):
    estimate = np.dot(A, x / sum(x))

    # squared error on each term (return this for least_squares)
    residual_term = [math.pow(estimate[i] - b[i], 2) for i in range(len(b))]

    # total error
    error = math.sqrt(sum(residual_term))
    return error

  def cosine(x):
    '''
      x is the current estimate of exposures
    '''
    estimate = np.dot(A, x)
    similarity = b.dot(estimate) / (np.linalg.norm(b) * np.linalg.norm(estimate))
    return -similarity

  def l1(x):
    '''
      l1 difference
    '''
    estimate = np.dot(A, x / sum(x))

    # squared error on each term (return this for least_squares)
    residual_term = [abs(estimate[i] - b[i]) for i in range(len(b))]

    # total error
    error = sum(residual_term)
    return error


  return distance

def grid_search(A, b, max_sigs=4):
  best = (None, 1e9)
  dist_fn = make_distance(A, b)
  for sigs in range(1, max_sigs + 1):
    logging.debug('sig length %i', sigs)
    for subsig in itertools.combinations(range(0, 30), sigs):
      for weights in itertools.product(range(1, len(subsig) + 1), repeat=len(subsig)):
        candidate = [0] * A.shape[1]
        for idx, sig in enumerate(subsig):
          candidate[sig] = weights[idx]
        #candidate = [x / sum(candidate) for x in candidate]
        distance = dist_fn(np.array(candidate))
        if distance < best[1]:
          logging.debug('new best %f: %s', distance, candidate)
          best = (candidate, distance)
  return best[0]

def decompose(signatures, counts, out, metric, seed, evaluate):
  logging.info('finding signatures {} in {}...'.format(signatures, counts))

  if seed is not None:
    logging.info('setting random seed to %i', seed)
    random.seed(seed)

  # make array of signatures (A = sigs x classes)
  names = []
  A = np.empty((0, 96), float)
  first = True
  for line in open(signatures, 'r'):
    fields = line.strip('\n').split('\t')
    if first:
      first = False
      if '>' in fields[1]:
        signature_classes = fields[1:]
      else:
        signature_classes = ['{}{}{}>{}'.format(f[0], f[1], f[3], f[2]) for f in fields[1:]]
      logging.debug('%i signature classes', len(signature_classes))
      continue
    names.append(fields[0])
    row = [float(x) for x in fields[1:]]
    A = np.vstack([A, row]) # 30x96

  A = np.transpose(A)

  
  # make target (b = observed_classes)
  b = [0] * len(signature_classes)


  first = True
  for idx, line in enumerate(open(counts, 'r')):
    if first:
      first = False
      continue
    fields = line.strip('\n').split('\t')
    signature_index = signature_classes.index(fields[0])
    #b.append(float(fields[1])) # count
    b[signature_index] = float(fields[2]) # percentage

  # find x for Ax = b, x > 0 x = exposure to signature
  b = np.array(b)

  if evaluate is None: # find a solution
    x0 = np.array([random.random() for _ in range(0, A.shape[1])])
    x0 = x0 / sum(x0) # normalize
  
    # solve with least squared (must use euclidean distance)
    #result = scipy.optimize.least_squares(make_distance(A, b), x0, bounds=(0.0, np.inf)).x
  
    # solve with basinhopping
    bounds=[(0.0, np.inf)] * len(x0)
    minimizer_kwargs = dict(method="L-BFGS-B", bounds=bounds)
    result = scipy.optimize.basinhopping(make_distance(A, b, metric), x0, minimizer_kwargs=minimizer_kwargs, stepsize=5, T=5).x

    # solve with grid
    #result = grid_search(A, b)

  else: # use provided solution
    result = [0] * len(names)
    for line in open(args.evaluate, 'r'):
      fields = line.strip('\n').split('\t')
      if fields[0] in names:
        result[names.index(fields[0])] = float(fields[1])
      else:
        logging.info('skipped %s', line.strip('\n'))
    result = np.array(result)

  # write signature exposure
  total = sum(result)

  #sorted_indices = sorted(range(len(result)), key=lambda k: result[k])
  #for i in reversed(sorted_indices):
  #  sys.stdout.write('{}\t{:.3f}\n'.format(names[i], result[i] / total))

  # sort by name
  for i in sorted(range(len(result)), key=lambda k: names[k]):
    sys.stdout.write('{}\t{:.3f}\n'.format(names[i], result[i] / total))

  # compare reconstruction
  for metric in ('euclidean', 'cosine', 'l1'):
    error = make_distance(A, b, metric)(result)
    logging.info('%s error:\t%.5f', metric, error)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='mutational signature finder')
  parser.add_argument('--signatures', required=True, help='signatures')
  parser.add_argument('--counts', required=True, help='counts')
  parser.add_argument('--metric', required=False, default='cosine', help='metric. cosine, euclidean, or l1')
  parser.add_argument('--seed', required=False, type=int, help='random number seed for reproducibility')
  parser.add_argument('--evaluate', required=False, help='evaluate a list of exposures')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  decompose(args.signatures, args.counts, sys.stdout, args.metric, args.seed, args.evaluate)