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

def make_distance(A, b, metric, with_composition=False):
  '''
    assess the current candidate - how close is Ax to b?
    A = signature definitions: mutation_types x signatures
    b = observed mutation frequencies: mutation_types x 1
    x = predicted exposures: signatures x 1
  '''

  def distance(x):
    return distance_and_composition(x)[0]

  def distance_and_composition(x):
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
    return (error, residual_term)

  def cosine(x):
    '''
      x is the current estimate of exposures
    '''
    estimate = np.dot(A, x)
    similarity = b.dot(estimate) / (np.linalg.norm(b) * np.linalg.norm(estimate))

    # contribution to error
    normalised_error = [abs(estimate[i] - b[i]) for i in range(len(b))]
    #normalised_error = normalised_error / np.linalg.norm(normalised_error)
    normalised_error = normalised_error / sum(normalised_error)

    return (-similarity, normalised_error)

  def l1(x):
    '''
      l1 difference
    '''
    estimate = np.dot(A, x / sum(x))

    # squared error on each term (return this for least_squares)
    residual_term = [abs(estimate[i] - b[i]) for i in range(len(b))]

    # total error
    error = sum(residual_term)
    return (error, residual_term)

  if with_composition:
    return distance_and_composition
  else:
    return distance

def context_difference(A, b, x):
  '''
    wrongness in context counts
    positive value means signature has overshot
  '''
  estimate = np.dot(A, x)
  return estimate - b

def grid_search(A, b, metric, max_sigs):
  if max_sigs is None:
    max_sigs = 3

  logging.info('grid_search with up to %i signatures', max_sigs)

  if max_sigs > 3:
    logging.warn('grid_search with more than 3 sigs can be slow')

  best = (None, 1e9)
  dist_fn = make_distance(A, b, metric)
  for sigs in range(1, max_sigs + 1):
    logging.info('finding best with sig length %i', sigs)
    for tried, subsig in enumerate(itertools.combinations(range(0, A.shape[1]), sigs)): # all combinations of the specified number of sigs
      for weights in itertools.product(range(1, len(subsig) + 1), repeat=len(subsig)):
        candidate = [0] * A.shape[1]
        for idx, sig in enumerate(subsig):
          candidate[sig] = weights[idx]
        #candidate = [x / sum(candidate) for x in candidate]
        distance = dist_fn(np.array(candidate))
        if distance < best[1]:
          logging.debug('new best %f: %s', distance, candidate)
          best = (candidate, distance)
      if tried % 1000 == 0:
        logging.debug('tested %i...', tried)
  return np.array(best[0])

def basin_hopping_solver(A, b, metric, max_sigs):

  if max_sigs is not None:
    logging.warn('Ignoring max_sigs parameter')

  x0 = np.array([random.random() for _ in range(0, A.shape[1])])
  x0 = x0 / sum(x0) # normalize
 
  # solve with least squared (must use euclidean distance)
  #result = scipy.optimize.least_squares(make_distance(A, b), x0, bounds=(0.0, np.inf)).x
 
  # solve with basinhopping
  bounds=[(0.0, np.inf)] * len(x0)
  minimizer_kwargs = dict(method="L-BFGS-B", bounds=bounds)
  result = scipy.optimize.basinhopping(make_distance(A, b, metric), x0, minimizer_kwargs=minimizer_kwargs, stepsize=5, T=5).x

  return result


def decompose(signatures, counts_fh, out, metric, seed, evaluate, solver, max_sigs, context_cutoff, error_contribution=False):
  logging.info('starting...')

  if seed is not None:
    logging.info('setting random seed to %i', seed)
    random.seed(seed)

  # make array of signatures (A = sigs x classes)
  names = []
  first = True
  exclude_map = collections.defaultdict(set)
  for line in open(signatures, 'r'):
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

    # list of signatures to potentially exclude
    row = [float(x) for x in fields[1:]]
    for row_id, row_val in enumerate(row):
      if row_val > context_cutoff:
        logging.debug('potentially excluding signature %i if %s is empty because contribution is %s', len(names) + 1, signature_classes[row_id], fields[1 + row_id])
        exclude_map[signature_classes[row_id]].add(len(names))

    names.append(fields[0])
    A = np.vstack([A, row]) # 30x96 -> each row is a set of context values

  A = np.transpose(A) # -> each row is a signature

  # make target (b = observed_classes)
  b = [0] * len(signature_classes) # i.e. 96

  # read the count percentages

  first = True
  total_count = 0
  excluded_signatures = set()

  for idx, line in enumerate(counts_fh):
    if first:
      first = False
      continue
    if line.startswith('#'):
      continue
    fields = line.strip('\n').split('\t')
    if fields[0] not in signature_classes:
      logging.info('context %s not in signature definitions', fields[0])
      continue
    signature_index = signature_classes.index(fields[0])
    b[signature_index] = float(fields[1]) # counts
    #b[signature_index] = float(fields[2]) # percentage
    total_count += int(fields[1])

    # check for excluded signatures
    if int(fields[1]) == 0:
      if len(exclude_map[fields[0]]) > 0:
        logging.debug('excluding %s because %s is not present and has proportion >%.2f', ' '.join([names[x] for x in exclude_map[fields[0]]]), fields[0], context_cutoff)
      [excluded_signatures.add(x) for x in exclude_map[fields[0]]]

  if len(excluded_signatures) > 0:
    logging.info('signatures to exclude: %s', ' '.join([names[x] for x in sorted(list(excluded_signatures))]))
    all_names = np.copy(names)
    names = np.delete(names, list(excluded_signatures), axis=0)
    A = np.delete(A, list(excluded_signatures), axis=1)
  else:
    all_names = names

  # find x for Ax = b, x > 0 x = exposure to signature
  b = np.array(b)

  if evaluate is None: # find a solution
    logging.info('finding signatures from %s...', signatures)
    if solver == 'basin':
      solver = basin_hopping_solver
    elif solver == 'grid':
      solver = grid_search

    result = solver(A, b, metric, max_sigs)


  else: # use provided solution
    logging.info('evaluating signatures from %s...', signatures)
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
  
  # expand result
  if len(all_names) > len(names):
    all_result = [0] * len(all_names)
    for i in range(len(names)):
      value = result[i]
      name = names[i]
      all_index = list(all_names).index(name)
      all_result[all_index] = value
  else:
    all_names = names
    all_result = result

  if out is not None:
    # sort by name
    for i in sorted(range(len(all_result)), key=lambda k: all_names[k]):
      out.write('{}\t{:.3f}\n'.format(all_names[i], all_result[i] / total))

    out.write('Mutations\t{}\n'.format(total_count))

  # compare reconstruction
  for m in ('euclidean', 'cosine', 'l1'):
    error = make_distance(A, b, m, with_composition=True)(result)
    logging.info('%s error:\t%.5f', m, error[0])
    if metric == m:
      if metric == 'cosine':
        total_error = 1.0 + error[0]
      else:
        total_error = error[0]
      if out is not None:
        out.write('Error\t{:.3f}\n'.format(total_error))
        if error_contribution:
          for context_name, error_contribution, difference in zip(signature_classes, error[1], context_difference(A, b, result)):
            # context \t error % \t number of mutations
            out.write('Error {}\t{:.3f}\t{:.1f}\n'.format(context_name, error_contribution, difference))
      target_error = (total_error, error[1])

  return {'signature_names': all_names, 'signature_values': all_result, 'total': total, 'error': target_error}

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='mutational signature finder')
  parser.add_argument('--signatures', required=True, help='mutational signatures e.g. cosmic')
  parser.add_argument('--counts', required=True, help='counts file')
  parser.add_argument('--metric', required=False, default='cosine', help='metric. cosine, euclidean, or l1')
  parser.add_argument('--solver', required=False, default='basin', help='solver. basin, grid')
  parser.add_argument('--max_sigs', required=False, type=int, help='maximum number of sigs to use')
  parser.add_argument('--context_cutoff', required=False, type=float, default=1e6, help='exclude signatures with contexts above this percent that are not represented in the sample') # deconstructSigs = 0.2
  parser.add_argument('--error_contribution', action='store_true', help='show contribution of each context to error')
  parser.add_argument('--seed', required=False, type=int, help='random number seed for reproducibility')
  parser.add_argument('--evaluate', required=False, help='evaluate a mutational profile instead of calculating')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  decompose(args.signatures, open(args.counts, 'r'), sys.stdout, args.metric, args.seed, args.evaluate, args.solver, args.max_sigs, args.context_cutoff, args.error_contribution)
