#!/usr/bin/env python
'''
  probability of signatures arising by chance
'''

import argparse
import collections
import csv
import io
import logging
import random
import sys

import numpy

import plotme.box

import mutational_signature.decompose

def main(signatures, signature_sum=None, variant_count=10, threshold=0.5, bootstrap_count=1000):
  logging.info('reading counts from stdin...')
  variants = 0
  counts = collections.defaultdict(int)
  all_contexts = set()
  probs = []
  total_prob = 0
  
  for row in csv.DictReader(sys.stdin, delimiter='\t'):
    context = row['Variation']
    count = int(row['Count'])
    if count > 0:
      all_contexts.add(context)
      counts[context] += count
      probs.append((context, total_prob))
      total_prob += count

  all_contexts_list = list(all_contexts)

  logging.info('generating point estimate from %i context types', len(counts))
  #def decompose(signatures, counts_fh, out, metric, seed, evaluate, solver, max_sigs, context_cutoff, error_contribution=False):
  #result = mutational_signature.decompose.decompose(signatures, out, None, 'cosine', None, None, 'basin', None, context_cutoff, False)
  dummy = io.StringIO()
  context_cutoff = 1e6
  counts_fh = ['Variation\tCount\n'] + ['{}\t{}\n'.format(k, counts[k]) for k in counts]
  point_estimate = mutational_signature.decompose.decompose(signatures=signatures, counts_fh=counts_fh, out=dummy, metric='cosine', seed=None, evaluate=None, solver='basin', max_sigs=None, context_cutoff=context_cutoff, error_contribution=False)

  point_error = point_estimate['error'][0]
  point_signatures = {k[0]: k[1] / point_estimate['total'] for k in zip(point_estimate['signature_names'], point_estimate['signature_values'])}
  if signature_sum is not None and len(signature_sum) > 0:
    point_signatures['+'.join(signature_sum)] = sum([point_signatures[name] for name in signature_sum])
  contribution = {k: '{:.3f}'.format(point_signatures[k]) for k in point_signatures.keys()}
  logging.info('point error is %f; signatures are %s...', point_error, point_signatures)

  out = csv.DictWriter(sys.stdout, ['Context', 'Error'] + list(point_signatures.keys()), delimiter='\t')
  out.writeheader()
  contribution['Context'] = 'Point'
  contribution['Error'] = '{:.3f}'.format(point_error)
  out.writerow(contribution) 

  # now bootstrap
  passed = 0
  for bc in range(bootstrap_count):
    logging.info('bootstrap %i of %i...', bc, bootstrap_count)
    # choose some variants
    contexts = collections.defaultdict(int)
    for v in range(variant_count):
      r = random.randint(0, total_prob)
      for cand in probs:
        if r < cand[1]:
          contexts[cand[0]] += 1
          break
    # now get sigs
    dummy = io.StringIO()
    context_cutoff = 1e6
    counts_fh = ['Variation\tCount\n'] + ['{}\t{}\n'.format(k, contexts[k]) for k in contexts]
    point_estimate = mutational_signature.decompose.decompose(signatures=signatures, counts_fh=counts_fh, out=dummy, metric='cosine', seed=None, evaluate=None, solver='basin', max_sigs=None, context_cutoff=context_cutoff, error_contribution=False)
    point_signatures = {k[0]: k[1] / point_estimate['total'] for k in zip(point_estimate['signature_names'], point_estimate['signature_values'])}
    if signature_sum is not None and len(signature_sum) > 0:
      value = sum([point_signatures[name] for name in signature_sum])
      point_signatures['+'.join(signature_sum)] =  value
      if  value > threshold:
        passed += 1
    contribution = {k: '{:.3f}'.format(point_signatures[k]) for k in point_signatures.keys()}
    contribution['Context'] = 'Bootstrap'
    contribution['Error'] = '{:.3f}'.format(point_estimate['error'][0])
    out.writerow(contribution) 
  logging.info('done. %i of %i passed (%.2f)', passed, bootstrap_count, passed / bootstrap_count * 100)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Bootstrap reads counts file from stdin')
  parser.add_argument('--signatures', required=True, help='signatures to decompose to')
  parser.add_argument('--signature_sum', required=False, nargs='+', help='Show the variation in a sum of signatures')
  parser.add_argument('--threshold', required=False, default=0.5, type=float, help='threshold of interest')
  parser.add_argument('--variant_count', required=False, default=10, type=int, help='number of variants to choose')
  parser.add_argument('--bootstrap_count', required=False, default=100, type=int, help='number of bootsrap to choose')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.signatures, signature_sum=args.signature_sum, variant_count=args.variant_count, threshold=args.threshold, bootstrap_count=args.bootstrap_count)

