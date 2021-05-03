#!/usr/bin/env python
'''
  generate the impact of the loss of 1 of each seen context
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

def main(signatures, signature_sum=None):
  logging.info('reading counts from stdin...')
  variants = 0
  counts = collections.defaultdict(int)
  all_contexts = set()
  for row in csv.DictReader(sys.stdin, delimiter='\t'):
    context = row['Variation']
    count = int(row['Count'])
    if count > 0:
      all_contexts.add(context)
      counts[context] += count

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

  out = csv.DictWriter(sys.stdout, ['Context'] + list(point_signatures.keys()), delimiter='\t')
  out.writeheader()
  contribution['Context'] = 'Point'
  out.writerow(contribution) 
  
  for i in all_contexts_list:
    logging.debug('context %s with count %i...', i, counts[i])
    counts[i] -= 1

    dummy = io.StringIO()
    context_cutoff = 1e6
    counts_fh = ['Variation\tCount\n'] + ['{}\t{}\n'.format(k, counts[k]) for k in counts]
    estimate = mutational_signature.decompose.decompose(signatures=signatures, counts_fh=counts_fh, out=dummy, metric='cosine', seed=None, evaluate=None, solver='basin', max_sigs=None, context_cutoff=context_cutoff, error_contribution=False)
    new_signatures = {k[0]: k[1] / estimate['total'] for k in zip(estimate['signature_names'], estimate['signature_values'])}
    if signature_sum is not None and len(signature_sum) > 0:
      new_signatures['+'.join(signature_sum)] = sum([new_signatures[name] for name in signature_sum])
    # difference
    contribution = {k: '{:.3f}'.format(point_signatures[k] - new_signatures[k]) for k in new_signatures.keys()}
    contribution['Context'] = i
    out.writerow(contribution) 
    
    counts[i] += 1
  logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Bootstrap reads counts file from stdin')
  parser.add_argument('--signatures', required=True, help='signatures to decompose to')
  parser.add_argument('--signature_sum', required=False, nargs='+', help='Show the variation in a sum of signatures')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.signatures, signature_sum=args.signature_sum)

