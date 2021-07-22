#!/usr/bin/env python
'''
  improve signatures by removing redundant signatures
'''

import argparse
import csv
import collections
import io
import logging
import sys

import mutational_signature
import mutational_signature.decompose

def remove_sig(sig, lines):
  result = []
  for line in lines:
    if line.split('\t')[0] != sig:
      result.append(line)
    else:
      pass #logging.debug('removed %s', sig)
  return result

def main(signatures, minimum_error=0.01, just_sbs=True, out=sys.stdout):
  logging.info('reading counts from stdin...')

  variants = 0
  samples = []
  counts = collections.defaultdict(int)
  all_contexts = set()
  for row in csv.DictReader(sys.stdin, delimiter='\t'):
    context = row['Variation']
    all_contexts.add(context)
    if just_sbs and len(context) != 5 and context[3] != '>':
      continue
    count = int(row['Count'])
    counts[context] += count

  counts_fh = ['Variation\tCount\n'] + ['{}\t{}\n'.format(k, counts[k]) for k in counts]

  signatures_fh = open(signatures).readlines()
  logging.debug('refining %i signatures...', len(signatures_fh) - 1)
  dummy = io.StringIO()
  estimate = mutational_signature.decompose.decompose(signatures_fh=signatures_fh, counts_fh=counts_fh, out=dummy, metric='cosine', seed=None, evaluate=None, solver='basin', max_sigs=None, context_cutoff=1e6, error_contribution=False)
  error = estimate['error'][0]
  sigs = {k[0]: k[1] / estimate['total'] for k in zip(estimate['signature_names'], estimate['signature_values'])}

  # remove empty
  logging.debug('removing empty sigs...')
  for sig in sigs:
    if sigs[sig] < 0.001:
      signatures_fh = remove_sig(sig, signatures_fh)

  logging.debug('removing unnecessary sigs...')
  refining = True
  while refining:
    refining = False
    # remove unhelpful if no empties removed
    for sig in sigs:
      dummy = io.StringIO()
      logging.debug('refining %i signatures...', len(signatures_fh) - 1)
      new_signatures_fh = remove_sig(sig, signatures_fh)
      new_estimate = mutational_signature.decompose.decompose(signatures_fh=new_signatures_fh, counts_fh=counts_fh, out=dummy, metric='cosine', seed=None, evaluate=None, solver='basin', max_sigs=None, context_cutoff=1e6, error_contribution=False)
      new_error = new_estimate['error'][0]
      sigs = {k[0]: k[1] / new_estimate['total'] for k in zip(new_estimate['signature_names'], new_estimate['signature_values'])}
      if new_error < error + minimum_error: # error didn't increase
        logging.info('removing %s: error increase was only %.4f -> %.4f', sig, error, new_error)
        signatures_fh = new_signatures_fh
        error = new_error
        refining = True
      else:
        logging.info('keeping %s: error increase was %.4f -> %.4f', sig, error, new_error)

  # final result
  sys.stdout.write(dummy.getvalue())
  
  logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Refine')
  parser.add_argument('--signatures', required=True, help='signatures to decompose to')
  parser.add_argument('--minimum_error', required=False, type=float, default=0.01, help='minimum improvement in cosine similarity required')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.signatures, args.minimum_error)
