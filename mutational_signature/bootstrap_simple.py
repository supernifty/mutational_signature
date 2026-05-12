#!/usr/bin/env python
'''
  decompose signatures using a bootstrap to generate confidence intervals
'''

import argparse
import collections
import csv
import io
import logging
import random
import sys

import numpy

import mutational_signature.decompose

MAX_FRACTION=0.5
MIN_MUTS=5

def cosine(a, b):
  similarity = b.dot(a) / (numpy.linalg.norm(b) * numpy.linalg.norm(a))
  return similarity

def main(ofh, definition, count_files, replicates, downsample_style, downsample_points, aggregate, start, count):
  odr = csv.DictReader(open(definition, 'rt'), delimiter='\t')
  defs = odr.fieldnames[1:]
  sigs = []
  for r in odr:
    logging.info(r)
    sigs.append(r['Sig'])

  odw = csv.DictWriter(ofh, delimiter='\t', fieldnames=['sample', 'description', 'similarity', 'mutations', 'error'] + sigs)
  odw.writeheader()
  is_sbs = all([len(x) == 4 for x in defs])
  logging.info('processing %i files...', len(count_files))
  for i, f in enumerate(count_files):
    logging.info('%i: %s', i + 1, f)
    if i < start or i >= start + count:
      continue
    # read in contexts
    ctx = {}
    muts = 0
    for row in csv.DictReader(open(f, 'rt'), delimiter='\t'):
      if row['Variation'] in defs or is_sbs and len(row['Variation']) == 5 and row['Variation'][3] == '>':
        ctx[row['Variation']] = int(row['Count'])
        muts += int(row['Count'])
    if muts < 1:
      logging.info('skipping %i mutation tumour')
      continue
    logging.info('%i mutations total', muts)
    
    # first do reference
    dummy = io.StringIO()
    context_cutoff = 1e6
    counts_fh = ['Variation\tCount\n'] + ['{}\t{}\n'.format(k, ctx[k]) for k in ctx]
    point_estimate = mutational_signature.decompose.decompose(signatures_fh=open(definition, 'rt'), counts_fh=counts_fh, out=dummy, metric='cosine', seed=None, evaluate=None, solver='basin', max_sigs=None, context_cutoff=context_cutoff, error_contribution=False)
    point_error = point_estimate['error'][0]
    orow = {k[0]: '{:.6f}'.format(k[1] / point_estimate['total']) for k in zip(point_estimate['signature_names'], point_estimate['signature_values'])}
    orow.update({'sample': f, 'description': 'reference', 'mutations': muts, 'error': point_error, 'similarity': 'n/a'})
    odw.writerow(orow)

    ctxs = []
    for k in ctx:
      ctxs.extend([k] * ctx[k])
    for pnt in downsample_points:
      agg = []
      for rep in range(replicates):
        if downsample_style == 'percent':
          mutation_count = max([1, int(pnt * muts / 100)])
        else:
          mutation_count = int(pnt)
          if mutation_count >= muts * MAX_FRACTION or muts < MIN_MUTS: # skip if mutation count is too close to original
            continue
        logging.info('replicate %i: point pnt %f: %i mutations', rep, pnt, mutation_count)
        random.shuffle(ctxs)
        counts = collections.Counter(ctxs[:mutation_count])
        # shuffle the contexts and choose the first mutation count
        dummy = io.StringIO()
        context_cutoff = 1e6
        counts_fh = ['Variation\tCount\n'] + ['{}\t{}\n'.format(k, counts[k]) for k in counts]
        bootstrap = mutational_signature.decompose.decompose(signatures_fh=open(definition, 'rt'), counts_fh=counts_fh, out=dummy, metric='cosine', seed=None, evaluate=None, solver='basin', max_sigs=None, context_cutoff=context_cutoff, error_contribution=False)
        point_error = bootstrap['error'][0]
        orow = {k[0]: k[1] / bootstrap['total'] for k in zip(bootstrap['signature_names'], bootstrap['signature_values'])}
        similarity = cosine(point_estimate['signature_values'], bootstrap['signature_values'])
        if aggregate is None:
          orow.update({'sample': f, 'description': 'replicate {} downsample {}'.format(rep, pnt), 'mutations': mutation_count, 'error': point_error, 'similarity': similarity})
          odw.writerow(orow)
        else:
          orow.update({'mutations': mutation_count, 'error': point_error, 'similarity': similarity})
          agg.append(orow) 
      if aggregate == 'median' and len(agg) > 0:
        logging.info(agg)
        median_row = {k: '{:.6f}'.format(numpy.median([agg[i][k] for i in range(len(agg))])) for k in agg[0]}
        median_row.update({'sample': f, 'description': 'median downsample {}: {}'.format(pnt, ','.join(['{:.6f}'.format(agg[i]['similarity']) for i in range(len(agg))])), 'mutations': mutation_count})
        odw.writerow(median_row)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Bootstrap reads counts file from stdin')
  parser.add_argument('--definition', required=True, help='signatures to decompose to')
  parser.add_argument('--replicates', required=False, default=10, type=int, help='how many times to run the bootstrap')
  parser.add_argument('--count_files', required=True, nargs='+', help='how many times to run the bootstrap')
  parser.add_argument('--downsample_style', required=False, default='percent', help='how many times to run the bootstrap')
  parser.add_argument('--downsample_points', nargs='+', type=float, help='how many times to run the bootstrap')
  parser.add_argument('--aggregate', help='set to median optionally to merge replicates')
  parser.add_argument('--start', type=int, help='set to median optionally to merge replicates')
  parser.add_argument('--count', type=int, help='set to median optionally to merge replicates')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(sys.stdout, args.definition, args.count_files, args.replicates, args.downsample_style, args.downsample_points, args.aggregate, args.start, args.count)
