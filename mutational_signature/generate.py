#!/usr/bin/env python
'''
  generate signatures from count files
'''

import argparse
import collections
import csv
import logging
import sys

import numpy as np
import sklearn.decomposition

import count_maf

def main(counts, components, use_proportion, stats):
  logging.info('processing %i count files to generate %s components...', len(counts), components)

  contexts = set()
  samples = collections.defaultdict(dict)
  for count_file in counts:
    sample = count_file.split('/')[-1].split('.')[0]
    for row in csv.DictReader(open(count_file, 'r'), delimiter='\t'):
      # expct Variation Count Probability
      contexts.add(row['Variation'])
      if use_proportion:
        samples[sample][row['Variation']] = float(row['Probability'])
      else:
        samples[sample][row['Variation']] = int(row['Count'])
  logging.info('found %i contexts in %i samples', len(contexts), len(counts))

  X = []
  contexts_list = sorted(list(contexts))
  for sample in sorted(samples.keys()): # each relevant sample
    features = [samples[sample].get(context, 0) for context in contexts_list]
    X.append(features)

  X = np.array(X)
  logging.info('decomposing %i samples to %s components...', len(X), components)

  if stats:
    stats_fh = open(stats, 'w')
    stats_fh.write('Components\tFrobenius\tRSS\tExplainedVariance\n')
    
  #for components in range(3,16):
  for component in components:
    model = sklearn.decomposition.NMF(n_components=component, init="nndsvd")
    W = model.fit_transform(X)
    H = model.components_

    # stats
    estimate = np.dot(W, H)
    difference = X - estimate
    rss = (difference ** 2).sum() # frobenius is just the sqrt of this
    explained_variance = sum(np.var(W, axis=0)) # components
    logging.info('signature count %i with %i features and %i samples gives frobenius error: %.2f rss: %.2f explained_variance: %.2f', component, len(contexts), len(X), model.reconstruction_err_, rss, explained_variance)
    if stats:
      stats_fh.write('{}\t{:.1f}\t{:.1f}\t{:.1f}\n'.format(component, model.reconstruction_err_, rss, explained_variance))

    logging.info('writing signatures...')
    sys.stdout.write('Name\t{}\n'.format('\t'.join(contexts_list)))
    for sig_num, sig in enumerate(H):
      total = sum(sig)
      sys.stdout.write('Signature.{}\t{}\n'.format(sig_num + 1, '\t'.join(['{:.6f}'.format(x / total) for x in sig])))

  logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Find mutational signatures')
  parser.add_argument('--counts', required=True, nargs='+', help='context counts for each sample')
  parser.add_argument('--components', required=True, nargs='+', type=int, help='number of decomposition components')
  parser.add_argument('--use_proportion', action='store_true', help='use variation proportion instead of count')
  parser.add_argument('--stats', required=False, help='write stats to file')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.counts, args.components, args.use_proportion, args.stats)
