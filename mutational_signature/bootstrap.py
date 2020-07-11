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

import plotme.box

import mutational_signature.decompose

def main(signatures, bootstraps, confidence, just_sbs=True, all_contexts_possible=True, out=sys.stdout, subsample=1.0, subsample_count=None, plot=None, plot_title=None, error_probability=0.01):
  confidence=[float(c) for c in confidence]
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
    if all_contexts_possible:
      count = max(1, int(row['Count'])) # every context is possible TODO but need a better calculation
    else:
      count = max(0, int(row['Count'])) 
    samples = samples + [context] * count
    counts[context] += count

  all_contexts_list = list(all_contexts)

  logging.info('generating point estimate from %i contexts', len(counts))
  #def decompose(signatures, counts_fh, out, metric, seed, evaluate, solver, max_sigs, context_cutoff, error_contribution=False):
  #result = mutational_signature.decompose.decompose(signatures, out, None, 'cosine', None, None, 'basin', None, context_cutoff, False)
  dummy = io.StringIO()
  context_cutoff = 1e6
  counts_fh = ['Variation\tCount\n'] + ['{}\t{}\n'.format(k, counts[k]) for k in counts]
  point_estimate = mutational_signature.decompose.decompose(signatures=signatures, counts_fh=counts_fh, out=dummy, metric='cosine', seed=None, evaluate=None, solver='basin', max_sigs=None, context_cutoff=context_cutoff, error_contribution=False)

  point_error = point_estimate['error'][0]
  point_signatures = {k[0]: k[1] / point_estimate['total'] for k in zip(point_estimate['signature_names'], point_estimate['signature_values'])}
  logging.info('point error is %f; signatures are %s...', point_error, point_signatures)
  
  logging.info('bootstrapping from %i variants...', len(samples))

  bootstrap_error = []
  bootstrap_signatures = collections.defaultdict(list)
  
  for i in range(bootstraps):
    # make bootstrap
    counts = collections.defaultdict(int)
    if subsample_count is not None:
      logging.info('selecting %i contexts...', subsample_count)
      for _ in range(subsample_count):
        if random.random() < error_probability:
          # take any random context at all
          counts[all_contexts_list[random.randrange(len(all_contexts_list))]] += 1 
        else:
          counts[samples[random.randrange(len(samples))]] += 1
      subsampled = subsample_count
    else:
      subsampled = int(subsample * len(samples))
      # this way adds variability to the number of mutations
      logging.info('selecting approximately %i contexts...', subsampled)
      subsampled = 0
      for sample in samples: # every mutation
        if random.random() < subsample: # we'll take it
          if random.random() < error_probability:
            # take any random context at all
            counts[all_contexts_list[random.randrange(len(all_contexts_list))]] += 1
          else:
            # take an existing context
            counts[samples[random.randrange(len(samples))]] += 1
          subsampled += 1

    counts_fh = ['Variation\tCount\n'] + ['{}\t{}\n'.format(k, counts[k]) for k in counts]
    logging.info('selected %i contexts...%s', subsampled, counts_fh)
    bootstrap_estimate = mutational_signature.decompose.decompose(signatures=signatures, counts_fh=counts_fh, out=dummy, metric='cosine', seed=None, evaluate=None, solver='basin', max_sigs=None, context_cutoff=context_cutoff, error_contribution=False)
    for name, value in zip(bootstrap_estimate['signature_names'], bootstrap_estimate['signature_values']):
      bootstrap_signatures[name].append(value / bootstrap_estimate['total'])
    bootstrap_error.append(bootstrap_estimate['error'][0])
    
    logging.info('bootstrap %i: done', i + 1)

  # point is the value of the non-bootstrapped data
  out.write('Name\tPoint\tMean\tSD\t{}\tMin\tMax\tValues\n'.format('\t'.join([str(c) for c in confidence])))
  values = [numpy.percentile(bootstrap_error, c * 100) for c in confidence]
  out.write('Error\t{:.3f}\t{:.3f}\t{:.3f}\t{}\t{:.3f}\t{:.3f}\t{}\n'.format(point_error, numpy.mean(bootstrap_error), numpy.std(bootstrap_error), '\t'.join(['{:.3f}'.format(v) for v in values]), min(bootstrap_error), max(bootstrap_error), ','.join(['{:.2f}'.format(x) for x in bootstrap_error])))
  for sig in sorted(bootstrap_signatures.keys()):
    values = [numpy.percentile(bootstrap_signatures[sig], c * 100) for c in confidence]
    out.write('{}\t{:.3f}\t{:.3f}\t{:.3f}\t{}\t{:.3f}\t{:.3f}\t{}\n'.format(sig, point_signatures[sig], numpy.mean(bootstrap_signatures[sig]), numpy.std(bootstrap_signatures[sig]), '\t'.join(['{:.3f}'.format(v) for v in values]), min(bootstrap_signatures[sig]), max(bootstrap_signatures[sig]), ','.join(['{:.2f}'.format(x) for x in bootstrap_signatures[sig]])))

  if plot is not None:
    logging.info('writing boxplots to %s...', plot)
    fh = io.StringIO()
    fh.write('x\ty\tz\n')
    # write original values
    for sig in sorted(bootstrap_signatures.keys()):
      for value in bootstrap_signatures[sig]:
        fh.write('{}\t{}\t{}\n'.format(sig, 'Pct', value))

    for value in bootstrap_error:
      fh.write('{}\t{}\t{}\n'.format('Error', 'Pct', value))

    fh.seek(0)
    plotme.box.plot_box(fh, plot, 'x', 'y', 'z', plot_title or 'Variability in bootstrapped signatures', x_label=None, y_label=None, x_order=None, y_order=None, fig_width=16, fig_height=8, fontsize=8, significance=None, significance_nobar=False, separator=False, include_zero=True, x_label_rotation='vertical', annotate=None, include_other=1)
    logging.info('writing boxplots to %s: done', plot)
  
  logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Bootstrap reads counts file from stdin')
  parser.add_argument('--signatures', required=True, help='signatures to decompose to')
  parser.add_argument('--count', required=False, default=10, type=int, help='how many times to run the bootstrap')
  parser.add_argument('--confidence', required=False, nargs='+', default=[0.05, 0.5, 0.95], help='what confidence intervals to include in output')
  parser.add_argument('--subsample', required=False, type=float, default=1.0, help='reduce context count by this percentage')
  parser.add_argument('--subsample_count', required=False, type=int, help='reduce context to this number of variants')
  parser.add_argument('--plot', required=False, help='filename for boxplot')
  parser.add_argument('--plot_title', required=False, help='title for boxplot')
  parser.add_argument('--all_contexts_possible', action='store_true', help='set empty contexts to have a count of 1')
  parser.add_argument('--error_probability', required=False, default=0.01, help='probability of a draw being a random context')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.signatures, args.count, args.confidence, subsample=args.subsample, subsample_count=args.subsample_count, all_contexts_possible=args.all_contexts_possible, plot=args.plot, plot_title=args.plot_title, error_probability=args.error_probability)

