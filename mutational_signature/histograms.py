#!/usr/bin/env python
'''
  make a panel of histograms from signatures
'''

import argparse
import collections
import csv
import logging
import math
import matplotlib.pyplot as plt
import sys

import numpy as np

COLORS=['red', 'violet']

def main(exposures, combined, cols, target, dpi=300, width=None, height=None, exclude_below=None, include=None, exclude=None, percentiles=None):
  logging.info('reading %i files...', len(exposures))
  ds = collections.defaultdict(list)
  skip = set()
  for fn in exposures:
    with open(fn, 'rt') as fh:
      if combined: # Sample Sigs...
        for r in csv.DictReader(fh, delimiter='\t'):
          for name in r:
            if name == 'Sample':
              continue
            if include is not None and name not in include:
              skip.add(name)
              continue
            if exclude is not None and name in exclude:
              skip.add(name)
              continue
            value = float(r[name])
            if value > 1:
              skip.add(name)
            ds[name].append(value)
          
      else:
        for line in fh:
          name, value = line.strip().split('\t')
          value = float(value)
          if include is not None and name not in include:
            skip.add(name)
          if exclude is not None and name in exclude:
            skip.add(name)
          if value > 1:
            skip.add(name)
          if math.isnan(value):
            continue
          ds[name].append(value)

  for s in ds:
    if len(ds[s]) == 0:
      skip.add(s)
    if max(ds[s]) < exclude_below:
      skip.add(s)

  # remove skips
  for s in skip:
    ds.pop(s, None)
    logging.info('skipping %s', s)

  plot(ds, cols, target, dpi, width, height, percentiles)

def plot(ds, cols, target, dpi=300, width=None, height=None, percentiles=None):
        
  rows = math.ceil(len(ds) / cols)
  if width is None: 
    width = 3 * cols
  if height is None:
    height = 2 * rows

  logging.info('plotting %i signatures on a %ix%i grid with size %ix%i...', len(ds), cols, rows, width, height)
  fig, axs = plt.subplots(rows, cols, sharex=True, sharey=True)
  fig.set_size_inches(width, height, forward=True)
  plt.xlabel('Tumours')
  plt.ylabel('Signature Proportion')
  for sig, ax in zip(ds, axs.flat):
    logging.info('plotting %s', sig)
    ax.hist(ds[sig])
    ax.set_title(sig)
    if percentiles is not None:
      for idx, p in enumerate(percentiles):
        position = np.percentile(ds[sig], p)
        logging.debug('percentile %i is at %.2f', p, position)
        ax.axvline(position, color=COLORS[idx % len(COLORS)], label='{} percentile'.format(p)) 
        ax.text(position + 0.01, 1, '{}'.format(int(position * 100)), color=COLORS[idx % len(COLORS)], va='bottom')
  logging.info('saving to %s', target)
  plt.tight_layout()
  plt.savefig(target, dpi=dpi, transparent=False)

  logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Generate ?')
  parser.add_argument('--exposures', required=True, nargs='+', help='calculated exposures')
  parser.add_argument('--combined', action='store_true', help='input is a combined sig file')
  parser.add_argument('--cols', type=int, default=4, help='plots per row')
  parser.add_argument('--percentiles', nargs='*', type=int, default=[95], required=False, help='percentiles to show')
  parser.add_argument('--figure_width', type=int, required=False, help='width')
  parser.add_argument('--figure_height', type=int, required=False, help='height')
  parser.add_argument('--exclude_below', type=float, default=0, required=False, help='exclude max below')
  parser.add_argument('--include', nargs='+', required=False, help='explicit list of sigs to include')
  parser.add_argument('--exclude', nargs='+', required=False, help='explicit list of sigs to exclude')
  parser.add_argument('--target', default='plot.png', help='output')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.exposures, args.combined, args.cols, args.target, width=args.figure_width, height=args.figure_height, exclude_below=args.exclude_below, include=args.include, exclude=args.exclude, percentiles=args.percentiles)
