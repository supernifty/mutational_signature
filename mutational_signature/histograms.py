#!/usr/bin/env python
'''
  make a panel of histograms from signatures
'''

import argparse
import collections
import logging
import math
import matplotlib.pyplot as plt
import sys

def main(exposures, cols, target, dpi=300, width=None, height=None, exclude_below=None):
  logging.info('reading %i files...', len(exposures))
  ds = collections.defaultdict(list)
  skip = set()
  for fn in exposures:
    with open(fn, 'rt') as fh:
      for line in fh:
        name, value = line.strip().split('\t')
        value = float(value)
        if value > 1:
          skip.add(name)
        ds[name].append(value)

  for s in ds:
    if max(ds[s]) < exclude_below:
      skip.add(s)

  # remove skips
  for s in skip:
    ds.pop(s, None)
    logging.info('skipping %s', s)

        
  rows = math.ceil(len(ds) / cols)
  if width is None: 
    width = 3 * cols
  if height is None:
    height = 2 * rows

  logging.info('plotting %i signatures on a %ix%i grid with size %ix%i...', len(ds), cols, rows, width, height)
  fig, axs = plt.subplots(rows, cols, sharex=True, sharey=True)
  fig.set_size_inches(width, height, forward=True)
  for sig, ax in zip(ds, axs.flat):
    logging.info('plotting %s', sig)
    ax.hist(ds[sig])
    ax.set_title(sig)
    
  logging.info('saving to %s', target)
  plt.tight_layout()
  plt.savefig(target, dpi=dpi, transparent=False)

  logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Assess MSI')
  parser.add_argument('--exposures', required=True, nargs='+', help='calculated exposures')
  parser.add_argument('--cols', type=int, default=4, help='plots per row')
  parser.add_argument('--figure_width', type=int, required=False, help='width')
  parser.add_argument('--figure_height', type=int, required=False, help='height')
  parser.add_argument('--exclude_below', type=float, default=0, required=False, help='exclude max below')
  parser.add_argument('--target', default='plot.png', help='output')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.exposures, args.cols, args.target, width=args.figure_width, height=args.figure_height, exclude_below=args.exclude_below)
