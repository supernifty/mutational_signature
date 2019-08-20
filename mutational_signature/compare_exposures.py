#!/usr/bin/env python
'''
  how similar are the input signatures?
'''

import argparse
import collections
import logging
import sys

import numpy as np

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt

def cosine(x, y):
  similarity = x.dot(y) / (np.linalg.norm(x) * np.linalg.norm(y))
  return similarity

def compare(x_exposures, y_exposures, out, plot, x_label, y_label, title, log):

  exclude = ('Error', 'Mutations')
  seen = set()

  logging.info('processing %s...', x_exposures)
  x = collections.defaultdict(float)
  for line in open(x_exposures, 'r'):
    k, v = line.strip('\n').split('\t')
    if k not in exclude:
      x[k] = float(v)
      seen.add(k)

  logging.info('processing %s...', y_exposures)
  y = collections.defaultdict(float)
  for line in open(y_exposures, 'r'):
    k, v = line.strip('\n').split('\t')
    if k not in exclude:
      y[k] = float(v)
      seen.add(k)

  logging.info('%i signatures seen', len(seen))
  
  # calculate cosine
  ls = sorted(seen)
  xs = [x[k] + 0.001 for k in ls]
  ys = [y[k] + 0.001 for k in ls]
  similarity = cosine(np.array(xs), np.array(ys))
  out.write('{}\t{}\t{}\n'.format(x_label, y_label, similarity))

  # scatter
  if plot is not None:
    import matplotlib.style
    matplotlib.style.use('seaborn')
    fig, ax = plt.subplots()
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    if log:
      ax.set_yscale('log')
      ax.set_xscale('log')

    min_axis = min(ax.get_xlim()[1], ax.get_ylim()[1])
    ax.plot([0.001, min_axis], [0.001, min_axis], c=".3", alpha=0.5)

    ax.scatter(xs, ys)
    ax.text(0.0, 0.95, 'Correlation {}%'.format(int(similarity * 100)), transform=ax.transAxes)
    if title is not None:
      ax.set_title(title)
    for i, txt in enumerate(ls):
      ax.annotate(txt, (xs[i], ys[i]))
    plt.tight_layout()
    plt.savefig(plot)
    matplotlib.pyplot.close('all')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='compare signatures')
  parser.add_argument('--x_exposures', required=True, help='exposures')
  parser.add_argument('--y_exposures', required=True, help='second group of exposures')
  parser.add_argument('--x_label', required=True, help='exposures')
  parser.add_argument('--y_label', required=True, help='second group of exposures')
  parser.add_argument('--title', required=False, help='title')
  parser.add_argument('--plot', required=False, help='correlation plot')
  parser.add_argument('--log', action='store_true', help='log scale')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)
  compare(args.x_exposures, args.y_exposures, sys.stdout, args.plot, args.x_label, args.y_label, args.title, args.log)

