#!/usr/bin/env python
'''
  plots the 96 snv contexts given a count input
'''

import logging
import sys

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt

from pylab import rcParams
rcParams['figure.figsize'] = 16, 10

def plot(counts, target):
  first = True
  vals = {}
  for line in counts:
    if first:
      first = False
      continue
    if line.count('\t') > 2:
      v, count, probability, _ = line.strip('\n').split('\t', 3)
    elif line.count('\t') > 1:
      v, count, probability = line.strip('\n').split('\t')
    else:
      v, probability = line.strip('\n').split('\t')
    if len(v) == 5 and v[3] == '>':
      variation = '{}>{} {}'.format(v[1], v[4], v[0:3])
      vals[variation] = float(probability)
    else:
      sys.stderr.write('skipped {}\n'.format(v))

  sys.stderr.write('{} contexts\n'.format(len(vals)))
  plot_signature(vals, target)

def plot_signature(vals, target):
    xs = sorted(vals.keys())
    ys = list([vals[x] for x in xs])

    color = ((0.2,0.7,0.9),)*16 + ((0.1,0.1,0.1),)*16 + ((0.8,0.2,0.2),)*16 + ((0.8,0.8,0.8),)*16 + ((0.6,0.8,0.4),)*16 + ((0.9,0.8,0.7),)*16
    color = list(color)

    ylim = max(ys)
    width = len(xs)
    x = range(width)
    f,ax = plt.subplots(1)
    bars = ax.bar(x, ys)
    for h in range(len(x)):
        bars[h].set_color(color[h])
    ax.set_xticks(x)
    ax.set_xticklabels([x.split(' ')[1] for x in xs], minor=False, rotation=90)
    plt.ylim(0, ylim)
    plt.xlim(0, width)

    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks([ width * (x/6.0 + 1/12.0) for x in range(6) ])
    ax2.set_xticklabels(['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G'])

    plt.savefig(target)

if __name__ == '__main__':
  plot(sys.stdin, sys.argv[1])
