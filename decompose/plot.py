#!/usr/bin/env python

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
    v, count, probability = line.strip('\n').split('\t')
    variation = '{}>{} {}'.format(v[1], v[4], v[0:3])
    vals[variation] = float(probability)

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
