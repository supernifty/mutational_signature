#!/usr/bin/env python
'''
  plots the breakdown of sigs
'''

import argparse
import csv
import logging
import sys

import numpy as np

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt

from pylab import rcParams
# rcParams['figure.figsize'] = 16, 10
#FIGSIZE = (16, 8)

HEIGHT_MULTIPLIER = 2
rcParams.update({'font.size': 16})

def plot(sigs, threshold, order, target, show_name, descriptions, description_threshold, highlight):
  logging.info('reading from stdin with threshold %f and order %s...', threshold, order)

  header = next(sigs)

  if order is None:
    order = header[1:]

  logging.info('plotting up to %i sigs...', len(order))

  samples = []
  data = []
  seen = set()
  seen_index = set()

  for row in sigs:
    xs = []
    samples.append(row[0])
    for i, o in enumerate(order):
      val = float(row[header.index(o)])
      xs.append(val * 100)
      if val > threshold:
        seen_index.add(i)
    data.append(xs)

  logging.info('saw %i sigs...', len(seen_index))

  # filter on seen
  order = [ o for i, o in enumerate(order) if i in seen_index ]
  if descriptions is not None:
    descriptions = [ o for i, o in enumerate(descriptions) if i in seen_index ]

  filtered = []
  for row in data:
    new_row = [ x for idx, x in enumerate(row) if idx in seen_index ]
    filtered.append(new_row)
  data = filtered

  data = list(reversed(data))
  samples = list(reversed(samples))

  # now plot
  #colors = [ '#f6597B', '#5cc45b', '#fff139', '#6373e8', '#f59251', '#e16e44', '#62e4f4', '#f052f6', '#dfff55', '#fadece', '#66a9a0', '#f6deff', '#aA7344', '#fffae8', '#901020', '#caffd3', '#909010', '#ffe8d1' ]
  colors = [
    '#cc5ac5',
    '#61c350',
    '#9358ca',
    '#9bb932',
    '#626cd4',
    '#d8ab39',
    '#5e73b8',
    '#dc842f',
    '#60a4da',
    '#cd4b33', # 10
    '#4ec584',
    '#db3666',
    '#3d9037',
    '#d5509a',
    '#69862a',
    '#cb91d4',
    '#b4ad4a',
    '#98518b',
    '#9db76f',
    '#a44657', # 20
    '#42c3c2',
    '#e07a8a',
    '#2b7650',
    '#e38f6b',
    '#5aad86',
    '#995a2b',
    '#5d8547',
    '#c39e64',
    '#676828',
    '#977b1f' # 30
  ]
  
  
  fig = plt.figure(figsize=(24, HEIGHT_MULTIPLIER * len(samples)))
  ax = fig.add_subplot(111)
  patch_handles = []
  left = np.zeros(len(samples))
  sample_id = np.arange(len(samples))

  for i in range(len(order)): # each signature
    vals = [row[i] for row in data] # all values for that signature
    if highlight is None or order[i] in highlight:
      logging.info('highlighting %s', order[i])
      alpha = 0.9
    else:
      logging.info('not highlighting %s', order[i])
      alpha = 0.5
    if show_name and descriptions is not None and descriptions[i] != '':
      patch_handles.append(ax.barh(sample_id, vals, color=colors[(int(order[i]) - 1) % len(colors)], alpha=alpha, align='center', left=left, label='{} - {}'.format(order[i], descriptions[i])))
    else:
      patch_handles.append(ax.barh(sample_id, vals, color=colors[(int(order[i]) - 1) % len(colors)], alpha=alpha, align='center', left=left, label=order[i]))
    # accumulate the left-hand offsets
    left += vals

  # go through all of the bar segments and annotate
  for j in range(len(patch_handles)):
    for i, patch in enumerate(patch_handles[j].get_children()):
        if data[i][j] >= 0.01:
          bl = patch.get_xy()
          x = 0.5 * patch.get_width() + bl[0]
          y = 0.5 * patch.get_height() + bl[1] - 0.2
          if show_name and data[i][j] > description_threshold:
            if descriptions is not None and descriptions[j] != '':
              y = 0.5 * patch.get_height() + bl[1] - 0.3
              ax.text(x,y, '%s\n%s\n%d%%' % (order[j], descriptions[j], data[i][j]), ha='center')
            else:
              ax.text(x,y, '%s\n%d%%' % (order[j], data[i][j]), ha='center')
          elif data[i][j] > 5:
            ax.text(x,y, "%d%%" % (data[i][j]), ha='center')

  ax.set_yticks(sample_id)
  ax.set_yticklabels(samples)
  ax.set_xlabel('Contribution of signatures to somatic mutations')
  ax.set_ylabel('Sample')
  ax.set_title('Somatic mutational signatures per sample')

  # place legend at right based on https://stackoverflow.com/questions/10101700/moving-matplotlib-legend-outside-of-the-axis-makes-it-cutoff-by-the-figure-box/10154763#10154763
  handles, labels = ax.get_legend_handles_labels()
  lgd = ax.legend(handles, labels, loc='upper left', bbox_to_anchor=(1.01,1.0), borderaxespad=0)
  lgd.get_frame().set_edgecolor('#000000')
  fig.savefig(target, bbox_extra_artists=(lgd,), bbox_inches='tight')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Plot signature breakdown')
  parser.add_argument('--threshold', default=0.05, type=float, help='ignore sigs below this')
  parser.add_argument('--description_threshold', default=10, type=int, help='show description if above this value')
  parser.add_argument('--order', nargs='+', required=False, help='order of signatures')
  parser.add_argument('--descriptions', nargs='+', required=False, help='ignore sigs below this')
  parser.add_argument('--target', default='sigs.png', required=False, help='output file')
  parser.add_argument('--show_signature', action='store_true', help='more logging')
  parser.add_argument('--highlight', nargs='*', required=False, help='signatures to highlight')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  plot(csv.reader(sys.stdin, delimiter='\t'), args.threshold, args.order, args.target, args.show_signature, args.descriptions, args.description_threshold, args.highlight)
