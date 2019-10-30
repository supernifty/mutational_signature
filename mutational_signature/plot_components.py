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

HEIGHT_MULTIPLIER = 1.8
WIDTH_MULTIPLIER = 0.6
DPI=300

def plot(sigs, threshold, order, target, show_name, descriptions, description_threshold, highlight, xlabel=None, ylabel=None, title=None, vertical=False, figure_height=None, figure_width=None, legend_height=None, legend_width=None, legend_y_offset=None, fontsize=12, legend_cols=1, denormalize=False):
  logging.info('reading from stdin with threshold %f and order %s...', threshold, order)
  rcParams.update({'font.size': fontsize})

  header = next(sigs)
  logging.info('header is %s', header)

  if order is None:
    order = header[1:]

  logging.info('plotting up to %i sigs...', len(order))

  samples = []
  data = []
  variants = []
  seen = set()
  seen_index = set()

  for row in sigs: # each line
    xs = []
    samples.append(row[0]) # assume sample name first column
    for i, o in enumerate(order): # each signature
      val = float(row[header.index(o)])
      xs.append(val)
      if val > threshold:
        seen_index.add(i)
    data.append(xs)
    # denormalize
    if denormalize:
      variants.append(float(row[header.index('Variants')]))
    else:
      variants.append(100.0) # do not consider variant count

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

  colors = {"SBS1": "#c0c0c0", "SBS2": "#41ac2f", "SBS3": "#7951d0", "SBS4": "#73d053", "SBS5": "#b969e9", "SBS6": "#91ba2c", "SBS7a": "#b4b42f", "SBS7b": "#5276ec", "SBS7c": "#daae36", "SBS7d": "#9e40b5", "SBS8": "#43c673", "SBS9": "#dd4cb0", "SBS10a": "#3d9332", "SBS10b": "#de77dd", "SBS11": "#7bad47", "SBS12": "#9479e8", "SBS13": "#487b21", "SBS14": "#a83292", "SBS15": "#83c67d", "SBS16": "#664db1", "SBS17a": "#e18d28", "SBS17b": "#588de5", "SBS18": "#e2672a", "SBS19": "#34c7dd", "SBS20": "#cf402b", "SBS21": "#5acdaf", "SBS22": "#d74587", "SBS23": "#428647", "SBS24": "#7b51a7", "SBS25": "#b4ba64", "SBS26": "#646cc1", "SBS27": "#a27f1f", "SBS28": "#3b63ac", "SBS29": "#dca653", "SBS30": "#505099", "SBS31": "#7d8529", "SBS32": "#bf8ade", "SBS33": "#516615", "SBS34": "#b65da7", "SBS35": "#57a87a", "SBS36": "#c84249", "SBS37": "#37b5b1", "SBS38": "#a14622", "SBS39": "#58b5e1", "SBS40": "#ba6e2f", "SBS41": "#589ed8", "SBS42": "#e98261", "SBS43": "#3176ae", "SBS44": "#656413", "SBS45": "#a19fe2", "SBS46": "#756121", "SBS47": "#7e4a8d", "SBS48": "#326a38", "SBS49": "#dd8abf", "SBS50": "#1a6447", "SBS51": "#e78492", "SBS52": "#30876c", "SBS53": "#9d4d7c", "SBS54": "#919d5b", "SBS55": "#9d70ac", "SBS56": "#5b6f34", "SBS57": "#65659c", "SBS58": "#c9a865", "SBS59": "#a1455d", "SBS60": "#5e622c", "SBS84": "#b66057", "SBS85": "#dca173", "DBS1": "#855524", "DBS2": "#9f7846", "DBS3": "#7951d0", "DBS4": "#73d053", "DBS5": "#b969e9", "DBS6": "#91ba2c", "DBS7": "#3656ca", "DBS8": "#b4b42f", "DBS9": "#5276ec", "DBS10": "#daae36", "DBS11": "#9e40b5", "ID1": "#de3860", "ID2": "#41ac2f", "ID3": "#7951d0", "ID4": "#73d053", "ID5": "#b969e9", "ID6": "#91ba2c", "ID7": "#9e40b5", "ID8": "#43c673", "ID9": "#dd4cb0", "ID10": "#3d9332", "ID11": "#de77dd", "ID12": "#7bad47", "ID13": "#9479e8", "ID14": "#487b21", "ID15": "#a83292", "ID16": "#83c67d", "ID17": "#664db1", "1": "#b66057", "2": "#dca173", "3": "#855524", "4": "#9f7846", "5": "#7951d0", "6": "#73d053", "7": "#b969e9", "8": "#91ba2c", "9": "#3656ca", "10": "#b4b42f", "11": "#5276ec", "12": "#daae36", "13": "#9e40b5", "14": "#de3860", "15": "#41ac2f", "16": "#7951d0", "17": "#73d053", "18": "#b969e9", "19": "#91ba2c", "20": "#9e40b5", "21": "#43c673", "22": "#dd4cb0", "23": "#3d9332", "24": "#de77dd", "25": "#7bad47", "26": "#9479e8", "27": "#487b21", "28": "#a83292", "29": "#83c67d", "30": "#664db1"}
  
  if vertical:
    import matplotlib.style
    matplotlib.style.use('seaborn') 

    figure_height = figure_height or 12
    if figure_width is None:
      fig = plt.figure(figsize=(WIDTH_MULTIPLIER * len(samples), figure_height))
    else:
      figure_width = figure_width or 24
      fig = plt.figure(figsize=(figure_width, figure_height))
    ax = fig.add_subplot(111)
    patch_handles = []

    sample_id = np.arange(len(samples))
    width=0.8
    bottom=[0] * len(samples)
    for i in range(len(order)): # each signature
      vals = [row[i] * variant for row, variant in zip(data, variants)] # all values for that signature
      if highlight is None or len(highlight) == 0 or order[i] in highlight:
        logging.info('highlighting %s', order[i])
        alpha = 0.9
      else:
        logging.info('not highlighting %s', order[i])
        alpha = 0.5
  
      # choose a new color
      color = colors[order[i]]
      logging.debug('sample_id %s vals %s bottom %s', sample_id, vals, bottom)

      if show_name and descriptions is not None and descriptions[i] != '':
        patch_handles.append(ax.bar(sample_id, vals, width, bottom=bottom, color=color, alpha=alpha, label='{} - {}'.format(order[i], descriptions[i])))
      else:
        patch_handles.append(ax.bar(sample_id, vals, width, bottom=bottom, color=color, alpha=alpha, label=order[i]))

      bottom = [bottom[i] + vals[i] for i in range(len(bottom))]
 
    ax.set_xticks(sample_id)
    ax.set_xticklabels(samples, rotation=90)

    if denormalize:
      ax.set_ylabel(ylabel or 'Contribution of signatures to somatic mutations (variants)')
    else:
      ax.set_ylabel(ylabel or 'Contribution of signatures to somatic mutations (%)')
    ax.set_xlabel(xlabel or 'Sample')

    if not denormalize:
      ax.set_ylim(0, 100)

    ax.set_title(title or 'Somatic mutational signatures per sample')
  
    # place legend at right based on https://stackoverflow.com/questions/10101700/moving-matplotlib-legend-outside-of-the-axis-makes-it-cutoff-by-the-figure-box/10154763#10154763
    handles, labels = ax.get_legend_handles_labels()

    if legend_height is not None and legend_width is not None and legend_y_offset is not None:
      lgd = ax.legend(handles, labels, loc='upper left', bbox_to_anchor=(1.01, legend_y_offset, legend_width, legend_height), borderaxespad=0, ncol=legend_cols)
    else:
      lgd = ax.legend(handles, labels, loc='upper left', bbox_to_anchor=(1.01, 1.0), borderaxespad=0, ncol=legend_cols)
    lgd.get_frame().set_edgecolor('#000000')
    #fig.savefig(target, transparent=True, dpi=DPI, bbox_extra_artists=(lgd,), bbox_inches='tight')
    fig.savefig(target, transparent=False, dpi=DPI, bbox_extra_artists=(lgd,), bbox_inches='tight')

  else: # horizontal
    figure_width = figure_width or 24
    if figure_height is None:
      fig = plt.figure(figsize=(figure_width, HEIGHT_MULTIPLIER * len(samples)))
    else:
      figure_height = figure_height or 12
      fig = plt.figure(figsize=(figure_width, figure_width))
    ax = fig.add_subplot(111)

    data = list(reversed(data))
    samples = list(reversed(samples))

    patch_handles = []
    left = np.zeros(len(samples))
    sample_id = np.arange(len(samples))
  
    for i in range(len(order)): # each signature
      vals = [row[i] for row in data] # all values for that signature
      if highlight is None or len(highlight) == 0 or order[i] in highlight:
        logging.info('highlighting %s', order[i])
        alpha = 0.9
      else:
        logging.info('not highlighting %s', order[i])
        alpha = 0.5
  
      # choose a new color
      color = colors[order[i]]
  
      if show_name and descriptions is not None and descriptions[i] != '':
        patch_handles.append(ax.barh(sample_id, vals, color=color, alpha=alpha, align='center', left=left, label='{} - {}'.format(order[i], descriptions[i])))
      else:
        patch_handles.append(ax.barh(sample_id, vals, color=color, alpha=alpha, align='center', left=left, label=order[i]))
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
                # signature, description, percentage
                y = 0.5 * patch.get_height() + bl[1] - 0.3
                ax.text(x,y, '%s\n%s\n%d%%' % (order[j], descriptions[j], data[i][j]), ha='center')
              else:
                # signature, percentage
                ax.text(x,y, '%s\n%d%%' % (order[j], data[i][j]), ha='center')
            elif data[i][j] > 5:
              # signature, percentage
              ax.text(x,y, '%s\n%d%%' % (order[j], data[i][j]), ha='center')
            # else nothing
  
    ax.set_yticks(sample_id)
    ax.set_yticklabels(samples)
    ax.set_xlabel(xlabel or 'Contribution of signatures to somatic mutations (%)')
    ax.set_ylabel(ylabel or 'Sample')
    ax.set_title(title or 'Somatic mutational signatures per sample')
  
    # place legend at right based on https://stackoverflow.com/questions/10101700/moving-matplotlib-legend-outside-of-the-axis-makes-it-cutoff-by-the-figure-box/10154763#10154763
    handles, labels = ax.get_legend_handles_labels()
    lgd = ax.legend(handles, labels, loc='upper left', bbox_to_anchor=(1.01,1.0), borderaxespad=0)
    lgd.get_frame().set_edgecolor('#000000')
    #fig.savefig(target, transparent=True, dpi=DPI, bbox_extra_artists=(lgd,), bbox_inches='tight')
    fig.savefig(target, transparent=False, dpi=DPI, bbox_extra_artists=(lgd,), bbox_inches='tight')
  plt.close('all')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Plot signature breakdown')
  parser.add_argument('--threshold', default=0.05, type=float, help='ignore sigs below this')
  parser.add_argument('--description_threshold', default=10, type=int, help='show description if above this value')
  parser.add_argument('--order', nargs='+', required=False, help='order of signatures')
  parser.add_argument('--descriptions', nargs='+', required=False, help='ignore sigs below this')
  parser.add_argument('--target', default='sigs.png', required=False, help='output file')
  parser.add_argument('--show_signature', action='store_true', help='more logging')
  parser.add_argument('--highlight', nargs='*', required=False, help='signatures to highlight')
  parser.add_argument('--xlabel', required=False, help='x axis label')
  parser.add_argument('--ylabel', required=False, help='y axis label')
  parser.add_argument('--height', required=False, type=float, help='height of vertical plot')
  parser.add_argument('--width', required=False, type=float, help='width of vertical plot')
  parser.add_argument('--legend_height', required=False, type=float, help='legend height of vertical plot')
  parser.add_argument('--legend_width', required=False, type=float, help='legend width of vertical plot')
  parser.add_argument('--legend_y_offset', required=False, type=float, help='legend y offset')
  parser.add_argument('--legend_cols', required=False, type=int, default=1, help='legend cols')
  parser.add_argument('--title', required=False, help='title')
  parser.add_argument('--fontsize', required=False, default=12, type=int, help='plot font size')
  parser.add_argument('--vertical', action='store_true', help='plot vertically')
  parser.add_argument('--denormalize', action='store_true', help='do not constrain to 100%')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  plot(csv.reader(sys.stdin, delimiter='\t'), args.threshold, args.order, args.target, args.show_signature, args.descriptions, args.description_threshold, args.highlight, args.xlabel, args.ylabel, args.title, args.vertical, args.height, args.width, args.legend_height, args.legend_width, args.legend_y_offset, args.fontsize, args.legend_cols, args.denormalize)
