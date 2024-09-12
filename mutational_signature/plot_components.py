#!/usr/bin/env python
'''
  plots the breakdown of sigs with the expectation that it contains headers SBS1, ...
'''

import argparse
import collections
import csv
import logging
import sys

import numpy as np

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt

import matplotlib.font_manager

from pylab import rcParams
# rcParams['figure.figsize'] = 16, 10
#FIGSIZE = (16, 8)

HEIGHT_MULTIPLIER = 1.8
WIDTH_MULTIPLIER = 0.6
DPI=300

def formatter_container(vals):
  def formatter(val, loc):
    return vals[val]

  return formatter

def transpose(m):
  '''
    not really a full transpose
    [[x1], [x2], ...] -> [[x1, x2, ...]]
  '''
  return [[x[0] for x in m]]

def plot(sigs, threshold, order, target, show_name, descriptions, description_threshold, highlight, xlabel=None, ylabel=None, title=None, vertical=False, figure_height=None, figure_width=None, legend_height=None, legend_width=None, legend_y_offset=None, fontsize=12, legend_fontsize=None, legend_cols=1, denormalize=False, transparent=False, custom_colors=None, fontfamily=None, indicators=None, indicator_cmaps=None, indicator_cat=None, auto_max=False, xaxis_fontsize=None, yaxis_fontsize=None, dpi=300, linewidth=None):
  logging.info('reading from stdin with threshold %f and order %s...', threshold, order)

  import matplotlib # again!
  if fontfamily is not None:
    matplotlib.rcParams.update({'font.family': fontfamily})

  if linewidth is not None:
    matplotlib.rcParams.update({'axes.linewidth': linewidth})
    matplotlib.rcParams.update({'patch.linewidth': linewidth})
    matplotlib.rcParams.update({'xtick.major.width': linewidth})
    matplotlib.rcParams.update({'ytick.major.width': linewidth})

  header = next(sigs)
  logging.info('header is %s', header)

  if indicators is None:
    indicators = []
  if indicator_cat is None:
    indicator_cat = []
  else:
    indicator_vals = collections.defaultdict(list)

  if order is None:
    order = [x for x in header[1:] if x not in indicators]

  logging.info('plotting up to %i sig: %s...', len(order), order)

  samples = []
  data = []
  variants = []
  data_ind = collections.defaultdict(list)
  seen = set()
  seen_index = set()

  for idx, row in enumerate(sigs): # each line
    logging.debug('processing line %i: %s...', idx, row)
    xs = []
    samples.append(row[0]) # assume sample name first column
    for i, o in enumerate(order): # each signature
      logging.debug('processing line %i: extracting %s...', idx, o)
      if o in ('Error', 'Mutations'):
        continue
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

    for i in indicators:
      val = row[header.index(i)]
      if i in indicator_cat:
        if val not in indicator_vals[i]:
          indicator_vals[i].append(val)
        #data_ind[i].append([indicator_vals[i].index(val)])
        data_ind[i].append(val)
      else:
        if val == '':
          logging.debug('empty value for %s on row %i', i, idx)
          data_ind[i].append([0.0])
        else:
          data_ind[i].append([float(row[header.index(i)])])

  # sort the categorical indicators
  for i in indicator_cat:
    indicator_vals[i] = sorted(indicator_vals[i])
    data_ind[i] = [ [indicator_vals[i].index(x)] for x in data_ind[i] ]

  logging.info('saw %i sigs and data_ind: %s...', len(seen_index), data_ind)

  # filter on seen
  order = [ o for i, o in enumerate(order) if i in seen_index ]
  if descriptions is not None:
    descriptions = [ o for i, o in enumerate(descriptions) if i in seen_index ]

  filtered = []
  for row in data:
    new_row = [ x for idx, x in enumerate(row) if idx in seen_index ]
    filtered.append(new_row)
  data = filtered

  colors = {"SBS1": "#c0c0c0", "SBS2": "#41ac2f", "SBS3": "#7951d0", "SBS4": "#73d053", "SBS5": "#b969e9", "SBS6": "#3176ae", "SBS7a": "#b4b42f", "SBS7b": "#5276ec", "SBS7c": "#daae36", "SBS7d": "#9e40b5", "SBS8": "#43c673", "SBS9": "#dd4cb0", "SBS10a": "#3d9332", "SBS10b": "#de77dd", "SBS10c": "#9e77dd", "SBS10d": "#7e77dd", "SBS11": "#7bad47", "SBS12": "#9479e8", "SBS13": "#487b21", "SBS14": "#a83292", "SBS15": "#664db1", "SBS16": "#83c67d", "SBS17a": "#e18d28", "SBS17b": "#588de5", "SBS18": "#e2672a", "SBS19": "#34c7dd", "SBS20": "#646cc1", "SBS21": "#3b63ac", "SBS22": "#d74587", "SBS23": "#428647", "SBS24": "#7b51a7", "SBS25": "#505099", "SBS26": "#58b5e1", "SBS27": "#a27f1f", "SBS28": "#5acdaf", "SBS29": "#dca653", "SBS30": "#b4ba64", "SBS31": "#7d8529", "SBS32": "#bf8ade", "SBS33": "#516615", "SBS34": "#b65da7", "SBS35": "#57a87a", "SBS36": "#c84249", "SBS37": "#37b5b1", "SBS38": "#a14622", "SBS39": "#df503b", "SBS40": "#ba6e2f", "SBS41": "#589ed8", "SBS42": "#e98261", "SBS43": "#91ba2c", "SBS44": "#656413", "SBS45": "#a19fe2", "SBS46": "#756121", "SBS47": "#7e4a8d", "SBS48": "#326a38", "SBS49": "#dd8abf", "SBS50": "#1a6447", "SBS51": "#e78492", "SBS52": "#30876c", "SBS53": "#9d4d7c", "SBS54": "#919d5b", "SBS55": "#9d70ac", "SBS56": "#5b6f34", "SBS57": "#65659c", "SBS58": "#c9a865", "SBS59": "#a1455d", "SBS60": "#5e622c", "SBS84": "#b66057", "SBS85": "#dca173", 'SBS86': '#ffa600', 'SBS87': '#ffac73', 'SBS88': '#f95d6a', 'SBS89': '#d45087', 'SBS90': '#a05195', 'SBS91': '#665191', 'SBS92': '#2f4b7c', 'SBS93': '#003f5c', 'SBS94': '#123456', 'SBS96': '#64CE74', "DBS1": "#855524", "DBS2": "#9f7846", "DBS3": "#7951d0", "DBS4": "#73d053", "DBS5": "#b969e9", "DBS6": "#91ba2c", "DBS7": "#3656ca", "DBS8": "#b4b42f", "DBS9": "#5276ec", "DBS10": "#daae36", "DBS11": "#9e40b5", "ID1": "#c0c0c0", "ID2": "#7951d0", "ID3": "#41ac2f", "ID4": "#73d053", "ID5": "#b969e9", "ID6": "#91ba2c", "ID7": "#9e40b5", "ID8": "#43c673", "ID9": "#dd4cb0", "ID10": "#3d9332", "ID11": "#de77dd", "ID12": "#9479e8", "ID13": "#7bad47", "ID14": "#487b21", "ID15": "#a83292", "ID16": "#83c67d", "ID17": "#664db1", "ID18": "#964db1", "1": "#b66057", "2": "#dca173", "3": "#855524", "4": "#9f7846", "5": "#7951d0", "6": "#73d053", "7": "#b969e9", "8": "#91ba2c", "9": "#3656ca", "10": "#b4b42f", "11": "#5276ec", "12": "#daae36", "13": "#9e40b5", "14": "#de3860", "15": "#41ac2f", "16": "#7951d0", "17": "#73d053", "18": "#b969e9", "19": "#91ba2c", "20": "#9e40b5", "21": "#43c673", "22": "#dd4cb0", "23": "#3d9332", "24": "#de77dd", "25": "#7bad47", "26": "#9479e8", "27": "#487b21", "28": "#a83292", "29": "#83c67d", "30": "#664db1"}
  if custom_colors is not None:
    for sc in custom_colors:
      s, c = sc.split('=')
      colors[s] = c
  
  if vertical: ##### vertical view
    logging.info('going vertical...')
    import matplotlib.style
    matplotlib.style.use('seaborn') 

    if len(indicators) > 0:
      logging.warning('indicators beta on vertical view')

    if linewidth is not None:
      matplotlib.rcParams.update({'axes.linewidth': linewidth})

    figure_height = figure_height or 12
    if figure_width is None:
      fig = plt.figure(figsize=(WIDTH_MULTIPLIER * len(samples), figure_height))
    else:
      figure_width = figure_width or 24
      fig = plt.figure(figsize=(figure_width, figure_height))

    if len(indicators) > 0:
      logging.info('%i indicators', len(indicators))
      grid = plt.GridSpec(nrows=40, ncols=1, hspace=0.2, wspace=0.2) # number of rows, number of columns
      ax = fig.add_subplot(grid[len(indicators):, :])
      axis = []
      for i in range(len(indicators)):
        axis.append(fig.add_subplot(grid[i, :], sharex=ax)) # first rows
        axis[i].axes.get_xaxis().set_visible(False) # and no x-axis
        axis[i].axes.get_yaxis().set_visible(False) # and no y-axis
      logging.info('axis is %s', axis)
    else:
      logging.info('no indicators')
      ax = fig.add_subplot(111)

    if xaxis_fontsize is not None:
      ax.tick_params(axis='x', labelsize=xaxis_fontsize)
    if yaxis_fontsize is not None:
      ax.tick_params(axis='y', labelsize=yaxis_fontsize)

    ax.xaxis.label.set_size(fontsize)
    ax.yaxis.label.set_size(fontsize)
    if len(indicators) > 0:
      axis[0].title.set_fontsize(fontsize)
    else:
      ax.title.set_fontsize(fontsize)

    if legend_fontsize is not None:
      plt.rc('legend',fontsize=legend_fontsize)

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

    if not denormalize and not auto_max:
      ax.set_ylim(0, 100)

    if len(indicators) > 0:
      axis[0].set_title(title or 'Somatic mutational signatures per sample', fontsize=fontsize)
    else:
      ax.set_title(title or 'Somatic mutational signatures per sample', fontsize=fontsize)
  
    # place legend at right based on https://stackoverflow.com/questions/10101700/moving-matplotlib-legend-outside-of-the-axis-makes-it-cutoff-by-the-figure-box/10154763#10154763
    handles, labels = ax.get_legend_handles_labels()

    if len(indicators) > 0:
      logging.info('adding indicators: %s', data_ind)
      for i in range(len(indicators)):
        # and the legend
        cbaxes = fig.add_axes([0.91, 0.11 + 0.11 * i, 0.01, 0.1]) # left bottom width height - position of legend
        if indicators[i] in indicator_cat:
          vals = indicator_vals[indicators[i]].copy()
          tx = transpose(data_ind[indicators[i]])
          im = axis[i].imshow(tx, cmap=plt.cm.get_cmap(indicator_cmaps[i], len(vals)), aspect="auto")
          logging.debug('indicator_vals: %s', vals)
          #formatter = plt.FuncFormatter(lambda val, loc: vals[val])
          formatter = plt.FuncFormatter(formatter_container(vals))
          cbar = axis[i].figure.colorbar(im, ax=axis[i], cax=cbaxes, ticks=range(len(vals)), format=formatter) #fraction=0.4, pad=0.2, anchor=(0.1, 0.4)) #, fraction=0.04, pad=0.01, shrink=0.5)
          im.set_clim(-0.5, len(vals) - 0.5)
        else:
          tx = transpose(data_ind[indicators[i]])
          im = axis[i].imshow(tx, cmap=indicator_cmaps[i], interpolation='nearest', aspect="auto")
          cbar = axis[i].figure.colorbar(im, ax=axis[i], cax=cbaxes) #fraction=0.4, pad=0.2, anchor=(0.1, 0.4)) #, fraction=0.04, pad=0.01, shrink=0.5) # this is the legend colorbar
        cbar.set_label(indicators[i]) # this is the legend text for the indicators
        axis[i].set_yticklabels([''] + [indicators[i]], rotation=90)
        axis[i].margins(x=0)
      #cb = plt.colorbar(axi, cax=cbaxes)  
      lgd = axis[-len(indicators)].legend(handles, labels, loc='upper left', bbox_to_anchor=(1.1,1.0), borderaxespad=0.1) # these are the signature labels
    else:
      logging.info('no indicators to add')
      #lgd = ax.legend(handles, labels, loc='upper left', bbox_to_anchor=(1.01,1.0), borderaxespad=0, ncol=legend_cols)
      if legend_height is not None and legend_width is not None and legend_y_offset is not None:
        lgd = ax.legend(handles, labels, loc='upper left', bbox_to_anchor=(1.01, legend_y_offset, legend_width, legend_height), borderaxespad=0, ncol=legend_cols)
      else:
        lgd = ax.legend(handles, labels, loc='upper left', bbox_to_anchor=(1.01, 1.0), borderaxespad=0, ncol=legend_cols)

    lgd.get_frame().set_edgecolor('#000000')

    #fig.savefig(target, transparent=True, dpi=DPI, bbox_extra_artists=(lgd,), bbox_inches='tight')
    logging.info('writing to %s with dpi %i', target, dpi)
    ax.margins(x=0)
    fig.savefig(target, transparent=transparent, dpi=dpi, bbox_extra_artists=(lgd,), bbox_inches='tight')

  else: ###### horizontal
    logging.info('going horizontal...')
    figure_width = figure_width or 24
    if figure_height is None:
      fig = plt.figure(figsize=(figure_width, HEIGHT_MULTIPLIER * len(samples)))
    else:
      figure_height = figure_height or 12
      logging.info('figure size is %i, %i', figure_width, figure_height)
      fig = plt.figure(figsize=(figure_width, figure_height))

    if len(indicators) > 0:
      logging.info('%i indicators', len(indicators))
      grid = plt.GridSpec(1, 40, hspace=0, wspace=0.2)
      ax = fig.add_subplot(grid[:, 0:-len(indicators)])
      axis = []
      for i in range(len(indicators)):
        axis.append(fig.add_subplot(grid[:, -i-1], sharey=ax))
        axis[-1].axes.get_yaxis().set_visible(False)
    else:
      logging.info('no indicators')
      ax = fig.add_subplot(111)

    if xaxis_fontsize is not None:
      ax.tick_params(axis='x', labelsize=xaxis_fontsize)
    if yaxis_fontsize is not None:
      ax.tick_params(axis='y', labelsize=yaxis_fontsize)

    ax.xaxis.label.set_size(fontsize)
    ax.yaxis.label.set_size(fontsize)
    ax.title.set_fontsize(fontsize)

    if legend_fontsize is not None:
      plt.rc('legend',fontsize=legend_fontsize)

    data = list(reversed(data))
    samples = list(reversed(samples))

    patch_handles = []
    left = np.zeros(len(samples))
    sample_id = np.arange(len(samples))
  
    for i in range(len(order)): # each signature
      vals = [100 * row[i] for row in data] # all values for that signature
      if highlight is None or len(highlight) == 0 or order[i] in highlight:
        logging.debug('highlighting %s', order[i])
        alpha = 0.9
      else:
        logging.debug('not highlighting %s', order[i])
        alpha = 0.5
  
      # choose a new color
      color = colors[order[i]]
  
      if show_name and descriptions is not None and descriptions[i] != '':
        patch_handles.append(ax.barh(sample_id, vals, color=color, alpha=alpha, align='center', left=left, label='{} - {}'.format(order[i], descriptions[i])))
      else:
        patch_handles.append(ax.barh(sample_id, vals, color=color, alpha=alpha, align='center', left=left, label=order[i]))
      # accumulate the left-hand offsets
      left += vals

    if not auto_max:
      ax.set_xlim([0, 100])
  
    #PATCH_OFFSET=-0.3
    PATCH_OFFSET=-0.2
    PATCH_OFFSET_LONG=-0.3

    # go through all of the bar segments and annotate
    for j in range(len(patch_handles)):
      for i, patch in enumerate(patch_handles[j].get_children()):
          if data[i][j] >= 0.01:
            bl = patch.get_xy()
            x = 0.5 * patch.get_width() + bl[0]
            y = 0.5 * patch.get_height() + bl[1] - PATCH_OFFSET
            if show_name and data[i][j] > description_threshold:
              if descriptions is not None and descriptions[j] != '':
                # signature, description, percentage
                y = 0.5 * patch.get_height() + bl[1] - PATCH_OFFSET_LONG
                ax.annotate('%s\n%s\n%d%%' % (order[j], descriptions[j], 100 * data[i][j]), (x, y), ha='center')
              else:
                # signature
                logging.debug('placing text %s at %.2f %.2f with patch at %.2f %.2f width %.2f', order[j], x, y, bl[0], bl[1], patch.get_width())
                #ax.text(x,y, '%s' % (order[j],), ha='center')
                ax.annotate(order[j], (x, y), ha='center')
            elif data[i][j] > 5:
              # signature, percentage
              ax.text('%s\n%d%%' % (order[j], 100 * data[i][j]), (x, y), ha='center')
            # else nothing
  
    ax.set_yticks(sample_id)
    ax.set_yticklabels(samples)
    ax.set_xlabel(xlabel or 'Contribution of signatures to somatic mutations (%)')
    ax.set_ylabel(ylabel or 'Sample')
    ax.set_title(title or 'Somatic mutational signatures per sample', fontsize=fontsize)
  
    # place legend at right based on https://stackoverflow.com/questions/10101700/moving-matplotlib-legend-outside-of-the-axis-makes-it-cutoff-by-the-figure-box/10154763#10154763
    handles, labels = ax.get_legend_handles_labels()
   #fig.savefig(target, transparent=True, dpi=DPI, bbox_extra_artists=(lgd,), bbox_inches='tight')
    
    if len(indicators) > 0:
      logging.info('adding indicators: %s', data_ind)
      for i in range(len(indicators)):
        # and the legend
        cbaxes = fig.add_axes([0.91, 0.11 + 0.11 * i, 0.01, 0.1]) # left bottom width height - position of legend
        if indicators[i] in indicator_cat:
          vals = indicator_vals[indicators[i]].copy()
          im = axis[i].imshow(data_ind[indicators[i]], cmap=plt.cm.get_cmap(indicator_cmaps[i], len(vals)), aspect="auto")
          logging.debug('indicator_vals: %s', vals)
          #formatter = plt.FuncFormatter(lambda val, loc: vals[val])
          formatter = plt.FuncFormatter(formatter_container(vals))
          cbar = axis[i].figure.colorbar(im, ax=axis[i], cax=cbaxes, ticks=range(len(vals)), format=formatter) #fraction=0.4, pad=0.2, anchor=(0.1, 0.4)) #, fraction=0.04, pad=0.01, shrink=0.5)
          im.set_clim(-0.5, len(vals) - 0.5)
        else:
          im = axis[i].imshow(data_ind[indicators[i]], cmap=indicator_cmaps[i], aspect="auto")
          cbar = axis[i].figure.colorbar(im, ax=axis[i], cax=cbaxes) #fraction=0.4, pad=0.2, anchor=(0.1, 0.4)) #, fraction=0.04, pad=0.01, shrink=0.5)
        cbar.set_label(indicators[i])
        axis[i].set_xticklabels([''] + [indicators[i]], rotation=90)
      #cb = plt.colorbar(axi, cax=cbaxes)  
      lgd = axis[-len(indicators)].legend(handles, labels, loc='upper left', bbox_to_anchor=(1.1,1.0), borderaxespad=0.1)
    else:
      logging.info('no indicators to add')
      lgd = ax.legend(handles, labels, loc='upper left', bbox_to_anchor=(1.01,1.0), borderaxespad=0, ncol=legend_cols)
    lgd.get_frame().set_edgecolor('#000000')
 
    logging.info('writing to %s with dpi %i', target, dpi)
    ax.margins(x=0)
    fig.savefig(target, transparent=transparent, dpi=dpi, bbox_extra_artists=(lgd,), bbox_inches='tight')
  plt.close('all')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Plot signature breakdown')
  parser.add_argument('--threshold', default=0.05, type=float, help='ignore sigs below this')
  parser.add_argument('--description_threshold', default=0.1, type=float, help='show description if above this value')
  parser.add_argument('--order', nargs='+', required=False, help='order of signatures')
  parser.add_argument('--colors', nargs='+', required=False, help='custom signature colors of the form Sig=Color...')
  parser.add_argument('--descriptions', nargs='+', required=False, help='ignore sigs below this')
  parser.add_argument('--target', default='sigs.png', required=False, help='output file')
  parser.add_argument('--show_signature', action='store_true', help='show signature name on graph')
  parser.add_argument('--highlight', nargs='*', required=False, help='signatures to highlight')
  parser.add_argument('--xlabel', required=False, help='x axis label')
  parser.add_argument('--ylabel', required=False, help='y axis label')
  parser.add_argument('--height', required=False, type=float, help='height of vertical plot')
  parser.add_argument('--width', required=False, type=float, help='width of vertical plot')
  parser.add_argument('--legend_height', required=False, type=float, help='legend height of vertical plot')
  parser.add_argument('--legend_width', required=False, type=float, help='legend width of vertical plot')
  parser.add_argument('--legend_y_offset', required=False, type=float, help='legend y offset')
  parser.add_argument('--legend_cols', required=False, type=int, default=1, help='legend cols')
  parser.add_argument('--linewidth', required=False, type=float, help='linewidth')
  parser.add_argument('--title', required=False, help='title')
  parser.add_argument('--fontsize', required=False, default=12, type=int, help='plot font size')
  parser.add_argument('--fontfamily', required=False, help='plot font family')
  parser.add_argument('--dpi', required=False, default=300, type=int, help='plot dpi')
  parser.add_argument('--legend_fontsize', required=False, type=int, help='legend font size')
  parser.add_argument('--xaxis_fontsize', required=False, type=int, help='xaxis font size')
  parser.add_argument('--yaxis_fontsize', required=False, type=int, help='yaxis font size')
  parser.add_argument('--vertical', action='store_true', help='plot vertically')
  parser.add_argument('--denormalize', action='store_true', help='do not constrain to 100 percent')
  parser.add_argument('--auto_max', action='store_true', help='do not set max to 100 percent')
  parser.add_argument('--transparent', action='store_true', help='transparent')
  parser.add_argument('--indicators', nargs='*', required=False, help='columns to use as indicators')
  parser.add_argument('--indicator_cmaps', nargs='*', required=False, help='cmap to use for each indicator')
  parser.add_argument('--indicator_cat', nargs='*', required=False, help='categorical indicators')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  plot(csv.reader(sys.stdin, delimiter='\t'), args.threshold, args.order, args.target, args.show_signature, args.descriptions, args.description_threshold, args.highlight, args.xlabel, args.ylabel, args.title, args.vertical, args.height, args.width, args.legend_height, args.legend_width, args.legend_y_offset, args.fontsize, args.legend_fontsize, args.legend_cols, args.denormalize, args.transparent, args.colors, args.fontfamily, args.indicators, args.indicator_cmaps, args.indicator_cat, args.auto_max, args.xaxis_fontsize, args.yaxis_fontsize, args.dpi, args.linewidth)
