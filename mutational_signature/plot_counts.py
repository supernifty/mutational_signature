#!/usr/bin/env python
'''
  plots the 96 snv contexts given a count input
'''

import argparse
import logging
import sys

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.patches

from pylab import rcParams

WIDTH=18
HEIGHT=6

rcParams['figure.figsize'] = WIDTH, HEIGHT

def plot(counts, target, style='sbs', name=None, dpi=72):
  logging.info('reading from stdin...')
  first = True
  vals = {}
  for line in counts:
    if first:
      first = False
      continue
    if line.count('\t') > 2: # 4+ fields
      v, count, probability, _ = line.strip('\n').split('\t', 3)
    elif line.count('\t') > 1: # 3 fields
      v, count, probability = line.strip('\n').split('\t')
    else: # 2 fields
      v, probability = line.strip('\n').split('\t')
    if style == 'sbs':
      if len(v) == 5 and v[3] == '>':
        variation = '{}>{} {}'.format(v[1], v[4], v[0:3])
        vals[variation] = float(probability)
      else:
        logging.info('skipped {}\n'.format(v))
    else:
      vals[v] = float(probability)

  sys.stderr.write('{} contexts\n'.format(len(vals)))
  if style == 'sbs':
    plot_signature(vals, target, name=name, dpi=dpi)
  elif style == 'id':
    plot_signature_ids(vals, target, name=name, dpi=dpi)
  else:
    logging.warn('unrecognised plot type %s', style)

def plot_signature(vals, target, name=None, fontsize=14, figure_width=10, figure_height=4, dpi=72):
  #xs = sorted(vals.keys())
  xs = sorted(['{}>{} {}{}{}'.format(k[1], k[2], k[0], k[1], k[3]) for k in ['ACAA', 'ACAC', 'ACAG', 'ACAT', 'ACGA', 'ACGC', 'ACGG', 'ACGT', 'ACTA', 'ACTC', 'ACTG', 'ACTT', 'ATAA', 'ATAC', 'ATAG', 'ATAT', 'ATCA', 'ATCC', 'ATCG', 'ATCT', 'ATGA', 'ATGC', 'ATGG', 'ATGT', 'CCAA', 'CCAC', 'CCAG', 'CCAT', 'CCGA', 'CCGC', 'CCGG', 'CCGT', 'CCTA', 'CCTC', 'CCTG', 'CCTT', 'CTAA', 'CTAC', 'CTAG', 'CTAT', 'CTCA', 'CTCC', 'CTCG', 'CTCT', 'CTGA', 'CTGC', 'CTGG', 'CTGT', 'GCAA', 'GCAC', 'GCAG', 'GCAT', 'GCGA', 'GCGC', 'GCGG', 'GCGT', 'GCTA', 'GCTC', 'GCTG', 'GCTT', 'GTAA', 'GTAC', 'GTAG', 'GTAT', 'GTCA', 'GTCC', 'GTCG', 'GTCT', 'GTGA', 'GTGC', 'GTGG', 'GTGT', 'TCAA', 'TCAC', 'TCAG', 'TCAT', 'TCGA', 'TCGC', 'TCGG', 'TCGT', 'TCTA', 'TCTC', 'TCTG', 'TCTT', 'TTAA', 'TTAC', 'TTAG', 'TTAT', 'TTCA', 'TTCC', 'TTCG', 'TTCT', 'TTGA', 'TTGC', 'TTGG', 'TTGT']])
  ys = list([100 * vals.get(x, 0) for x in xs]) # convert to %

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
  plt.xlim(-0.5, width)

  ax2 = ax.twiny()
  ax2.set_xlim(ax.get_xlim())
  ax2.set_xticks([ width * (x/6.0 + 1/12.0) for x in range(6) ])
  ax2.set_xticklabels(['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G'], fontsize=fontsize)

  for x in range(1, 6):
    ax.axvline(x=x * 16 - 0.5, color='#e0e0e0')

  if name is not None:
    plt.annotate(name, xy=(0.01, 1-fontsize * 0.003), xycoords='axes fraction', fontsize=fontsize)

  ax.set_ylabel('Mutation probability (%)', fontsize=fontsize)
  ax.set_xlabel('Mutational Context', fontsize=fontsize)

  #plt.figure(figsize=(figure_width, figure_height))
  plt.savefig(target, bbox_inches='tight', dpi=dpi)
  plt.close()


def plot_signature_ids(vals, target, name=None, fontsize=14, figure_width=6, figure_height=2, dpi=72):
  xs = sorted(vals.keys())

  contexts = ('DEL_C_1_0', 'DEL_C_1_1', 'DEL_C_1_2', 'DEL_C_1_3', 'DEL_C_1_4', 'DEL_C_1_5+', 
    'DEL_T_1_0', 'DEL_T_1_1', 'DEL_T_1_2', 'DEL_T_1_3', 'DEL_T_1_4', 'DEL_T_1_5+', 
    'INS_C_1_0', 'INS_C_1_1', 'INS_C_1_2', 'INS_C_1_3', 'INS_C_1_4', 'INS_C_1_5+', 
    'INS_T_1_0', 'INS_T_1_1', 'INS_T_1_2', 'INS_T_1_3', 'INS_T_1_4', 'INS_T_1_5+', 
    'DEL_repeats_2_0', 'DEL_repeats_2_1', 'DEL_repeats_2_2', 'DEL_repeats_2_3', 'DEL_repeats_2_4', 'DEL_repeats_2_5+', 'DEL_repeats_3_0', 'DEL_repeats_3_1', 'DEL_repeats_3_2', 'DEL_repeats_3_3', 'DEL_repeats_3_4', 'DEL_repeats_3_5+', 'DEL_repeats_4_0', 'DEL_repeats_4_1', 'DEL_repeats_4_2', 'DEL_repeats_4_3', 'DEL_repeats_4_4', 'DEL_repeats_4_5+', 'DEL_repeats_5+_0', 'DEL_repeats_5+_1', 'DEL_repeats_5+_2', 'DEL_repeats_5+_3', 'DEL_repeats_5+_4', 'DEL_repeats_5+_5+', 
    'INS_repeats_2_0', 'INS_repeats_2_1', 'INS_repeats_2_2', 'INS_repeats_2_3', 'INS_repeats_2_4', 'INS_repeats_2_5+', 'INS_repeats_3_0', 'INS_repeats_3_1', 'INS_repeats_3_2', 'INS_repeats_3_3', 'INS_repeats_3_4', 'INS_repeats_3_5+', 'INS_repeats_4_0', 'INS_repeats_4_1', 'INS_repeats_4_2', 'INS_repeats_4_3', 'INS_repeats_4_4', 'INS_repeats_4_5+', 'INS_repeats_5+_0', 'INS_repeats_5+_1', 'INS_repeats_5+_2', 'INS_repeats_5+_3', 'INS_repeats_5+_4', 'INS_repeats_5+_5+',
    'DEL_MH_2_1', 'DEL_MH_3_1', 'DEL_MH_3_2', 'DEL_MH_4_1', 'DEL_MH_4_2', 'DEL_MH_4_3', 'DEL_MH_5+_1', 'DEL_MH_5+_2', 'DEL_MH_5+_3', 'DEL_MH_5+_4', 'DEL_MH_5+_5+') 

  ys = list([vals[x] for x in xs if x in contexts])

  display = ('1', '2', '3', '4', '5', '6+', 
    '1', '2', '3', '4', '5', '6+', 
    '0', '1', '2', '3', '4', '5+', 
    '0', '1', '2', '3', '4', '5+', 
    '1', '2', '3', '4', '5', '6+', 
    '1', '2', '3', '4', '5', '6+', 
    '1', '2', '3', '4', '5', '6+', 
    '1', '2', '3', '4', '5', '6+', 
    '0', '1', '2', '3', '4', '5+', 
    '0', '1', '2', '3', '4', '5+', 
    '0', '1', '2', '3', '4', '5+', 
    '0', '1', '2', '3', '4', '5+',
    '1', '1', '2', '1', '2', '3', '1', '2', '3', '4', '5') 


  color = ('#FECA82',) * 6 + ('#FF9400',) * 6 + ('#BBE19E',) * 6 + ('#3EAD3D',) * 6 + ('#FCD6C1',) * 6 + ('#FE9E7D',) * 6 + ('#F65B47',) * 6 + ('#C82F1D',) * 6 + ('#D8E8F7',) * 6 + ('#A5CFE5',) * 6 + ('#5AAAD3',) * 6 + ('#167AB7',) * 6 + ('#E8E8EE',) * 1 + ('#C2C3DE',) * 2 + ('#999ACE',) * 3 + ('#7758A8',) * 5
  color = list(color)

  ylim = max(ys)
  width = len(contexts) # 84
  x = range(width)
  fig, ax = plt.subplots(1)
  fig.subplots_adjust(bottom=0.10, top=0.85)

  bars = ax.bar(x, [vals[x] for x in contexts])
  for h in range(len(x)):
    bars[h].set_color(color[h])

  # bottom axes takes display list
  ax.set_xticks(x)
  ax.set_xticklabels([x for x in display], minor=False, rotation=90)
  plt.ylim(0, ylim)
  plt.xlim(-0.5, width)

  # top axes
  ax2 = ax.twiny()
  ax2.spines['top'].set_position(('outward', 20))
  ax2.spines['top'].set_visible(False)
  ax2.tick_params(axis='both', which='both', length=0)
  ax2.set_xlim(ax.get_xlim())
  ax2.set_xticks([ width * (0/84 + 6/84), width * (12/84 + 6/84), width * (24/84 + 12/84), width * (48/84 + 12/84), width * (72/84 + 6/84)])
  ax2.set_xticklabels(['1bp deletion', '1bp insertion', '>1bp deletion at repeat\n(deletion length)', '>1bp insertion at repeat\n(insertion length)', 'Deletion with microhomology\n(deletion length)'])

  # bottom axes
  ax3 = ax.twiny()
  ax3.xaxis.set_ticks_position('bottom')
  ax3.set_xlim(ax.get_xlim())
  ax3.spines['bottom'].set_position(('outward', 20))
  ax3.spines['bottom'].set_visible(False)
  ax3.tick_params(axis='both', which='both', length=0)
  ax3.set_xticks([ width * (0/84 + 6/84), width * (12/84 + 6/84), width * (24/84 + 12/84), width * (48/84 + 12/84), width * (72/84 + 6/84)])
  ax3.set_xticklabels(['Homopolymer length', 'Homopolymer length', 'Number of repeat units', 'Number of repeat units', 'Microhomology length'])

  # additional help
  squiggem = 0.01 * ylim
  patches = (('C', 3, '#feca82', 6), ('T', 9, '#FF9400', 6), 
    ('C', 15, '#BBE19E', 6), ('T', 21, '#3EAD3D', 6), 
    ('2', 27, '#FCD6C1', 6), ('3', 33, '#FE9E7D', 6), ('4', 39, '#F65B47', 6), ('5+', 45, '#C82F1D', 6), 
    ('2', 51, '#D8E8F7', 6), ('3', 57, '#A5CFE5', 6), ('4', 63, '#5AAAD3', 6), ('5+', 69, '#167AB7', 6),
    ('2', 72.5, '#E8E8EE', 1), ('3', 74, '#C2C3DE', 2), ('4', 76.5, '#5AAAD3', 3), ('5+', 80.5, '#7758A8', 5))
  for i, patch in enumerate(patches):
    ax.add_patch(matplotlib.patches.Rectangle(xy=(patch[1] - patch[3]/2 - 0.5, ylim + squiggem), width=patch[3], height=ylim/20, color=patch[2], clip_on=False))
    ax.text(s=patch[0], x=patch[1] - 0.75, y=ylim + 2 * squiggem, fontsize=10, ha="left", color='black')
    if i < len(patches) - 1:
      ax.axvline(x=patch[1] + patch[3]/2 - 0.5, color='#e0e0e0')

  if name is not None:
    plt.annotate(name, xy=(0.01, 1-fontsize * 0.003), xycoords='axes fraction', fontsize=fontsize)
  #plt.figure(figsize=(figure_width, figure_height))
  plt.savefig(target, bbox_inches='tight', dpi=dpi)
  plt.close()
 
if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Plot contexts')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  parser.add_argument('--type', required=False, default='sbs', help='sbs or id')
  parser.add_argument('--name', required=False, help='name of plot')
  parser.add_argument('--target', required=False, default='plot.png', help='output filename')
  parser.add_argument('--dpi', required=False, default=72, type=int, help='dpi')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  plot(sys.stdin, args.target, args.type, args.name, args.dpi)

