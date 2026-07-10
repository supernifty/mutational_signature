#!/usr/bin/env python
'''
  plot signature-stability point exposures with bootstrap intervals
'''

import argparse
import csv
import logging
import os
import re
import sys
from collections import defaultdict

import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


REQUIRED_COLUMNS = (
  'Sample',
  'Signature',
  'point_exposure_proportion',
  'detection_frequency',
)

DETECTION_CATEGORIES = (
  ('<50%', 0.0, 0.5, '#d73027'),
  ('50-75%', 0.5, 0.75, '#fc8d59'),
  ('75-95%', 0.75, 0.95, '#fee08b'),
  ('95-99%', 0.95, 0.99, '#91cf60'),
  ('>99%', 0.99, 1.0000001, '#1a9850'),
)


def safe_name(value):
  return re.sub(r'[^A-Za-z0-9_.-]+', '_', value)


def detect_interval_columns(fieldnames):
  proportion_columns = sorted([c for c in fieldnames if c.startswith('bootstrap_percentile_') and c.endswith('_proportion')])
  if len(proportion_columns) < 2:
    raise ValueError('could not find lower and upper bootstrap percentile proportion columns')
  return proportion_columns[0], proportion_columns[-1]


def detection_category(value):
  for label, low, high, colour in DETECTION_CATEGORIES:
    if low <= value < high:
      return label, colour
  return 'unknown', '#666666'


def read_rows(input_fh, include=None, exclude=None):
  reader = csv.DictReader(input_fh, delimiter='\t')
  if reader.fieldnames is None:
    raise ValueError('input has no header')
  missing = [c for c in REQUIRED_COLUMNS if c not in reader.fieldnames]
  if len(missing) > 0:
    raise ValueError('input is missing required columns: {}'.format(', '.join(missing)))
  lower_column, upper_column = detect_interval_columns(reader.fieldnames)

  by_signature = defaultdict(list)
  for row in reader:
    if include is not None and row['Signature'] not in include:
      continue
    if exclude is not None and row['Signature'] in exclude:
      continue
    try:
      row['_point'] = float(row['point_exposure_proportion'])
      row['_lower'] = float(row[lower_column])
      row['_upper'] = float(row[upper_column])
      row['_detect'] = float(row['detection_frequency'])
    except ValueError:
      logging.warning('skipping non-numeric row for %s %s', row.get('Sample'), row.get('Signature'))
      continue
    by_signature[row['Signature']].append(row)

  if len(by_signature) == 0:
    raise ValueError('no rows to plot')

  return by_signature, lower_column, upper_column


def plot_signature(signature, rows, output_dir, figure_width=None, figure_height=6, dpi=180, max_samples=None):
  rows = sorted(rows, key=lambda r: (-r['_point'], r['Sample']))
  if max_samples is not None:
    rows = rows[:max_samples]

  samples = [r['Sample'] for r in rows]
  point = np.array([r['_point'] for r in rows], dtype=float)
  lower = np.array([r['_lower'] for r in rows], dtype=float)
  upper = np.array([r['_upper'] for r in rows], dtype=float)
  detect = np.array([r['_detect'] for r in rows], dtype=float)
  yerr = np.vstack([np.maximum(0, point - lower), np.maximum(0, upper - point)])
  colours = [detection_category(v)[1] for v in detect]

  if figure_width is None:
    figure_width = max(10, min(28, 0.28 * len(samples) + 4))

  fig, ax = plt.subplots(figsize=(figure_width, figure_height), dpi=dpi)
  x = np.arange(len(samples))
  ax.bar(x, point, color=colours, edgecolor='black', linewidth=0.3)
  ax.errorbar(x, point, yerr=yerr, fmt='none', ecolor='black', elinewidth=0.8, capsize=2, capthick=0.8)
  ax.set_title('{} point exposure with bootstrap percentile interval'.format(signature))
  ax.set_ylabel('Exposure proportion')
  ax.set_xlabel('Samples ordered by point exposure')
  ax.set_xticks(x)
  ax.set_xticklabels(samples, rotation=90, fontsize=6)
  ax.set_ylim(0, min(1.0, max(0.02, float(np.max(upper)) * 1.12)))
  ax.grid(axis='y', alpha=0.25, linewidth=0.5)
  legend_handles = [mpatches.Patch(color=colour, label=label) for label, _, _, colour in DETECTION_CATEGORIES]
  ax.legend(handles=legend_handles, title='Detection frequency', loc='upper right', frameon=False, fontsize=8, title_fontsize=8)
  fig.tight_layout()

  target = os.path.join(output_dir, '{}.signature_stability.png'.format(safe_name(signature)))
  fig.savefig(target)
  plt.close(fig)
  return {
    'Signature': signature,
    'Plot': target,
    'Samples': len(samples),
    'MaxPointExposureProportion': '{:.6g}'.format(float(np.max(point))),
    'MaxUpperPercentileProportion': '{:.6g}'.format(float(np.max(upper))),
  }


def write_index(path, rows):
  fieldnames = ('Signature', 'Plot', 'Samples', 'MaxPointExposureProportion', 'MaxUpperPercentileProportion')
  with open(path, 'w', newline='') as fh:
    writer = csv.DictWriter(fh, delimiter='\t', fieldnames=fieldnames, lineterminator='\n')
    writer.writeheader()
    writer.writerows(rows)


def run(args):
  include = set(args.include) if args.include is not None else None
  exclude = set(args.exclude) if args.exclude is not None else None
  os.makedirs(args.output_dir, exist_ok=True)

  with open(args.input, 'r', newline='') as fh:
    by_signature, lower_column, upper_column = read_rows(fh, include=include, exclude=exclude)
  logging.info('using interval columns %s and %s', lower_column, upper_column)

  index_rows = []
  for signature in sorted(by_signature):
    index_rows.append(plot_signature(
      signature,
      by_signature[signature],
      args.output_dir,
      figure_width=args.figure_width,
      figure_height=args.figure_height,
      dpi=args.dpi,
      max_samples=args.max_samples,
    ))

  index_path = os.path.join(args.output_dir, 'index.tsv')
  write_index(index_path, index_rows)
  logging.info('wrote %i plots to %s', len(index_rows), args.output_dir)
  return index_rows


def build_parser():
  parser = argparse.ArgumentParser(description='Plot signature-stability point exposures with bootstrap intervals')
  parser.add_argument('--input', required=True, help='combined long signature-stability summary TSV')
  parser.add_argument('--output-dir', required=True, help='directory for per-signature plots')
  parser.add_argument('--include', nargs='+', help='only plot these signatures')
  parser.add_argument('--exclude', nargs='+', help='exclude these signatures')
  parser.add_argument('--max-samples', type=int, help='plot only the top N samples by point exposure for each signature')
  parser.add_argument('--figure-width', type=float, help='figure width in inches; default scales with sample count')
  parser.add_argument('--figure-height', type=float, default=6, help='figure height in inches')
  parser.add_argument('--dpi', type=int, default=180, help='output image DPI')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  return parser


if __name__ == '__main__':
  parser = build_parser()
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)
  try:
    run(args)
  except ValueError as exc:
    raise SystemExit(str(exc)) from exc
