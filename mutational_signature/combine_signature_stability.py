#!/usr/bin/env python
'''
  combine signature-stability summary files
'''

import argparse
import csv
import logging
import os
import sys


REQUIRED_COLUMNS = (
  'Sample',
  'Signature',
  'point_exposure_proportion',
  'detection_frequency',
)

DEFAULT_COLUMNS = (
  'Sample',
  'Signature',
  'point_exposure_proportion',
  'detection_frequency',
  'point_exposure_mutations',
  'total_mutations',
  'perturbation_method',
  'seed',
  'fitting_mode',
  'metric',
  'solver',
  'Source',
)


def detect_interval_columns(fieldnames):
  proportion = [c for c in fieldnames if c.startswith('bootstrap_percentile_') and c.endswith('_proportion')]
  mutations = [c for c in fieldnames if c.startswith('bootstrap_percentile_') and c.endswith('_mutations')]
  if len(proportion) < 2:
    raise ValueError('could not find lower and upper bootstrap percentile proportion columns')
  return sorted(mutations), sorted(proportion)


def read_summary(path):
  with open(path, 'r', newline='') as fh:
    reader = csv.DictReader(fh, delimiter='\t')
    if reader.fieldnames is None:
      raise ValueError('{} has no header'.format(path))
    missing = [c for c in REQUIRED_COLUMNS if c not in reader.fieldnames]
    if len(missing) > 0:
      raise ValueError('{} is missing required columns: {}'.format(path, ', '.join(missing)))
    detect_interval_columns(reader.fieldnames)
    rows = []
    samples = set()
    for row in reader:
      row['Source'] = path
      rows.append(row)
      samples.add(row['Sample'])
    if len(rows) == 0:
      raise ValueError('{} has no rows'.format(path))
    if len(samples) > 1:
      raise ValueError('{} has multiple Sample values: {}'.format(path, ', '.join(sorted(samples))))
    return reader.fieldnames + ['Source'], rows


def include_row(row, include, exclude):
  if include is not None and row['Signature'] not in include:
    return False
  if exclude is not None and row['Signature'] in exclude:
    return False
  return True


def load_rows(files, include=None, exclude=None):
  all_rows = []
  all_columns = []
  seen = set()
  for path in files:
    columns, rows = read_summary(path)
    for column in columns:
      if column not in all_columns:
        all_columns.append(column)
    for row in rows:
      if not include_row(row, include, exclude):
        continue
      key = (row['Sample'], row['Signature'])
      if key in seen:
        raise ValueError('duplicate Sample+Signature row: {} {}'.format(row['Sample'], row['Signature']))
      seen.add(key)
      all_rows.append(row)
  return all_columns, all_rows


def default_columns(all_columns):
  columns = []
  mutation_intervals, proportion_intervals = detect_interval_columns(all_columns)
  for column in DEFAULT_COLUMNS:
    if column in all_columns:
      columns.append(column)
    if column == 'point_exposure_proportion':
      for interval_column in proportion_intervals:
        if interval_column not in columns:
          columns.append(interval_column)
    if column == 'point_exposure_mutations':
      for interval_column in mutation_intervals:
        if interval_column not in columns:
          columns.append(interval_column)
  return columns


def write_long(rows, columns, out):
  writer = csv.DictWriter(out, delimiter='\t', fieldnames=columns, lineterminator='\n', extrasaction='ignore')
  writer.writeheader()
  for row in rows:
    writer.writerow({column: row.get(column, '') for column in columns})


def metric_columns(columns):
  return [c for c in columns if c not in ('Sample', 'Signature', 'Source')]


def write_wide(rows, columns, out):
  metrics = metric_columns(columns)
  samples = []
  signatures = []
  by_sample = {}
  for row in rows:
    sample = row['Sample']
    signature = row['Signature']
    if sample not in samples:
      samples.append(sample)
    if signature not in signatures:
      signatures.append(signature)
    by_sample.setdefault(sample, {})[signature] = row

  fieldnames = ['Sample']
  for signature in signatures:
    for metric in metrics:
      fieldnames.append('{}_{}'.format(signature, metric))

  writer = csv.DictWriter(out, delimiter='\t', fieldnames=fieldnames, lineterminator='\n')
  writer.writeheader()
  for sample in samples:
    out_row = {'Sample': sample}
    for signature in signatures:
      row = by_sample.get(sample, {}).get(signature, {})
      for metric in metrics:
        out_row['{}_{}'.format(signature, metric)] = row.get(metric, '')
    writer.writerow(out_row)


def run(args, out=None):
  if out is None:
    out = sys.stdout
  include = set(args.include) if args.include is not None else None
  exclude = set(args.exclude) if args.exclude is not None else None
  all_columns, rows = load_rows(args.files, include=include, exclude=exclude)
  if len(rows) == 0:
    raise ValueError('no rows to combine')
  columns = all_columns if args.all_columns else default_columns(all_columns)
  if args.wide:
    write_wide(rows, columns, out)
  else:
    write_long(rows, columns, out)


def build_parser():
  parser = argparse.ArgumentParser(description='Combine signature-stability summary TSV files')
  parser.add_argument('--files', required=True, nargs='+', help='signature-stability summary TSV files')
  parser.add_argument('--all-columns', action='store_true', help='include every input summary column')
  parser.add_argument('--wide', action='store_true', help='write one row per sample with signature_metric columns')
  parser.add_argument('--include', nargs='+', help='only include these signatures')
  parser.add_argument('--exclude', nargs='+', help='exclude these signatures')
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
  except BrokenPipeError:
    sys.stdout = open(os.devnull, 'w')
    raise SystemExit(0)
  except ValueError as exc:
    raise SystemExit(str(exc)) from exc
