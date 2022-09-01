#!/usr/bin/env python
'''
  given signatures, plot the 96 contexts - one file for each signature
'''

import argparse
import csv
import logging
import sys

import mutational_signature.plot_counts

def main(prefix, category, sig_col, fontsize, dpi, suffix):
  logging.info('reading from stdin...')
  for row in csv.DictReader(sys.stdin, delimiter='\t'):
    output = '{}{}.{}'.format(prefix, row[sig_col], suffix)
    logging.info('plotting %s...', output)
    if category == 'sbs':
      # ACAA -> C>A ACA
      vals = {'{}>{} {}{}{}'.format(k[1], k[2], k[0], k[1], k[3]): float(row[k]) for k in row if k != sig_col}
    
      logging.debug('plotting %s with %i vals i.e.. %s...', output, len(vals), vals)
      mutational_signature.plot_counts.plot_signature(vals, output, name=row[sig_col], fontsize=fontsize, dpi=dpi)
    elif category == 'ids':
      vals = {k: float(row[k]) for k in row if k != sig_col}
      mutational_signature.plot_counts.plot_signature_ids(vals, output, name=row[sig_col], fontsize=fontsize, dpi=dpi)

    logging.info('plotting %s: done', output)
  logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Plot signature contexts')
  parser.add_argument('--prefix', required=False, default='', help='prefix output files')
  parser.add_argument('--suffix', required=False, default='png', help='image type')
  parser.add_argument('--category', required=False, default='sbs', help='sbs or ids')
  parser.add_argument('--sig_col', required=False, default='Sig', help='name of signature column')
  parser.add_argument('--fontsize', required=False, default=14, type=int, help='signature name fontsize')
  parser.add_argument('--dpi', required=False, default=72, type=int, help='dpi')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.prefix, args.category, args.sig_col, args.fontsize, args.dpi, args.suffix)
