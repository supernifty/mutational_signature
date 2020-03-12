#!/usr/bin/env python
'''
  given signatures, plot the 96 contexts
'''

import argparse
import csv
import logging
import sys

import mutational_signature.plot_counts

def main(prefix):
  logging.info('reading from stdin...')
  for row in csv.DictReader(sys.stdin, delimiter='\t'):
    output = '{}{}.png'.format(prefix, row['Sig'])
    logging.info('plotting %s...', output)
    # ACAA -> C>A ACA
    vals = {'{}>{} {}{}{}'.format(k[1], k[2], k[0], k[1], k[3]): float(row[k]) for k in row if k != 'Sig'}
    
    logging.debug('plotting %s with %i vals i.e.. %s...', output, len(vals), vals)
    mutational_signature.plot_counts.plot_signature(vals, output)
    logging.info('plotting %s: done', output)
  logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Assess MSI')
  parser.add_argument('--prefix', required=True, default='', help='prefix output files')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.prefix)
