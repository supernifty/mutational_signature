#!/usr/bin/env python
'''
  adjust exposures based on given frequencies
'''

import argparse
import csv
import logging
import sys

def main(input, output, adjust_from_fh, adjust_to_fh):
  logging.info('starting...')
  adjust_from = {}
  for row in csv.DictReader(open(adjust_from_fh, 'r'), delimiter='\t'):
    adjust_from[row['Variation']] = int(row['Count'])
  adjust_to = {}
  for row in csv.DictReader(open(adjust_to_fh, 'r'), delimiter='\t'):
    adjust_to[row['Variation']] = int(row['Count'])

  fh = csv.DictWriter(sys.stdout, delimiter='\t', fieldnames=['Variation', 'Count', 'Probability'])
  fh.writeheader()
  updated = 0
  for row in csv.DictReader(sys.stdin, delimiter='\t'):
    if '>' in row['Variation']:
      context = row['Variation'].split('>')[0]
      if context in adjust_from and context in adjust_to:
        ratio = adjust_to[context] / adjust_from[context]
        logging.debug('context %s from variation %s adjusting by %f = %i / %i', context, row['Variation'], ratio, adjust_to[context], adjust_from[context])
        fh.writerow({'Variation': row['Variation'], 'Count': int(int(row['Count']) * ratio), 'Probability': float(row['Probability']) * ratio})
        updated += 1
      else:
        fh.writerow({'Variation': row['Variation'], 'Count': row['Count'], 'Probability': row['Probability']})

  logging.info('done. updated %i', updated)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Adjust exposures')
  parser.add_argument('--adjust_from', required=True, help='context of source of variants')
  parser.add_argument('--adjust_to', required=True, help='context of signatures')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(sys.stdin, sys.stdout, args.adjust_from, args.adjust_to)
