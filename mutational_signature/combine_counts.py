#!/usr/bin/env python
'''
  combine signature outputs
'''

import argparse
import csv
import logging
import sys

def main(files, use_probability):
  logging.info('starting...')

  rows = {}
  seen = set()
  for fn in files:
    logging.debug('reading %s...', fn)
    keys = {'Sample': fn.split('/')[-1]}
    for row in csv.DictReader(open(fn, 'rt'), delimiter='\t'):
      # variation count prob
      if use_probability:
        keys[row['Variation']] = row['Probability']
      else:
        keys[row['Variation']] = row['Count']
      seen.add(row['Variation'])
    # generate result
    rows[fn] = keys

  sys.stdout.write('Sample\t{}\n'.format('\t'.join(sorted(seen))))

  ofh = csv.DictWriter(sys.stdout, delimiter='\t', fieldnames=['Sample'] + sorted(seen))
  ofh.writeheader()
  
  logging.info('writing...')
  for fn in files:
    logging.debug('writing %s...', fn)
    ofh.writerow(rows[fn])
  logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Combine signatures')
  parser.add_argument('--files', required=True, nargs='+', help='count files')
  parser.add_argument('--use_probability', action='store_true', help='probability instead of count')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.files, args.use_probability)
