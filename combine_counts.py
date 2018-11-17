#!/usr/bin/env python
'''
  combine signature outputs
'''

import argparse
import logging
import sys

def main(files):
  logging.info('starting...')

  rows = {}
  seen = set()
  for file in files:
    logging.info('processing %s...', file)
    keys = {}
    first = True
    for line in open(file, 'r'):
      if first:
        first = False
        continue
      # variation count prob
      key, count, probability = line.strip('\n').split('\t')
      keys[key] = probability
      seen.add(key)
    # generate result
    rows[file] = keys

  sys.stdout.write('Sample\t{}\n'.format('\t'.join(sorted(seen))))
  
  for file in files:
    result = [file.split('/')[-1].split('.')[0]]
    
    for key in sorted(seen):
      result.append(rows[file][key])

    sys.stdout.write('{}\n'.format('\t'.join(result)))

  logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Combine signatures')
  parser.add_argument('--files', required=True, nargs='+', help='count files')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.files)
