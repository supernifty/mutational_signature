#!/usr/bin/env python
'''
  combine signature outputs
'''

import argparse
import logging
import sys

def main(signatures, files):
  logging.info('starting...')

  sigs = []
  first = True
  for line in open(signatures, 'r'):
    if first:
      first = False
      continue
    sigs.append(line.strip('\n').split('\t')[0])

  sys.stdout.write('Sample\t{}\n'.format('\t'.join(sigs)))
  
  for file in files:
    #result = [file.split('.')[0]]
    result = [file.split('/')[-1]]
    vals = {}
    for line in open(file, 'r'):
      sig, val = line.strip('\n').split('\t')
      vals[sig] = val
    # generate result
    for sig in sigs:
      result.append(vals[sig])
    sys.stdout.write('{}\n'.format('\t'.join(result)))

  logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Combine signatures')
  parser.add_argument('--signatures', required=True, help='signature files')
  parser.add_argument('--files', required=True, nargs='+', help='signature files')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.signatures, args.files)
