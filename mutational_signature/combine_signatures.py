#!/usr/bin/env python
'''
  combine signature outputs
'''

import argparse
import logging
import sys

ADDITIONAL_FIELDS = ('Error', 'Mutations')

def main(signatures, files):
  logging.info('starting...')

  sigs = []
  first = True
  for line in open(signatures, 'r'):
    if first:
      first = False
      continue
    sigs.append(line.strip('\n').split('\t')[0])
  sigs.extend(ADDITIONAL_FIELDS)
  #logging.debug('sigs: %s', sigs)

  sys.stdout.write('Sample\t{}\n'.format('\t'.join([x.replace('Signature.', '') for x in sigs])))
  
  for file in files:
    logging.debug('processing %s...', file)
    #result = [file.split('.')[0]]
    #result = [file.split('/')[-1].split('.')[0]]
    result = [file.split('/')[-1]]
    vals = {}
    for line in open(file, 'r'):
      sig, val = line.strip('\n').split('\t')[:2]
      #logging.debug('sig: %s, val: %s', sig, val)
      vals[sig] = val
    # generate result
    for sig in sigs:
      if sig in vals:
        result.append(vals[sig])
      else:
        result.append('0') # default
    sys.stdout.write('{}\n'.format('\t'.join(result)))
    logging.debug('processing %s: done', file)

  logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Combine signatures')
  parser.add_argument('--signatures', required=True, help='signature definition')
  parser.add_argument('--files', required=True, nargs='+', help='output from decompose')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.signatures, args.files)
