#!/usr/bin/env python
'''
  find the main context for each signature
'''

import argparse
import csv
import logging
import sys

import numpy as np

def convert(context):
  if len(context) == 4 and all([x in 'ACTG' for x in context]):
    return '{lhs}{ref}{rhs}>{lhs}{alt}{rhs}'.format(lhs=context[0], rhs=context[3], ref=context[1], alt=context[2])
  else:
    return context

def main(fh, out):
  logging.info('reading from stdin...')
  
  first = True
  out.write('Sig\tContext\n')
  for row in csv.reader(fh, delimiter='\t'):
    if first:
      header = row
      first = False
    else:
      context = header[np.argmax([float(x) for x in row[1:]]) + 1]
      out.write('{}\t{}\n'.format(row[0], convert(context)))
      
  logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Find context')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(sys.stdin, sys.stdout)
