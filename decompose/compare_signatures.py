#!/usr/bin/env python
'''
  how similar are the input signatures?
'''

import argparse
import collections
import logging
import sys

import numpy as np

def cosine(x, y):
  similarity = x.dot(y) / (np.linalg.norm(x) * np.linalg.norm(y))
  return similarity

def compare(signatures, out):
  logging.info('processing %s...', signatures)

  names = []
  rows = []
  first = True
  for line in open(signatures, 'r'):
    fields = line.strip('\n').split('\t')
    if first:
      first = False
      if '>' in fields[1]:
        signature_classes = fields[1:]
      else:
        signature_classes = ['{}{}{}>{}'.format(f[0], f[1], f[3], f[2]) for f in fields[1:]]
      logging.debug('%i signature classes', len(signature_classes))
      continue
    names.append(fields[0].replace('Signature.', ''))
    row = [float(x) for x in fields[1:]]
    rows.append(np.array(row))

  logging.debug('%i rows', len(rows))

  # pairwise distance
  distances = np.zeros((len(names),len(names)), dtype=float)
  for aid, arow in enumerate(rows):
    for bid, brow in enumerate(rows):
      distance = cosine(arow, brow)
      distances[aid][bid] = distance

  # print all
  sys.stdout.write('Sig\t{}\tSimilar Signatures\n'.format('\t'.join(names)))
  for row, name in enumerate(names):
    best = np.argsort(distances[row])[::-1]
    similar = '{} ({:.2f}) {} {}'.format(names[best[1]], distances[row][best[2]], names[best[2]], names[best[3]])
    sys.stdout.write('{}\t{}\t{}\n'.format(name, '\t'.join(['{:.2f}'.format(x) for x in distances[row]]), similar))

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='compare signatures')
  parser.add_argument('--signatures', required=True, help='signatures')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)
  compare(args.signatures, sys.stdout)
