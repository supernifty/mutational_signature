#!/usr/bin/env python
'''
  predict lynch pathway
'''

import argparse
import collections
import csv
import logging
import sys

import numpy as np

def cosine(x, y):
  x = np.array(x)
  y = np.array(y)
  similarity = x.dot(y) / (np.linalg.norm(x) * np.linalg.norm(y))
  return similarity

def compare(x, y):
  keys = sorted([a for a in x.keys() if a != 'Sig'])
  v1 = [x[z] for z in keys]
  v2 = [y[z] for z in keys]
  return cosine(v1, v2)

def merge_signatures(v1, v2, weights, merge_strategy):
  if merge_strategy == 'first':
    # just return the first
    return v1
  else:
    # average them
    return {k: (weights[0] * v1[k] + weights[1] * v2[k]) / sum(weights) for k in v1.keys()}

def most_similar_pair(vals):
  best = (None, None, 0)
  names = sorted(vals.keys())
  logging.info('%i signatures to compare...', len(names))
  for i1 in range(len(vals)):
    for i2 in range(i1, len(vals)):
      if i1 == i2:
        continue
      s = compare(vals[names[i1]], vals[names[i2]])
      if s > best[2]:
        best = (names[i1], names[i2], s)
  return best

def main(ifh, ofh, max_similarity, merge_strategy):
  logging.info('starting...')
  idr = csv.DictReader(ifh, delimiter='\t')
  vals = {}
  for r in idr:
    name = r['Sig']
    vals[name] = {x: float(r[x]) for x in r if x != 'Sig'}

  logging.info('starting with %i signatures', len(vals))

  while True:
    pp = most_similar_pair(vals) # (name, name, sim)
    if pp[2] > max_similarity:
      vals['{}-{}'.format(pp[0], pp[1])] = merge_signatures(vals[pp[0]], vals[pp[1]], weights=(1 + pp[0].count('-'), 1 + pp[1].count('-')), merge_strategy=merge_strategy)
      del vals[pp[0]]
      del vals[pp[1]]
      logging.info('merged %s and %s with similarity %.3f. signatures is now %i', pp[0], pp[1], pp[2], len(vals))
    else:
      break # done

  # write them out
  odw = csv.DictWriter(ofh, delimiter='\t', fieldnames=['Sig'] + list(vals[list(vals)[0]].keys()))
  odw.writeheader()
  for v in vals:
    r = vals[v]
    r['Sig'] = v
    odw.writerow(r)

  logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Assess MSI')
  parser.add_argument('--max_similarity', type=float, default=0.8, required=False, help='max allowed similarity in final set')
  parser.add_argument('--merge_strategy', required=False, default='average', help='how to merge. average or first')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(sys.stdin, sys.stdout, args.max_similarity, args.merge_strategy)
