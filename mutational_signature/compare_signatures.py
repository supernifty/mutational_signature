#!/usr/bin/env python
'''
  how similar are the input signatures?
'''

import argparse
import collections
import logging
import sys

import numpy as np
import scipy.stats

def cosine(x, y):
  similarity = x.dot(y) / (np.linalg.norm(x) * np.linalg.norm(y))
  return similarity

def pearson(x, y):
  similarity =  scipy.stats.pearsonr(x, y) # returns correlation, p-value
  return similarity[0]

def compare(signatures, signatures2, measure, out, no_convert):
  logging.info('processing %s...', signatures)

  names1 = []
  rows1 = []
  first = True
  for line in open(signatures, 'r'):
    fields = line.strip('\n').split('\t')
    if first:
      first = False
      if '>' in fields[1] or no_convert:
        signature_classes = fields[1:]
      else:
        signature_classes = ['{}{}{}>{}'.format(f[0], f[1], f[3], f[2]) for f in fields[1:]]
      logging.debug('%i signature classes in first dataset', len(signature_classes))
      continue
    names1.append(fields[0].replace('Signature.', ''))
    row = [float(x) for x in fields[1:]]
    rows1.append(np.array(row))
    logging.debug('row sum is %.2f', sum(rows1[-1]))

  logging.debug('%i signatures in first dataset', len(rows1))

  if signatures2 is None:
    names2 = names1
    rows2 = rows1
  else:
    names2 = []
    rows2 = []
    first = True
    for line in open(signatures2, 'r'):
      fields = line.strip('\n').split('\t')
      if first:
        first = False
        if '>' in fields[1] or no_convert:
          signature_classes_2 = fields[1:]
        else:
          signature_classes_2 = ['{}{}{}>{}'.format(f[0], f[1], f[3], f[2]) for f in fields[1:]]
        logging.debug('%i signature classes in second dataset', len(signature_classes_2))
        continue
      names2.append(fields[0].replace('Signature.', ''))
      #row = [float(x) for x in fields[1:]]
      #row = [float(fields[signature_classes_2.index(signature_classes[idx])]) for idx in range(1, len(fields))]
      row = []
      for idx in range(1, len(fields)):
        try:
          jdx = signature_classes_2.index(signature_classes[idx - 1]) # matching context
          row.append(float(fields[jdx + 1]))
        except IndexError:
          logging.debug('skipping idx %i from %i', idx-1, len(signature_classes))
        except ValueError:
          logging.debug('skipping idx %i from %i', idx-1, len(signature_classes))

      rows2.append(np.array(row))
      logging.debug('row sum is %.2f', sum(rows2[-1]))

    logging.debug('%i signatures in second dataset', len(rows2))

  # pairwise distance
  similarities = np.zeros((len(names1),len(names2)), dtype=float)
  for aid, arow in enumerate(rows1):
    for bid, brow in enumerate(rows2):
      logging.debug('%s vs %s: %i vs %i len', aid, bid, len(arow), len(brow))
      if measure == 'pearson':
        similarity = pearson(arow, brow)
      elif measure == 'cosine':
        similarity = cosine(arow, brow)
      else:
        logging.warn('unrecognised measure %s', measure)
        similarity = cosine(arow, brow)
      similarities[aid][bid] = similarity

  # print all
  sys.stdout.write('Sig\t{}\tSimilar Signatures\n'.format('\t'.join(names2)))
  for row, name in enumerate(names1):
    best = np.argsort(similarities[row])[::-1]
    similar = '{} ({:.2f}) {} {}'.format(names2[best[0]], similarities[row][best[0]], names2[best[1]], names2[best[2]])
    sys.stdout.write('{}\t{}\t{}\n'.format(name, '\t'.join(['{:.4f}'.format(x) for x in similarities[row]]), similar))

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='compare signatures')
  parser.add_argument('--signatures', required=True, help='signatures')
  parser.add_argument('--signatures2', required=False, help='second group of signatures')
  parser.add_argument('--measure', required=False, default='cosine', help='cosine or pearson')
  parser.add_argument('--no_convert', action='store_true', help='do not try to convert column names')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)
  compare(args.signatures, args.signatures2, args.measure, sys.stdout, args.no_convert)
