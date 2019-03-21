#!/usr/bin/env python

import argparse
import csv
import collections
import gzip
import logging
import sys

import cyvcf2
import intervaltree

def rotate(r):
  '''
    normalize repeats to start from earliest base e.g. GCAA -> AAGC
  '''
  best = 0
  for idx, el in enumerate(r):
    if r[idx] < r[best]:
      best = idx

  return r[best:] + r[:best]

def complement(r):
  rc = reverse_complement(r)
  if r <= rc:
    return r
  else:
    return rc

RC = {'a': 't', 'c': 'g', 't': 'a', 'g': 'c', 'n': 'n', 'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C', 'N': 'N'}
def reverse_complement(repeat):
  return ''.join([RC[x] for x in repeat][::-1])

def get_value(header, col, row):
  return row[header.index(col)]

def make_tree(bed):
  tree = {}

  for idx, line in enumerate(open(bed, 'r')):
    chrom, start, finish, annotation = line.strip('\n').split('\t')
    if chrom not in tree:
      tree[chrom] = intervaltree.IntervalTree()
    tree[chrom][int(start):int(finish)] = annotation
    if idx % 10000 == 0:
      logging.debug('parsing {}: {} lines parsed'.format(bed, idx))

  logging.info('parsing %s: done.', bed)
  return tree

def main(maf, sample, bed, comp, rot, out_fh, maxlen, tree=None):

  if tree is None:
    logging.info('parsing %s...', bed)
    tree = make_tree(bed)

  logging.info('parsing %s...', maf)
  reject = accept = 0

  counts = collections.defaultdict(int)
  header = None
  for line, row in enumerate(csv.reader(gzip.open(maf, 'rt'), delimiter='\t')):
    if row[0].startswith('#'):
      continue
    if header is None:
      header = row
      continue

    row_sample = get_value(header, "Tumor_Sample_Barcode", row)
    if sample is not None and row_sample != sample:
      continue

    chrom = get_value(header, "Chromosome", row)
    pos = int(get_value(header, "Start_Position", row))
    ref = get_value(header, "Reference_Allele", row).replace('-', '')
    alt = get_value(header, "Tumor_Seq_Allele2", row).replace('-', '')

    # does it overlap the bed?
    if chrom in tree:
      overlap = tree[chrom].at(pos)
      if len(overlap) == 0:
        reject += 1
      else:
        accept += 1
        annotations = {}
        for namevalue in list(overlap)[0].data.split(';'):
          name, value = namevalue.split('=')
          annotations[name] = value

        indel_len = len(alt) - len(ref)

        if indel_len > maxlen:
          indel_len = maxlen
        if indel_len < -maxlen:
          indel_len = -maxlen

        # if we want only indels can filter on len here
        context = annotations['repeat']
        if comp:
          context = complement(context)
        if rot:
          context = rotate(context)

        counts['{}{}'.format(context, indel_len)] += 1
        
    else:
      reject += 1

  logging.info('included %i variants. rejected %i variants.', accept, reject)

  # write results
  out_fh.write('Type\tCount\tProportion\n')
  for count in sorted(counts):
    out_fh.write('{}\t{}\t{:.3f}\n'.format(count, counts[count], counts[count] / accept))

  logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='mutational signature counter')
  parser.add_argument('--annotation', required=True, help='bed')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  parser.add_argument('--complement', action='store_true', help='ignore direction')
  parser.add_argument('--rotate', action='store_true', help='merge same repeat types')
  parser.add_argument('--maf', required=True, help='maf')
  parser.add_argument('--sample', required=True, help='sample')
  parser.add_argument('--maxlen', type=int, default=1000, help='max len indel')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.maf, args.sample, args.annotation, args.complement, args.rotate, sys.stdout, args.maxlen)
