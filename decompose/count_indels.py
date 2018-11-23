#!/usr/bin/env python

import argparse
import collections
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

def main(bed, comp, rot):
  logging.info('parsing {}...'.format(bed))

  tree = {}

  for idx, line in enumerate(open(bed, 'r')):
    chrom, start, finish, annotation = line.strip('\n').split('\t')
    if chrom not in tree:
      tree[chrom] = intervaltree.IntervalTree()
    tree[chrom][int(start):int(finish)] = annotation
    if idx % 10000 == 0:
      logging.debug('parsing {}: {} lines parsed'.format(bed, idx))
  logging.info('parsing %s: done.', bed)

  logging.info('parsing vcf from stdin...')
  vcf_in = cyvcf2.VCF('-')
  reject = accept = 0

  counts = collections.defaultdict(int)
  for variant in vcf_in:
    # does it overlap the bed?
    if variant.CHROM in tree:
      overlap = tree[variant.CHROM].search(variant.POS)
      if len(overlap) == 0:
        reject += 1
      else:
        accept += 1
        annotations = {}
        for namevalue in list(overlap)[0].data.split(';'):
          name, value = namevalue.split('=')
          annotations[name] = value

        indel_len = len(variant.ALT[0].replace('-', '')) - len(variant.REF.replace('-', '')) # maf uses -
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

  # write results - TODO strand conversion, repeat rotation
  sys.stdout.write('Type\tCount\tProportion\n')
  for count in sorted(counts):
    sys.stdout.write('{}\t{}\t{:.3f}\n'.format(count, counts[count], counts[count] / accept))

  logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='mutational signature counter')
  parser.add_argument('--annotation', required=True, help='bed')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  parser.add_argument('--complement', action='store_true', help='ignore direction')
  parser.add_argument('--rotate', action='store_true', help='merge same repeat types')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.annotation, args.complement, args.rotate)
