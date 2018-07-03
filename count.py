#!/usr/bin/env python

import argparse
import collections
import logging
import sys

import cyvcf2

COMP = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

def normalize(v):
  '''
    input GAT>G => ATC>C
  '''
  if v[1] in ('C', 'T'):
    return v # ok
  else:
    return ''.join([COMP[v[2]], COMP[v[1]], COMP[v[0]], '>', COMP[v[4]]])
    
def update_chroms(required, chroms, genome, next_chrom):
  '''
    pull the entire chromosome into memory
  '''
  seq = []
  for line in genome:
    line = line.strip('\n')
    if line.startswith('>'):
      if next_chrom is None: # first line of file
        next_chrom = line[1:].split(' ')[0]
        logging.debug('reading chrom %s from genome...', next_chrom)
      else:
        # remove previous chromosomes
        chroms[next_chrom] = ''.join(seq)
        seq = []
        logging.debug('reading chrom %s from genome. size is %i: done', next_chrom, len(chroms[next_chrom]))
        if required == next_chrom:
          next_chrom = line[1:].split(' ')[0]
          logging.debug('required chrom %s found', next_chrom)
          return next_chrom
        else:
          next_chrom = line[1:].split(' ')[0]
          logging.debug('reading chrom %s from genome...', next_chrom)
    else:
      seq.append(line)

  # end of file
  chroms[next_chrom] = ''.join(seq)
  logging.debug('reading chrom %s from genome. size is %i: done', next_chrom, len(chroms[next_chrom]))
  return None

def update_counts(counts, variant, chroms):
  if variant.POS == 1 or variant.POS > len(chroms[variant.CHROM]):
    logging.info('skipped edge variant at %s:%i', variant.CHROM, variant.POS)
    return
  if len(variant.REF) != 1 or len(variant.ALT[0]) != 1:
    logging.info('skipped indel at %s:%i', variant.CHROM, variant.POS)
    return

  # 0 1 2 -> my indexes
  # 1 2 3 -> vcf indexes
  fragment = chroms[variant.CHROM][variant.POS - 2:variant.POS + 1].upper() # TODO should we instead skip lower case
  if fragment[1] != variant.REF:
    logging.warn('skipping variant with position mismatch at %s:%i: VCF: %s genome: %s[%s]%s', variant.CHROM, variant.POS, variant.REF, fragment[0], fragment[1], fragment[2])
    return

  if any([x not in 'ACGT' for x in ''.join([fragment, variant.ALT[0]])]):
    logging.warn('skipping variant with ill-defined transition {}>{} at {}:{}'.format(fragment, variant.ALT[0], variant.CHROM, variant.POS))
    return
    
  v = '{}>{}'.format(fragment, variant.ALT[0]) # TODO multiallele
  v = normalize(v)
  counts[v] += 1

def count(genome_fh, vcf, out):
  logging.info('processing %s...', vcf)
  chroms = {}
  chroms_seen = set()
  counts = collections.defaultdict(int)
  next_chrom = None
  for line, variant in enumerate(cyvcf2.VCF(vcf)):

    if variant.CHROM not in chroms_seen:
      logging.debug('chrom %s seen in vcf', variant.CHROM)
      chroms = {} # wipe previous chromosomes
      next_chrom = update_chroms(variant.CHROM, chroms, genome_fh, next_chrom)
      chroms_seen.add(variant.CHROM)

    update_counts(counts, variant, chroms)

    if (line + 1) % 100000 == 0:
      logging.debug('processed %i lines. current counts: %s...', line + 1, ' '.join(['{}:{}'.format(k, counts[k]) for k in counts]))
    
  logging.info('processing %s: done', vcf)

  # write out results
  total_count = sum([counts[v] for v in counts])
  out.write('{}\t{}\t{}\n'.format('Variation', 'Count', 'Probability'))
  for k in sorted(counts):
    out.write('{}\t{}\t{:.3f}\n'.format(k, counts[k], counts[k] / total_count))
  # add zero results
  for ref in ('C', 'T'):
    for alt in ('A', 'C', 'G', 'T'):
      if ref == alt:
        continue
      for prefix in ('A', 'C', 'G', 'T'):
        for suffix in ('A', 'C', 'G', 'T'):
          count = '{}{}{}>{}'.format(prefix, ref, suffix, alt)
          if count not in counts:
            out.write('{}\t{}\t{:.3f}\n'.format(count, 0, 0))

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='mutational signature counter')
  parser.add_argument('--genome', required=True, help='reference genome')
  parser.add_argument('--vcf', required=True, help='vcf')
  args = parser.parse_args()
  logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  count(open(args.genome, 'r'), args.vcf, sys.stdout)
