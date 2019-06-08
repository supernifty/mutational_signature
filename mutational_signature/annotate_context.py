#!/usr/bin/env python
'''
  adds snv_context to a vcf
'''

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

def context(variant, chroms):
  if variant.POS == 1 or variant.POS > len(chroms[variant.CHROM]):
    logging.info('skipped edge variant at %s:%i', variant.CHROM, variant.POS)
    return None
  if len(variant.REF) != 1 or len(variant.ALT[0]) != 1:
    logging.debug('skipped indel at %s:%i', variant.CHROM, variant.POS)
    return None

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
  return v

def annotate(genome_fh, vcf, out=None, chroms=None, variant_filter=None):
  logging.info('processing %s...', vcf)

  if chroms is None:
    chroms = {}
    chroms_seen = set()
  else:
    chroms_seen = set(chroms.keys())
    logging.debug('using existing genome with %i chromosomes', len(chroms_seen))

  next_chrom = None
  filtered = 0

  vcf_in = cyvcf2.VCF(vcf)
  vcf_in.add_info_to_header({'ID': 'snv_context', 'Description': 'mutational signature trinucleotide context', 'Type':'Character', 'Number': '1'})
  sys.stdout.write(vcf_in.raw_header)

  for line, variant in enumerate(vcf_in):
    if variant.CHROM not in chroms_seen:
      logging.debug('chrom %s seen in vcf', variant.CHROM)
      next_chrom = update_chroms(variant.CHROM, chroms, genome_fh, next_chrom)
      chroms_seen.add(variant.CHROM)

    if variant_filter is not None and not variant_filter(vcf_in, variant):
      filtered += 1
      continue

    snv_context = context(variant, chroms)
    if snv_context is not None:
      variant.INFO["snv_context"] = snv_context
    sys.stdout.write(str(variant))

    if (line + 1) % 100000 == 0:
      logging.debug('processed %i lines...', line + 1)
    
  logging.info('processing %s: filtered %i. done', vcf, filtered)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='mutational signature counter')
  parser.add_argument('--genome', required=True, help='reference genome')
  parser.add_argument('--vcf', required=True, help='vcf')
  args = parser.parse_args()
  logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  annotate(open(args.genome, 'r'), args.vcf, sys.stdout)
