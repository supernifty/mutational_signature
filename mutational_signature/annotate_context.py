#!/usr/bin/env python
'''
  adds snv_context to a vcf
'''

import argparse
import collections
import csv
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
        next_chrom = line[1:].split(' ')[0].replace('chr', '')
        logging.debug('reading chrom %s from genome...', next_chrom)
      else:
        # remove previous chromosomes
        chroms[next_chrom] = ''.join(seq)
        seq = []
        logging.debug('reading chrom %s from genome. size is %i: done', next_chrom, len(chroms[next_chrom]))
        if required == next_chrom:
          next_chrom = line[1:].split(' ')[0].replace('chr', '')
          logging.debug('required chrom %s found', next_chrom)
          return next_chrom
        else:
          next_chrom = line[1:].split(' ')[0].replace('chr', '')
          logging.debug('reading chrom %s from genome...', next_chrom)
    else:
      seq.append(line)

  # end of file
  chroms[next_chrom] = ''.join(seq)
  logging.debug('reading chrom %s from genome. size is %i: done', next_chrom, len(chroms[next_chrom]))
  return None

def context(variant, chroms):
  chrom = variant.CHROM.replace('chr', '')
  if chrom not in chroms:
    logging.info('skipping chromosome %s', chrom)
    return None
  if variant.POS == 1 or variant.POS > len(chroms[chrom]):
    logging.info('skipped edge variant at %s:%i', chrom, variant.POS)
    return None
  if len(variant.REF) != 1 or len(variant.ALT) == 0 or len(variant.ALT[0]) != 1:
    logging.debug('skipped indel at %s:%i', chrom, variant.POS)
    return None

  # 0 1 2 -> my indexes
  # 1 2 3 -> vcf indexes
  fragment = chroms[chrom][variant.POS - 2:variant.POS + 1].upper() # TODO should we instead skip lower case
  if fragment[1] != variant.REF:
    logging.warn('skipping variant with position mismatch at %s:%i: VCF: %s genome: %s[%s]%s', chrom, variant.POS, variant.REF, fragment[0], fragment[1], fragment[2])
    return

  if any([x not in 'ACGT' for x in ''.join([fragment, variant.ALT[0]])]):
    logging.warn('skipping variant with ill-defined transition {}>{} at {}:{}'.format(fragment, variant.ALT[0], chrom, variant.POS))
    return
    
  v = '{}>{}'.format(fragment, variant.ALT[0]) # TODO multiallele
  v = normalize(v)
  return v

def surrounding(variant, sequence, chroms):
  ''' note that we do not normalize anything here'''
  if sequence == 0:
    return None
  chrom = variant.CHROM.replace('chr', '')
  if chrom not in chroms:
    logging.info('skipping chromosome %s', chrom)
    return None
  if variant.POS < sequence or variant.POS > len(chroms[chrom]) - sequence:
    logging.info('skipped edge variant at %s:%i', chrom, variant.POS)
    return None

  # 0 1 2 -> my indexes
  # 1 2 3 -> vcf indexes
  fragment = chroms[chrom][variant.POS - sequence - 1:variant.POS + sequence + len(variant.REF) - 1].upper() # TODO should we instead skip lower case
  if fragment[sequence:sequence + len(variant.REF)] != variant.REF:
    logging.warn('skipping variant with position mismatch at %s:%i: VCF: %s genome: %s', chrom, variant.POS, variant.REF, fragment)
    return

  return fragment

def vcf_writer(out):
  def write_header(vcf_in):
    out.write(vcf_in.raw_header)

  def write_variant(variant, snv_context, surrounding_context):
    if snv_context is not None:
      variant.INFO["snv_context"] = snv_context
    if surrounding_context is not None:
      variant.INFO["surrounding"] = surrounding_context
    out.write(str(variant))

  def add_to_header(d):
    vcf_in.add_info_to_header(d)

  return {'write_header': write_header, 'write_variant': write_variant, 'add_to_header': add_to_header}

def annotate(genome_fh, vcf_in, out=None, chroms=None, variant_filter=None, sequence=0):
  logging.info('processing...')

  if chroms is None:
    chroms = {}
    chroms_seen = set()
  else:
    chroms_seen = set(chroms.keys())
    logging.debug('using existing genome with %i chromosomes', len(chroms_seen))

  next_chrom = None
  filtered = 0

  #vcf_in = cyvcf2.VCF(vcf)
  out['add_to_header']({'ID': 'snv_context', 'Description': 'mutational signature trinucleotide context', 'Type':'Character', 'Number': '1'})
  if sequence > 0:
    out['add_to_header']({'ID': 'surrounding', 'Description': 'reference sequence surrounding variant', 'Type':'Character', 'Number': '1'})
  out['write_header'](vcf_in)

  line = 0
  for line, variant in enumerate(vcf_in):
    chrom = variant.CHROM.replace('chr', '')
    if chrom not in chroms_seen:
      logging.debug('chrom %s seen in vcf', chrom)
      next_chrom = update_chroms(chrom, chroms, genome_fh, next_chrom)
      chroms_seen.add(chrom)

    if variant_filter is not None and not variant_filter(vcf_in, variant):
      filtered += 1
      continue

    snv_context = context(variant, chroms)
    surrounding_context = surrounding(variant, sequence, chroms)

    out['write_variant'](variant, snv_context, surrounding_context)

    if (line + 1) % 10000 == 0:
      logging.debug('processed %i lines...', line + 1)
    
  logging.info('processing: filtered %i of %i. done', filtered, line)

def maf_to_vcf(maf, chrom_col, pos_col, ref_col, alt_col):

  def maf_reader(reader):
    Variant = collections.namedtuple('Variant', 'CHROM POS REF ALT row')
  
    for line, row in enumerate(reader):
      if line % 1000 == 0:
        logging.debug('processed %i lines of %s...', line, maf)
  
  
      chrom = row[chrom_col]
      pos = int(row[pos_col])
      ref = row[ref_col]
      if ref == '-':
        pos += 1 # fix for TCGA mafs
      ref = ref.replace('-', '')
      alt = row[alt_col]
  
      yield Variant(chrom, pos, ref, (alt,), row)

  # enumeration a maf into a variant
  reader = csv.DictReader(open(maf, 'r'), delimiter='\t')
  return {'tsv_reader': reader, 'maf_reader': maf_reader(reader)}

def maf_writer(out): # DictWriter

  def write_header(vcf_in):
    out.writeheader()

  def write_variant(variant, snv_context, surrounding_context):
    if snv_context is not None:
      variant.row["snv_context"] = snv_context
    if surrounding_context is not None:
      variant.row["surrounding"] = surrounding_context
    out.writerow(variant.row)

  def dummy(d):
    pass

  return {'write_header': write_header, 'write_variant': write_variant, 'add_to_header': dummy}


if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='mutational signature counter')
  parser.add_argument('--genome', required=True, help='reference genome')
  parser.add_argument('--sequence', required=False, default=0, type=int, help='surrounding sequence in each direction to annotate, 0 to skip')
  parser.add_argument('--vcf', required=True, help='vcf')
  parser.add_argument('--is_maf', action='store_true', help='vcf is a maf')
  parser.add_argument('--maf_chrom_column', required=False, default='Chromosome', help='maf chrom column name')
  parser.add_argument('--maf_pos_column', required=False, default='Start_Position', help='maf pos column name')
  parser.add_argument('--maf_ref_column', required=False, default='Reference_Allele', help='maf ref column name')
  parser.add_argument('--maf_alt_column', required=False, default='Tumor_Seq_Allele2', help='maf alt column name')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  if args.is_maf:
    vcf_in = maf_to_vcf(args.vcf, args.maf_chrom_column, args.maf_pos_column, args.maf_ref_column, args.maf_alt_column)
    writer = csv.DictWriter(sys.stdout, delimiter='\t', fieldnames=vcf_in['tsv_reader'].fieldnames + ['snv_context', 'surrounding'])
    vcf_out = maf_writer(writer)
    annotate(open(args.genome, 'r'), vcf_in['maf_reader'], vcf_out, sequence=args.sequence)
  else:
    vcf_in = cyvcf2.VCF(args.vcf)
    annotate(open(args.genome, 'r'), vcf_in, vcf_writer(sys.stdout), sequence=args.sequence)
