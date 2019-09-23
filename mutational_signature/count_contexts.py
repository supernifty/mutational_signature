#!/usr/bin/env python
'''
  count opportunities for a particular context to occur

  notes:
  * for indels, variants must be left normalized
'''

import argparse
import collections
import logging
import sys

import cyvcf2

COMP = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
INDEL_COMP = {'A': 'T', 'C': 'C', 'G': 'C', 'T': 'T'} # normalize to C or T


def normalize_sbs(v):
  '''
    input GAT>G => ATC>C
  '''
  if v[1] in ('C', 'T'):
    return v # ok
  else:
    return ''.join([COMP[v[2]], COMP[v[1]], COMP[v[0]], '>', COMP[v[4]]])

DOUBLETS = set(['AC>CA', 'AC>CG', 'AC>CT', 'AC>GA', 'AC>GG', 'AC>GT', 'AC>TA', 'AC>TG', 'AC>TT', 'AT>CA', 'AT>CC', 'AT>CG', 'AT>GA', 'AT>GC', 'AT>TA', 'CC>AA', 'CC>AG', 'CC>AT', 'CC>GA', 'CC>GG', 'CC>GT', 'CC>TA', 'CC>TG', 'CC>TT', 'CG>AT', 'CG>GC', 'CG>GT', 'CG>TA', 'CG>TC', 'CG>TT', 'CT>AA', 'CT>AC', 'CT>AG', 'CT>GA', 'CT>GC', 'CT>GG', 'CT>TA', 'CT>TC', 'CT>TG', 'GC>AA', 'GC>AG', 'GC>AT', 'GC>CA', 'GC>CG', 'GC>TA', 'TA>AT', 'TA>CG', 'TA>CT', 'TA>GC', 'TA>GG', 'TA>GT', 'TC>AA', 'TC>AG', 'TC>AT', 'TC>CA', 'TC>CG', 'TC>CT', 'TC>GA', 'TC>GG', 'TC>GT', 'TG>AA', 'TG>AC', 'TG>AT', 'TG>CA', 'TG>CC', 'TG>CT', 'TG>GA', 'TG>GC', 'TG>GT', 'TT>AA', 'TT>AC', 'TT>AG', 'TT>CA', 'TT>CC', 'TT>CG', 'TT>GA', 'TT>GC', 'TT>GG'])

INDELS = set(['DEL_C_1_0', 'DEL_C_1_1', 'DEL_C_1_2', 'DEL_C_1_3', 'DEL_C_1_4', 'DEL_C_1_5+', 'DEL_MH_2_1', 'DEL_MH_3_1', 'DEL_MH_3_2', 'DEL_MH_4_1', 'DEL_MH_4_2', 'DEL_MH_4_3', 'DEL_MH_5+_1', 'DEL_MH_5+_2', 'DEL_MH_5+_3', 'DEL_MH_5+_4', 'DEL_MH_5+_5+', 'DEL_T_1_0', 'DEL_T_1_1', 'DEL_T_1_2', 'DEL_T_1_3', 'DEL_T_1_4', 'DEL_T_1_5+', 'DEL_repeats_2_0', 'DEL_repeats_2_1', 'DEL_repeats_2_2', 'DEL_repeats_2_3', 'DEL_repeats_2_4', 'DEL_repeats_2_5+', 'DEL_repeats_3_0', 'DEL_repeats_3_1', 'DEL_repeats_3_2', 'DEL_repeats_3_3', 'DEL_repeats_3_4', 'DEL_repeats_3_5+', 'DEL_repeats_4_0', 'DEL_repeats_4_1', 'DEL_repeats_4_2', 'DEL_repeats_4_3', 'DEL_repeats_4_4', 'DEL_repeats_4_5+', 'DEL_repeats_5+_0', 'DEL_repeats_5+_1', 'DEL_repeats_5+_2', 'DEL_repeats_5+_3', 'DEL_repeats_5+_4', 'DEL_repeats_5+_5+', 'INS_C_1_0', 'INS_C_1_1', 'INS_C_1_2', 'INS_C_1_3', 'INS_C_1_4', 'INS_C_1_5+', 'INS_T_1_0', 'INS_T_1_1', 'INS_T_1_2', 'INS_T_1_3', 'INS_T_1_4', 'INS_T_1_5+', 'INS_repeats_2_0', 'INS_repeats_2_1', 'INS_repeats_2_2', 'INS_repeats_2_3', 'INS_repeats_2_4', 'INS_repeats_2_5+', 'INS_repeats_3_0', 'INS_repeats_3_1', 'INS_repeats_3_2', 'INS_repeats_3_3', 'INS_repeats_3_4', 'INS_repeats_3_5+', 'INS_repeats_4_0', 'INS_repeats_4_1', 'INS_repeats_4_2', 'INS_repeats_4_3', 'INS_repeats_4_4', 'INS_repeats_4_5+', 'INS_repeats_5+_0', 'INS_repeats_5+_1', 'INS_repeats_5+_2', 'INS_repeats_5+_3', 'INS_repeats_5+_4', 'INS_repeats_5+_5+'])

def normalize_doublet(v):
  if v not in DOUBLETS:
    # reverse it
    v = '{}{}>{}{}'.format(COMP[v[1]], COMP[v[0]], COMP[v[4]], COMP[v[3]])
    if v not in DOUBLETS:
      logging.warn('failed to solve doublet %s', v)
  return v

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
        logging.info('reading chrom %s from genome. size is %i: done', next_chrom, len(chroms[next_chrom]))
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
  logging.info('reading chrom %s from genome. size is %i: done', next_chrom, len(chroms[next_chrom]))
  return None

def update_counts(counts, chrom, pos, chroms, indels=False, just_indels=False, doublets=False):
  if pos == 1 or pos > len(chroms[chrom]):
    logging.info('skipped edge variant at %s:%i', chrom, pos)
    return

  if indels and len(variant.REF) != len(variant.ALT[0]):
    logging.warn('not supported')
    return

    category = {}

    del_length = len(variant.REF) - len(variant.ALT[0]) # -ve for insertions
    # indel can be INS or DEL
    if del_length > 0:
      category['category'] = 'DEL'
    else:
      category['category'] = 'INS'

    # length can be 1 2 3 4 5+
    if abs(del_length) >= 5:
      category['length'] = '5+'
    else:
      category['length'] = abs(del_length)

    # measure repeat length
    repeats = 0
    if del_length > 0: # deletion
      # deletions look like TAAA -> A POS is where the T is
      deleted_sequence = variant.REF[1:]
      # first look for ongoing repeated sequence. POS-1 puts us at the T
      current_pos = variant.POS - 1 + 1 + len(deleted_sequence)
      while chroms[variant.CHROM][current_pos:current_pos + len(deleted_sequence)] == deleted_sequence:
        current_pos += len(deleted_sequence)
        repeats += 1

      if del_length == 1:
        category['repeats'] = '5+' if repeats >= 5 else repeats
        category['content'] = INDEL_COMP[deleted_sequence]
      else: # del length > 1
        if repeats == 0:
          # measure microhomology if deletion is not repeat
          microhomology = 0
          # first look to the right
          for pos in range(0, len(deleted_sequence)):
            if chroms[variant.CHROM][variant.POS - 1 + 1 + len(deleted_sequence) + pos] == deleted_sequence[pos]:
              microhomology += 1
            else:
              break
          # now look to the left
          for pos in range(len(deleted_sequence) - 1, -1, -1):
            if chroms[variant.CHROM][variant.POS - 1 + 1 - len(deleted_sequence) + pos] == deleted_sequence[pos]:
              microhomology += 1
            else:
              break

          if microhomology >= del_length:
            logging.warn('Microhomology calculation went wrong at %s:%i %i', variant.CHROM, variant.POS, microhomology)
          
          if microhomology > 0:
            category['repeats'] = '5+' if microhomology >= 5 else microhomology
            category['content'] = 'MH'
          else:
            category['repeats'] = 0
            category['content'] = 'repeats'
        else:
          category['repeats'] = '5+' if repeats >= 5 else repeats
          category['content'] = 'repeats'

    else: # insertion
      # insertions look like G -> GCA. POS is where the G is
      inserted_sequence = variant.ALT[0][1:]
      current_pos = variant.POS - 1 + 1
      while chroms[variant.CHROM][current_pos:current_pos + len(inserted_sequence)] == inserted_sequence:
        current_pos += len(inserted_sequence)
        repeats += 1
      category['repeats'] = '5+' if repeats >= 5 else repeats
      if del_length == -1:
        category['content'] = INDEL_COMP[inserted_sequence]
      else: # insertion length > 1
        category['content'] = 'repeats'

    indel_category = '{category}_{content}_{length}_{repeats}'.format(**category)
    counts[indel_category] += 1
    #logging.debug('indel at %s:%s %s -> %s classified as %s', variant.CHROM, variant.POS, variant.REF, variant.ALT[0], indel_category)
    if indel_category not in INDELS:
      logging.warn('unexpected indel category %s', indel_category)

  # no need to look at this indel any more
  #if len(variant.REF) != 1 or len(variant.ALT[0]) != 1:
  #  if not indels:
  #    logging.debug('skipped indel at %s:%i', variant.CHROM, variant.POS)
  #  return

  if just_indels:
    return

  # --- sbs
  # 0 1 2 -> my indexes
  # 1 2 3 -> vcf indexes
  # pulling in -1 0 +1
  fragment = chroms[chrom][pos - 2:pos + 1].upper() # potentially could instead skip lower case
  if any([x not in 'ACGTacgt' for x in fragment]):
    return
  v = '{}>{}'.format(fragment, 'A') # G is dummy
  v = normalize_sbs(v)[:3] # just the 1st three
  counts[v] += 1

  # --- doublets
  if doublets:
    if last_variant is not None and last_variant.CHROM == variant.CHROM and last_variant.POS == variant.POS - 1:
      doublet = '{}{}>{}{}'.format(last_variant.REF, variant.REF, last_variant.ALT[0], variant.ALT[0])
      if len(doublet) != 5:
        logging.warn('skipping doublet %s at %s:%i', doublet, variant.CHROM, variant.POS)
      else:
        doublet = normalize_doublet(doublet)
        counts[doublet] += 1
        logging.debug('doublet found at %s:%s: %s', variant.CHROM, variant.POS, doublet)

def count(genome_fh, bed, out=None, chroms=None, just_indels=False, doublets=False, indels=False):
  logging.info('processing %s...', bed)

  if chroms is None:
    chroms = {}
    chroms_seen = set()
  else:
    chroms_seen = set(chroms.keys())
    logging.info('using existing genome with %i chromosomes', len(chroms_seen))

  counts = collections.defaultdict(int)
  next_chrom = None
  filtered = 0

  for idx, line in enumerate(open(bed, 'r')):
    fields = line.strip('\n').split('\t')
    if len(fields) < 3:
      continue
    chrom, start, finish = fields[:3]
    start = int(start)
    finish = int(finish)
    if chrom not in chroms_seen:
      logging.info('chrom %s seen in bed', chrom)
      next_chrom = update_chroms(chrom, chroms, genome_fh, next_chrom)
      chroms_seen.add(chrom)

    for pos in range(start, finish):
      update_counts(counts, chrom, pos, chroms)

    if (idx + 1) % 100 == 0:
      logging.debug('processed %i lines. current counts: %s...', idx + 1, ' '.join(['{}:{}'.format(k, counts[k]) for k in counts]))

  logging.info('processing %s: filtered %i. included %i. done', bed, filtered, sum([counts[k] for k in counts]))

  # write out results
  total_count = sum([counts[v] for v in counts])
  if out is not None:
    out.write('{}\t{}\t{}\n'.format('Variation', 'Count', 'Probability'))
    for k in sorted(counts):
      out.write('{}\t{}\t{:.3f}\n'.format(k, counts[k], counts[k] / total_count))

    # add zero results for SBS
    if not just_indels:
      for ref in ('C', 'T'):
        for prefix in ('A', 'C', 'G', 'T'):
          for suffix in ('A', 'C', 'G', 'T'):
            count = '{}{}{}'.format(prefix, ref, suffix)
            if count not in counts:
              out.write('{}\t{}\t{:.3f}\n'.format(count, 0, 0))

    # add zero results for doublets
    if doublets:
      for doublet in DOUBLETS:
        if doublet not in counts:
          out.write('{}\t{}\t{:.3f}\n'.format(doublet, 0, 0))

    if indels:
      for indel in INDELS:
        if indel not in counts:
          out.write('{}\t{}\t{:.3f}\n'.format(indel, 0, 0))


  return {'chroms': chroms, 'counts': counts, 'total': total_count}

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='mutational signature counter')
  parser.add_argument('--genome', required=True, help='reference genome')
  parser.add_argument('--bed', required=True, help='regions to consider')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)
  count(genome_fh=open(args.genome, 'r'), bed=args.bed, out=sys.stdout)
