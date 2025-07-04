#!/usr/bin/env python
'''
  count variants for the purpose of calculating mutational signatures

  notes:
  * for indels, variants must be left normalized
'''

import argparse
import csv
import collections
import gzip
import itertools
import logging
import sys

import cyvcf2
import intervaltree
import scipy.stats

COMP = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
COMP_TX = {'+': '-', '-': '+', None: None}
INDEL_COMP = {'A': 'T', 'C': 'C', 'G': 'C', 'T': 'T'} # normalize to C or T

EXCLUDE_UTR=True # transcription bias
SKIP_ALT_CHROM=True
SKIP_M=True
SKIP_M=False

def normalize_sbs(v, strand_tx, strand_exon):
  '''
    input GAT>G => ATC>C
  '''
  mid = int((len(v) - 3) / 2)
  if v[mid] in ('C', 'T'):
    return v, strand_tx, strand_exon # ok
  else:
    #return ''.join([COMP[v[2]], COMP[v[1]], COMP[v[0]], '>', COMP[v[-1]]]), COMP_TX[strand_tx], COMP_TX[strand_exon]
    return ''.join([''.join([COMP[x] for x in v[:-2]][::-1]), '>', COMP[v[-1]]]), COMP_TX[strand_tx], COMP_TX[strand_exon]

DOUBLETS = set(['AC>CA', 'AC>CG', 'AC>CT', 'AC>GA', 'AC>GG', 'AC>GT', 'AC>TA', 'AC>TG', 'AC>TT', 'AT>CA', 'AT>CC', 'AT>CG', 'AT>GA', 'AT>GC', 'AT>TA', 'CC>AA', 'CC>AG', 'CC>AT', 'CC>GA', 'CC>GG', 'CC>GT', 'CC>TA', 'CC>TG', 'CC>TT', 'CG>AT', 'CG>GC', 'CG>GT', 'CG>TA', 'CG>TC', 'CG>TT', 'CT>AA', 'CT>AC', 'CT>AG', 'CT>GA', 'CT>GC', 'CT>GG', 'CT>TA', 'CT>TC', 'CT>TG', 'GC>AA', 'GC>AG', 'GC>AT', 'GC>CA', 'GC>CG', 'GC>TA', 'TA>AT', 'TA>CG', 'TA>CT', 'TA>GC', 'TA>GG', 'TA>GT', 'TC>AA', 'TC>AG', 'TC>AT', 'TC>CA', 'TC>CG', 'TC>CT', 'TC>GA', 'TC>GG', 'TC>GT', 'TG>AA', 'TG>AC', 'TG>AT', 'TG>CA', 'TG>CC', 'TG>CT', 'TG>GA', 'TG>GC', 'TG>GT', 'TT>AA', 'TT>AC', 'TT>AG', 'TT>CA', 'TT>CC', 'TT>CG', 'TT>GA', 'TT>GC', 'TT>GG'])

INDELS_83 = set(['DEL_C_1_0', 'DEL_C_1_1', 'DEL_C_1_2', 'DEL_C_1_3', 'DEL_C_1_4', 'DEL_C_1_5+', 'DEL_MH_2_1', 'DEL_MH_3_1', 'DEL_MH_3_2', 'DEL_MH_4_1', 'DEL_MH_4_2', 'DEL_MH_4_3', 'DEL_MH_5+_1', 'DEL_MH_5+_2', 'DEL_MH_5+_3', 'DEL_MH_5+_4', 'DEL_MH_5+_5+', 'DEL_T_1_0', 'DEL_T_1_1', 'DEL_T_1_2', 'DEL_T_1_3', 'DEL_T_1_4', 'DEL_T_1_5+', 'DEL_repeats_2_0', 'DEL_repeats_2_1', 'DEL_repeats_2_2', 'DEL_repeats_2_3', 'DEL_repeats_2_4', 'DEL_repeats_2_5+', 'DEL_repeats_3_0', 'DEL_repeats_3_1', 'DEL_repeats_3_2', 'DEL_repeats_3_3', 'DEL_repeats_3_4', 'DEL_repeats_3_5+', 'DEL_repeats_4_0', 'DEL_repeats_4_1', 'DEL_repeats_4_2', 'DEL_repeats_4_3', 'DEL_repeats_4_4', 'DEL_repeats_4_5+', 'DEL_repeats_5+_0', 'DEL_repeats_5+_1', 'DEL_repeats_5+_2', 'DEL_repeats_5+_3', 'DEL_repeats_5+_4', 'DEL_repeats_5+_5+', 'INS_C_1_0', 'INS_C_1_1', 'INS_C_1_2', 'INS_C_1_3', 'INS_C_1_4', 'INS_C_1_5+', 'INS_T_1_0', 'INS_T_1_1', 'INS_T_1_2', 'INS_T_1_3', 'INS_T_1_4', 'INS_T_1_5+', 'INS_repeats_2_0', 'INS_repeats_2_1', 'INS_repeats_2_2', 'INS_repeats_2_3', 'INS_repeats_2_4', 'INS_repeats_2_5+', 'INS_repeats_3_0', 'INS_repeats_3_1', 'INS_repeats_3_2', 'INS_repeats_3_3', 'INS_repeats_3_4', 'INS_repeats_3_5+', 'INS_repeats_4_0', 'INS_repeats_4_1', 'INS_repeats_4_2', 'INS_repeats_4_3', 'INS_repeats_4_4', 'INS_repeats_4_5+', 'INS_repeats_5+_0', 'INS_repeats_5+_1', 'INS_repeats_5+_2', 'INS_repeats_5+_3', 'INS_repeats_5+_4', 'INS_repeats_5+_5+'])

# Koh 2025 https://www.nature.com/articles/s41588-025-02152-y see supp table 8
INDELS_89 = set([
  'A[Ins(C):R0]A',
  'A[Ins(C):R0]T',
  'Ins(C):R(0,3)',
  'Ins(C):R(4,6)',
  'Ins(C):R(7,9)',
  'A[Ins(T):R(0,4)]A',
  'A[Ins(T):R(0,4)]C',
  'A[Ins(T):R(0,4)]G',
  'C[Ins(T):R(0,4)]A',
  'C[Ins(T):R(0,4)]C',
  'C[Ins(T):R(0,4)]G',
  'G[Ins(T):R(0,4)]A',
  'G[Ins(T):R(0,4)]C',
  'G[Ins(T):R(0,4)]G',
  'A[Ins(T):R(5,7)]A',
  'A[Ins(T):R(5,7)]C',
  'A[Ins(T):R(5,7)]G',
  'C[Ins(T):R(5,7)]A',
  'C[Ins(T):R(5,7)]C',
  'C[Ins(T):R(5,7)]G',
  'G[Ins(T):R(5,7)]A',
  'G[Ins(T):R(5,7)]C',
  'G[Ins(T):R(5,7)]G',
  'A[Ins(T):R(8,9)]A',
  'A[Ins(T):R(8,9)]C',
  'A[Ins(T):R(8,9)]G',
  'C[Ins(T):R(8,9)]A',
  'C[Ins(T):R(8,9)]C',
  'C[Ins(T):R(8,9)]G',
  'G[Ins(T):R(8,9)]A',
  'G[Ins(T):R(8,9)]C',
  'G[Ins(T):R(8,9)]G',
  'Ins(2,4):R0',
  'Ins(5,):R0',
  'Ins(2,4):R1',
  'Ins(5,):R1',
  'Ins(2,):R(2,4)',
  'Ins(2,):R(5,9)',
  '[Del(C):R1]A',
  '[Del(C):R1]T',
  '[Del(C):R2]A',
  '[Del(C):R2]T',
  '[Del(C):R3]A',
  '[Del(C):R3]T',
  '[Del(C):R(4,5)]A',
  '[Del(C):R(4,5)]T',
  '[Del(C):R(1,5)]G',
  'Del(C):R(6,9)',
  'A[Del(T):R(1,4)]A',
  'A[Del(T):R(1,4)]C',
  'A[Del(T):R(1,4)]G',
  'C[Del(T):R(1,4)]A',
  'C[Del(T):R(1,4)]C',
  'C[Del(T):R(1,4)]G',
  'G[Del(T):R(1,4)]A',
  'G[Del(T):R(1,4)]C',
  'G[Del(T):R(1,4)]G',
  'A[Del(T):R(5,7)]A',
  'A[Del(T):R(5,7)]C',
  'A[Del(T):R(5,7)]G',
  'C[Del(T):R(5,7)]A',
  'C[Del(T):R(5,7)]C',
  'C[Del(T):R(5,7)]G',
  'G[Del(T):R(5,7)]A',
  'G[Del(T):R(5,7)]C',
  'G[Del(T):R(5,7)]G',
  'A[Del(T):R(8,9)]A',
  'A[Del(T):R(8,9)]C',
  'A[Del(T):R(8,9)]G',
  'C[Del(T):R(8,9)]A',
  'C[Del(T):R(8,9)]C',
  'C[Del(T):R(8,9)]G',
  'G[Del(T):R(8,9)]A',
  'G[Del(T):R(8,9)]C',
  'G[Del(T):R(8,9)]G',
  'Del(2,4):R1',
  'Del(5,):R1',
  'Del(2,8):U(1,2):R(2,4)',
  'Del(2,):U(1,2):R(5,9)',
  'Del(3,):U(3,):R2',
  'Del(3,):U(3,):R(3,9)',
  'Del(2,5):M1',
  'Del(3,5):M2',
  'Del(4,5):M(3,4)',
  'Del(6,):M1',
  'Del(6,):M2',
  'Del(6,):M3',
  'Del(6,):M(4,)',
  'Complex'
])

def normalize_doublet(v):
  if v not in DOUBLETS:
    # reverse it
    try:
      v = '{}{}>{}{}'.format(COMP[v[1]], COMP[v[0]], COMP[v[4]], COMP[v[3]])
    except:
      logging.warning('something weird with doublet %s', v)

    if v not in DOUBLETS:
      logging.warning('failed to solve doublet %s', v)
  return v

def update_chroms(required, chroms, genome, next_chrom):
  '''
    pull the entire chromosome into memory
  '''
  seq = []
  for linenum, line in enumerate(genome):
    line = line.strip('\n')
    if line.startswith('>'):
      if next_chrom is None: # first line of file
        next_chrom = no_chr(line[1:].split(' ')[0])
        logging.debug('reading chrom %s from genome...', next_chrom)
      else:
        # remove previous chromosomes
        chroms[next_chrom] = ''.join(seq)
        seq = []
        logging.debug('reading chrom %s from genome. size is %i: done', next_chrom, len(chroms[next_chrom]))
        if required == next_chrom:
          next_chrom = no_chr(line[1:].split(' ')[0])
          logging.debug('required chrom %s found', next_chrom)
          return next_chrom
        else:
          next_chrom = no_chr(line[1:].split(' ')[0])
          logging.debug('reading chrom %s from genome...', next_chrom)
    else:
      seq.append(line)
    if linenum % 1000000 == 0:
      logging.debug('processed %i lines of genome...', linenum)

  # end of file
  chroms[next_chrom] = ''.join(seq)
  logging.info('reading chrom %s from genome. size is %i: done', next_chrom, len(chroms[next_chrom]))
  return None

def no_chr(chrom):
  if chrom == 'MT':
    return 'M'
  # deals with chrUn_gl000220.1
  return chrom.split('.')[0].split('_', maxsplit=1)[-1].replace('chr', '').upper()

def get_indel_category(variant, chroms):
  category = {}
  # tcga/intogen put - in the alt
  alt = variant.ALT[0]
  if alt == '-':
    alt = ''
  del_length = len(variant.REF) - len(alt) # -ve for insertions
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
    if len(alt) == 0:
      # deletions look like TA -> - e.g. in TCGA
      deleted_sequence = variant.REF
    else:  
      # deletions look like TAAA -> T POS is where the T is
      deleted_sequence = variant.REF[len(alt):]
    logging.debug('deleted sequence is %s from %s to %s', deleted_sequence, variant.REF, alt)
    # first look for ongoing repeated sequence. POS-1 puts us at the T
    current_pos = variant.POS + len(alt) - 1 + len(deleted_sequence)
    while chroms[no_chr(variant.CHROM)][current_pos:current_pos + len(deleted_sequence)] == deleted_sequence:
      current_pos += len(deleted_sequence)
      repeats += 1
      if repeats > 10000:
        break

    if del_length == 1:
      category['repeats'] = '5+' if repeats >= 5 else repeats
      category['content'] = INDEL_COMP[deleted_sequence]
    else: # del length > 1
      if repeats == 0:
        # measure microhomology if deletion is not repeat
        microhomology = 0
        # first look to the right
        for pos in range(0, len(deleted_sequence)):
          if chroms[no_chr(variant.CHROM)][variant.POS + len(alt) - 1 + len(deleted_sequence) + pos] == deleted_sequence[pos]:
            microhomology += 1
          else:
            break
        # now look to the left
        for pos in range(len(deleted_sequence) - 1, -1, -1):
          if chroms[no_chr(variant.CHROM)][variant.POS + len(alt) - 1 - len(deleted_sequence) + pos] == deleted_sequence[pos]:
            microhomology += 1
          else:
            break

        if microhomology >= del_length:
          logging.warning('Microhomology calculation went wrong at %s:%i %i del length %i del sequence %s', no_chr(variant.CHROM), variant.POS, microhomology, del_length, deleted_sequence)
        
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
    if len(variant.REF) == 0:
      # ref used to be -
      inserted_sequence = alt
    else:
      # insertions look like G -> GCA. POS is where the G is
      inserted_sequence = alt[len(variant.REF):]
    current_pos = variant.POS + len(variant.REF) - 1
    while chroms[no_chr(variant.CHROM)][current_pos:current_pos + len(inserted_sequence)] == inserted_sequence:
      current_pos += len(inserted_sequence)
      repeats += 1
      if repeats > 10000:
        break
    category['repeats'] = '5+' if repeats >= 5 else repeats
    if del_length == -1:
      if inserted_sequence not in INDEL_COMP:
        logging.warn('unexpected inserted sequence %s', inserted_sequence)
        category['content'] = 'unknown'
      else:
        category['content'] = INDEL_COMP[inserted_sequence]
    else: # insertion length > 1
      category['content'] = 'repeats'

  return category

def update_counts(counts, variant, last_variant, chroms, doublets, indels, just_indels, transcripts=None, exons=None, tx_counts=None, mer=3, weight=None):
  if no_chr(variant.CHROM) not in chroms:
    logging.warning('chromosome %s not found in chroms: %s', no_chr(variant.CHROM), chroms.keys())
    return

  if variant.POS == 1 or variant.POS > len(chroms[no_chr(variant.CHROM)]):
    logging.debug('skipped edge variant at %s:%i', no_chr(variant.CHROM), variant.POS)
    return

  # tcga/intogen put - in the alt
  alt = variant.ALT[0]
  if alt == '-':
    alt = ''

  exon_strand = None
  tx_strand = None

  # determine strand if relevant
  if exons is not None and no_chr(variant.CHROM) in exons:
    intersections = exons[no_chr(variant.CHROM)][variant.POS]
    if len(intersections) == 0:
      logging.debug('%s:%i is not coding', no_chr(variant.CHROM), variant.POS)
    elif len(intersections) > 1:
      logging.debug('%i intersections at %s:%i', len(intersections), no_chr(variant.CHROM), variant.POS) # can't decide
      # pick the first
      exon_strand = list(intersections)[0][2]
      logging.debug('%i intersections at %s:%i: first is %s', len(intersections), no_chr(variant.CHROM), variant.POS, exon_strand) # can't decide
    else:
      exon_strand = list(intersections)[0][2]

  # determine transcript if relevant
  if transcripts is not None and no_chr(variant.CHROM) in transcripts:
    intersections = transcripts[no_chr(variant.CHROM)][variant.POS]
    if len(intersections) == 0:
      logging.debug('%s:%i is not coding', no_chr(variant.CHROM), variant.POS)
    elif len(intersections) > 1:
      logging.debug('%i multiple intersections at %s:%i', len(intersections), no_chr(variant.CHROM), variant.POS) # can't decide
      tx_strand = list(intersections)[0][2]
      logging.debug('%i multiple intersections at %s:%i: first is %s', len(intersections), no_chr(variant.CHROM), variant.POS, tx_strand) # can't decide
    else:
      tx_strand = list(intersections)[0][2]

  # determine weight if relevant
  if weight is None:
    count_weight = 1
  else:
    count_weight = variant.INFO[weight] # let's hope for a float
  logging.debug('weight is %s', count_weight)

  if indels and len(variant.REF) != len(alt):
    category = get_indel_category(variant, chroms)
    indel_category = '{category}_{content}_{length}_{repeats}'.format(**category)
    counts[indel_category] += count_weight
    #logging.debug('indel at %s:%s %s -> %s classified as %s', no_chr(variant.CHROM), variant.POS, variant.REF, alt, indel_category)
    if indel_category not in INDELS_83:
      logging.warning('unexpected indel category %s', indel_category)

  # no need to look at this indel any more
  if len(variant.REF) != 1 or len(alt) != 1:
    if not indels:
      logging.debug('skipped indel at %s:%i', no_chr(variant.CHROM), variant.POS)
    return

  if just_indels:
    return

  # --- sbs
  # 0 1 2 -> my indexes
  # 1 2 3 -> vcf indexes
  # pulling in -1 0 +1
  context_extension = int((mer - 1) / 2) # 3 -> 1, 5 -> 2, etc
  fragment = chroms[no_chr(variant.CHROM)][variant.POS - 1 - context_extension:variant.POS + context_extension].upper() # potentially could instead skip lower case
  if fragment[context_extension] != variant.REF:
    logging.warning('skipping variant with position mismatch at %s:%i: VCF: %s genome: %s[%s]%s', no_chr(variant.CHROM), variant.POS, variant.REF, fragment[0], fragment[1], fragment[2])
    return

  if any([x not in 'ACGT' for x in ''.join([fragment, alt])]):
    logging.warning('skipping variant with ill-defined transition {}>{} at {}:{}'.format(fragment, alt, no_chr(variant.CHROM), variant.POS))
    return
    
  v = '{}>{}'.format(fragment, alt) # TODO multiallele
  v, tx_strand, exon_strand = normalize_sbs(v, tx_strand, exon_strand)
  counts[v] += count_weight
  if tx_strand is not None:
    tx_counts['{}/{}'.format(v, tx_strand)] += count_weight
  if exon_strand is not None:
    tx_counts['{}|{}'.format(v, exon_strand)] += count_weight

  # --- doublets
  if doublets:
    if last_variant is not None and no_chr(last_variant.CHROM) == no_chr(variant.CHROM) and last_variant.POS == variant.POS - 1:
      doublet = '{}{}>{}{}'.format(last_variant.REF, variant.REF, last_variant.ALT[0], alt) # TODO ALT[0] could be weird if '-' ?
      if len(doublet) != 5:
        logging.warning('skipping doublet %s at %s:%i', doublet, no_chr(variant.CHROM), variant.POS)
      else:
        doublet = normalize_doublet(doublet)
        counts[doublet] += count_weight
        logging.debug('doublet found at %s:%s: %s', no_chr(variant.CHROM), variant.POS, doublet)

def multi_count(genome_fh, vcf_in, outs=None, chroms=None, variant_filters=None, doublets=False, indels=False, just_indels=False, mer=3, weight=None):
  logging.info('multi_count with %i filters...', len(variant_filters))

  if chroms is None:
    chroms = {}
    chroms_seen = set()
  else:
    chroms_seen = set(chroms.keys())
    logging.debug('using existing genome with %i chromosomes', len(chroms_seen))

  all_counts = [collections.defaultdict(int) for _ in range(len(variant_filters))]

  next_chrom = None

  last_variant = None
  #vcf_in = cyvcf2.VCF(vcf)
  for line, variant in enumerate(vcf_in):
    if no_chr(variant.CHROM) not in chroms_seen:
      logging.debug('chrom %s seen in vcf', no_chr(variant.CHROM))
      next_chrom = update_chroms(no_chr(variant.CHROM), chroms, genome_fh, next_chrom)
      chroms_seen.add(no_chr(variant.CHROM))

    # keep track of lots of counts
    for idx, variant_filter in enumerate(variant_filters):
      if variant_filter(vcf_in, variant):
        update_counts(all_counts[idx], variant, last_variant, chroms, doublets, indels, just_indels, mer=mer, weight=weight)

    last_variant = variant

    if (line + 1) % 100000 == 0:
      logging.debug('processed %i lines', line + 1)
    
  logging.info('multi_count: done')

  # write out results
  all_total_counts = []
  for idx, counts in enumerate(all_counts):
    out = outs[idx]
    total_count = sum([counts[v] for v in counts])
    all_total_counts.append(total_count)
    if out is not None:
      out.write('{}\t{}\t{}\n'.format('Variation', 'Count', 'Probability'))
      for k in sorted(counts):
        out.write('{}\t{}\t{:.6f}\n'.format(k, counts[k], counts[k] / max(1, total_count)))

      # add zero results for SBS
      if not just_indels:
        for ref in ('C', 'T'):
          for alt in ('A', 'C', 'G', 'T'):
            if ref == alt:
              continue
            for prefix in ('A', 'C', 'G', 'T'):
              for suffix in ('A', 'C', 'G', 'T'):
                count = '{}{}{}>{}'.format(prefix, ref, suffix, alt)
                if count not in counts:
                  out.write('{}\t{}\t{:.6f}\n'.format(count, 0, 0))

      # add zero results for doublets
      if doublets:
        for doublet in DOUBLETS:
          if doublet not in counts:
            out.write('{}\t{}\t{:.6f}\n'.format(doublet, 0, 0))

      if indels:
        for indel in INDELS_83:
          if indel not in counts:
            out.write('{}\t{}\t{:.6f}\n'.format(indel, 0, 0))

  return {'chroms': chroms, 'all_counts': all_counts, 'all_totals': all_total_counts}

def read_transcripts(transcripts):
  logging.info('reading %s...', transcripts)
  tree = {}
  txtree = {}
  bases = 0
  txbases = 0
  #bin    name    chrom   strand  txStart txEnd   cdsStart        cdsEnd  exonCount       exonStarts      exonEnds        score   name2   cdsStartStat    cdsEndStat      exonFrames
  for rowline, row in enumerate(csv.DictReader(gzip.open(transcripts, 'rt'), delimiter='\t')):
    chrom = no_chr(row['chrom'])
    if chrom not in tree:
      logging.debug('added %s to transcripts. %i exon bases and %i transcript bases so far', chrom, bases, txbases)
      tree[chrom] = intervaltree.IntervalTree()
      txtree[chrom] = intervaltree.IntervalTree()
    txtree[chrom][int(row['txStart']):int(row['txEnd'])] = row['strand']
    txbases += int(row['txEnd']) - int(row['txStart'])
    for start, end in zip(row['exonStarts'].split(','), row['exonEnds'].split(',')):
      if start == '' or end == '':
        continue
      if EXCLUDE_UTR:
        start = max(int(start), int(row['cdsStart']))
        end = min(int(end), int(row['cdsEnd']))
        if start < end:
          tree[chrom][start:end] = row['strand'] # + or -
          bases += end - start
      else:
        tree[chrom][start:end] = row['strand']
        bases += end - start
    if rowline % 1000 == 0:
      logging.debug('processed %i lines...', rowline)
  logging.info('reading %s: done with %i exon bases and %i tx bases', transcripts, bases, txbases)
  return txtree, tree

def count(genome_fh, vcf_in, out=None, chroms=None, variant_filter=None, doublets=False, indels=False, just_indels=False, transcripts_fn=None, mer=3, weight=None, extended=None):
  '''
    extended not used
  '''
  logging.info('processing vcf...')

  if chroms is None:
    chroms = {}
    chroms_seen = set()
  else:
    chroms_seen = set(chroms.keys())
    logging.debug('using existing genome with %i chromosomes', len(chroms_seen))

  if transcripts_fn is not None:
    transcripts, exons = read_transcripts(transcripts_fn)
  else:
    transcripts, exons = (None, None)

  counts = collections.defaultdict(int)
  tx_counts = collections.defaultdict(int)
  next_chrom = None
  filtered = 0

  last_variant = None
  for line, variant in enumerate(vcf_in):
    if SKIP_ALT_CHROM and '_' in variant.CHROM:
      logging.info('skipping %s', variant.CHROM)
      continue

    if SKIP_M and (variant.CHROM == 'M' or variant.CHROM == 'MT'):
      logging.info('skipping %s', variant.CHROM)
      continue

    if no_chr(variant.CHROM) not in chroms_seen:
      logging.debug('chrom %s seen in vcf but not in %s', no_chr(variant.CHROM), chroms_seen)
      next_chrom = update_chroms(no_chr(variant.CHROM), chroms, genome_fh, next_chrom)
      chroms_seen.add(no_chr(variant.CHROM))

    if variant_filter is not None and not variant_filter(vcf_in, variant):
      filtered += 1
      continue

    update_counts(counts, variant, last_variant, chroms, doublets, indels, just_indels, transcripts, exons, tx_counts, mer, weight=weight)

    last_variant = variant

    if (line + 1) % 100000 == 0:
      logging.debug('processed %i lines. current counts: %s...', line + 1, ' '.join(['{}:{}'.format(k, counts[k]) for k in counts]))

  logging.info('processing vcf: filtered %i. included %i. done', filtered, sum([counts[k] for k in counts]))

  # write out results
  total_count = sum([counts[v] for v in counts])
  if out is not None:
    out.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format('Variation', 'Count', 'Probability', 'CodingTx', 'NonCodingTx', 'TxP', 'CodingExon', 'NonCodingExon', 'ExonP'))

    for k in sorted(counts):
      codingTx = tx_counts['{}/{}'.format(k, '+')]
      nonCodingTx = tx_counts['{}/{}'.format(k, '-')]
      codingExon = tx_counts['{}|{}'.format(k, '+')]
      nonCodingExon = tx_counts['{}|{}'.format(k, '-')]

      out.write('{}\t{}\t{:.6f}\t{}\t{}\t{:.6f}\t{}\t{}\t{:.6f}\n'.format(k, counts[k], counts[k] / total_count, 
        codingTx, nonCodingTx, scipy.stats.binomtest(k=codingTx, n=max(1, codingTx + nonCodingTx)).pvalue,
        codingExon, nonCodingExon, scipy.stats.binomtest(k=codingExon, n=max(1, codingExon + nonCodingExon)).pvalue))

    # add zero results for SBS
    context_extension = int((mer - 1) / 2) # 3 -> 1, 5 -> 2, etc
    if not just_indels:
      for ref in ('C', 'T'):
        for alt in ('A', 'C', 'G', 'T'):
          if ref == alt:
            continue
          for prefix in itertools.product('ACGT', repeat=context_extension):
            for suffix in itertools.product('ACGT', repeat=context_extension):
              count = '{}{}{}>{}'.format(''.join(prefix), ref, ''.join(suffix), alt)
              if count not in counts:
                out.write('{}\t{}\t{:.6f}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(count, 0, 0, 0, 0, 1, 0, 0, 1))

    # add zero results for doublets
    if doublets:
      for doublet in DOUBLETS:
        if doublet not in counts:
          out.write('{}\t{}\t{:.6f}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(doublet, 0, 0, 0, 0, 1, 0, 0, 1))

    if indels:
      for indel in INDELS_83:
        if indel not in counts:
          out.write('{}\t{}\t{:.6f}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(indel, 0, 0, 0, 0, 1, 0, 0, 1))

  return {'chroms': chroms, 'counts': counts, 'total': total_count}

def get_value(header, col, row):
  return row[header.index(col)]

def get_info(header, row):
  return {header[i]: row[i] for i in range(len(header))}

def open_file(fn, is_gzipped):
  if is_gzipped:
    return gzip.open(fn, 'rt')
  else:
    return open(fn, 'rt')

def maf_to_vcf(maf, sample, sample_col, chrom_col, pos_col, ref_col, alt_col, is_not_zipped, maf_filter_column):

  Variant = collections.namedtuple('Variant', 'CHROM POS REF ALT FILTER INFO')

  # enumeration a maf into a variant
  header = None
  for line, row in enumerate(csv.reader(open_file(maf, not is_not_zipped), delimiter='\t')):
    if line % 1000 == 0:
      logging.debug('processed %i lines of %s...', line, maf)

    if row[0].startswith('#'):
      continue
    if header is None:
      header = row
      continue

    #Hugo_Symbol     Entrez_Gene_Id  Center  NCBI_Build      Chromosome      Start_Position  End_Position    Strand  Variant_Classification  Variant_Type    Reference_Allele        Tumor_Seq_Allele1       Tumor_Seq_Allele2       dbSNP_RS        dbSNP_Val_Status        Tumor_Sample_Barcode    Matched_Norm_Sample_Barcode     Match_Norm_Seq_Allele1  Match_Norm_Seq_Allele2  Tumor_Validation_Allele1        Tumor_Validation_Allele2        Match_Norm_Validation_Allele1   Match_Norm_Validation_Allele2   Verification_Status     Validation_Status       Mutation_Status Sequencing_Phase        Sequence_Source Validation_Method       Score   BAM_File        Sequencer       Tumor_Sample_UUID       Matched_Norm_Sample_UUID        HGVSc   HGVSp   HGVSp_Short     Transcript_ID   Exon_Number     t_depth t_ref_count     t_alt_count     n_depth n_ref_count     n_alt_count     all_effects     Allele  Gene    Feature Feature_type    One_Consequence Consequence     cDNA_position   CDS_position    Protein_position        Amino_acids Codons  Existing_variation      ALLELE_NUM      DISTANCE        TRANSCRIPT_STRAND       SYMBOL  SYMBOL_SOURCE   HGNC_ID BIOTYPE CANONICAL       CCDS    ENSP    SWISSPROT       TREMBL  UNIPARC RefSeq  SIFT    PolyPhen        EXON    INTRON  DOMAINS GMAF    AFR_MAF AMR_MAF ASN_MAF EAS_MAF EUR_MAF SAS_MAF AA_MAF  EA_MAF  CLIN_SIG        SOMATIC PUBMED  MOTIF_NAME      MOTIF_POS       HIGH_INF_POS    MOTIF_SCORE_CHANGE      IMPACT  PICK    VARIANT_CLASS   TSL     HGVS_OFFSET     PHENO   MINIMISED       ExAC_AF ExAC_AF_Adj     ExAC_AF_AFR     ExAC_AF_AMR     ExAC_AF_EAS     ExAC_AF_FIN     ExAC_AF_NFE     ExAC_AF_OTH     ExAC_AF_SAS     GENE_PHENO      FILTER  CONTEXT src_vcf_id      tumor_bam_uuid  normal_bam_uuid case_id GDC_FILTER      COSMIC  MC3_Overlap     GDC_Validation_Status

    row_sample = get_value(header, sample_col, row)
    if sample is not None and row_sample != sample:
      continue

    chrom = no_chr(get_value(header, chrom_col, row))
    pos = int(get_value(header, pos_col, row))
    ref = get_value(header, ref_col, row)
    if ref == '-':
      pos += 1 # fix for TCGA mafs
    ref = ref.replace('-', '')
    alt = get_value(header, alt_col, row).replace('-', '')
    flter = get_value(header, maf_filter_column, row)
    info = get_info(header, row)

    yield Variant(chrom, pos, ref, (alt,), flter, info)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='mutational signature counter')
  parser.add_argument('--genome', required=True, help='reference genome')
  parser.add_argument('--vcf', required=True, nargs='+', help='vcfs or mafs')
  parser.add_argument('--out', required=False, nargs='+', help='output files or stdout by default')
  parser.add_argument('--vcf_not_zipped', action='store_true', help='do not try to unzip (only matters for maf)')
  parser.add_argument('--maf_sample', required=False, help='vcf is actually a maf with this sample of interest')
  parser.add_argument('--maf_sample_column', required=False, default='Tumor_Sample_Barcode', help='maf chrom column name')
  parser.add_argument('--maf_chrom_column', required=False, default='Chromosome', help='maf chrom column name')
  parser.add_argument('--maf_pos_column', required=False, default='Start_Position', help='maf pos column name')
  parser.add_argument('--maf_ref_column', required=False, default='Reference_Allele', help='maf ref column name')
  parser.add_argument('--maf_alt_column', required=False, default='Tumor_Seq_Allele2', help='maf alt column name')
  parser.add_argument('--maf_filter_column', required=False, default='FILTER', help='maf alt column name')
  parser.add_argument('--doublets', action='store_true', help='count doublets')
  parser.add_argument('--indels', action='store_true', help='count indels')
  parser.add_argument('--just_indels', action='store_true', help='count only indels')
  parser.add_argument('--transcripts', required=False, help='refseq transcript file')
  parser.add_argument('--extended', required=False, nargs='+', help='the 5 prime base count e.g. at -3=A -4=A')
  parser.add_argument('--mer', required=False, default=3, type=int, help='context length to consider for sbs')
  parser.add_argument('--weight_field', required=False, help='weight mutations by vcf field')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  chroms = None
  for idx, v in enumerate(args.vcf):
    if args.maf_sample is not None:
      vcf_in = maf_to_vcf(v, args.maf_sample, args.maf_sample_column, args.maf_chrom_column, args.maf_pos_column, args.maf_ref_column, args.maf_alt_column, args.vcf_not_zipped, args.maf_filter_column)
    else:
      vcf_in = cyvcf2.VCF(v)
    if args.out is None:
      out = sys.stdout
      out_fn = 'stdout'
    else:
      out = open(args.out[idx], 'w')
      out_fn = args.out[idx]
    logging.info('processing %i of %i: %s -> %s...', idx, len(args.vcf), v, out_fn)
    result = count(genome_fh=open(args.genome, 'r'), vcf_in=vcf_in, out=out, chroms=chroms, doublets=args.doublets, indels=args.indels, just_indels=args.just_indels, transcripts_fn=args.transcripts, mer=args.mer, weight=args.weight_field, extended=args.extended)
    chroms = result['chroms']
    logging.info('chroms is %s', chroms.keys())
    logging.info('processing %i of %i: %s -> %s: done', idx, len(args.vcf), v, out_fn)
