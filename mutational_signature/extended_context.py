#!/usr/bin/env python
'''
  look for mutations in extended contexts
  example rules
  substitutions
  T>*,-3=A - any substitution from T on 5' to 3' with A three bases earlier
  A>*,+3=T - reverse strand version of the same thing

  deletions
  T>,-1=A,+1=T,+2=T - deletion with previous A and following TT
  A>,+1=T,-1=A,-2=A - revserse strand version of the same thing
'''

import argparse
import csv
import collections
import gzip
import logging
import sys

import cyvcf2
import intervaltree

EXCLUDE_UTR=True

def no_chr(chrom):
  if chrom == 'MT':
    return 'M'
  # deals with chrUn_gl000220.1
  return chrom.split('.')[0].split('_', maxsplit=1)[-1].replace('chr', '').upper()

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

def update_chroms(required, chroms, genome, next_chrom):
  '''
    pull the entire chromosome into memory
  '''
  seq = []
  for linenum, line in enumerate(genome):
    line = line.strip('\n')
    if line.startswith('>'):
      if next_chrom is None: # first line of file
        next_chrom = no_chr(line[1:].split(' ')[0].replace('chr', ''))
        logging.debug('reading chrom %s from genome...', next_chrom)
      else:
        # remove previous chromosomes
        chroms[next_chrom] = ''.join(seq)
        seq = []
        logging.info('reading chrom %s from genome. size is %i: done', next_chrom, len(chroms[next_chrom]))
        if required == next_chrom:
          next_chrom = no_chr(line[1:].split(' ')[0].replace('chr', ''))
          logging.debug('required chrom %s found', next_chrom)
          return next_chrom
        else:
          next_chrom = line[1:].split(' ')[0].replace('chr', '')
          logging.debug('reading chrom %s from genome...', next_chrom)
    else:
      seq.append(line)
    if linenum % 1000000 == 0:
      logging.debug('processed %i lines of genome...', linenum)

  # end of file
  chroms[next_chrom] = ''.join(seq)
  logging.info('reading chrom %s from genome. size is %i: done', next_chrom, len(chroms[next_chrom]))
  return None

def assess(genome_fh, vcf_in, out, transcripts_fn, rules, output):
  logging.info('processing vcf...')

  chroms = {}
  chroms_seen = set()

  # read transcripts if provided
  if transcripts_fn is not None:
    transcripts, exons = read_transcripts(transcripts_fn)
  else:
    transcripts, exons = (None, None)

  counts = collections.defaultdict(int)
  tx_counts = collections.defaultdict(int)
  next_chrom = None
  filtered = 0

  last_variant = None
  rules = [rule.split(',') for rule in rules] # [[T>*, -3=A], [A>*, +3=T]]
  counts = [[0, 0] for _ in rules] # any mutation, specific mutation
  tx_counts = [[0, 0, 0, 0] for _ in rules] # any mutation pos, specific mutation pos, any mutation neg, specific mutation neg
  tx_strand = None
  considered = indels = 0

  if output is not None:
    logging.info('writing variant annotation to %s', output)
    output_fh = csv.DictWriter(open(output, 'w'), delimiter='\t', fieldnames=['chrom', 'pos', 'ref', 'alt', 'result', 'detail'])
    output_fh.writeheader()
  else:
    output_fh = None

  for line, variant in enumerate(vcf_in):
    variant_passed = False
    chrom = no_chr(variant.CHROM)
    if chrom not in chroms_seen:
      logging.debug('chrom %s seen in vcf', chrom)
      next_chrom = update_chroms(chrom, chroms, genome_fh, next_chrom)
      chroms_seen.add(chrom)

    if False:
      filtered += 1
      continue

    if len(variant.REF) != len(variant.ALT[0]):
      indels += 1

    considered += 1

    # check transcript
    if transcripts is not None and chrom in transcripts:
      intersections = transcripts[chrom][variant.POS]
      if len(intersections) == 0:
        logging.debug('%s:%i is not coding', chrom, variant.POS)
        tx_strand = None
      elif len(intersections) > 1:
        logging.debug('%i multiple intersections at %s:%i', len(intersections), chrom, variant.POS) # can't decide
        tx_strand = list(intersections)[0][2]
        logging.debug('%i multiple intersections at %s:%i: first is %s', len(intersections), chrom, variant.POS, tx_strand) # can't decide
      else:
        tx_strand = list(intersections)[0][2]

      if tx_strand == '+':
        tx_strand_idx = 0
      elif tx_strand == '-':
        tx_strand_idx = 2
      else:
        pass
        # logging.warn('surprising tx strand: %s', tx_strand) # actually not surprising

    # assess each rule
    for idx, rule in enumerate(rules):
      # does it pass the variant requirement
      ref, alt = rule[0].split('>')
      rel_pos = 0
      variant_ok = len(variant.REF) == len(variant.ALT[0]) and (ref == '*' or ref == variant.REF) and (alt == '*' or alt == variant.ALT[0])
      if not variant_ok and len(variant.REF) != len(variant.ALT[0]):
        # see if it passes as an indel
        indel_size = len(variant.ALT[0]) - len(variant.REF)
        logging.debug('variant at %s:%i %s>%s: indel_size is %i. checking against %s>%s', chrom, variant.POS, variant.REF, variant.ALT[0], indel_size, ref, alt)
        if ref == variant.REF and alt == variant.ALT[0]: # T> (TCGA mafs)
          variant_ok = True
        elif indel_size < 0 and variant.REF[indel_size:] == ref and alt == '': # TA>T => T> # compare the deletion to ref
          variant_ok = True
          rel_pos = len(variant.ALT[0]) # e.g. +1
        elif indel_size > 0 and variant.REF == ref and variant.ALT[0][1:] == alt: # T>TAA TODO untested
          variant_ok = True
          rel_pos = len(variant.REF) # e.g. +1
        else:
          variant_ok = False
      if variant_ok:
        logging.debug('variant at %s:%i %s>%s will be evaluated against %s', chrom, variant.POS, variant.REF, variant.ALT[0], rule)
        passes = False
        counts[idx][0] += 1
        if tx_strand is not None:
          tx_counts[idx][0 + tx_strand_idx] += 1

        # does it pass the context requirement(s)
        current_pos = 0 # iterate along the repeated region
        while chroms[chrom][variant.POS - 1 + rel_pos + current_pos] == ref[0]: # should continue to be true
          passes = True # ok until something goes wrong
          logging.debug('trying rule %s at relative position %i', rule, current_pos)
          for context in rule[1:]:
            logging.debug('evaluating context %s...', context)
            if '=' in context:
              pos, expected_base = context.split('=')
              pos = variant.POS + rel_pos + int(pos) + current_pos
              if pos < 1 or pos > len(chroms[chrom]):
                logging.warn('unable to assess variant at %s:%i', chrom, variant.POS)
                passes = False
                break
              genome_base = chroms[chrom][pos - 1]
              logging.debug('genome base is %s at %s:%i, required base is %s', genome_base, chrom, pos, expected_base)
              if genome_base != expected_base:
                logging.debug('genome base is %s at %s:%i, required base is %s: failed', genome_base, chrom, pos, expected_base)
                passes = False
                break
            elif '~' in context:
              pos, expected_bases = context.split('~')
              pos = variant.POS + rel_pos + int(pos) + current_pos
              if pos < 1 or pos > len(chroms[chrom]):
                logging.warn('unable to assess variant at %s:%i', chrom, variant.POS)
                passes = False
                break
              genome_base = chroms[chrom][pos - 1]
              logging.debug('genome base is %s at %s:%i, required base is any of %s', genome_base, chrom, pos, expected_bases)
              if genome_base not in expected_bases:
                logging.debug('genome base is %s at %s:%i, required base is any of %s: failed', genome_base, chrom, pos, expected_bases)
                passes = False
                break
            elif '!' in context:
              pos, unexpected_base = context.split('!')
              pos = variant.POS + rel_pos + int(pos) + current_pos
              if pos < 1 or pos > len(chroms[chrom]):
                logging.warn('unable to assess variant at %s:%i', chrom, variant.POS)
                passes = False
                break
              genome_base = chroms[chrom][pos - 1]
              logging.debug('genome base is %s at %s:%i, base cannot be %s', genome_base, chrom, pos, unexpected_base)
              if genome_base == unexpected_base:
                logging.debug('genome base is %s at %s:%i, base cannot be %s: failed', genome_base, chrom, pos, unexpected_base)
                passes = False
                break
                
            else:
              logging.warn('invalid context %s', context)

          logging.debug('finished trying rules at %i, passes: %s', current_pos, passes)
          if passes: # found a position that succeeds for all contexts
            logging.debug('rule %s succeeded at relative position %i', rule, current_pos)
            break # while
          current_pos += 1

        if passes: # this rule passed
          logging.debug('variant at %s:%i %s>%s passed rule %s', chrom, variant.POS, variant.REF, variant.ALT[0], rule)
          counts[idx][1] += 1
          if tx_strand is not None:
            tx_counts[idx][1 + tx_strand_idx] += 1
          variant_passed = True
          variant_rule = rule
        else: # this rule failed
          pass

    if variant_passed:
      if output_fh is not None:
        output_fh.writerow({'chrom': chrom, 'pos': variant.POS, 'ref': variant.REF, 'alt': variant.ALT[0], 'result': 'PASS', 'detail': ','.join(variant_rule)}) # the last rule passed
    else:
      if output_fh is not None:
        output_fh.writerow({'chrom': chrom, 'pos': variant.POS, 'ref': variant.REF, 'alt': variant.ALT[0], 'result': 'FAIL', 'detail': ''})

    if line % 1000 == 0:
      logging.info('processed %i variants with interim results %s', line, counts)

  # now write out results
  out.write('Rule\tMutationCount\tPassCount\tPctPass\n')
  for rule, count in zip(rules, counts):
    if count[0] > 0: 
      out.write('{}\t{}\t{}\t{:.3f}\n'.format(','.join(rule), count[0], count[1], count[1] / count[0])) 
    else: # deal with zero denominator
      out.write('{}\t{}\t{}\t{}\n'.format(','.join(rule), count[0], count[1], 0))

  if transcripts is not None:
    for rule, count in zip(rules, tx_counts):
      # count contains +mut +pass -mut -pass
      if count[0] > 0: 
        out.write('{}_tx\t{}\t{}\t{:.3f}\n'.format(','.join(rule), count[0], count[1], count[1] / count[0])) 
      else: # deal with zero denominator
        out.write('{}_tx\t{}\t{}\t{}\n'.format(','.join(rule), count[0], count[1], 0))
      if count[2] > 0: 
        out.write('{}_notx\t{}\t{}\t{:.3f}\n'.format(','.join(rule), count[2], count[3], count[3] / count[2])) 
      else: # deal with zero denominator
        out.write('{}_notx\t{}\t{}\t{}\n'.format(','.join(rule), count[2], count[3], 0))

  # and write total
  counts_0 = sum([count[0] for count in counts])
  counts_1 = sum([count[1] for count in counts])
  if counts_0 > 0:
    out.write('{}\t{}\t{}\t{:.3f}\n'.format('TOTAL', counts_0, counts_1, counts_1 / counts_0)) 
  else:
    out.write('{}\t{}\t{}\t{}\n'.format('TOTAL', counts_0, counts_1, 0)) 

  # all mutations
  out.write('{}\t{}\t{}\t{}\n'.format('Indels', indels, 0, 0)) 
  out.write('{}\t{}\t{}\t{}\n'.format('Considered', considered, 0, 0)) 

def get_value(header, col, row):
  return row[header.index(col)]

def open_file(fn, is_gzipped):
  if is_gzipped:
    return gzip.open(fn, 'rt')
  else:
    return open(fn, 'rt')

def maf_to_vcf(maf, sample, sample_col, chrom_col, pos_col, ref_col, alt_col, is_not_zipped):

  Variant = collections.namedtuple('Variant', 'CHROM POS REF ALT')

  # enumeration a maf into a variant
  header = None
  for line, row in enumerate(csv.reader(open_file(maf, not is_not_zipped), delimiter='\t')):
    if line % 1000 == 0:
      logging.debug('processed %i lines of %s...', line, maf)

    if row[0].startswith('#'):
      continue
    if header is None:
      header = row
      logging.debug('header is %s', header)
      continue

    #Hugo_Symbol     Entrez_Gene_Id  Center  NCBI_Build      Chromosome      Start_Position  End_Position    Strand  Variant_Classification  Variant_Type    Reference_Allele        Tumor_Seq_Allele1       Tumor_Seq_Allele2       dbSNP_RS        dbSNP_Val_Status        Tumor_Sample_Barcode    Matched_Norm_Sample_Barcode     Match_Norm_Seq_Allele1  Match_Norm_Seq_Allele2  Tumor_Validation_Allele1        Tumor_Validation_Allele2        Match_Norm_Validation_Allele1   Match_Norm_Validation_Allele2   Verification_Status     Validation_Status       Mutation_Status Sequencing_Phase        Sequence_Source Validation_Method       Score   BAM_File        Sequencer       Tumor_Sample_UUID       Matched_Norm_Sample_UUID        HGVSc   HGVSp   HGVSp_Short     Transcript_ID   Exon_Number     t_depth t_ref_count     t_alt_count     n_depth n_ref_count     n_alt_count     all_effects     Allele  Gene    Feature Feature_type    One_Consequence Consequence     cDNA_position   CDS_position    Protein_position        Amino_acids Codons  Existing_variation      ALLELE_NUM      DISTANCE        TRANSCRIPT_STRAND       SYMBOL  SYMBOL_SOURCE   HGNC_ID BIOTYPE CANONICAL       CCDS    ENSP    SWISSPROT       TREMBL  UNIPARC RefSeq  SIFT    PolyPhen        EXON    INTRON  DOMAINS GMAF    AFR_MAF AMR_MAF ASN_MAF EAS_MAF EUR_MAF SAS_MAF AA_MAF  EA_MAF  CLIN_SIG        SOMATIC PUBMED  MOTIF_NAME      MOTIF_POS       HIGH_INF_POS    MOTIF_SCORE_CHANGE      IMPACT  PICK    VARIANT_CLASS   TSL     HGVS_OFFSET     PHENO   MINIMISED       ExAC_AF ExAC_AF_Adj     ExAC_AF_AFR     ExAC_AF_AMR     ExAC_AF_EAS     ExAC_AF_FIN     ExAC_AF_NFE     ExAC_AF_OTH     ExAC_AF_SAS     GENE_PHENO      FILTER  CONTEXT src_vcf_id      tumor_bam_uuid  normal_bam_uuid case_id GDC_FILTER      COSMIC  MC3_Overlap     GDC_Validation_Status

    row_sample = get_value(header, sample_col, row)
    if sample is not None and row_sample != sample:
      continue

    chrom = get_value(header, chrom_col, row).replace('chr', '')
    pos = int(get_value(header, pos_col, row))
    ref = get_value(header, ref_col, row)
    if ref == '-':
      pos += 1 # fix for TCGA mafs
    ref = ref.replace('-', '')
    alt = get_value(header, alt_col, row).replace('-', '')

    yield Variant(chrom, pos, ref, (alt,))


if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='mutational signature counter')
  parser.add_argument('--genome', required=True, help='reference genome')
  parser.add_argument('--transcripts', required=False, help='transcript from ucsc')
  parser.add_argument('--vcf', required=True, help='vcf or maf')
  parser.add_argument('--rules', nargs='+', required=True, help='extended contexts to look for. of the form e.g. "T>*,-3=A A>*,+3=T"')
  parser.add_argument('--vcf_not_zipped', action='store_true', help='do not try to unzip (only matters for maf)')
  parser.add_argument('--maf_sample', required=False, help='vcf is actually a maf with this sample of interest')
  parser.add_argument('--maf_sample_column', required=False, default='Tumor_Sample_Barcode', help='maf chrom column name')
  parser.add_argument('--maf_chrom_column', required=False, default='Chromosome', help='maf chrom column name')
  parser.add_argument('--maf_pos_column', required=False, default='Start_Position', help='maf pos column name')
  parser.add_argument('--maf_ref_column', required=False, default='Reference_Allele', help='maf ref column name')
  parser.add_argument('--maf_alt_column', required=False, default='Tumor_Seq_Allele2', help='maf alt column name')
  parser.add_argument('--output', required=False, help='write individual variant annotations to file')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)
 
  if args.maf_sample is not None:
    vcf_in = maf_to_vcf(args.vcf, args.maf_sample, args.maf_sample_column, args.maf_chrom_column, args.maf_pos_column, args.maf_ref_column, args.maf_alt_column, args.vcf_not_zipped)
  else:
    vcf_in = cyvcf2.VCF(args.vcf)

  assess(genome_fh=open(args.genome, 'r'), vcf_in=vcf_in, out=sys.stdout, transcripts_fn=args.transcripts, rules=args.rules, output=args.output)
