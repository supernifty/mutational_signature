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

INDELS = set(['DEL_C_1_0', 'DEL_C_1_1', 'DEL_C_1_2', 'DEL_C_1_3', 'DEL_C_1_4', 'DEL_C_1_5+', 'DEL_MH_2_1', 'DEL_MH_3_1', 'DEL_MH_3_2', 'DEL_MH_4_1', 'DEL_MH_4_2', 'DEL_MH_4_3', 'DEL_MH_5+_1', 'DEL_MH_5+_2', 'DEL_MH_5+_3', 'DEL_MH_5+_4', 'DEL_MH_5+_5+', 'DEL_T_1_0', 'DEL_T_1_1', 'DEL_T_1_2', 'DEL_T_1_3', 'DEL_T_1_4', 'DEL_T_1_5+', 'DEL_repeats_2_0', 'DEL_repeats_2_1', 'DEL_repeats_2_2', 'DEL_repeats_2_3', 'DEL_repeats_2_4', 'DEL_repeats_2_5+', 'DEL_repeats_3_0', 'DEL_repeats_3_1', 'DEL_repeats_3_2', 'DEL_repeats_3_3', 'DEL_repeats_3_4', 'DEL_repeats_3_5+', 'DEL_repeats_4_0', 'DEL_repeats_4_1', 'DEL_repeats_4_2', 'DEL_repeats_4_3', 'DEL_repeats_4_4', 'DEL_repeats_4_5+', 'DEL_repeats_5+_0', 'DEL_repeats_5+_1', 'DEL_repeats_5+_2', 'DEL_repeats_5+_3', 'DEL_repeats_5+_4', 'DEL_repeats_5+_5+', 'INS_C_1_0', 'INS_C_1_1', 'INS_C_1_2', 'INS_C_1_3', 'INS_C_1_4', 'INS_C_1_5+', 'INS_T_1_0', 'INS_T_1_1', 'INS_T_1_2', 'INS_T_1_3', 'INS_T_1_4', 'INS_T_1_5+', 'INS_repeats_2_0', 'INS_repeats_2_1', 'INS_repeats_2_2', 'INS_repeats_2_3', 'INS_repeats_2_4', 'INS_repeats_2_5+', 'INS_repeats_3_0', 'INS_repeats_3_1', 'INS_repeats_3_2', 'INS_repeats_3_3', 'INS_repeats_3_4', 'INS_repeats_3_5+', 'INS_repeats_4_0', 'INS_repeats_4_1', 'INS_repeats_4_2', 'INS_repeats_4_3', 'INS_repeats_4_4', 'INS_repeats_4_5+', 'INS_repeats_5+_0', 'INS_repeats_5+_1', 'INS_repeats_5+_2', 'INS_repeats_5+_3', 'INS_repeats_5+_4', 'INS_repeats_5+_5+'])

def normalize_doublet(v):
  if v not in DOUBLETS:
    # reverse it
    v = '{}{}>{}{}'.format(COMP[v[1]], COMP[v[0]], COMP[v[4]], COMP[v[3]])
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
        logging.info('reading chrom %s from genome. size is %i: done', next_chrom, len(chroms[next_chrom]))
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

def update_counts(variant, last_variant, chroms, doublets, indels, just_indels, mer=3):
  if no_chr(variant.CHROM) not in chroms:
    logging.warning('chromosome %s not found in chroms: %s', no_chr(variant.CHROM), chroms.keys())
    return

  if variant.POS == 1 or variant.POS > len(chroms[no_chr(variant.CHROM)]):
    logging.info('skipped edge variant at %s:%i', no_chr(variant.CHROM), variant.POS)
    return

  exon_strand = None
  tx_strand = None

  if indels and len(variant.REF) != len(variant.ALT[0]):
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
      if len(variant.ALT[0]) == 0:
        # deletions look like TA -> - e.g. in TCGA
        deleted_sequence = variant.REF
      else:  
        # deletions look like TAAA -> T POS is where the T is
        deleted_sequence = variant.REF[len(variant.ALT[0]):]
      logging.debug('deleted sequence is %s from %s to %s', deleted_sequence, variant.REF, variant.ALT[0])
      # first look for ongoing repeated sequence. POS-1 puts us at the T
      current_pos = variant.POS + len(variant.ALT[0]) - 1 + len(deleted_sequence)
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
            if chroms[no_chr(variant.CHROM)][variant.POS + len(variant.ALT[0]) - 1 + len(deleted_sequence) + pos] == deleted_sequence[pos]:
              microhomology += 1
            else:
              break
          # now look to the left
          for pos in range(len(deleted_sequence) - 1, -1, -1):
            if chroms[no_chr(variant.CHROM)][variant.POS + len(variant.ALT[0]) - 1 - len(deleted_sequence) + pos] == deleted_sequence[pos]:
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
        inserted_sequence = variant.ALT[0]
      else:
        # insertions look like G -> GCA. POS is where the G is
        inserted_sequence = variant.ALT[0][len(variant.REF):]
      current_pos = variant.POS + len(variant.REF) - 1
      while chroms[no_chr(variant.CHROM)][current_pos:current_pos + len(inserted_sequence)] == inserted_sequence:
        current_pos += len(inserted_sequence)
        repeats += 1
        if repeats > 10000:
          break
      category['repeats'] = '5+' if repeats >= 5 else repeats
      if del_length == -1:
        category['content'] = INDEL_COMP[inserted_sequence]
      else: # insertion length > 1
        category['content'] = 'repeats'

    indel_category = '{category}_{content}_{length}_{repeats}'.format(**category)
    return indel_category
    #logging.debug('indel at %s:%s %s -> %s classified as %s', no_chr(variant.CHROM), variant.POS, variant.REF, variant.ALT[0], indel_category)
    if indel_category not in INDELS:
      logging.warning('unexpected indel category %s', indel_category)

  # no need to look at this indel any more
  if len(variant.REF) != 1 or len(variant.ALT[0]) != 1:
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

  if any([x not in 'ACGT' for x in ''.join([fragment, variant.ALT[0]])]):
    logging.warning('skipping variant with ill-defined transition {}>{} at {}:{}'.format(fragment, variant.ALT[0], no_chr(variant.CHROM), variant.POS))
    return
    
  v = '{}>{}'.format(fragment, variant.ALT[0]) # TODO multiallele
  v, tx_strand, exon_strand = normalize_sbs(v, tx_strand, exon_strand)
  return v
  if tx_strand is not None:
    tx_counts['{}/{}'.format(v, tx_strand)] += 1
  if exon_strand is not None:
    tx_counts['{}|{}'.format(v, exon_strand)] += 1

  # --- doublets
  if doublets:
    if last_variant is not None and no_chr(last_variant.CHROM) == no_chr(variant.CHROM) and last_variant.POS == variant.POS - 1:
      doublet = '{}{}>{}{}'.format(last_variant.REF, variant.REF, last_variant.ALT[0], variant.ALT[0])
      if len(doublet) != 5:
        logging.warning('skipping doublet %s at %s:%i', doublet, no_chr(variant.CHROM), variant.POS)
      else:
        doublet = normalize_doublet(doublet)
        return doublet
        logging.debug('doublet found at %s:%s: %s', no_chr(variant.CHROM), variant.POS, doublet)

def count(genome_fh, vcf_in, out=None, chroms=None, variant_filter=None, doublets=False, indels=False, just_indels=False, mer=3):
  logging.info('processing vcf...')

  if chroms is None:
    chroms = {}
    chroms_seen = set()
  else:
    chroms_seen = set(chroms.keys())
    logging.info('using existing genome with %i chromosomes', len(chroms_seen))

  counts = collections.defaultdict(int)
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

    context = update_counts(variant, last_variant, chroms, doublets, indels, just_indels, mer)

    # write annotated row
    out.write('{}\t{}\n'.format('\t'.join(variant.row), context))

    last_variant = variant

def get_value(header, col, row):
  return row[header.index(col)]

def open_file(fn, is_gzipped):
  if is_gzipped:
    return gzip.open(fn, 'rt')
  else:
    return open(fn, 'rt')

def maf_to_vcf(maf, chrom_col, pos_col, ref_col, alt_col, is_not_zipped, out):

  Variant = collections.namedtuple('Variant', 'CHROM POS REF ALT row')

  # enumeration a maf into a variant
  header = None
  for line, row in enumerate(csv.reader(open_file(maf, not is_not_zipped), delimiter='\t')):
    if line % 1000 == 0:
      logging.debug('processed %i lines of %s...', line, maf)

    if row[0].startswith('#'):
      continue
    if header is None:
      header = row
      header.append('context')
      out.write('{}\n'.format('\t'.join(header)))
      continue

    #Hugo_Symbol     Entrez_Gene_Id  Center  NCBI_Build      Chromosome      Start_Position  End_Position    Strand  Variant_Classification  Variant_Type    Reference_Allele        Tumor_Seq_Allele1       Tumor_Seq_Allele2       dbSNP_RS        dbSNP_Val_Status        Tumor_Sample_Barcode    Matched_Norm_Sample_Barcode     Match_Norm_Seq_Allele1  Match_Norm_Seq_Allele2  Tumor_Validation_Allele1        Tumor_Validation_Allele2        Match_Norm_Validation_Allele1   Match_Norm_Validation_Allele2   Verification_Status     Validation_Status       Mutation_Status Sequencing_Phase        Sequence_Source Validation_Method       Score   BAM_File        Sequencer       Tumor_Sample_UUID       Matched_Norm_Sample_UUID        HGVSc   HGVSp   HGVSp_Short     Transcript_ID   Exon_Number     t_depth t_ref_count     t_alt_count     n_depth n_ref_count     n_alt_count     all_effects     Allele  Gene    Feature Feature_type    One_Consequence Consequence     cDNA_position   CDS_position    Protein_position        Amino_acids Codons  Existing_variation      ALLELE_NUM      DISTANCE        TRANSCRIPT_STRAND       SYMBOL  SYMBOL_SOURCE   HGNC_ID BIOTYPE CANONICAL       CCDS    ENSP    SWISSPROT       TREMBL  UNIPARC RefSeq  SIFT    PolyPhen        EXON    INTRON  DOMAINS GMAF    AFR_MAF AMR_MAF ASN_MAF EAS_MAF EUR_MAF SAS_MAF AA_MAF  EA_MAF  CLIN_SIG        SOMATIC PUBMED  MOTIF_NAME      MOTIF_POS       HIGH_INF_POS    MOTIF_SCORE_CHANGE      IMPACT  PICK    VARIANT_CLASS   TSL     HGVS_OFFSET     PHENO   MINIMISED       ExAC_AF ExAC_AF_Adj     ExAC_AF_AFR     ExAC_AF_AMR     ExAC_AF_EAS     ExAC_AF_FIN     ExAC_AF_NFE     ExAC_AF_OTH     ExAC_AF_SAS     GENE_PHENO      FILTER  CONTEXT src_vcf_id      tumor_bam_uuid  normal_bam_uuid case_id GDC_FILTER      COSMIC  MC3_Overlap     GDC_Validation_Status

    chrom = no_chr(get_value(header, chrom_col, row))
    pos = int(get_value(header, pos_col, row))
    ref = get_value(header, ref_col, row)
    if ref == '-':
      pos += 1 # fix for TCGA mafs
    ref = ref.replace('-', '')
    alt = get_value(header, alt_col, row).replace('-', '')

    yield Variant(chrom, pos, ref, (alt,), row)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='mutational signature counter')
  parser.add_argument('--genome', required=True, help='reference genome')
  parser.add_argument('--maf', required=True, help='maf file')
  parser.add_argument('--maf_chrom_column', required=False, default='Chromosome', help='maf chrom column name')
  parser.add_argument('--maf_pos_column', required=False, default='Start_Position', help='maf pos column name')
  parser.add_argument('--maf_ref_column', required=False, default='Reference_Allele', help='maf ref column name')
  parser.add_argument('--maf_alt_column', required=False, default='Tumor_Seq_Allele2', help='maf alt column name')
  parser.add_argument('--mer', type=int, default=3, help='context width')
  parser.add_argument('--doublets', action='store_true', help='count doublets')
  parser.add_argument('--indels', action='store_true', help='count indels')
  parser.add_argument('--just_indels', action='store_true', help='count only indels')
  parser.add_argument('--not_zipped', action='store_true', help='do not try to unzip')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  chroms = None
  
  out = sys.stdout
  out_fn = 'stdout'
  vcf_in = maf_to_vcf(args.maf, args.maf_chrom_column, args.maf_pos_column, args.maf_ref_column, args.maf_alt_column, args.not_zipped, out)
  logging.info('starting...')
  result = count(genome_fh=open(args.genome, 'r'), vcf_in=vcf_in, out=out, chroms=chroms, doublets=args.doublets, indels=args.indels, just_indels=args.just_indels, mer=args.mer)
  logging.info('done')

