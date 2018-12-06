#!/usr/bin/env python

import argparse
import collections
import csv
import gzip
import logging
import sys

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
        # don't overwrite previous chromosomes
        if next_chrom in chroms:
          logging.debug('already have chrom %s. not overwriting', next_chrom)
        else:
          chroms[next_chrom] = ''.join(seq)
          logging.debug('reading chrom %s from genome. size is %i: done', next_chrom, len(chroms[next_chrom]))
        seq = []
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

def update_counts(counts, pos, chrom, ref, alt, chroms):
  if pos == 1 or pos > len(chroms[chrom]):
    logging.info('skipped edge variant at %s:%i < %i', chrom, pos, len(chroms[chrom]))
    return
  if len(ref.replace('-', '')) != 1 or len(alt.replace('-', '')) != 1:
    logging.info('skipped indel at %s:%i', chrom, pos)
    return

  # 0 1 2 -> my indexes
  # 1 2 3 -> vcf indexes
  fragment = chroms[chrom][pos - 2:pos + 1].upper() # TODO should we instead skip lower case
  if fragment[1] != ref:
    logging.warn('skipping variant with position mismatch at %s:%i: VCF: %s genome: %s[%s]%s', chrom, pos, ref, fragment[0], fragment[1], fragment[2])
    return

  if any([x not in 'ACGT' for x in ''.join([fragment, alt])]):
    logging.warn('skipping variant with ill-defined transition {}>{} at {}:{}'.format(fragment, alt, chrom, pos))
    return
    
  v = '{}>{}'.format(fragment, alt) # TODO multiallele
  v = normalize(v)
  counts[v] += 1

def get_value(header, col, row):
  return row[header.index(col)]

def count_contexts(genome_fh, maf, sample=None, chroms=None):
  logging.info('processing %s...', maf)
  if chroms is None:
    chroms = {}
  chroms_seen = set()
  counts = {}
  next_chrom = None

  header = None
  for line, row in enumerate(csv.reader(gzip.open(maf, 'rt'), delimiter='\t')):
    if row[0].startswith('#'):
      continue
    if header is None:
      header = row
      continue

    #Hugo_Symbol     Entrez_Gene_Id  Center  NCBI_Build      Chromosome      Start_Position  End_Position    Strand  Variant_Classification  Variant_Type    Reference_Allele        Tumor_Seq_Allele1       Tumor_Seq_Allele2       dbSNP_RS        dbSNP_Val_Status        Tumor_Sample_Barcode    Matched_Norm_Sample_Barcode     Match_Norm_Seq_Allele1  Match_Norm_Seq_Allele2  Tumor_Validation_Allele1        Tumor_Validation_Allele2        Match_Norm_Validation_Allele1   Match_Norm_Validation_Allele2   Verification_Status     Validation_Status       Mutation_Status Sequencing_Phase        Sequence_Source Validation_Method       Score   BAM_File        Sequencer       Tumor_Sample_UUID       Matched_Norm_Sample_UUID        HGVSc   HGVSp   HGVSp_Short     Transcript_ID   Exon_Number     t_depth t_ref_count     t_alt_count     n_depth n_ref_count     n_alt_count     all_effects     Allele  Gene    Feature Feature_type    One_Consequence Consequence     cDNA_position   CDS_position    Protein_position        Amino_acids Codons  Existing_variation      ALLELE_NUM      DISTANCE        TRANSCRIPT_STRAND       SYMBOL  SYMBOL_SOURCE   HGNC_ID BIOTYPE CANONICAL       CCDS    ENSP    SWISSPROT       TREMBL  UNIPARC RefSeq  SIFT    PolyPhen        EXON    INTRON  DOMAINS GMAF    AFR_MAF AMR_MAF ASN_MAF EAS_MAF EUR_MAF SAS_MAF AA_MAF  EA_MAF  CLIN_SIG        SOMATIC PUBMED  MOTIF_NAME      MOTIF_POS       HIGH_INF_POS    MOTIF_SCORE_CHANGE      IMPACT  PICK    VARIANT_CLASS   TSL     HGVS_OFFSET     PHENO   MINIMISED       ExAC_AF ExAC_AF_Adj     ExAC_AF_AFR     ExAC_AF_AMR     ExAC_AF_EAS     ExAC_AF_FIN     ExAC_AF_NFE     ExAC_AF_OTH     ExAC_AF_SAS     GENE_PHENO      FILTER  CONTEXT src_vcf_id      tumor_bam_uuid  normal_bam_uuid case_id GDC_FILTER      COSMIC  MC3_Overlap     GDC_Validation_Status

    row_sample = get_value(header, "Tumor_Sample_Barcode", row)
    if sample is not None and row_sample != sample:
      continue

    chrom = get_value(header, "Chromosome", row)

    if chrom not in chroms_seen:
      logging.debug('chrom %s seen in maf', chrom)
      #chroms = {} # wipe previous chromosomes (TODO if maf sorted)
      next_chrom = update_chroms(chrom, chroms, genome_fh, next_chrom)
      chroms_seen.add(chrom)

    pos = int(get_value(header, "Start_Position", row))
    ref = get_value(header, "Reference_Allele", row).replace('-', '')
    alt = get_value(header, "Tumor_Seq_Allele2", row).replace('-', '')

    if row_sample not in counts:
      counts[row_sample] = collections.defaultdict(int)

    update_counts(counts[row_sample], pos, chrom, ref, alt, chroms)

    if (line + 1) % 100000 == 0:
      logging.debug('processed %i lines. %i samples seen...', line + 1, len(counts))
    
  logging.info('processing %s: done', maf)

  return (counts, chroms)

def write_counts(genome_fh, maf, sample, out):
  counts, chroms = count_contexts(genome_fh, maf, sample)

  # write out results
  total_count = sum([counts[sample][v] for v in counts[sample]])
  out.write('{}\t{}\t{}\n'.format('Variation', 'Count', 'Probability'))
  for k in sorted(counts[sample]):
    out.write('{}\t{}\t{:.3f}\n'.format(k, counts[sample][k], counts[sample][k] / total_count))
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
  parser.add_argument('--maf', required=True, help='maf')
  parser.add_argument('--sample', required=True, help='sample')
  args = parser.parse_args()
  logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  write_counts(open(args.genome, 'r'), args.maf, args.sample, sys.stdout)
