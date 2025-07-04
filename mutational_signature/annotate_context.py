#!/usr/bin/env python
'''
  adds sig_context to a vcf
'''

import argparse
import collections
import csv
import gzip
import logging
import sys

import cyvcf2
import intervaltree

import mutational_signature.count

COMP = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
EXCLUDE_UTR=True # transcription bias

def normalize(v):
  '''
    input GAT>G => ATC>C
  '''
  if v[1] in ('C', 'T'):
    return v # ok
  else:
    return ''.join([COMP[v[2]], COMP[v[1]], COMP[v[0]], '>', COMP[v[4]]])

RC = {'a': 't', 'c': 'g', 't': 'a', 'g': 'c', 'n': 'n', 'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C', 'N': 'N'}
def reverse_complement(r):
  return ''.join([RC[x] for x in r][::-1])
    
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

def doublet_context(variant, chroms, last_variant):
  '''
    look for a doublet context!
  '''
  if last_variant is not None and no_chr(last_variant.CHROM) == no_chr(variant.CHROM) and last_variant.POS == variant.POS - 1:
    doublet = '{}{}>{}{}'.format(last_variant.REF, variant.REF, last_variant.ALT[0], variant.ALT[0]) # TODO ALT[0] could be weird if '-' ?
    if len(doublet) != 5:
      logging.warning('skipping doublet %s at %s:%i', doublet, no_chr(variant.CHROM), variant.POS)
    else:
      doublet = mutational_signature.count.normalize_doublet(doublet)
      logging.debug('doublet found at %s:%s: %s', no_chr(variant.CHROM), variant.POS, doublet)
      return doublet
  return None

def context(variant, chroms):
  chrom = variant.CHROM.replace('chr', '')
  if chrom not in chroms:
    logging.info('context: skipping chromosome %s', chrom)
    return None
  if variant.POS == 1 or variant.POS > len(chroms[chrom]):
    logging.info('context: skipped edge variant at %s:%i with chrom length %i', chrom, variant.POS, len(chroms[chrom]))
    return None
  if len(variant.REF) != 1 or len(variant.ALT) == 0 or len(variant.ALT[0]) != 1 or variant.ALT[0] == '-':
    # indel context!
    logging.debug('context: indel at %s:%i', chrom, variant.POS)
    category = mutational_signature.count.get_indel_category(variant, chroms)
    indel_category = '{category}_{content}_{length}_{repeats}'.format(**category)
    return indel_category

  # 0 1 2 -> my indexes
  # 1 2 3 -> vcf indexes
  fragment = chroms[chrom][variant.POS - 2:variant.POS + 1].upper() # TODO should we instead skip lower case
  if fragment[1] != variant.REF:
    logging.warn('context: skipping variant with position mismatch at %s:%i: VCF: %s genome: %s[%s]%s', chrom, variant.POS, variant.REF, fragment[0], fragment[1], fragment[2])
    return

  if any([x not in 'ACGT' for x in ''.join([fragment, variant.ALT[0]])]):
    logging.warn('context: skipping variant with ill-defined transition {}>{} at {}:{}'.format(fragment, variant.ALT[0], chrom, variant.POS))
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
    logging.info('surrounding: skipping chromosome %s', chrom)
    return None
  if variant.POS < sequence or variant.POS > len(chroms[chrom]) - sequence:
    logging.info('surrounding: skipped edge variant at %s:%i from chromosome length %i', chrom, variant.POS, len(chroms[chrom]))
    return None

  # 0 1 2 -> my indexes
  # 1 2 3 -> vcf indexes
  fragment = chroms[chrom][variant.POS - sequence - 1:variant.POS + sequence + len(variant.REF) - 1].upper() # TODO should we instead skip lower case
  if fragment[sequence:sequence + len(variant.REF)] != variant.REF:
    logging.warn('surrounding: skipping variant with position mismatch at %s:%i: VCF: %s genome: %s', chrom, variant.POS, variant.REF, fragment)
    return

  return fragment

def vcf_writer(out):
  def write_header(vcf_in):
    out.write(vcf_in.raw_header)

  def write_variant(variant, sig_context, surrounding_context, tx_strand, sequence_name, best_sig, ctx_totals):
    if sig_context is not None:
      variant.INFO["sig_context"] = sig_context
      if best_sig is not None and sig_context in best_sig:
        variant.INFO["best_sig"] = best_sig[sig_context][0]
        variant.INFO["best_sig_likelihood"] = '{:.3f}'.format(best_sig[sig_context][1] / ctx_totals[sig_context])
    if surrounding_context is not None:
      variant.INFO[sequence_name] = surrounding_context
    if tx_strand is not None:
      variant.INFO["tx_strand"] = tx_strand
    out.write(str(variant))

  def add_to_header(d):
    vcf_in.add_info_to_header(d)

  return {'write_header': write_header, 'write_variant': write_variant, 'add_to_header': add_to_header}

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

def annotate(genome_fh, vcf_in, out=None, chroms=None, variant_filter=None, sequence=0, plot=None, transcripts_fn=None, sequence_name=None, signatures=None, signatures_as_paths=False):
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
  out['add_to_header']({'ID': 'sig_context', 'Description': 'mutational signature trinucleotide context', 'Type':'Character', 'Number': '1'})
  out['add_to_header']({'ID': 'best_sig', 'Description': 'signature most likely associated with this context', 'Type':'Character', 'Number': '1'})
  out['add_to_header']({'ID': 'best_sig_likelihood', 'Description': 'likelihood of best signature associated with this context', 'Type':'Character', 'Number': '1'})
  if sequence > 0:
    out['add_to_header']({'ID': sequence_name, 'Description': 'reference sequence surrounding variant', 'Type':'Character', 'Number': '1'})
  out['write_header'](vcf_in)

  if transcripts_fn is not None:
    transcripts, exons = read_transcripts(transcripts_fn)
    out['add_to_header']({'ID': 'tx_strand', 'Description': 'which transcript if any is the variant on', 'Type':'Character', 'Number': '1'})
  else:
    transcripts, exons = (None, None)

  best_sig = {} # ctx -> (signame prop) (with highest prop)
  ctx_totals = collections.defaultdict(float)
  if signatures is not None:
    if not signatures_as_paths:
      from importlib.resources import files
      signatures = [files("mutational_signature").joinpath("data").joinpath(s) for s in signatures]
    for s in signatures:
      for r in csv.DictReader(open(s, 'rt'), delimiter='\t'):
        signame = r['Sig']
        for ctx in r:
          if ctx == 'Sig':
            continue
          prop = float(r[ctx])
          if len(ctx) == 4: # ACGT -> ACT>G
            ctx = '{}{}{}>{}'.format(ctx[0], ctx[1], ctx[3], ctx[2])
          ctx_totals[ctx] += prop
          if ctx not in best_sig or best_sig[ctx][1] < prop:
            best_sig[ctx] = (signame, prop)

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

    sig_context = context(variant, chroms)
    surrounding_context = surrounding(variant, sequence, chroms)

    # determine transcript if relevant
    if transcripts is not None and no_chr(variant.CHROM) in transcripts:
      intersections = transcripts[no_chr(variant.CHROM)][variant.POS]
      if len(intersections) == 0:
        logging.debug('%s:%i is not coding', no_chr(variant.CHROM), variant.POS)
        tx_strand = None
      elif len(intersections) > 1:
        logging.debug('%i multiple intersections at %s:%i', len(intersections), no_chr(variant.CHROM), variant.POS) # can't decide
        tx_strand = list(intersections)[0][2]
        logging.debug('%i multiple intersections at %s:%i: first is %s', len(intersections), no_chr(variant.CHROM), variant.POS, tx_strand) # can't decide
      else:
        tx_strand = list(intersections)[0][2]
    else:
      tx_strand = None

    out['write_variant'](variant, sig_context, surrounding_context, tx_strand, sequence_name, best_sig, ctx_totals)

    if (line + 1) % 10000 == 0:
      logging.debug('processed %i lines...', line + 1)
    
  logging.info('processing: filtered %i of %i. done', filtered, line)

# accept colnames of the form Variant/1
def maf_value(colname, row):
  if '/' in colname:
    c, n = colname.split('/')
    return row[c].split('/')[int(n)]
  elif '-' in colname:
    c, n = colname.split('-')
    return row[c].split('-')[int(n)]
  else:
    return row[colname]

def maf_to_vcf(maf, chrom_col, pos_col, ref_col, alt_col):

  def maf_reader(reader):
    Variant = collections.namedtuple('Variant', 'CHROM POS REF ALT row')
  
    for line, row in enumerate(reader):
      if line % 1000 == 0:
        logging.debug('processed %i lines of %s...', line, maf)
  
      if '/' in chrom_col or '-' in chrom_col or ':' in chrom_col:
        chrom = row 
      chrom = maf_value(chrom_col, row)
      pos = int(maf_value(pos_col, row))
      ref = maf_value(ref_col, row)
      if ref == '-':
        pos += 1 # fix for TCGA mafs
      ref = ref.replace('-', '')
      alt = maf_value(alt_col, row)
  
      yield Variant(chrom, pos, ref, (alt,), row)

  # enumeration a maf into a variant
  reader = csv.DictReader(open(maf, 'r'), delimiter='\t')
  return {'tsv_reader': reader, 'maf_reader': maf_reader(reader)}

def maf_writer(out): # DictWriter

  def write_header(vcf_in):
    out.writeheader()

  def write_variant(variant, sig_context, surrounding_context, tx_strand, sequence_name, best_sig, ctx_totals):
    if sig_context is not None:
      variant.row["sig_context"] = sig_context
      if best_sig is not None and sig_context in best_sig:
        #logging.info('writing ctx %s with %s and %s', sig_context, best_sig[sig_context], ctx_totals)
        variant.row["best_sig"] = best_sig[sig_context][0]
        variant.row["best_sig_likelihood"] = '{:.3f}'.format(best_sig[sig_context][1] / ctx_totals[sig_context])
    if surrounding_context is not None:
      variant.row[sequence_name] = surrounding_context
    if tx_strand is not None:
      variant.row["tx_strand"] = tx_strand
    logging.debug('%s', variant.row)
    out.writerow(variant.row)

  def dummy(d):
    pass

  return {'write_header': write_header, 'write_variant': write_variant, 'add_to_header': dummy}


if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='mutational signature counter')
  parser.add_argument('--genome', required=True, help='reference genome')
  parser.add_argument('--sequence', required=False, default=4, type=int, help='surrounding sequence in each direction to annotate, 0 to skip')
  parser.add_argument('--sequence_name', required=False, default='surrounding', help='surrounding sequence field name')
  parser.add_argument('--vcf', required=True, help='vcf')
  parser.add_argument('--is_maf', action='store_true', help='vcf is a maf')
  parser.add_argument('--maf_chrom_column', required=False, default='Chromosome', help='maf chrom column name')
  parser.add_argument('--maf_pos_column', required=False, default='Start_Position', help='maf pos column name')
  parser.add_argument('--maf_ref_column', required=False, default='Reference_Allele', help='maf ref column name')
  parser.add_argument('--maf_alt_column', required=False, default='Tumor_Seq_Allele2', help='maf alt column name')
  parser.add_argument('--plot', required=False, help='plot context breakdowns')
  parser.add_argument('--transcripts', required=False, help='refseq transcript file')
  parser.add_argument('--signatures', required=False, nargs='+', help='signature files to add a representative signature')
  parser.add_argument('--signatures_as_paths', required=False, default=False, action='store_true', help='signature files are paths, rather than filenames in the data directory')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  if args.is_maf:
    vcf_in = maf_to_vcf(args.vcf, args.maf_chrom_column, args.maf_pos_column, args.maf_ref_column, args.maf_alt_column)
    writer = csv.DictWriter(sys.stdout, delimiter='\t', fieldnames=vcf_in['tsv_reader'].fieldnames + ['sig_context', 'best_sig', 'best_sig_likelihood', args.sequence_name, 'tx_strand'])
    vcf_out = maf_writer(writer)
    annotate(open(args.genome, 'r'), vcf_in['maf_reader'], vcf_out, sequence=args.sequence, plot=args.plot, transcripts_fn=args.transcripts, sequence_name=args.sequence_name, signatures=args.signatures, signatures_as_paths=args.signatures_as_paths)
  else:
    vcf_in = cyvcf2.VCF(args.vcf)
    annotate(open(args.genome, 'r'), vcf_in, vcf_writer(sys.stdout), sequence=args.sequence, plot=args.plot, transcripts_fn=args.transcripts, sequence_name=args.sequence_name, signatures=args.signatures, signatures_as_paths=args.signatures_as_paths)
