#!/usr/bin/env python
'''
  summarize surrounding sequence
'''

import argparse
import collections
import logging
import sys

import cyvcf2

RC = {'a': 't', 'c': 'g', 't': 'a', 'g': 'c', 'n': 'n', 'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C', 'N': 'N'}
def reverse_complement(repeat):
  return ''.join([RC[x] for x in repeat][::-1])

# TODO handle indels
def normalized_sequence(variant, normalize):
  ctx = variant.INFO['surrounding']
  logging.debug('considering %s with %s>%s', ctx, variant.REF, variant.ALT[0])
  if not normalize or variant.REF[0] in ('C', 'T'):
    return ctx
  else:
    return reverse_complement(ctx)

def normalized_sequence_with_ref_alt(sub, variant, normalize):
  ctx = variant.INFO['surrounding']
  if not normalize or variant.REF[0] in ('C', 'T'):
    return '{}|{}>{}'.format(sub, variant.REF, variant.ALT[0])
  else:
    # sub is already normalised
    return '{}|{}>{}'.format(sub, reverse_complement(variant.REF), reverse_complement(variant.ALT[0]))

def main(vcfs, sequence, lengths, normalize, snvs, proportions, smooth):
  logging.info('starting...')
  ctxs = set()
  result = collections.defaultdict(dict)
  for vcf in vcfs:
    logging.info('processing %s...', vcf)
    result[vcf]['TOTAL'] = 0
    for line, variant in enumerate(cyvcf2.VCF(vcf)):
      if snvs and len(variant.INFO['surrounding']) != (2 * sequence) + 1:
        continue # skip non-snv
      ctx = normalized_sequence(variant, normalize)
      result[vcf]['TOTAL'] += 1
      for length in lengths:
        # e.g. sequence=6 length=3
        #sub = ctx[sequence-length:-(sequence-length)]
        sub = ctx[sequence-length:sequence+length+1] # will include some of ref
        logging.debug('adding context %s', sub)
        ctxs.add(sub)
        if sub not in result[vcf]:
          result[vcf][sub] = 0
        result[vcf][sub] += 1
        with_alt = normalized_sequence_with_ref_alt(sub, variant, normalize) #'{}>{}'.format(sub, variant.ALT[0])
        logging.debug('adding context %s', with_alt)
        ctxs.add(with_alt)
        if with_alt not in result[vcf]:
          result[vcf][with_alt] = 0
        result[vcf][with_alt] += 1

  logging.info('writing results to stdout')
  # write each vcf as a row
  #  cols = sorted(ctxs)
  #  sys.stdout.write('vcf\tTOTAL\t{}\n'.format('\t'.join(cols)))
  #  for vcf in vcfs:
  #    sys.stdout.write('{}\t{}\t{}\n'.format(vcf, result[vcf]['TOTAL'], '\t'.join([str(result[vcf].get(col, 0)) for col in cols])))
  #  sys.stdout.write('vcf\tTOTAL\t{}\n'.format('\t'.join(cols)))

  # write each ctx as a row
  sys.stdout.write('ctx\t{}\n'.format('\t'.join(vcfs)))
  for ctx in sorted(ctxs):
    const = 1 if smooth else 0
    if proportions:
      sys.stdout.write('{}\t{}\n'.format(ctx, '\t'.join([str((result[vcf].get(ctx, 0) + const) / result[vcf].get('TOTAL', 1)) for vcf in vcfs])))
    else:
      sys.stdout.write('{}\t{}\n'.format(ctx, '\t'.join([str(result[vcf].get(ctx, 0) + const) for vcf in vcfs])))
  sys.stdout.write('TOTAL\t{}\n'.format('\t'.join([str(result[vcf].get('TOTAL', 0)) for vcf in vcfs])))

  logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Assess MSI')
  parser.add_argument('--vcfs', nargs='+', required=True, help='vcfs to summarise')
  parser.add_argument('--sequence', required=True, type=int, help='sequence length provided')
  parser.add_argument('--lengths', required=True, nargs='+', type=int, help='sequence requested')
  parser.add_argument('--snvs', action='store_true', help='only snvs')
  parser.add_argument('--normalize', action='store_true', help='normalize variant sequence')
  parser.add_argument('--proportions', action='store_true', help='write counts as proportions')
  parser.add_argument('--smooth', action='store_true', help='add one to counts')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.vcfs, args.sequence, args.lengths, args.normalize, args.snvs, args.proportions, args.smooth)
