#!/usr/bin/env python
'''
  summarize surrounding sequence
'''

import argparse
import collections
import logging
import sys

import cyvcf2

def main(vcfs, sequence, lengths):
  logging.info('starting...')
  ctxs = set()
  result = collections.defaultdict(dict)
  for vcf in vcfs:
    logging.info('processing %s...', vcf)
    result[vcf]['TOTAL'] = 0
    for line, variant in enumerate(cyvcf2.VCF(vcf)):
      ctx = variant.INFO["surrounding"]
      result[vcf]['TOTAL'] += 1
      for length in lengths:
        sub = ctx[sequence-length:-(sequence-length)]
        ctxs.add(sub)
        if sub not in result[vcf]:
          result[vcf][sub] = 0
        result[vcf][sub] += 1
        with_alt = '{}>{}'.format(sub, variant.ALT[0])
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
    sys.stdout.write('{}\t{}\n'.format(ctx, '\t'.join([str(result[vcf].get(ctx, 0)) for vcf in vcfs])))
  sys.stdout.write('TOTAL\t{}\n'.format('\t'.join([str(result[vcf].get('TOTAL', 0)) for vcf in vcfs])))

  logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Assess MSI')
  parser.add_argument('--vcfs', nargs='+', required=True, help='vcfs to summarise')
  parser.add_argument('--sequence', required=True, type=int, help='sequence length provided')
  parser.add_argument('--lengths', required=True, nargs='+', type=int, help='sequence requested')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.vcfs, args.sequence, args.lengths)
