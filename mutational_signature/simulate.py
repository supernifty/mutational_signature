#!/usr/bin/env python
'''
  simulate contexts
'''

import argparse
import collections
import csv
import logging
import sys

import numpy.random

def main(n, output_template, defs, exposures, injections, noise, answers, mutation_counts, seed):
  # read defs
  logging.info('reading definitions from %s...', defs)
  # Sig     ACAA    ACAC    ACAG    ACAT    ACGA    ACGC    ACGG    ACGT    ACTA    ACTC    ACTG    ACTT    ATAA    ATAC    ATAG    ATAT    ATCA    ATCC    ATCG    ATCT    ATGA    ATGC    ATGG    ATGT    CCAA    CCAC CCAG     CCAT    CCGA    CCGC    CCGG    CCGT    CCTA    CCTC    CCTG    CCTT    CTAA    CTAC    CTAG    CTAT    CTCA    CTCC    CTCG    CTCT    CTGA    CTGC    CTGG    CTGT    GCAA    GCAC    GCAG    GCAT    GCGA GCGC     GCGG    GCGT    GCTA    GCTC    GCTG    GCTT    GTAA    GTAC    GTAG    GTAT    GTCA    GTCC    GTCG    GTCT    GTGA    GTGC    GTGG    GTGT    TCAA    TCAC    TCAG    TCAT    TCGA    TCGC    TCGG    TCGT TCTA     TCTC    TCTG    TCTT    TTAA    TTAC    TTAG    TTAT    TTCA    TTCC    TTCG    TTCT    TTGA    TTGC    TTGG    TTGT
  #SBS1    0.000886157230877471    0.00228040461219034     0.000177031410683197    0.00128022715070335     0.0018603300783658      0.00122021650301413     0.000115020408071004    0.00114020230609517     0.0250044365371747    0.00632112155659776     0.365064773442751       0.00958170008104535     0.00080014196918959     0.00223039573911599     0.00114020230609517     0.000183032475452119    0.00109019343302082     0.00304053948292045   0.000106018810917621    0.00574101862893531
  sigs = {}
  sfh = csv.DictReader(open(defs, 'rt'), delimiter='\t')
  contexts = [k for k in sfh.fieldnames if k != 'Sig']
  for r in sfh:
    sigs[r['Sig']] = [float(r[k]) for k in contexts]
  logging.info('%i sigs', len(sigs))

  if seed is not None:
    numpy.random.seed(seed)

  if answers is not None:
    adw = csv.DictWriter(open(answers, 'wt'), delimiter='\t', fieldnames=['iteration', 'signal', 'noise', 'total'] + sorted([b.split(',')[0] for b in exposures]))
    adw.writeheader()

  if mutation_counts is None:
    mutation_counts = [None]
 
  x = 0
  for injection in injections:
    for mutation_count in mutation_counts:
      # counts of contexts
      logging.info('starting %i iterations for %s...', n, injection)
      for _ in range(n):
        counts = collections.defaultdict(int) # empty counts
        # figure out total number of mutations from exposures
        total_muts = 0
        background = {}
        for exposure in exposures:
          name, mean, sd = exposure.split(',')
          mean = float(mean)
          sd = float(sd)
          muts = max(0, int(numpy.random.normal(loc=mean, scale=sd)))
          background[name] = muts #.append('{}={}'.format(name, muts))
          total_muts += muts
          logging.debug('%i muts for %s from %.3f±%.3f', muts, exposure, mean, sd)
          # pick each context with probability from signature definitions
          distn = numpy.random.multinomial(n=muts, pvals=sigs[name])
          for m in zip(contexts, distn):
            counts[m[0]] += m[1]
    
        # add injected signature
        injected_sig, injected_prop = injection.split(',')
        injected_prop = float(injected_prop) / (1 - float(injected_prop)) # so that final proportion is as expected
        injected_muts = int(injected_prop * total_muts + 0.5)
        logging.debug('injected %i from %s added to %i...', injected_muts, injected_sig, total_muts)
        distn = numpy.random.multinomial(n=injected_muts, pvals=sigs[injected_sig])
        for m in zip(contexts, distn):
          counts[m[0]] += m[1]
    
        # add noise
        mean, sd = noise.split(',')
        mean = float(mean)
        sd = float(sd)
        muts = max(0, int(numpy.random.normal(loc=mean, scale=sd)))
        if muts > 0:
          total_muts += muts
          distn = numpy.random.multinomial(n=muts, pvals=[1/len(contexts)]*len(contexts))
          for m in zip(contexts, distn):
            counts[m[0]] += m[1]
          logging.debug('%i noise muts', muts)
    
        # fixed mutation count
        total = sum([counts[v] for v in contexts])
        if mutation_count is not None:
          logging.debug('downsampling from %i to %i', total, mutation_count)
          distn = numpy.random.multinomial(n=mutation_count, pvals=[counts[v] / total for v in contexts])
          for m in zip(contexts, distn):
            counts[m[0]] = m[1]
     
        new_total = sum([counts[v] for v in contexts])
        target = output_template.replace('NUM', str(x))
        logging.debug('writing to %s...', target)
        odw = csv.DictWriter(open(target, 'wt'), delimiter='\t', fieldnames=['Variation', 'Count'])
        odw.writeheader()
        for v in contexts:
          count = counts[v]
          if len(v) == 4:
            v = '{}{}{}>{}'.format(v[0], v[1], v[3], v[2])
          odw.writerow({'Variation': v, 'Count': count})
        if answers is not None:
          #adw.writerow({'iteration': x, 'signal': injected_muts, 'noise': muts, 'total': sum([counts[v] for v in contexts]), 'background': ','.join(background)})
          # proportion
          result = {'iteration': x, 'signal': '{:.3f}'.format(injected_muts / total), 'noise': '{:.3f}'.format(muts / total), 'total': new_total}
          for sig in background:
            result[sig] = '{:.3f}'.format(background[sig] / total)
          adw.writerow(result)
        x += 1 # iteration

  logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='generate context counts')
  parser.add_argument('--n', type=int, default=10, help='number of simulations')
  parser.add_argument('--output_template', default='counts.NUM', help='output filenames NUM becomes simulation number')
  parser.add_argument('--sigdefs', required=True, help='signature definitions')
  parser.add_argument('--exposures', required=True, nargs='+', help='exposures of the form signame,meanmuts,sdmuts')
  parser.add_argument('--injections', required=True, nargs='+', help='exact exposure of the form signame,proportion each added separately not combined')
  parser.add_argument('--noise', required=False, default='0,0', help='noise term added to contexts mean,sd')
  parser.add_argument('--answers', required=False, help='write input parameters for each iteration to this file')
  parser.add_argument('--mutation_counts', required=False, nargs='+', type=int, help='fixed mutation count to downsample to')
  parser.add_argument('--seed', required=False, type=int, help='random seed for reproducibility')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.n, args.output_template, args.sigdefs, args.exposures, args.injections, args.noise, args.answers, args.mutation_counts, args.seed)
