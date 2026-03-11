#!/usr/bin/env python
'''
  given tumour and normal vcf pairs, explore msi status
'''

import argparse
import csv
import logging
import random
import sys

import numpy as np
import scipy.optimize
import scipy.spatial.distance

def make_distance(A, b):
  def cosine(x):
    estimate = np.dot(A, x)
    similarity = b.dot(estimate) / (np.linalg.norm(b) * np.linalg.norm(estimate))
    return -similarity

  return cosine

def main(signatures, signatures2, components=5):
  logging.info('starting...')

  #Sig     ACAA    ACAC    ACAG    ACAT    ACGA    ACGC    ACGG    ACGT    ACTA    ACTC    ACTG    ACTT    ATAA    ATAC    ATAG    ATAT    ATCA    ATCC    ATCG    ATCT    ATGA    ATGC    ATGG    ATGT    CCAA    CCAC    CCAG    CCAT    CCGA    CCGC    CCGG    CCGT      CCTA    CCTC    CCTG    CCTT    CTAA    CTAC    CTAG    CTAT    CTCA    CTCC    CTCG    CTCT    CTGA    CTGC    CTGG    CTGT    GCAA    GCAC    GCAG    GCAT    GCGA    GCGC    GCGG    GCGT    GCTA    GCTC    GCTG    GCTT    GTAA    GTAC    GTAG    GTAT      GTCA    GTCC    GTCG    GTCT    GTGA    GTGC    GTGG    GTGT    TCAA    TCAC    TCAG    TCAT    TCGA    TCGC    TCGG    TCGT    TCTA    TCTC    TCTG    TCTT    TTAA    TTAC    TTAG    TTAT    TTCA    TTCC    TTCG    TTCT    TTGA    TTGC    TTGG    TTGT

  defs = {}
  deps = {}
  for r in csv.DictReader(open(signatures, 'r'), delimiter='\t'):
    n = r['Sig']
    v = np.array([float(r[x]) for x in sorted(r) if x != 'Sig'])
    defs[n] = v

  if signatures2 is not None:
    for r in csv.DictReader(open(signatures2, 'r'), delimiter='\t'):
      n = r['Sig']
      v = np.array([float(r[x]) for x in sorted(r) if x != 'Sig'])
      deps[n] = v
  else:
    deps = defs


  # approximate linear combination
  sys.stdout.write('Sig\tError\t{}\n'.format('\t'.join(['Component {}'.format(x) for x in range(1, components)])))
  for s in defs:
    logging.info('solving for %s', s)
    A = np.array([deps[n] for n in sorted(deps) if n != s])
    A = np.swapaxes(A, 0, 1)
    #x = ?
    b = defs[s]
    #b = np.swapaxes(b, 0, 1)
    
    x0 = np.array([random.random() for _ in range(0, A.shape[1])])
    x0 = x0 / sum(x0) # normalize
    bounds=[(0.0, np.inf)] * len(x0)
    minimizer_kwargs = dict(method="L-BFGS-B", bounds=bounds)

    result = scipy.optimize.basinhopping(make_distance(A, b), x0, minimizer_kwargs=minimizer_kwargs, stepsize=5, T=5)

    error = 1 + result.fun
    formula = []
    for sr, vr in zip([x for x in sorted(deps) if x != s], result.x):
      if vr != 0:
        formula.append((vr / sum(result.x), sr))

    sys.stdout.write('{}\t{:.3f}\t{}\n'.format(s, error, '\t'.join(['{:.3f} * {}'.format(f[0], f[1]) for f in sorted(formula, key=lambda x: -x[0])[:components]])))
    #sys.exit(0)
    
    #sys.stdout.write('{} Formula\t{}\n'.format(s, ' + '.join(['{:.3f} * {}'.format(f[0], f[1]) for f in sorted(formula, key=lambda x: -x[0])])))

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Find linearly dependent combinations')
  parser.add_argument('--signatures', required=True, help='signature definitions')
  parser.add_argument('--signatures2', required=False, help='if specified, use these only as possible dependends')
  parser.add_argument('--components', required=False, type=int, default=5, help='how many components to use')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.signatures, args.signatures2, args.components)
