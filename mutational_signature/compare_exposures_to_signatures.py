#!/usr/bin/env python
'''
  for each signature, how similar are exposures?
'''

import argparse
import csv
import collections
import logging
import sys

import numpy as np
import scipy.stats

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt

def cosine(x, y):
  similarity = x.dot(y) / (np.linalg.norm(x) * np.linalg.norm(y))
  return similarity

def pearson(x, y):
  similarity =  scipy.stats.pearsonr(x, y) # returns correlation, p-value
  return similarity[0]

def compare(exposures, signatures, measure, out):

  exclude = ('Error', 'Mutations')

  logging.info('processing %s...', exposures)
  # Variation       Count   Probability
  # ACG>T   1.0102588553396328      0.14432254929725394
  # ATT>G   1.5636287170029137      0.22337530762488528
  x = {}
  for r in csv.DictReader(open(exposures, 'rt'), delimiter='\t'):
    x[r['Variation']] = float(r['Count'])

  logging.debug('contexts: %s', x.keys())

  odw = csv.DictWriter(out, delimiter='\t', fieldnames=['Sig', 'cosine_similarity'])
  odw.writeheader()
  logging.info('processing %s...', signatures)
  # Sig     ACAA    ACAC    ACAG    ACAT    ACGA    ACGC    ACGG    ACGT    ACTA    ACTC    ACTG    ACTT    ATAA    ATAC    ATAG    ATAT    ATCA    ATCC    ATCG    ATCT    ATGA    ATGC    ATGG    ATGT    CCAA    CCAC    CCAG    CCAT    CCGA    CCGC    CCGG    CCGT    CCTA    CCTC    CCTG    CCTT    CTAA    CTAC    CTAG    CTAT    CTCA    CTCC    CTCG    CTCT    CTGA    CTGC    CTGG    CTGT    GCAA    GCAC    GCAG    GCAT    GCGA    GCGC    GCGG    GCGT    GCTA    GCTC    GCTG    GCTT    GTAA    GTAC    GTAG    GTAT    GTCA    GTCC    GTCG    GTCT    GTGA    GTGC    GTGG    GTGT    TCAA    TCAC    TCAG    TCAT    TCGA    TCGC    TCGG    TCGT    TCTA    TCTC    TCTG    TCTT    TTAA    TTAC    TTAG    TTAT    TTCA    TTCC    TTCG    TTCT    TTGA    TTGC    TTGG    TTGT
  # SBS1    0.000886157230877471    0.00228040461219034     0.000177031410683197    0.00128022715070335     0.0018603300783658      0.00122021650301413     0.000115020408071004    0.00114020230609517     0.0250044365371747      0.00632112155659776     0.365064773442751       0.00958170008104535     0.00080014196918959     0.00223039573911599     0.00114020230609517     0.000183032475452119    0.00109019343302082     0.00304053948292045     0.000106018810917621    0.00574101862893531     0.000172030523375762    0.000207036734527807    0.000268047559678513    0.000112019875686543    0.000312055367983941    0.00179031765606171     9.32165394105873e-05    2.23039573911599e-16    2.41042768218364e-05    7.68136290422007e-05    0.00035206246644342     2.23039573911599e-16    0.00200035492297398     0.000270047914601487    0.19603478245145        0.00019603478245145     4.30076308439405e-05    0.000393069742364387    0.000324057497521784    0.000260046139986617    2.23039573911599e-16    0.00250044365371747     0.000360063886135316    4.26075598593457e-05    3.55062998827881e-05    0.000212037621835242    0.000128022715070335    0.000171030345914275    0.00158028038914944     0.000339060159444089    0.000587104169892862    2.23039573911599e-16    9.59170185566022e-06    0.000164029103683866    0.00016602945860684     2.23039573911599e-16    0.00444078792900222     9.28164684259925e-05    0.218038686604164       3.84068145211004e-05    8.12144098727434e-05    0.000129022892531822    0.000246043655525799    0.000258045785063643    0.00105018633456134     0.00190033717682528     0.00117020762993978     7.13126530040222e-05    2.23039573911599e-16    2.23039573911599e-16    0.000348061756597472    1.460259093771e-05      6.58116769658438e-05    0.00253044897756208     2.23039573911599e-16    5.88104347354349e-06    2.23039573911599e-16    0.000203036024681859    6.9412315827197e-05     0.0015602768399197      0.00111019698225056     3.73066193134647e-05    0.110019520763569       2.23039573911599e-16    0.00672119254119256     2.23039573911599e-16    2.80049689216357e-05    0.00225039928834573     0.000255045252679182    0.00339060159444089     0.000416073823978587    0.00433076840823866     2.23039573911599e-16    5.5109778127933e-05     0.000583103460046914    2.23039573911599e-16
  for r in csv.DictReader(open(signatures, 'rt'), delimiter='\t'):
    y = {}
    for c in r: # columns
      if c == 'Sig':
        continue
      context = '{}{}{}>{}'.format(c[0], c[1], c[3], c[2])
      if context in x:
        y[context] = float(r[c])
      else:
        logging.debug('skipping context %s', context)

    logging.debug('contexts: %s', y.keys())

    # now compare
    common = sorted(list(set(x.keys()).intersection(set(y.keys()))))
    xr = [x[c] for c in common]
    yr = [y[c] for c in common]
    logging.debug('comparing %s to %s', xr, yr)
    similarity = cosine(np.array(xr), np.array(yr))
    odw.writerow({'Sig': r['Sig'], 'cosine_similarity': similarity})

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='compare signatures')
  parser.add_argument('--exposures', required=True, help='exposures')
  parser.add_argument('--signatures', required=True, help='second group of exposures')
  parser.add_argument('--measure', required=False, default='cosine', help='cosine or pearson')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)
  compare(args.exposures, args.signatures, args.measure, sys.stdout)

