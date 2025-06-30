#!/usr/bin/env python
'''
  spreadsheet of signature probability by context
'''

import argparse
import collections
import csv
import logging
import sys

import matplotlib.pyplot as plt

def main(ifh, ofh):
  logging.info('starting...')
  
  ctxs = collections.defaultdict(dict)
  sigs = set()
  for r in csv.DictReader(ifh, delimiter='\t'):
    s = r['Sig']
    sigs.add(s)
    cands = {k: float(r[k]) for k in r if k != 'Sig'} # list of context: val for this sig
    for c in cands:
      ctxs[c][s] = cands[c]

  odw = csv.DictWriter(ofh, delimiter='\t', fieldnames=['context', 'best'] + sorted(list(sigs)))
  odw.writeheader()
  for c in ctxs:
    total = sum([ctxs[c][s] for s in ctxs[c]])
    cands = ctxs[c] # {'sig': prop...}
    bests = sorted(cands, key=cands.get)[::-1]
    if len(c) == 4 and all([x in 'ACGT' for x in c]):
      r = {'context': '{}{}{}>{}'.format(c[0], c[1], c[3], c[2])}
    else:
      r = {'context': c}
    t = 0
    for s in sigs:
      r[s] = '{:.1f}'.format(100 * cands[s] / total)
      t += 100 * cands[s] / total
    bests = ';'.join(['{} {}'.format(x, r[x]) for x in sorted(cands, key=cands.get)[::-1][:3]])
    r['best'] = bests
    odw.writerow(r)
  logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Likelihood by context')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(sys.stdin, sys.stdout)


