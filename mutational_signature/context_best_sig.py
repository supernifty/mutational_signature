#!/usr/bin/env python
'''
  highest signature by context
'''

import argparse
import collections
import csv
import logging
import operator
import sys

def main(ifh, ofh):
  logging.info('starting...')
  
  ctxs = collections.defaultdict(dict)
  for r in csv.DictReader(ifh, delimiter='\t'):
    s = r['Sig']
    cands = {k: float(r[k]) for k in r if k != 'Sig'}
    for c in cands:
      ctxs[c][s] = cands[c]
    
  # figure highest for each ctx
  odw = csv.DictWriter(ofh, delimiter='\t', fieldnames=['context', 'sig', 'proportion'])
  odw.writeheader()
  for c in ctxs:
    cands = ctxs[c]
    best = max(cands.items(), key=operator.itemgetter(1))
    odw.writerow({'context': c, 'sig': best[0], 'proportion': best[1]})

  logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Assess MSI')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(sys.stdin, sys.stdout)

