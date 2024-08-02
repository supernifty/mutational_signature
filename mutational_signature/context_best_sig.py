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

import matplotlib.pyplot as plt

def main(ifh, ofh, ranks, pie):
  logging.info('starting...')
  
  ctxs = collections.defaultdict(dict)
  for r in csv.DictReader(ifh, delimiter='\t'):
    s = r['Sig']
    cands = {k: float(r[k]) for k in r if k != 'Sig'}
    for c in cands:
      ctxs[c][s] = cands[c]
    
  # figure highest for each ctx
  fieldnames = ['context']
  for i in range(ranks):
    fieldnames.extend(['sig{}'.format(i+1), 'proportion{}'.format(i+1)])
  odw = csv.DictWriter(ofh, delimiter='\t', fieldnames=fieldnames)
  odw.writeheader()
  for c in ctxs:
    total = sum([ctxs[c][s] for s in ctxs[c]])
    logging.debug('total is %.3f', total)
    cands = ctxs[c] # {'sig': prop...}
    bests = sorted(cands, key=cands.get)[::-1]
    o = {'context': c}
    for i in range(ranks):
      o['sig{}'.format(i+1)] = bests[i]
      o['proportion{}'.format(i+1)] = '{:.6f}'.format(cands[bests[i]] / total)
    odw.writerow(o)
    if pie is not None:
       target = pie.replace('CTX', c)
       fig = plt.figure(figsize=(12, 6))
       ax = fig.add_subplot(111) 
       labels = bests[:ranks]
       values = [cands[bests[i]] / total for i in range(ranks)]
       labels.append('Other')
       values.append(1-sum(values))
       ax.pie(values, labels=labels, autopct='%.0f%%', labeldistance=None, textprops={'fontsize': 16})
       ax.legend(labels=labels, loc='best', bbox_to_anchor=(-0.1, 1.), fontsize=14)
       plt.savefig(target, bbox_inches="tight")
       logging.info('wrote %s', target)

  logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Assess MSI')
  parser.add_argument('--ranks', type=int, required=False, default=1, help='how many to show')
  parser.add_argument('--pie', required=False, help='template for pie file x.CTX.png')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(sys.stdin, sys.stdout, args.ranks, args.pie)

