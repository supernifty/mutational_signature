#!/usr/bin/env python
'''
  Excludes specified contexts and renormalises remaining contexts
  Usage example:
    exclude_contexts TCAA TCAG CCAA CCAG TCAT TCAC CCAC CCAT < signatures30.txt > signatures.reduced.txt
  
  note: outer two bases are flanking, inner two are variant
'''

import sys

def exclude(contexts):
  contexts = set(contexts)
  sys.stderr.write('excluding {}\n'.format(contexts))
  header = sys.stdin.readline().strip('\n').split('\t')
  sys.stdout.write('{}\n'.format('\t'.join([x for x in header if x not in contexts])))
  for line in sys.stdin:
    cols = line.strip('\n').split('\t')
    sig = cols[0]
    vals = [float(x) for i, x in enumerate(cols[1:]) if header[i + 1] not in contexts]
    sys.stderr.write('sig {} now has {}\n'.format(sig, len(vals)))
    new_sum = sum(vals) # normalise
    vals = [x / new_sum for x in vals]
    sys.stdout.write('{}\t{}\n'.format(sig, '\t'.join([str(x) for x in vals])))

if __name__ == '__main__':
  exclude(sys.argv[1:])
