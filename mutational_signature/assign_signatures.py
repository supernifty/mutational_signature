#!/usr/bin/env python
'''
  assign signature probabilities to variants
'''

import argparse
import collections
import logging
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pylab import rcParams
rcParams['figure.figsize'] = 16, 10

import cyvcf2

colors = {"SBS1": "#de3860", "SBS2": "#41ac2f", "SBS3": "#7951d0", "SBS4": "#73d053", "SBS5": "#b969e9", "SBS6": "#91ba2c", "SBS7a": "#b4b42f", "SBS7b": "#5276ec", "SBS7c": "#daae36", "SBS7d": "#9e40b5", "SBS8": "#43c673", "SBS9": "#dd4cb0", "SBS10a": "#3d9332", "SBS10b": "#de77dd", "SBS11": "#7bad47", "SBS12": "#9479e8", "SBS13": "#487b21", "SBS14": "#a83292", "SBS15": "#83c67d", "SBS16": "#664db1", "SBS17a": "#e18d28", "SBS17b": "#588de5", "SBS18": "#e2672a", "SBS19": "#34c7dd", "SBS20": "#cf402b", "SBS21": "#5acdaf", "SBS22": "#d74587", "SBS23": "#428647", "SBS24": "#7b51a7", "SBS25": "#b4ba64", "SBS26": "#646cc1", "SBS27": "#a27f1f", "SBS28": "#3b63ac", "SBS29": "#dca653", "SBS30": "#505099", "SBS31": "#7d8529", "SBS32": "#bf8ade", "SBS33": "#516615", "SBS34": "#b65da7", "SBS35": "#57a87a", "SBS36": "#c84249", "SBS37": "#37b5b1", "SBS38": "#a14622", "SBS39": "#58b5e1", "SBS40": "#ba6e2f", "SBS41": "#589ed8", "SBS42": "#e98261", "SBS43": "#3176ae", "SBS44": "#656413", "SBS45": "#a19fe2", "SBS46": "#756121", "SBS47": "#7e4a8d", "SBS48": "#326a38", "SBS49": "#dd8abf", "SBS50": "#1a6447", "SBS51": "#e78492", "SBS52": "#30876c", "SBS53": "#9d4d7c", "SBS54": "#919d5b", "SBS55": "#9d70ac", "SBS56": "#5b6f34", "SBS57": "#65659c", "SBS58": "#c9a865", "SBS59": "#a1455d", "SBS60": "#5e622c", "SBS84": "#b66057", "SBS85": "#dca173", "DBS1": "#855524", "DBS2": "#9f7846", "DBS3": "#7951d0", "DBS4": "#73d053", "DBS5": "#b969e9", "DBS6": "#91ba2c", "DBS7": "#3656ca", "DBS8": "#b4b42f", "DBS9": "#5276ec", "DBS10": "#daae36", "DBS11": "#9e40b5", "ID1": "#de3860", "ID2": "#41ac2f", "ID3": "#7951d0", "ID4": "#73d053", "ID5": "#b969e9", "ID6": "#91ba2c", "ID7": "#9e40b5", "ID8": "#43c673", "ID9": "#dd4cb0", "ID10": "#3d9332", "ID11": "#de77dd", "ID12": "#7bad47", "ID13": "#9479e8", "ID14": "#487b21", "ID15": "#a83292", "ID16": "#83c67d", "ID17": "#664db1"}

def plot_sbs_signature(vals, target, contexts, sigs):
  '''
    vals is a dictionary of context and values
    contexts: dict ACA>A => (prob, sig)...
    sigs is the list of signames to plot
  '''
  # ACG>A => CA ACG
  xs = sorted(['{}{} {}'.format(x[1], x[4], x[:3]) for x in vals.keys()])
  real_xs = ['{}>{}'.format(x[3:6], x[1]) for x in xs]
  ys = list([vals[x] for x in real_xs])

  color = ((0.2,0.7,0.9),)*16 + ((0.1,0.1,0.1),)*16 + ((0.8,0.2,0.2),)*16 + ((0.8,0.8,0.8),)*16 + ((0.6,0.8,0.4),)*16 + ((0.9,0.8,0.7),)*16
  color = list(color)

  ylim = max(ys)
  width = len(xs)
  x = range(width)
  f,ax = plt.subplots(1)

  if contexts is None: # no sigs
    bars = ax.bar(x, ys)
    for h in range(len(x)):
      bars[h].set_color(color[h])
  else:
    bottom = [0] * len(xs)
    for sig_idx, sig in enumerate(sigs):
      new_ys = []
      for context_idx, context in enumerate(real_xs):
        for item in contexts[context]:
          if item[1] == sig:
            new_ys.append(item[0] * ys[context_idx])
            break
      bars = ax.bar(x, new_ys, label=sig, bottom=bottom, color=colors[sig])
      bottom = [x + y for x,y in zip(bottom, new_ys)]
      for h in range(len(x)):
        bars[h].set_color(colors[sig])

  ax.set_xticks(x)
  ax.set_xticklabels([x.split(' ')[1] for x in xs], minor=False, rotation=90)
  ax.legend()
  plt.ylim(0, ylim)
  plt.xlim(0, width)

  ax2 = ax.twiny()
  ax2.set_xlim(ax.get_xlim())
  ax2.set_xticks([ width * (x/6.0 + 1/12.0) for x in range(6) ])
  ax2.set_xticklabels(['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G'])

  plt.savefig(target)

def main(vcf, signatures, definition, threshold, plot):
  logging.info('starting...')
  # signatures:
  # SBS1    0.234
  all_sigs = {}
  for line in open(signatures, 'r'):
    name, value = line.strip('\n').split('\t')
    if float(value) > threshold:
      all_sigs[name] = float(value)
    
  # definition:
  # Sig     ACAA     ACAC     ACAG     ACAT     ACGA    
  # SBS1    8.86E-04 2.28E-03 1.77E-04 1.28E-03 1.86E-03
  contexts = collections.defaultdict(list)
  sigs = {}

  with open(definition, 'r') as def_fh:
    header = def_fh.readline().strip('\n').split('\t')

    for line in def_fh:
      fields = line.strip('\n').split('\t')
      if fields[0] in all_sigs:
        sigs[fields[0]] = all_sigs[fields[0]]
        # keep a map of context -> sigs
        for ctx, context in enumerate(header[1:]):
          if len(context) == 4 and '>' not in context:
            context = '{}{}{}>{}'.format(context[0], context[1], context[3], context[2])
          context_val = float(fields[ctx + 1])
          # signature percent * context percent
          contexts[context].append((sigs[fields[0]] * context_val, fields[0])) # value, signature

  logging.info('%i signatures with value above %.2f', len(sigs), threshold)

  # normalize on contexts
  normalized_context = {}
  for context in contexts:
    total = sum([likelihood[0] for likelihood in contexts[context]])
    normalized_context[context] = [(contexts[context][idx][0] / total, contexts[context][idx][1]) for idx in range(0, len(contexts[context]))]
    logging.debug('normalized %s -> %s', context, normalized_context[context])

  contexts = normalized_context

  # normalize over total context
  annotation = {}
  for context in contexts:
    annotation[context] = ','.join(['{}/{:.3f}'.format(likelihood[1], likelihood[0]) for likelihood in contexts[context]])
    # find max
    best = None
    for item in contexts[context]:
      if best is None or item[0] > best[0]:
        best = item
    logging.debug('%s: %s most likely: %s', context, annotation[context], best)

  # vcf
  if vcf is not None:
    logging.info('annotating vcf...')
    
    vcf_in = cyvcf2.VCF(vcf)
    vcf_in.add_info_to_header({'ID': 'signature_likelihood', 'Description': 'signature likelihood', 'Type':'Character', 'Number': '1'})
    sys.stdout.write(vcf_in.raw_header)

    counts = collections.defaultdict(int)
    skipped_no_context = skipped_wrong_context = added = 0
    for line, variant in enumerate(vcf_in):
      try:
        context = variant.INFO["snv_context"]
      except:
        skipped_no_context += 1
        continue

      if context in contexts:
        variant.INFO["signature_likelihood"] = annotation[context]
        counts[context] += 1
        added += 1
      else:
        skipped_wrong_context += 1
      sys.stdout.write(str(variant))

    # reconstruction
    logging.debug('counts: %s', counts)
    expected = collections.defaultdict(float)
    for context in contexts:
      for prob, sig in contexts[context]:
        # how many variants are a result of this signature?
        expected[sig] += counts[context] * prob
        
    total = sum([expected[sig] for sig in expected])
    logging.debug('reconstruction: {}'.format(['{}: {:.2f}'.format(sig, expected[sig] / total) for sig in expected]))

    if plot is not None:
      logging.info('generating plot...')
      plot_sbs_signature(counts, plot, contexts, sigs.keys())

    logging.info('done: skipped no context %i skipped wrong context %i added %i', skipped_no_context, skipped_wrong_context, added)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Assign signature probbilities to variants')
  parser.add_argument('--vcf', required=False, help='annotated vcf file')
  parser.add_argument('--plot', required=False, help='plot context breakdowns')
  parser.add_argument('--signatures', required=True, help='calculated signatures')
  parser.add_argument('--definition', required=True, help='signature definition')
  parser.add_argument('--threshold', required=False, default=0, type=float, help='minimum value signature')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.vcf, args.signatures, args.definition, args.threshold, args.plot)
