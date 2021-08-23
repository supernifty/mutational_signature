#!/usr/bin/env python
'''
  generate de novo signatures from count files
  TODO heatmap as well as pca
'''

import argparse
import collections
import csv
import logging
import random
import sys

import matplotlib
import matplotlib.pyplot as plt

import numpy as np
import sklearn.cluster
import sklearn.decomposition
import sklearn.metrics

colors = ["#c0c0c0", "#41ac2f", "#7951d0", "#73d053", "#b969e9", "#91ba2c", "#b4b42f", "#5276ec", "#daae36", "#9e40b5", "#43c673", "#dd4cb0", "#3d9332", "#de77dd", "#7bad47", "#9479e8", "#487b21", "#a83292", "#83c67d", "#664db1", "#e18d28", "#588de5", "#e2672a", "#34c7dd", "#cf402b", "#5acdaf", "#d74587", "#428647", "#7b51a7", "#b4ba64", "#646cc1", "#a27f1f", "#3b63ac", "#dca653", "#505099", "#7d8529", "#bf8ade", "#516615", "#b65da7", "#57a87a", "#c84249", "#37b5b1", "#a14622", "#58b5e1", "#ba6e2f", "#589ed8", "#e98261", "#3176ae", "#656413", "#a19fe2", "#756121", "#7e4a8d", "#326a38", "#dd8abf", "#1a6447", "#e78492", "#30876c", "#9d4d7c", "#919d5b", "#9d70ac", "#5b6f34", "#65659c", "#c9a865", "#a1455d", "#5e622c", "#b66057", "#dca173", "#855524", "#9f7846", "#7951d0", "#73d053", "#b969e9", "#91ba2c", "#3656ca", "#b4b42f", "#5276ec", "#daae36", "#9e40b5", "#de3860", "#41ac2f", "#7951d0", "#73d053", "#b969e9", "#91ba2c", "#9e40b5", "#43c673", "#dd4cb0", "#3d9332", "#de77dd", "#7bad47", "#9479e8", "#487b21", "#a83292", "#83c67d", "#664db1", "#b66057", "#dca173", "#855524", "#9f7846", "#7951d0", "#73d053", "#b969e9", "#91ba2c", "#3656ca", "#b4b42f", "#5276ec", "#daae36", "#9e40b5", "#de3860", "#41ac2f", "#7951d0", "#73d053", "#b969e9", "#91ba2c", "#9e40b5", "#43c673", "#dd4cb0", "#3d9332", "#de77dd", "#7bad47", "#9479e8", "#487b21", "#a83292", "#83c67d", "#664db1"]

SUFFIXES = ['', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j']

def cosine(a, b):
  # make vectors
  keys = set(list(a.keys()) + list(b.keys()))
  av = np.array([a.get(k, 0) for k in sorted(keys)])
  bv = np.array([b.get(k, 0) for k in sorted(keys)])
  return av.dot(bv) / (np.linalg.norm(av) * np.linalg.norm(bv))

def converted_context(context):
  if '>' in context:
    return context
  else:
    return '{}{}{}>{}'.format(context[0], context[1], context[3], context[2]) # ACAA -> ACA>A

def prepare_bootstrap(X):
  X_bs = []
  for row in X:
    row_bs = [0] * len(row)
    mutation_count = sum(row)
    X_bs.append([a for s in [[i] * x for i, x in enumerate(row)] for a in s])
  return X_bs

def bootstrap(X, all_mutations):
  # for each row, keep overall count the same but probability of a mutation each time is that of the count / total count
  X_bs = []
  for row, mutations in zip(X, all_mutations):
    row_bs = [0] * len(row)
    for _ in range(len(mutations)):
      row_bs[mutations[random.randint(0, len(mutations) - 1)]] += 1
    X_bs.append(row_bs)
  return X_bs

def convert_context(c):
  '''
    CGT>A -> CGAT
  '''
  return '{}{}{}{}'.format(c[0], c[1], c[4], c[2])

def generate(counts, components, use_proportion, stats, min_mutation_proportion, compare_to_fn, compare_to_matches, bootstraps, bootstrap_plot, component_plot, signature_file):
  logging.info('processing %i count files to generate %s components...', len(counts), components)

  contexts = set()
  samples = collections.defaultdict(dict)
  for count_file in counts:
    logging.info('processing %s...', count_file)
    #sample = count_file.split('/')[-1].split('.')[0]
    sample = count_file.split('/')[-1]
    for row in csv.DictReader(open(count_file, 'r'), delimiter='\t'):
      # TODO filter on sbs (for now)
      if len(row['Variation']) != 5 or row['Variation'][3] != '>': # ACC>T
        logging.debug('skipping context %s', row['Variation'])
        continue
      # expct Variation Count Probability
      contexts.add(row['Variation'])
      if use_proportion:
        samples[sample][row['Variation']] = float(row['Probability'])
      else:
        samples[sample][row['Variation']] = int(row['Count'])
  logging.info('found %i contexts in %i samples: %s', len(contexts), len(counts), contexts)

  X = []
  contexts_list = sorted(list(contexts))
  for sample in sorted(samples.keys()): # each relevant sample
    features = [samples[sample].get(context, 0) for context in contexts_list]
    X.append(features)
  
  # samples are rows
  # contexts are columns

  X = np.array(X)

  # remove rows with too few mutations
  total_mutation_count = np.sum(X)
  logging.info('%i mutations in total', total_mutation_count)
  if min_mutation_proportion > 0:
    threshold = total_mutation_count * min_mutation_proportion
    logging.info('removing mutations types with less than %i mutations in total', threshold)
    so_far = 0
    to_delete = []
    sums = X.sum(axis=0)
    while True:
      # find the minimum sum not already in the list
      candidates = []
      indexes = []
      for i, x in enumerate(sums):
        if i not in to_delete:
          candidates.append(x)
          indexes.append(i)

      candidate = min(candidates)
      if candidate < threshold - so_far:
        so_far += candidate
        to_delete.append(indexes[candidates.index(candidate)])
      else:
        break
    exclude = set([contexts_list[idx] for idx in to_delete])
    logging.info('%i contexts to remove: %s', len(to_delete), ' '.join(sorted(list(exclude))))
    X = np.delete(X, to_delete, axis=1) 
    
    # all the contexts are messed up
    [contexts_list.remove(x) for x in exclude]

  logging.info('decomposing %i samples with %i contexts to %s components...', X.shape[0], X.shape[1], components)

  if stats:
    stats_fh = open(stats, 'w')
    stats_fh.write('Components\tFrobenius\tRSS\tExplainedVariance\n')

  if compare_to_fn:
    compare_to = {}
    for row in csv.DictReader(open(compare_to_fn, 'r'), delimiter='\t'):
      vals = {converted_context(key): float(row[key]) for key in row if key != 'Sig'}
      compare_to[row['Sig']] = vals
    logging.info('%i signatures as comparison', len(compare_to))

  logging.info('generating %i bootstraps...', bootstraps)
  Xs = []
  if bootstraps > 0:
    bootstrap_array = prepare_bootstrap(X)
    for i in range(bootstraps):
      if i % 10 == 0:
        logging.info('%i bootstraps done...', i)
      Xs.append(bootstrap(X, bootstrap_array))
  else:
    Xs.append(X)

  if signature_file is not None:
    # header
    logging.info('writing signatures to %s...', signature_file)
    sig_fh = open(signature_file, 'w')
    sig_fh.write('Name\t{}\n'.format('\t'.join([convert_context(c) for c in contexts_list])))

  #for components in range(3,16): # components is a list of signature-number to try
  component_results = []
  for component in components:
    logging.info('finding %i signatures', component)
    bootstrap_results = []
    reconstruction_errors = []
    for ix, X in enumerate(Xs):
      model = sklearn.decomposition.NMF(n_components=component, init="nndsvd")
      W = model.fit_transform(X)
      H = model.components_
      for sig in H:
        bootstrap_results.append(sig / sum(sig)) # normalized

      # stats
      estimate = np.dot(W, H)
      difference = X - estimate
      rss = (difference ** 2).sum() # frobenius is just the sqrt of this
      explained_variance = sum(np.var(W, axis=0)) # components
      reconstruction_errors.append(model.reconstruction_err_) # TODO normalize?
      logging.debug('signature count %i with %i features and %i samples gives frobenius error: %.2f rss: %.2f explained_variance: %.2f', component, len(contexts), len(X), model.reconstruction_err_, rss, explained_variance)
      if stats:
        stats_fh.write('{}\t{:.1f}\t{:.1f}\t{:.1f}\n'.format(component, model.reconstruction_err_, rss, explained_variance))
     
    if bootstraps > 0:
      logging.info('processing bootstrapped results: clustering %i signatures...', len(bootstrap_results))
      bootstrap_results = np.array(bootstrap_results)
      # perform clustering by cosine similarity and measure stability, etc
      clusterer = sklearn.cluster.KMeans(n_clusters=component)
      labels = clusterer.fit_predict(bootstrap_results)
      H_result = clusterer.cluster_centers_
      # measure reproducibility
      score = sklearn.metrics.silhouette_score(bootstrap_results, labels=labels, metric='cosine')
      component_results.append({'silhouette': score, 'frobenius': np.mean(model.reconstruction_err_)})
      logging.info('signature count %i clustered gives %.2f', component, score)
    else:
      H_result = H

    # normalize
    H_norm = []
    for row in H_result:
      H_norm.append(row / sum(row))

    H_result = H_norm

    if compare_to_fn is not None:
      logging.info('comparing to %s...', compare_to_fn)
      sys.stdout.write('\n*Comparison*\nGeneratedSignature\tProvidedSignature\tSimilarity\n')
      result = collections.defaultdict(dict)
      for sig_num, sig in enumerate(H_result):
        total = sum(sig)
        our_sig = {context: value for context, value in zip(contexts_list, sig)}
        for comp in compare_to:
          result['Signature.{}'.format(sig_num + 1)][comp] = cosine(our_sig, compare_to[comp])
      logging.info('comparing to %s: done', compare_to_fn)
    
      # write top n only
      highlight = set()
      for our_sig in result:
        for key, value in sorted(result[our_sig].items(), key=lambda item: -item[1])[:compare_to_matches]:
          sys.stdout.write('{}\t{}\t{:.2f}\n'.format(our_sig, key, value))
          highlight.add(key)

    if signature_file is not None:
      logging.info('writing component %i signatures to %s...', component, signature_file)
      closest_seen = set()
      for sig_num, sig in enumerate(H_result):
        total = sum(sig)
        if compare_to_fn is None:
          sig_fh.write('SigN{}.{}\t{}\n'.format(component, sig_num + 1, '\t'.join(['{:.6f}'.format(x / total) for x in sig])))
        else:
          closest = sorted(result['Signature.{}'.format(sig_num + 1)].items(), key=lambda item: -item[1])[0]
          for suffix in SUFFIXES:
            closest_cand = ''.join([closest[0], suffix]) 
            if closest_cand not in closest_seen:
              closest_seen.add(closest_cand)
              break
          sig_fh.write('SigN{}.{}\t{}\n'.format(component, closest_cand, '\t'.join(['{:.6f}'.format(x / total) for x in sig])))

    if bootstrap_plot is not None and bootstraps > 0:
      bootstrap_fn = bootstrap_plot.replace('.png', '_{}.png'.format(component))
      logging.info('plotting to %s...', bootstrap_fn)
      pca = sklearn.decomposition.PCA(n_components=2).fit(bootstrap_results)
      xy = pca.transform(bootstrap_results)
      matplotlib.style.use('seaborn')
      fig = plt.figure()
      ax = fig.add_subplot(111)
      ax.scatter([x[0] for x in xy], [x[1] for x in xy], c=[colors[x] for x in clusterer.labels_])
      # label centers
      centers = pca.transform(H_result)
      ax.scatter([x[0] for x in centers], [x[1] for x in centers], marker='+', c='k')
      for i, center in enumerate(centers):
        ax.annotate(i + 1, (center[0], center[1]))
      if compare_to_fn is not None:
        # add known sigs
        reduced = []
        for sig in compare_to:
          features = [compare_to[sig][feature] for feature in contexts_list]
          transformed = pca.transform([features])
          reduced.append(transformed[0])
          if sig in highlight:
            c = 'k'
            a = 1.0
          else:
            c = '#a0a0a0'
            a = 0.5
          ax.annotate(sig, (reduced[-1][0], reduced[-1][1]), color=c, alpha=a)
        ax.scatter([x[0] for x in reduced], [x[1] for x in reduced], marker='.', c='k')
      plt.tight_layout()
      plt.savefig(bootstrap_fn, transparent=False)
    logging.info('processing bootstrapped results: done')

  # plot component reconstruction
  if component_plot is not None:
    matplotlib.style.use('seaborn')
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.grid(axis='y')
    ax1.plot(components, [x['silhouette'] for x in component_results], c='red', label='Average Silhouette Width', linestyle='-', marker='o') # -1 to 1
    ax1.set_ylabel('Average Silhouette')
    ax1.set_xlabel('Number of Signatures')
    ax2 = ax1.twinx()
    ax2.grid(False)
    ax2.plot(components, [x['frobenius'] for x in component_results], c='blue', label='Average Reconstruction Error', linestyle='-', marker='v')
    ax2.set_ylabel('Reconstruction Error')

    lines, labels = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax2.legend(lines + lines2, labels + labels2, loc=0)

    plt.tight_layout()
    plt.savefig(component_plot, transparent=False)

  logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Find mutational signatures')
  parser.add_argument('--counts', required=True, nargs='+', help='context counts for each sample')
  parser.add_argument('--components', required=True, nargs='+', type=int, help='number of decomposition components')
  parser.add_argument('--use_proportion', action='store_true', help='use variation proportion instead of count')
  parser.add_argument('--stats', required=False, help='write stats to file')
  parser.add_argument('--compare_to', required=False, help='find cosine similarity to existing signature definitions')
  parser.add_argument('--compare_to_matches', required=False, default=3, type=int, help='show the top matches for each signature')
  # data massaging based on "Deciphering signatures of mutational processesoperative in human cancer"
  parser.add_argument('--min_mutation_proportion', required=False, default=0, type=float, help='ignore any mutation type with less than this proportion across all genome inputs')
  parser.add_argument('--bootstraps', required=False, default=0, type=int, help='number of bootstraps to perform')
  parser.add_argument('--bootstrap_plot', required=False, help='pca plot of bootstraps')
  parser.add_argument('--component_plot', required=False, help='plot of component quality')
  parser.add_argument('--signature_file', required=False, help='write signatures to this file')
  # logging
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  generate(args.counts, args.components, args.use_proportion, args.stats, args.min_mutation_proportion, args.compare_to, args.compare_to_matches, args.bootstraps, args.bootstrap_plot, args.component_plot, args.signature_file)
