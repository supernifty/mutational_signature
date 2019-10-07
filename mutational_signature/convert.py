#!/usr/bin/env python
'''
  converts cosmic signatures to our format
'''

import argparse
import collections
import logging
import sys

def main(conversion):
  logging.info('reading from stdin...')
  if conversion is None or conversion == 'sbs':
    # v4 to v3
    # source
    #Type,SubType,SBS1,SBS2,SBS3,SBS4,SBS5,SBS6,SBS7a,SBS7b,SBS7c,SBS7d,SBS8,SBS9,SBS10a,SBS10b,SBS11,SBS12,SBS13,SBS14,SBS15,SBS16,SBS17a,SBS17b,SBS18,SBS19,SBS20,SBS21,SBS22,SBS23,SBS24,SBS25,SBS26,SBS27,SBS28,SBS29,SBS30,SBS31,SBS32,SBS33,SBS34,SBS35,SBS36,SBS37,SBS38,SBS39,SBS40,SBS41,SBS42,SBS43,SBS44,SBS45,SBS46,SBS47,SBS48,SBS49,SBS50,SBS51,SBS52,SBS53,SBS54,SBS55,SBS56,SBS57,SBS58,SBS59,SBS60
    #C>A,ACA,0.000371161,5.06e-07,0.01737271,0.030941737,0.009810659,0.000212145,5.04e-05,0.001645611,0.004996391,4.03e-05,0.03758363,0.000539102,0.001973717,0.000112869,0.000109198,0.004525118,0.001605603,0.00087715,0.00048418,0.016685316,0.001777461,0.000584795,0.03816761,0.000947208,0.000454604,0.000146319,0.005746905,0.000572688,0.023347643,0.008321667,0.000830632,0.005158122,0.00086572,0.042561895,0.001344467,0.007371027,0.017291868,0.002755619,0.005540585,0.006795994,0.019958708,0.003933069,0.009135994,0.009048769,0.024336679,0.002082149,0.00086859,0.024404535,5.78e-18,0.006233354,0.003772972,0.069141862,0.000392435,0.011053526,0.097406364,0.124232109,0.011443862,0.003424626,0.001684535,0.004964496,0.011521043,0.011670733,0.054043866,0.002574643,0.005235164
  
    # target
    # Sig     ACAA    ACAC    ACAG    ACAT    CCAA    CCAC    CCAG    CCAT    GCAA    GCAC    GCAG    GCAT    TCAA    TCAC    TCAG    TCAT    ACGA    ACGC    ACGG    ACGT    CCGA    CCGC    CCGG    CCGT    GCGA    GCGC    GCGG    GCGT    TCGA    TCGC    TCGG    TCGT   ACTA     ACTC    ACTG    ACTT    CCTA    CCTC    CCTG    CCTT    GCTA    GCTC    GCTG    GCTT    TCTA    TCTC    TCTG    TCTT    ATAA    ATAC    ATAG    ATAT    CTAA    CTAC    CTAG    CTAT    GTAA    GTAC    GTAG    GTAT    TTAA    TTAC    TTAG    TTAT    ATCA   ATCC     ATCG    ATCT    CTCA    CTCC    CTCG    CTCT    GTCA    GTCC    GTCG    GTCT    TTCA    TTCC    TTCG    TTCT    ATGA    ATGC    ATGG    ATGT    CTGA    CTGC    CTGG    CTGT    GTGA    GTGC    GTGG    GTGT    TTGA    TTGC    TTGG    TTGT
    #Signature.1     0.011098326166  0.009149340734  0.001490070468  0.006233885236  0.006595870109  0.007342367815  0.00089284037   0.007186581642  0.008232603989  0.00575802141   0.000616335232  0.004459080311  0.012250063666  0.011162229329  0.002275495686  0.015259102491  0.001801068369  0.00258090852   0.000592548022  0.002963986287  0.001284983446  0.000702134818  0.000506289594  0.001381542722  0.000602122711  0.002393352209  2.48534e-07     0.000890080731  0.001874853199  0.002067418791  0.000304897004  0.00315157446  0.029514532745   0.014322747041  0.171646931305  0.012623763151  0.020896446965  0.01850170477   0.095577217268  0.017113307642  0.024943814154  0.027161494035  0.103570762296  0.017689854381  0.014492099634  0.017680775357  0.076002221712  0.013761704021  0.004021520333  0.002371144163  0.002810909959  0.008360909345  0.001182587416  0.001903166857  0.00148796063   0.002179344412  0.000689289439  0.000552409528  0.001200228847  0.002107136837  0.005600155423  0.00199907926   0.001090065693  0.003981022761  0.01391577303  0.0062749606     0.010137636154  0.009256316389  0.004176674882  0.005252593331  0.00701322531   0.006713813119  0.011247835116  0.006999724257  0.004977592617  0.010667406133  0.008073616351  0.004857381178  0.008325454207  0.006257105605  0.001587636423  0.001784091288  0.001385830552  0.003158539312  0.000302691186  0.00209850244   0.0015995485    0.002758537619  9.9045003e-05   0.000202365646  0.001188353185  0.000800723342  0.001397553749  0.001291736985  0.00203107688   0.00403012816
    header = sys.stdin.readline().strip('\r\n').split(',')
  
    result = collections.defaultdict(dict)
    contexts = set()
    for line in sys.stdin:
      fields = line.strip('\r\n').split(',')
      context = '{}{}{}{}'.format(fields[1][0], fields[1][1], fields[0][2], fields[1][2]) # C>A,ACA => ACAA
      contexts.add(context)
      for idx, value in enumerate(fields[2:]):
        result[header[idx + 2]][context] = value
  
    logging.info('read %i contexts', len(contexts))
    
    # write in new format
    contexts_list = sorted(list(contexts))
    sys.stdout.write('{}\t{}\n'.format('Sig', '\t'.join(contexts_list)))
    for sig in sorted(result.keys()):
      sys.stdout.write('{}\t{}\n'.format(sig, '\t'.join([result[sig][context] for context in contexts_list])))

  elif conversion == 'dbs':
    # Mutation Type,DBS1,DBS2,DBS3,DBS4,DBS5,DBS6,DBS7,DBS8,DBS9,DBS10,DBS11
    # AC>CA,4.97E-05,3.74E-04,6.35E-03,4.14E-03,2.14E-03,2.66E-03,4.59E-02,2.15E-01,2.46E-02,9.63E-03,1.37E-03
    header = sys.stdin.readline().strip('\r\n').split(',')
    result = collections.defaultdict(dict)
    contexts = set()
    for line in sys.stdin:
      fields = line.strip('\r\n').split(',') # 
      contexts.add(fields[0])
      for idx, value in enumerate(fields[1:]):
        result[header[idx + 1]][fields[0]] = value
  
    logging.info('read %i contexts', len(contexts))
    
    # write in new format
    contexts_list = sorted(list(contexts))
    sys.stdout.write('{}\t{}\n'.format('Sig', '\t'.join(contexts_list)))
    for sig in sorted(result.keys()):
      sys.stdout.write('{}\t{}\n'.format(sig, '\t'.join([result[sig][context] for context in contexts_list])))


  elif conversion == 'id':
    # Mutation Type,ID1,ID2,ID3,ID4,ID5,ID6,ID7,ID8,ID9,ID10,ID11,ID12,ID13,ID14,ID15,ID16,ID17
    # DEL_C_1_0,0.000159889,0.004824116,0.124727109,0.007249717,0.022202108,0.030506799,0.000466561,0.039827558,0.334939881,0.025105414,0.001576565,0.052787741,0.065548371,0.012016009,0.026527369337078000000000000000000000,0.000952768885842137000000000000000000,0.006117598535198470000000000000000000
    header = sys.stdin.readline().strip('\r\n').split(',')
    result = collections.defaultdict(dict)
    contexts = set()
    for line in sys.stdin:
      fields = line.strip('\r\n').split(',') # 
      contexts.add(fields[0])
      for idx, value in enumerate(fields[1:]):
        result[header[idx + 1]][fields[0]] = value
  
    logging.info('read %i contexts', len(contexts))
    
    # write in new format
    contexts_list = sorted(list(contexts))
    sys.stdout.write('{}\t{}\n'.format('Sig', '\t'.join(contexts_list)))
    for sig in sorted(result.keys()):
      sys.stdout.write('{}\t{}\n'.format(sig, '\t'.join([result[sig][context] for context in contexts_list])))

  elif conversion == 'sbstx':
    # Strand,Type,Subtype,SBS1,SBS2,SBS3,SBS4,SBS5,SBS6,SBS7a,SBS7b,SBS7c,SBS7d,SBS8,SBS9,SBS10a,SBS10b,SBS11,SBS12,SBS13,SBS14,SBS15,SBS16,SBS17a,SBS17b,SBS18,SBS19,SBS20,SBS21,SBS22,SBS23-E,SBS24,SBS25-E,SBS26,SBS27-E,SBS28,SBS29-E,SBS30,SBS31,SBS32,SBS33,SBS34,SBS35,SBS36,SBS37,SBS38,SBS39,SBS40,SBS41,SBS42-E,SBS43,SBS44,SBS45-E,SBS46-E,SBS47,SBS48,SBS49,SBS50,SBS51,SBS52,SBS53,SBS54,SBS55,SBS56,SBS57,SBS58,SBS59-E,SBS60
    # U,C>A,ACA,0.002428978,1.52E-05,0.0084202,0.013198985,0.005920439,4.31E-05,0.000112008,0.000739051,0.001299453,0.000542596,0.014899119,0.000338971,0.000649837,0.00011401,0.0002262,0.003165023,0.00070077,0.000223899,0,0.005121873,0.002103173,0.000979044,0.023890022,0.00326593,0.002635085,0.000134902,0.002757404,0.001646329,0.008658196,0.002518386,0.00140495,0.000665657,0.000698814,0.015190625,0.004706682,0.004562547,0.004733878,0.000736242,0.001664096,0.009406096,0.011051821,0.002345625,0.001423934,0.009894456,0.010829197,0.000695961,0.002380018,0.007570297,0.001835964,0.001735395,0.001315841,0.047842405,0.009351833,0.00788947,0.054581454,0.037635487,0.001461495,0.007468025,0.002856412,0.005105252,0.008890439,0.008272066,0.024372478,0.006436883,0.001256248
    # 
    # target
    # Sig     ACAAU    ACACU    ...
    #Signature.1     0.011098326166  0.009149340734  0.001490070468  ...
    header = sys.stdin.readline().strip('\r\n').split(',')
  
    result = collections.defaultdict(dict)
    contexts = set()
    for line in sys.stdin:
      fields = line.strip('\r\n').split(',')
      context = '{}{}{}{}{}'.format(fields[2][0], fields[2][1], fields[1][2], fields[2][2], fields[0]) # C>A,ACA => ACAA
      contexts.add(context)
      for idx, value in enumerate(fields[3:]):
        result[header[idx + 3]][context] = value
  
    logging.info('read %i contexts', len(contexts))
    
    # write in new format
    contexts_list = sorted(list(contexts))
    sys.stdout.write('{}\t{}\n'.format('Sig', '\t'.join(contexts_list)))
    for sig in sorted(result.keys()):
      sys.stdout.write('{}\t{}\n'.format(sig, '\t'.join([result[sig][context] for context in contexts_list])))

  else:
    logging.warn('unrecognized conversion')
  
  logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Assess MSI')
  parser.add_argument('--conversion', help='specify conversion type sbs db id sbstx')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.conversion)
