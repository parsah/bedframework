'''
'''

import argparse
import itertools
import numpy
from src import parser
from src.model import BEDFileFactory
from src.config import TISSUE_SPEC, UBIQUITOUS
from pandas import concat
from scipy.stats import ranksums


def wilcox_test(x, y):
    pval = ranksums(x, y)[-1].astype('float64')
    if numpy.isnan(pval):
        pval = 1.0
    if pval == 0.0:
        pval = numpy.finfo(numpy.float64).tiny
    return pval


def compute_bed_similarity(beds):
    df = concat([bed.get_data() for bed in beds])  # add to main data-frame
    lengths = sorted([int(i) for i in set(df['Length']) if i % 100 == 0])
    tissues = [i for i in set(df['Tissue'])]
    ts = df[df['Class'] == TISSUE_SPEC]
    ub = df[df['Class'] == UBIQUITOUS]
    print('Tissue,Ubiquitous,Tissue-Specific,Pvalue')
    for tissue in tissues:
        for len_ts in lengths:  # filter the respective data-frames
            for len_ub in lengths:  # filter the respective data-frames
                x = ts[(ts['Length'] == len_ts) & (ts['Tissue'] == tissue)]
                x = list(itertools.chain.from_iterable(list(x['Vectors'])))
                y = ub[(ub['Length'] == len_ub) & (ub['Tissue'] == tissue)]
                y = list(itertools.chain.from_iterable(list(y['Vectors'])))
                print(tissue + ',' + str(len_ts) + ',' + str(len_ub) + ',' +
                      str(wilcox_test(x, y)))

if __name__ == '__main__':
    try:
        argsparser = argparse.ArgumentParser()
        argsparser.add_argument('-in', metavar='XML', required=True,
                            help='XML configuration file [req]')
        args = vars(argsparser.parse_args())  # parse arguments
        elems = parser.parse_config(xml=args['in'])  # parse config
        beds = [BEDFileFactory(elem).build() for elem in elems]
        compute_bed_similarity(beds)
    except OSError as e:
        print(e)
    except KeyboardInterrupt:
        print()
