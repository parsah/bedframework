'''
Computes a measure of similarity given BED entries of tissues and lengths. The
resultant p-value is therefore used to contrast whether vectors of differing
lengths are statistically significant.
'''

import argparse
import numpy
import sys
from itertools import product, chain
from src.ioutils import parse_config
from src.model import BEDFileFactory
from src.config import TISSUE_SPEC, UBIQUITOUS
from pandas import concat, DataFrame
from scipy.stats import ranksums


def wilcox_test(x, y):
    '''
    Performs the Wilcoxon-Ranked Sums Test.
    @param x: collection of numerics.
    @param y: collection of numerics.
    @return: float referencing the test p-value.
    '''

    pval = ranksums(x, y)[-1].astype('float64')
    if numpy.isnan(pval):  # replace NaN with a poor p-value.
        pval = 1.0
    if pval == 0.0:  # sanity checks to ensure all values are non-zero.
        pval = numpy.finfo(numpy.float64).tiny.astype('float')
    return pval


def main(args):
    '''
    Given a list of BEDFile objects, compute their length distribution. Such
    length and the respective tissue and class are sent to standard-output.
    @param args: dictionary of command-line arguments.
    '''

    d = DataFrame()
    df = concat([BEDFileFactory(elem).build().get_data()
                 for elem in parse_config(args['in'])])
    ts, ub = df[df['Class'] == TISSUE_SPEC], df[df['Class'] == UBIQUITOUS]
    combs = product(*[df['Tissue'].unique(),
                                df['Length'].unique(), df['Length'].unique()])
    for t, l_ts, l_ub in list(combs):  # iterate over tissue, lengths
        x = list(chain.from_iterable(
            list(ts[(ts['Length'] == l_ts) & (ts['Tissue'] == t)]['Vectors'])))
        y = list(chain.from_iterable(
            list(ub[(ub['Length'] == l_ub) & (ub['Tissue'] == t)]['Vectors'])))
        d = d.append({'Tissue': t, 'PValue': wilcox_test(x, y),
                      'Length_TS': l_ts, 'Length_UB': l_ub}, ignore_index=True)
    d.to_csv(sys.stdout, index=False)

if __name__ == '__main__':
    try:
        parser = argparse.ArgumentParser()
        parser.add_argument('-in', metavar='XML', required=True,
                            help='XML configuration file [req]')
        args = vars(parser.parse_args())  # parse arguments
        main(args)
    except OSError as e:
        print(e)
    except KeyboardInterrupt:
        print()
