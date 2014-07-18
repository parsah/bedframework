'''
Computes a measure of similarity given BED entries of tissues and lengths. The
resultant p-value is therefore used to contrast whether vectors of differing
lengths are statistically significant.
'''

import argparse
import itertools
import numpy
import pandas
from src.ioutils import parse_config
from src.model import BEDFileFactory
from src.config import TISSUE_SPEC, UBIQUITOUS
from pandas import concat
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
        pval = numpy.finfo(numpy.float64).tiny
    return pval


def main(args):
    '''
    Given a list of BEDFile objects, compute their length distribution. Such
    length and the respective tissue and class are sent to standard-output.
    @param args: dictionary of command-line arguments.
    '''

    df_out = pandas.DataFrame()  # create empty data-frame to save results into
    beds = [BEDFileFactory(elem).build() for elem in parse_config(args['in'])]
    df = concat([bed.get_data() for bed in beds])  # add to main data-frame
    ts = df[df['Class'] == TISSUE_SPEC]
    ub = df[df['Class'] == UBIQUITOUS]
    combs = itertools.product(*[df['Tissue'].unique(),
                                df['Length'].unique(), df['Length'].unique()])
    for tissue, len_ts, len_ub in list(combs):  # iterate over tissue, lengths
        x = ts[(ts['Length'] == len_ts) & (ts['Tissue'] == tissue)]
        x = list(itertools.chain.from_iterable(list(x['Vectors'])))  # flatten
        y = ub[(ub['Length'] == len_ub) & (ub['Tissue'] == tissue)]
        y = list(itertools.chain.from_iterable(list(y['Vectors'])))  # flatten
        df_out = df_out.append({'Tissue': bed.get_tissue(), 'Len_Ubiq': len_ub,
                                'PValue': wilcox_test(x, y),
                                'Len_Tiss': len_ts}, ignore_index=True)
    df_out.to_csv(args['out'])  # save output

if __name__ == '__main__':
    try:
        parser = argparse.ArgumentParser()
        parser.add_argument('-in', metavar='XML', required=True,
                            help='XML configuration file [req]')
        parser.add_argument('-out', metavar='CSV', default='./out.csv',
                            help='Output filename [./out.csv]')
        args = vars(parser.parse_args())  # parse arguments
        main(args)
    except OSError as e:
        print(e)
    except KeyboardInterrupt:
        print()
