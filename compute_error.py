'''
A vectorized BED file is such that the very last column references a vector
of numerical values (i.e. conservation scores, or gene-expression signals).
This module enumerates each element in this respective vector such that
element x maps to integer x; for easy plotting on the x-y coordinate plane.
'''

import argparse
import numpy as np
import scipy.stats
import itertools
from pandas import concat
from src.ioutils import parse_config, as_delim
from src.model import BEDFileFactory


def confidence_interval(x, ci=0.95):
    '''
    Computes a confidence interval given a numerical vector.
    @param x: numerical vector.
    @param ci: Confidence interval; default = 0.95
    @return: tuple given mean, upper, and lower bounds given the vector.
    '''

    low_per = 100 * (1 - ci) / 2.
    high_per = 100 * ci + low_per
    mu = x.mean()
    lwr, upr = scipy.stats.scoreatpercentile(x, [low_per, high_per])
    return mu, upr, lwr  # return the mean and its respective bounds.


def main(args):
    '''
    Given a list of BEDFile objects, compute their length distribution. Such
    length and the respective tissue and class are sent to standard-output.
    @param args: dictionary of command-line arguments.
    '''

    print(as_delim('Tissue', 'Class', 'Length', 'Index', 'Mu', 'Upr', 'Lwr'))
    beds = [BEDFileFactory(elm).build() for elm in parse_config(args['in'])]
    df = concat([b.get_data() for b in beds])  # merge BEDs into one structure
    combs = itertools.product(*[df['Tissue'].unique(),
                                df['Length'].unique(), df['Class'].unique()])
    for t, le, c in list(combs):  # loop combinations of tissue, length, class
        mat = np.matrix(df[(df['Length'] == le) & (df['Tissue'] == t) &
                           (df['Class'] == c)]['Vectors'].tolist())
        for num in range(mat.shape[1]):  # iterate over columns, get interval
            mu, upr, lwr = confidence_interval(mat[:, num], args['conf'])
            print(as_delim(t, c, le, num + 1, mu, upr, lwr))

if __name__ == '__main__':
    try:
        argsparser = argparse.ArgumentParser()
        argsparser.add_argument('-in', metavar='XML', required=True,
                            help='XML configuration file [req]')
        argsparser.add_argument('-conf', metavar='FLOAT', default=0.95,
                            type=float, help='Confidence interval [0.95]')
        args = vars(argsparser.parse_args())  # parse arguments
        main(args)
    except KeyboardInterrupt:
        print()
