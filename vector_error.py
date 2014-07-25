'''
A vectorized BED file is such that the very last column references a vector
of numerical values (i.e. conservation scores, or gene-expression signals).
This module enumerates each element in this respective vector such that
element x maps to integer x; for easy plotting on the x-y coordinate plane.
'''

import argparse
import numpy
import scipy.stats
import sys
from itertools import product
from pandas import DataFrame, concat
from src.ioutils import parse_config
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
    BED vectors are essentially an i x j matrix, whereby you have i enhancers,
    each being j elements (bases) long. Thus, a matrix of dimensions i x j
    can be easily used to model the mean, as well as a corresponding
    confidence interval for each jth column.
    @param args: dictionary of command-line arguments.
    '''

    df = concat([BEDFileFactory(e).build().get_data()
                 for e in parse_config(args['in'])])
    combs = product(*[set(df['Tissue']), set(df['Length']), set(df['Class'])])
    data = DataFrame()
    for t, le, c in list(combs):  # loop combinations of tissue, length, class
        mat = numpy.matrix(df[(df['Length'] == le) & (df['Tissue'] == t) &
                           (df['Class'] == c)]['Vectors'].tolist())
        for num in range(mat.shape[1]):  # iterate over columns, get interval
            mu, upr, lwr = confidence_interval(mat[:, num], args['conf'])
            data = data.append({'Tissue': t, 'Class': c, 'Length': le,
                                'Index': num + 1, 'Mu': mu, 'Upr': upr,
                                'Lwr': lwr}, ignore_index=True)
    data.to_csv(sys.stdout, index=False)

if __name__ == '__main__':
    try:
        parser = argparse.ArgumentParser()
        parser.add_argument('-in', metavar='XML', required=True,
                            help='XML configuration file [req]')
        parser.add_argument('-conf', metavar='FLOAT', default=0.95,
                            type=float, help='Confidence interval [0.95]')
        args = vars(parser.parse_args())  # parse arguments
        main(args)
    except KeyboardInterrupt:
        print()
