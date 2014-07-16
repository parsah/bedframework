'''
A vectorized BED file is such that the very last column references a vector
of numerical values (i.e. conservation scores, or gene-expression signals).
This module enumerates each element in this respective vector such that
element x maps to integer x; for easy plotting on the x-y coordinate plane.
'''

import argparse
import numpy as np
import scipy.stats
from src import parser
from src.model import BEDFileFactory
from pandas import concat


def confidence_interval(x, ci=0.95):
    low_per = 100 * (1 - ci) / 2.
    high_per = 100 * ci + low_per
    mu = x.mean()
    lwr, upr = scipy.stats.scoreatpercentile(x, [low_per, high_per])
    return mu, upr, lwr


def compute(l):
    '''
    Given a list of BEDFile objects, compute their length distribution. Such
    length and the respective tissue and class are sent to standard-output.
    @param beds: collection of BEDFile objects.
    '''
    df = concat([b.get_data() for b in l])
    print('Tissue,Class,Length,Index,Mean,Upper,Lower')
    for name in df['Tissue'].unique():
        for l in df['Length'].unique():  # pull all vectors matching the length
            for tiss_class in df['Class'].unique():
                mat = np.matrix(df[(df['Length'] == l) &
                                   (df['Tissue'] == name) &
                                   (df['Class'] == tiss_class)]['Vectors'].
                                tolist())
                for num in range(mat.shape[1]):  # iterate over each column
                    a_column = mat[:, num]  # derive confidence intervals
                    mu, upr, lwr = confidence_interval(a_column, 0.90)
                    print(name + ',' + tiss_class + ',' + str(l) + ',' +
                            str(num + 1) + ',' + str(mu) + ',' + str(upr) +
                            ',' + str(lwr))

if __name__ == '__main__':
    try:
        argsparser = argparse.ArgumentParser()
        argsparser.add_argument('-in', metavar='XML', required=True,
                            help='XML configuration file [req]')
        argsparser.add_argument('-conf', metavar='FLOAT', default=0.95,
                            type=float, help='Confidence interval [0.95]')
        args = vars(argsparser.parse_args())  # parse arguments
        elems = parser.parse_config(xml=args['in'])  # parse config
        beds = [BEDFileFactory(elem).build() for elem in elems]
        compute(beds)
    except KeyboardInterrupt:
        print()
