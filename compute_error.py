'''
A vectorized BED file is such that the very last column references a vector
of numerical values (i.e. conservation scores, or gene-expression signals).
This module enumerates each element in this respective vector such that
element x maps to integer x; for easy plotting on the x-y coordinate plane.
'''

import argparse
from src import parser
from src.model import BEDFileFactory
from pandas import concat
from src.config import TISSUE_SPEC, UBIQUITOUS


def compute(l):
    '''
    Given a list of BEDFile objects, compute their length distribution. Such
    length and the respective tissue and class are sent to standard-output.
    @param beds: collection of BEDFile objects.
    '''
    df_ts = concat([b.get_data() for b in l if b.get_class() == TISSUE_SPEC])
    df_ub = concat([b.get_data() for b in l if b.get_class() == UBIQUITOUS])
    for df in [df_ub]:
        lens = df['Length']  # get all the lengths for respective data-frame
        for length in lens.unique():  # pull all vectors matching the length
            print(length)


if __name__ == '__main__':
    try:
        argsparser = argparse.ArgumentParser()
        argsparser.add_argument('-in', metavar='XML', required=True,
                            help='XML configuration file [req]')
        args = vars(argsparser.parse_args())  # parse arguments
        elems = parser.parse_config(xml=args['in'])  # parse config
        beds = [BEDFileFactory(elem).build() for elem in elems]
        compute(beds)
    except KeyboardInterrupt:
        print()
