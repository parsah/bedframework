'''
A vectorized BED file is such that the very last column references a vector
of numerical values (i.e. conservation scores, or gene-expression signals).
This module enumerates each element in this respective vector such that
element x maps to integer x; for easy plotting on the x-y coordinate plane.
'''

import argparse
from src import parser


def vectorize(beds):
    '''
    Given a list of BEDFile objects, compute their length distribution. Such
    length and the respective tissue and class are sent to standard-output.
    @param beds: collection of BEDFile objects.
    '''
    for bed in beds:
        df = parser.parse_vectorized_bed(bed.get_filename())
        # TODO enumerate through each df vector

if __name__ == '__main__':
    try:
        argsparser = argparse.ArgumentParser()
        argsparser.add_argument('-in', metavar='XML', required=True,
                            help='XML configuration file [req]')
        args = vars(argsparser.parse_args())  # parse arguments
        bed_list = parser.parse_config(xml=args['in'])  # parse config
        vectorize(beds=bed_list)
    except KeyboardInterrupt:
        print()
