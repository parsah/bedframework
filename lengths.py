'''
Enables derivation of BED length distributions given their tissue-name and
tissue-class.
'''

import argparse
import sys
from src.model import BEDFileFactory
from src.ioutils import parse_config
from pandas import concat


def main(args):
    '''
    Given a list of BEDFile objects, compute their length distribution. Such
    length and the respective tissue and class are sent to standard-output.
    @param args: dictionary of command-line arguments.
    '''

    df = concat([BEDFileFactory(e).build().get_data()
                 for e in parse_config(args['in'])])
    df[['Length', 'Class', 'Tissue']].to_csv(sys.stdout, index=False)


if __name__ == '__main__':
    try:
        parser = argparse.ArgumentParser()
        parser.add_argument('-in', metavar='XML', required=True,
                            help='XML configuration file [req]')
        args = vars(parser.parse_args())  # parse arguments
        main(args)
    except KeyboardInterrupt:
        print()
