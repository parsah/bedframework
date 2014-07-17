'''
Enables derivation of BED length distributions given their tissue-name and
tissue-class.
'''

import argparse
from src.ioutils import parse_config, as_delim
from src.model import BEDFileFactory


def main(args):
    '''
    Given a list of BEDFile objects, compute their length distribution. Such
    length and the respective tissue and class are sent to standard-output.
    @param args: dictionary of command-line arguments.
    '''

    print(as_delim('Length', 'Class', 'Tissue'))  # header
    beds = [BEDFileFactory(elem).build() for elem in parse_config(args['in'])]
    for bed in beds:
        for le in list(bed.get_data()['Length']):  # iterate over each length
            print(as_delim(le, bed.get_class(), bed.get_tissue()))

if __name__ == '__main__':
    try:
        argsparser = argparse.ArgumentParser()
        argsparser.add_argument('-in', metavar='XML', required=True,
                            help='XML configuration file [req]')
        args = vars(argsparser.parse_args())  # parse arguments
        main(args)
    except KeyboardInterrupt:
        print()
