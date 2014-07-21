'''
Enables derivation of BED length distributions given their tissue-name and
tissue-class.
'''

import argparse
from src.ioutils import parse_config
from src.model import BEDFileFactory


def main(args):
    '''
    Given a list of BEDFile objects, compute their length distribution. Such
    length and the respective tissue and class are sent to standard-output.
    @param args: dictionary of command-line arguments.
    '''

    beds = [BEDFileFactory(elem).build() for elem in parse_config(args['in'])]
    print('Length,Class,Tissue')
    for bed in beds:
        for le in list(bed.get_data()['Length']):  # iterate over each length
            print(str(le) + ',' + bed.get_class() + ',' + bed.get_tissue())

if __name__ == '__main__':
    try:
        parser = argparse.ArgumentParser()
        parser.add_argument('-in', metavar='XML', required=True,
                            help='XML configuration file [req]')
        args = vars(parser.parse_args())  # parse arguments
        main(args)
    except KeyboardInterrupt:
        print()
