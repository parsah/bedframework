'''
Enables derivation of BED length distributions given their tissue-name and
tissue-class.
'''

import argparse
from src import parser
from src.model import BEDFileFactory


def compute_distributions(beds):
    '''
    Given a list of BEDFile objects, compute their length distribution. Such
    length and the respective tissue and class are sent to standard-output.
    @param beds: collection of BEDFile objects.
    '''
    print('Length,Class,Tissue')  # header
    for bed in beds:
        lengths = bed.get_data()['Length']
        for l in list(lengths):  # iterate over each length, write properties
            print(str(l) + ',' + bed.get_class() + ',' + bed.get_tissue())

if __name__ == '__main__':
    try:
        argsparser = argparse.ArgumentParser()
        argsparser.add_argument('-in', metavar='XML', required=True,
                            help='XML configuration file [req]')
        args = vars(argsparser.parse_args())  # parse arguments
        elems = parser.parse_config(xml=args['in'])  # parse config
        beds = [BEDFileFactory(elem).build() for elem in elems]
        compute_distributions(beds)
    except KeyboardInterrupt:
        print()
