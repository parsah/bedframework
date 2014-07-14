'''
Enables derivation of BED length distributions given their tissue-name and
tissue-class.
'''

import argparse
from src import parser


def compute_distributions(beds, is_divisible):
    '''
    Given a list of BEDFile objects, compute their length distribution. Such
    length and the respective tissue and class are sent to standard-output.
    @param beds: collection of BEDFile objects.
    '''
    print('Length , Class , Tissue')  # header
    for bed in beds:
        bed_df = parser.parse_abstract_bed(bed.get_filename())
        lengths = bed_df[2] - bed_df[1]
        for l in list(lengths):  # iterate over each length, write properties
            print(l, ',', bed.get_tissue_class(), ',', bed.get_tissue_name())

if __name__ == '__main__':
    try:
        argsparser = argparse.ArgumentParser()
        argsparser.add_argument('-in', metavar='XML', required=True,
                            help='XML configuration file [req]')
        args = vars(argsparser.parse_args())  # parse arguments
        bed_list = parser.parse_config(xml=args['in'])  # parse config
        compute_distributions(beds=bed_list, is_divisible=args['divisible'])
    except KeyboardInterrupt:
        print()
