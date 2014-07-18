'''
Enables derivation of BED length distributions given their tissue-name and
tissue-class.
'''

import argparse
import pandas
from src.ioutils import parse_config
from src.model import BEDFileFactory


def main(args):
    '''
    Given a list of BEDFile objects, compute their length distribution. Such
    length and the respective tissue and class are sent to standard-output.
    @param args: dictionary of command-line arguments.
    '''

    df_out = pandas.DataFrame()  # create empty data-frame to save results into
    beds = [BEDFileFactory(elem).build() for elem in parse_config(args['in'])]
    for bed in beds:
        for le in list(bed.get_data()['Length']):  # iterate over each length
            df_out = df_out.append({'Length': le, 'Class': bed.get_class(),
                            'Tissue': bed.get_tissue()}, ignore_index=True)
    df_out.to_csv(args['out'])  # save output

if __name__ == '__main__':
    try:
        parser = argparse.ArgumentParser()
        parser.add_argument('-in', metavar='XML', required=True,
                            help='XML configuration file [req]')
        parser.add_argument('-out', metavar='CSV', default='./out.csv',
                            help='Output filename [./out.csv]')
        args = vars(parser.parse_args())  # parse arguments
        main(args)
    except KeyboardInterrupt:
        print()
