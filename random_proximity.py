'''
'''

import argparse
import sys
from src.model import BEDFileFactory
from src.ioutils import parse_config, build_random_bed
from annotated_proximity import map_features


def main(args):
    '''
    @param args: dictionary of command-line arguments.
    '''

    beds = [BEDFileFactory(e).build() for e in parse_config(args['in'])]
    if len(beds) != 2:
        raise IOError('2x BED files needed; 1x tissue-specific, 1x ubiquitous')
    bed_ts = [b for b in beds if b.get_class() == 'Tissue-Specific'][0]
    bed_ub = [b for b in beds if b.get_class() == 'Ubiquitous'][0]
    map_features(bed_ub, bed_ts.get_file())  # map TS BEDs onto ubiquitous
    randbed = build_random_bed(args['fasta'], bed_ts.get_data().shape[0])
    map_features(bed_ub, randbed)
    bed_ub.get_data().to_csv(sys.stdout)

if __name__ == '__main__':
    try:
        parser = argparse.ArgumentParser()
        parser.add_argument('-in', metavar='XML', required=True,
                            help='XML configuration file [req]')
        parser.add_argument('-fasta', metavar='FASTA', required=True,
                            help='FASTA for creating random BED entries [req]')
        args = vars(parser.parse_args())  # parse arguments
        main(args)
    except KeyboardInterrupt:
        print()
