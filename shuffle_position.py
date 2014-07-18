'''
Given BED files(s) and intervals in-which their positions are located within,
this script randomly shuffles each BED entry within the respective interval.
This shuffling can help shed light on the positional landscape of BED entries,
and whether specific BED entries are truly located at specific positions within
an interval
'''

import argparse
import bisect
from src.model import BEDFileFactory
from src.ioutils import parse_config


def shuffle(bed, ranges):
    start_pos = list(bed.get_data()['Start'])  # map start position to interval
    intervals = [bisect.bisect_left(ranges, pos) for pos in start_pos]
    bed.get_data()['Interval'] = intervals
    print(bed.get_data())


def main(args):
    ranges = list(range(args['start'], args['stop'], args['step']))[1:]
    beds = [BEDFileFactory(elem).build() for elem in parse_config(args['in'])]
    for bed in beds:
        shuffle(bed, ranges)

    print(args['start'], args['stop'], args['step'])
    print(ranges)
#     for num in [22, 14, 555, 2, 0, 1, 120000]:
#         print(bisect.bisect_left(ranges, num))

if __name__ == '__main__':
    try:
        parser = argparse.ArgumentParser()
        parser.add_argument('-in', metavar='XML', required=True,
                            help='XML configuration file [req]')
        parser.add_argument('-dir', metavar='DIR', default='./out/',
                            help='Output folder to save BED files [./out/]')
        parser.add_argument('-start', metavar='INT', default=0,
                            help='Start base-pair [0]')
        parser.add_argument('-stop', metavar='INT', default=int(1e7),
                            help='End base-pair [1,000,000]')
        parser.add_argument('-step', metavar='INT', default=int(0.5e6),
                            help='Base-pair interval size [10,000,000]')
        args = vars(parser.parse_args())  # parse arguments
        main(args)
    except KeyboardInterrupt:
        print()
