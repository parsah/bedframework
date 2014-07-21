'''
Given BED files(s) and intervals in-which their positions are located within,
this script randomly shuffles each BED entry within the respective interval.
This shuffling can help shed light on the positional landscape of BED entries,
and whether specific BED entries are truly located at specific positions within
an interval
'''

import argparse
import numpy
from os.path import basename
from bisect import bisect_left
from numpy.random import randint
from src.model import BEDFileFactory
from src.ioutils import parse_config, mkdir


def shuffle(bed, nums):
    '''
    Fits BED start indices given a fixed number of intervals, and subsequently
    generates a new start index per BED entry within each interval. Such logic
    serves the purpose of ensuring that BED entries truly exhibit
    position-specific events.
    @param bed: BEDFile object.
    @param nums: Number of intervals to create of start indices.
    '''

    start = bed.get_data()['Start']  # get all Start indices
    ranges = numpy.linspace(min(start), max(start), nums)  # fit to intervals
    bed.get_data()['Interval'] = [bisect_left(ranges, pos) for pos in start]
    new_df = ''
    for each_int in sorted(bed.get_data()['Interval'].unique()):
        data = bed.get_data()[bed.get_data()['Interval'] == each_int]
        for _, row in data.iterrows():  # per row, compute new start index
            new_start = randint(min(data['Start']) - 1, max(data['Start']))
            new_df += str(row['Chr']) + '\t' + str(new_start) + '\t' +\
                str(new_start + int(row['Length'])) + '\n'
    return new_df  # return new BED data since it is now shuffled


def main(args):
    mkdir(args['dir'])
    beds = [BEDFileFactory(elem).build() for elem in parse_config(args['in'])]
    for b in beds:
        new_bed_data = shuffle(b, args['num'])
        loc = args['dir'] + '/' + b.get_tissue() + '-' + basename(b.get_file())
        out_handle = open(loc, 'w')
        out_handle.write(new_bed_data)
        out_handle.flush()
        out_handle.close()


if __name__ == '__main__':
    try:
        parser = argparse.ArgumentParser()
        parser.add_argument('-in', metavar='XML', required=True,
                            help='XML configuration file [req]')
        parser.add_argument('-dir', metavar='DIR', default='./out/',
                            help='Output folder to save BED files [./out/]')
        parser.add_argument('-num', metavar='INT', default=10, type=int,
                            help='# intervals to fit BED start indices [10]')
        args = vars(parser.parse_args())  # parse arguments
        main(args)
    except KeyboardInterrupt:
        print()
    except OSError as e:
        print(e)