'''
Given a BED file and a corresponding set of all genes, this module helps
map each BED entry to its most proximal gene. The bulk of such logic is made
possible using the BEDtools package.
'''

import argparse
import numpy as np
from bisect import bisect_left
from subprocess import Popen, PIPE
from src.ioutils import parse_config
from src.model import BEDFileFactory


def get_gtf_positions(l):
    return l[3], l[4]


def get_bed_positions(l):
    return l[1], l[2]


def map_features(b, gtf):
    '''
    Find the most proximal GTF entry given a BEDFile object. Each BED entry
    is iteratively queried to find the most proximal GTF entry, resulting in
    a distance (base-pair; integer) separating this entry and its respective
    GTF entry. This GTF file, however, must have a 'gene_id' attribute so that
    the most proximal GTF entry can explicitly be identified.
    @param b: BEDFile object.
    @param gtf: GTF filename
    @return: list referencing most proximal GTF entry for each BED entry.
    '''

    hit_start, hit_end, distance = [], [], []
    fname = b.get_file()
    n_cols = len(open(fname).readline().split('\t'))  # get # columns of BED
    cmd = 'bedtools closest -io -d -a ' + fname + ' -b ' + gtf + ' -t first'
    out, err = Popen(cmd, shell=True, stdout=PIPE).communicate()
    if err:  # stdout, stderr are outputs; stderr must be None
        raise IOError('Error found following mapping features to BED')
    for l in out.splitlines():
        l = l.decode('utf-8').split('\t')
        if gtf.endswith('bed'):
            start, end = get_bed_positions(l[n_cols:])
        elif gtf.endswith('gtf'):
            start, end = get_gtf_positions(l[n_cols:])
        hit_start.append(start)  # append the start and end bases
        hit_end.append(end)
        distance.append(int(l[-1]))  # add most proximal distance
    b.get_data()['Hit_Start'] = hit_start  # add new columns to data-structure
    b.get_data()['Hit_End'] = hit_end  # add new columns to data-structure
    b.get_data()['Distance'] = distance
    ints = list(range(0, 200000, 10000))  # fit distances to intervals
    b.get_data()['Interval'] = [bisect_left(ints, pos) for pos in distance]


def map_signals(b):
    '''
    '''
    signal = []
    for _, row in b.get_data().iterrows():  # row-number, row, respectively
        reps = []  # store bigwig replicate-averages across bigwig files
        for bigwig in b.get_bigwigs():  # compute BED average given bigwig
            cmd = 'echo ' + str(row['Chr']) + ' ' + str(row['Hit_Start']) +\
                ' ' + str(row['Hit_End']) + ' 1 | bigWigAverageOverBed ' +\
                bigwig + ' stdin stdout'
            out, err = Popen(cmd, shell=True, stdout=PIPE).communicate()
            if err:  # stdout, stderr are outputs; stderr must be None
                raise IOError('Error found following mapping features to BED')
            reps.append(float(out.decode('utf-8').strip().split('\t')[-1]))
        signal.append(np.mean(reps))  # append the mean
    b.get_data()['Signal'] = signal


def main(args):
    beds = [BEDFileFactory(elm).build() for elm in parse_config(args['in'])]
    bed_objs = [b for b in beds]  # merge BEDs into one structure
    print('Chr,Len,Start,End,Tissue,Class,Dist,Inter,Hit_Start,Hit_End,Signal')
    for b in bed_objs:
        map_features(b, args['annot'])
        if args['signals']:  # map signals (BigWig files), if need-be
            map_signals(b)
        else:
            b.get_data()['Signal'] = 'None'  # add invalid values if no signal
        data = b.get_data()
        for _, row in data.iterrows():
            if row['Distance'] >= 0:  # BED entry lacks proximal, i.e. chrX
                print(row['Chr'] + ',' + str(row['Length']) + ',' +
                      str(row['Start']) + ',' + str(row['End']) + ',' +
                      str(row['Tissue']) + ',' + str(row['Class']) + ',' +
                      str(row['Distance']) + ',' + str(row['Interval']) + ',' +
                      str(row['Hit_Start']) + ',' + str(row['Hit_End']) + ',' +
                      str(row['Signal']))

if __name__ == '__main__':
    try:
        parser = argparse.ArgumentParser()
        parser.add_argument('-in', metavar='XML', required=True,
                            help='XML configuration file [req]')
        parser.add_argument('-annot', metavar='FILE', required=True,
                            help='BED/GTF annotations file [req]')
        parser.add_argument('--signals', action='store_true', default=False,
                            help='Map BED bigwigs signals to features [false]')
        args = vars(parser.parse_args())
        main(args)
    except KeyboardInterrupt:
        print()
