'''
Given a BED file and a corresponding set of all genes, this module helps
map each BED entry to its most proximal gene. The bulk of such logic is made
possible using the BEDtools package.
'''

import argparse
import numpy as np
import sys
import os
from pandas import DataFrame
from bisect import bisect_left as bleft
from src.ioutils import parse_config, exec_closest_app, exec_average_app
from src.model import BEDFileFactory


def get_gtf_positions(l):
    '''
    List-representation of a strong produced from a GTF file.
    @param l: List object
    @return: tuple containing start and end indices of the provided list.
    '''

    return l[3], l[4]


def get_bed_positions(l):
    '''
    List-representation of a strong produced from a BED file.
    @param l: List object
    @return: tuple containing start and end indices of the provided list.
    '''

    return l[1], l[2]


def map_features(b, annot):
    '''
    Find the most proximal GTF entry given a BEDFile object. Each BED entry
    is iteratively queried to find the most proximal GTF entry, resulting in
    a distance (base-pair; integer) separating this entry and its respective
    GTF entry. This GTF file, however, must have a 'gene_id' attribute so that
    the most proximal GTF entry can explicitly be identified.
    @param b: BEDFile object.
    @param annot: BED/GTF filename
    @return: list referencing most proximal GTF entry for each BED entry.
    '''

    hit_start, hit_end, ds = [], [], []
    n_cols = len(open(b.get_file()).readline().split('\t'))  # columns of BED
    out = exec_closest_app(b.get_file(), annot)
    for l in out.splitlines():
        l = l.decode('utf-8').split('\t')
        if annot.endswith('bed'):
            start, end = get_bed_positions(l[n_cols:])
        elif annot.endswith('gtf'):
            start, end = get_gtf_positions(l[n_cols:])
        hit_start.append(start)  # append the start and end bases
        hit_end.append(end)
        ds.append(int(l[-1]))  # add most proximal distance
    b.get_data()['Hit_Start-' + os.path.basename(annot)] = hit_start
    b.get_data()['Hit_End-' + os.path.basename(annot)] = hit_end
    b.get_data()['Distance-' + os.path.basename(annot)] = ds
    b.get_data()['Intv-' + os.path.basename(annot)] =\
                [bleft(range(0, int(2e5), int(1e4)), p) for p in ds]


def map_signals(b, annot):
    '''
    Maps BigWig files onto the user-provided BED file. In many cases, multiple
    BigWig files are provided; each is iteratively mapped to each BED entry
    and the mean bigwig value is derived.
    @param b: BEDFile object
    '''

    signal = []
    for _, row in b.get_data().iterrows():  # row-number, row, respectively
        reps = []  # store bigwig replicate-averages across bigwig files
        for bw in b.get_bigwigs():  # compute BED average given bigwig
            out = exec_average_app(row['Chr'],
                                   row['Hit_Start-' + os.path.basename(annot)],
                                   row['Hit_End-' + os.path.basename(annot)],
                                   bw)
            reps.append(float(out.decode('utf-8').strip().split('\t')[-1]))
        signal.append(np.mean(reps))  # append the mean
    b.get_data()['Signal'] = signal


def main(args):
    '''
    Given a list of BEDFile objects, compute the distance of each BEDFile
    object-entry to its most proximal feature as defined by -annot. If any
    BEDFile object references a set of BigWig files, subsequently map the
    proximal features start and end onto these bigwig files. Doing so
    facilitates derivation of a mean value that helps infer a signal for the
    proximal annotation entry.
    @param args: dictionary of command-line arguments.
    '''

    beds = [BEDFileFactory(elm).build() for elm in parse_config(args['in'])]
    df = DataFrame()
    for b in beds:  # otherwise, work with all other BED files.
        map_features(b, args['annot'])
        if args['signals']:  # map signals (BigWig files), if need-be
            map_signals(b, args['annot'])
        df = df.append(b.get_data(), ignore_index=True)
    df.to_csv(sys.stdout, index=False)

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
    except IOError as e:
        print(e)
