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

    genes, distance = [], []
    fname = b.get_file()
    cmd = 'bedtools closest -io -d -a ' + fname + ' -b ' + gtf + ' -t first'
    out, err = Popen(cmd, shell=True, stdout=PIPE).communicate()
    if err:  # stdout, stderr are outputs; stderr must be None
        raise IOError('Error found following mapping features to BED')
    for l in out.splitlines():
        l = l.decode('utf-8').split('\t')
        gene = l[-2].split(';')[0].replace('gene_id', '').replace('"', "")
        genes.append(gene.strip())  # add the entrie's most proximal gene
        distance.append(int(l[-1]))  # add most proximal distance
    b.get_data()['Gene_ID'] = genes  # add new columns to the data-structure
    b.get_data()['Distance'] = distance
    ints = list(range(0, 200000, 10000))  # fit distances to intervals
    b.get_data()['Interval'] = [bisect_left(ints, pos) for pos in distance]


def map_signals(b):
    '''
    '''
    assert 'Gene_ID' in b.get_data().columns  # a Gene_ID columns is required
    signal = []
    for _, row in b.get_data().iterrows():  # row-number, row, respectively
        reps = []  # store bigwig replicate-averages across bigwig files
        for bigwig in b.get_bigwigs():  # compute BED average given bigwig
            cmd = 'echo ' + str(row['Chr']) + ' ' + str(row['Start']) + ' ' +\
                str(row['End']) + ' 1 | bigWigAverageOverBed ' + bigwig +\
                 ' stdin stdout'
            out, err = Popen(cmd, shell=True, stdout=PIPE).communicate()
            if err:  # stdout, stderr are outputs; stderr must be None
                raise IOError('Error found following mapping features to BED')
            reps.append(float(out.decode('utf-8').strip().split('\t')[-1]))
        signal.append(np.mean(reps))  # append the mean
    b.get_data()['Signal'] = signal


def main(args):
    beds = [BEDFileFactory(elm).build() for elm in parse_config(args['in'])]
    bed_objs = [b for b in beds]  # merge BEDs into one structure
    print('Chr,Length,Tissue,Class,GeneID,Distance,Interval,Signal')  # header
    for b in bed_objs:
        map_features(b, args['gtf'])
        if args['signals']:  # map signals (BigWig files), if need-be
            map_signals(b)
        else:
            b.get_data()['Signal'] = 'None'  # add invalid values if no signal
        data = b.get_data()
        for _, row in data.iterrows():
            if row['Distance'] >= 0:  # BED entry lacks proximal, i.e. chrX
                print(row['Chr'] + ',' + str(row['Length']) + ',' +
                      str(row['Tissue']) + ',' + str(row['Class']) + ',' +
                      str(row['Gene_ID']) + ',' + str(row['Distance']) + ',' +
                      str(row['Interval']) + ',' + str(row['Signal']))

if __name__ == '__main__':
    try:
        parser = argparse.ArgumentParser()
        parser.add_argument('-in', metavar='XML', required=True,
                            help='XML configuration file [req]')
        parser.add_argument('-gtf', metavar='FILE', required=True,
                            help='GTF with gene_id as first attribute [req]')
        parser.add_argument('--signals', action='store_true', default=False,
                            help='Map BED bigwigs signals to features [false]')
        args = vars(parser.parse_args())
        main(args)
    except KeyboardInterrupt:
        print()
