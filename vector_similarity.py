'''
'''

import argparse
import itertools
from src import parser
from pandas import DataFrame
from scipy.stats import ranksums


def compute_bed_similarity(beds):
    df = DataFrame()
    for bed in beds:
        bed_df = parser.parse_vectorized_bed(bed.get_filename())
        d = {'Length': bed_df[2] - bed_df[1], 'Tissue': bed.get_tissue_name(),
             'Class': bed.get_tissue_class(),
             'Vectors': bed_df[bed_df.shape[1] - 1]}
        df = df.append(DataFrame(d))  # add object to centralized data-frame

    lengths = sorted([int(i) for i in set(df['Length']) if i % 100 == 0])
    tissues = [i for i in set(df['Tissue'])]
    ts = df[df['Class'] == 'Tissue-Specific']
    ub = df[df['Class'] == 'Ubiquitous']
    print('Tissue , Ubiquitous , Tissue.Specific , PValue')
    for tissue in tissues:
        for len_ts in lengths:  # filter the respective data-frames
            for len_ub in lengths:  # filter the respective data-frames
                x = ts[(ts['Length'] == len_ts) & (ts['Tissue'] == tissue)]
                x = list(itertools.chain.from_iterable(list(x['Vectors'])))
                y = ub[(ub['Length'] == len_ub) & (ub['Tissue'] == tissue)]
                y = list(itertools.chain.from_iterable(list(y['Vectors'])))
                pval = ranksums(x, y)[-1]
                print(tissue, ',', len_ts, ',', len_ub, ',', pval)

if __name__ == '__main__':
    try:
        argsparser = argparse.ArgumentParser()
        argsparser.add_argument('-in', metavar='XML', required=True,
                            help='XML configuration file [req]')
        args = vars(argsparser.parse_args())  # parse arguments
        bed_list = parser.parse_config(xml=args['in'])  # parse config
        compute_bed_similarity(beds=bed_list)
    except OSError as e:
        print(e)
    except KeyboardInterrupt:
        print()
