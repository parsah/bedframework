'''
Enables parsing of all user-provided input files, namely the configuration
file and all accompanying BED files.
'''

import os
import random
from xml.etree import ElementTree
from pandas import read_table
from subprocess import Popen, PIPE
from Bio import SeqIO


def exec_closest_app(bed, annot):
    cmd = 'bedtools closest -io -d -a ' + bed + ' -b ' + annot + ' -t first'
    out, err = Popen(cmd, shell=True, stdout=PIPE).communicate()
    if err:
        raise IOError('Exception raised running bedtools')
    return out


def exec_average_app(chrom, start, end, bigwig):
    if not os.path.exists(bigwig):  # signaling is sought but file is invalid
        raise IOError(bigwig + ' is an invalid BigWig file.')
    cmd = 'echo ' + chrom + ' ' + str(start) + ' ' + str(end) +\
             ' 1 | bigWigAverageOverBed ' + bigwig + ' stdin stdout'
    out, err = Popen(cmd, shell=True, stdout=PIPE).communicate()
    if err:  # stdout, stderr are outputs; stderr must be None
        raise IOError('Error found following mapping features to BED')
    return out


def build_random_bed(fasta, n, l=400):
    '''
    Create a new BED file that is build solely given randomly-selected
    genomic coordinates from a FASTA file. The corresponding BED filename
    will be bedfile.bed, and will be saved in the current working directory.
    @param fasta: FASTA filename.
    @param n: Number of randomly-selected BED entries to make.
    @param l: Length of each randomly-selected genomic sequence.
    '''

    fname = './bedfile.bed'
    out = open(fname, 'w')
    records = list(SeqIO.parse(open(fasta), 'fasta'))
    for _ in range(n):
        record = random.choice(records)
        chrm = record.description
        start = random.randint(0, len(record.seq))
        out.write(chrm + '\t' + str(start) + '\t' + str(start + l) + '\n')
    return fname


def mkdir(folder):
    try:
        os.makedirs(folder)
    except FileExistsError as e:
        raise OSError(e)


def parse_config(xml):
    '''
    Parses the user-provided configuration XML file.
    @param f: BED file.
    @return: list of XML objects referencing BEDFile elements.
    '''

    tree = ElementTree.parse(xml)
    return list(tree.iter('bed'))


def parse_abstract_bed(f):
    '''
    Parses a BED file only chromosome, start index, and end index are saved.
    Even if the BED file explicitly is not abstract, only these first three
    columns are saved.
    @param f: BED file.
    '''

    df = read_table(f, header=None, sep='\t')  # BED files have no header
    df.insert(3, None, df[2] - df[1])  # insert BED entry length
    colnames = ['Chr', 'Start', 'End', 'Length']
    names = colnames + list(range(df.shape[1] - len(colnames)))
    df.columns = names
    return df


def parse_vectorized_bed(f):
    '''
    Parses a BED file whereby the last column is exclusively dedicated to
    referencing a vector of values. Such a vector may reference conservation
    scores or gene-expression signals per base of the BED  Thus, due to the
    granularity of this vector, addition computation is required.
    @param f: BED file.
    '''

    df = parse_abstract_bed(f)  # vectors are the last column (other)
    vector_data = df[df.columns[-1]].astype('str')  # last column is vector
    all_vectors = []
    for vector in vector_data:
        vector = [float(i) if i != 'n/a' else 0.0 for i in vector.split(',')]
        vector[vector == 'nan'] = 0
        all_vectors.append(vector)
    df['Vectors'] = all_vectors  # add vectors to the actual data-frame
    df['Vector_Length'] = [len(vector) for vector in df['Vectors']]
    df = df[df['Length'] == df['Vector_Length']]  # length must match vector
    df = df[df['Length'] % 100 == 0]  # only save divisible entries
    return df
