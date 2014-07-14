'''
Enables parsing of all user-provided input files, namely the configuration
file and all accompanying BED files.
'''

import numpy
from xml.etree import ElementTree
from src.model import BEDFile
from pandas import read_table


def parse_config(xml):
    '''
    Parses the user-provided configuration XML file.
    @param f: BED file.
    @return: list of objects of type BEDfile.
    '''
    global GENOMIC_BW
    tree = ElementTree.parse(xml)
    GENOMIC_BW = tree.find('genomebw').text  # genomic bigwig file
    bed_files = []  # store all BED files
    for elem in tree.iter('bed'):
        bf = BEDFile()
        bf.set_fasta(elem.find('fasta').text)
        bf.set_filename(elem.find('file').text)
        bf.set_tissue_class(elem.find('class').text)
        bf.set_tissue_name(elem.find('tissue').text)
        bf.set_bigwigs([i.text for i in elem.iter('bw')])  # set bigwigs
        bed_files.append(bf)
    return bed_files


def parse_abstract_bed(f):
    '''
    Parses a BED file only chromosome, start index, and end index are saved.
    Even if the BED file explicitly is not abstract, only these first three
    columns are saved.
    @param f: BED file.
    '''
    df = read_table(f, header=None, sep='\t')  # BED files have no header
    return df


def parse_vectorized_bed(f):
    '''
    Parses a BED file whereby the last column is exclusively dedicated to
    referencing a vector of values. Such a vector may reference conservation
    scores or gene-expression signals per base of the BED  Thus, due to the
    granularity of this vector, addition computation is required.
    @param f: BED file.
    '''
    df = parse_abstract_bed(f)  # vectorized BEDs extend abstract BEDs
    vector_data = df[df.shape[1] - 1].astype('str')  # get last column
    all_vectors = []
    for vector in vector_data:
        #bed_len = int(df.iloc[rownum][2] - df.iloc[rownum][1])
        vector = [float(i) if i != 'n/a' else 0 for i in vector.split(',')]
        vector[vector == 'nan'] = 0
        all_vectors.append(vector)
    df = df.drop(df.shape[1] - 1, axis=1)
    df[df.shape[1]] = all_vectors  # add vectors as final column
    return df
