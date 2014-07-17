'''
Enables parsing of all user-provided input files, namely the configuration
file and all accompanying BED files.
'''

from xml.etree import ElementTree
from pandas import read_table


def as_delim(*args, delim=','):
    '''
    Helpful function that takes a collection and outputs such elements as a
    string. Each element is delimited by a user-provided delimiter.
    @param x: collection populated with objects intent on being delimited.
    @param delim: character delimited.
    '''

    return delim.join([str(i) for i in args])


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
    data = df[[0, 1, 2]]
    data.columns = ['Chr', 'Start', 'End']
    data['Length'] = data['End'] - data['Start']  # get only core data
    if df.shape[1] > 3:  # if there are more than 3x columns, pull them out
        other = df[list(range(3, df.shape[1]))]
    else:
        other = None  # otherwise, return nothing
    return data, other


def parse_vectorized_bed(f):
    '''
    Parses a BED file whereby the last column is exclusively dedicated to
    referencing a vector of values. Such a vector may reference conservation
    scores or gene-expression signals per base of the BED  Thus, due to the
    granularity of this vector, addition computation is required.
    @param f: BED file.
    '''

    df, other = parse_abstract_bed(f)  # vectors are the last column (other)
    vector_data = other[other.columns[-1]].astype('str')  # get last column
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
