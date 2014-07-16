'''
Enables parsing of all user-provided input files, namely the configuration
file and all accompanying BED files.
'''

from xml.etree import ElementTree
from pandas import DataFrame, read_table


def parse_config(xml):
    '''
    Parses the user-provided configuration XML file.
    @param f: BED file.
    @return: list of XML objects referencing BEDFile elements.
    '''
    global GENOMIC_BW, GENOMIC_BED, IS_SCALAR
    tree = ElementTree.parse(xml)
    GENOMIC_BW = tree.find('genomebw').text  # genomic bigwig file
    GENOMIC_BED = tree.find('genomebed').text  # genomic features as BED file
    IS_SCALAR = True if tree.find('isscalar').text == 'true' else False
    return list(tree.iter('bed'))


def parse_abstract_bed(f):
    '''
    Parses a BED file only chromosome, start index, and end index are saved.
    Even if the BED file explicitly is not abstract, only these first three
    columns are saved.
    @param f: BED file.
    '''
    df = read_table(f, header=None, sep='\t')  # BED files have no header
    data = DataFrame(data=df.ix[:, 0: 2])  # annotate the data component
    data.columns = ['Chr', 'Start', 'End']
    data['Length'] = data['End'] - data['Start']
    other = df.ix[:, 3: df.shape[1] - 1]  # all other columns, if available
    return (data, other)


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
