'''
Enables parsing of all user-provided input files, namely the configuration
file and all accompanying BED files.
'''


def parse_config(xml):
    '''
    Parses the user-provided configuration XML file.
    @param f: BED file.
    '''


def parse_abstract_bed(f):
    '''
    Parses a BED file only chromosome, start index, and end index are saved.
    Even if the BED file explicitly is not abstract, only these first three
    columns are saved.
    @param f: BED file.
    '''


def parse_vectorized_bed(f):
    '''
    Parses a BED file whereby the last column is exclusively dedicated to
    referencing a vector of values. Such a vector may reference conservation
    scores or gene-expression signals per base of the BED entry.
    @param f: BED file.
    '''
    pass
