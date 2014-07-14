'''
Enables parsing of all user-provided input files, namely the configuration
file and all accompanying BED files.
'''

from xml.etree import ElementTree
from src.model import BEDFile
from pandas import DataFrame
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


def parse_vectorized_bed(f):
    '''
    Parses a BED file whereby the last column is exclusively dedicated to
    referencing a vector of values. Such a vector may reference conservation
    scores or gene-expression signals per base of the BED entry.
    @param f: BED file.
    '''
    pass
