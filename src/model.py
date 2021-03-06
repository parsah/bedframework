'''
Represents empirical classes that are used throughout application runtime.
Such classes allow for interaction between data-structures its underlying
BED contents.
'''

import os
from ast import literal_eval  # for translating 'True' to boolean True
from xml.etree.ElementTree import Element
from src import ioutils
from pandas import concat


class BEDFileFactory():
    '''
    Constructs an object of type BEDFile given its respective XML element from
    the user-provided configuration file. Such an element is parsed and
    mapped to the respective state, allowing for construction of a custom
    BEDFile object.
    '''
    def __init__(self, elem):
        assert isinstance(elem, Element)  # element must be an XML object.
        self._element = elem

    def build(self):
        '''
        Construct a BEDFile object given its own respective XML element.
        @return: object of type BEDFile.
        '''
        bf = BEDFile()
        bf.set_fasta(self.element().find('fasta').text)
        bf.set_file(self.element().find('file').text)
        bf.set_class(self.element().find('class').text)
        bf.set_tissue(self.element().find('tissue').text)
        bf.set_bigwigs([i.text for i in self.element().iter('bw')])
        bf.set_is_scalar(literal_eval(self.element().find('is_scalar').text))
        if bf.is_scalar():  # only save actual BED details; nothing else
            bf.set_data(ioutils.parse_abstract_bed(bf.get_file()))
        else:
            bf.set_data(ioutils.parse_vectorized_bed(bf.get_file()))
        bf.get_data()['Tissue'] = bf.get_tissue()  # add information to BED
        bf.get_data()['Class'] = bf.get_class()
        return bf

    def element(self):
        return self._element

    @staticmethod
    def combine(x):
        return concat([i.get_data() for i in x], ignore_index=True)


class BEDFile():
    '''
    Encapsulates various properties of a BED file, features such as a
    corresponding FASTA file, the respective tissue of the BED file, and its
    degree of tissue-specificity. Accompanying BigWig files may also be
    present for the BED file; in-cases whereby assays were performed.
    '''
    def __init__(self):
        self._bedfile = None  # input filename
        self._fasta = None  # corresponding FASTA sequences
        self._data = None  # parsed BED contents
        self._tissue_name = None  # tissue BED file references
        self._tissue_class = None  # magnitude of tissue-specificity
        self._bigwigs = []  # BED graph files useful in expression analysis
        self._is_scalar = False  # whether the BED file is scalar (default)

    def get_file(self):
        return self._bedfile

    def set_file(self, f):
        self._bedfile = f

    def get_fasta(self):
        return self._fasta

    def set_fasta(self, f):
        self._fasta = f

    def get_data(self):
        return self._data

    def set_data(self, x):
        self._data = x

    def get_tissue(self):
        return self._tissue_name

    def set_tissue(self, x):
        self._tissue_name = x

    def get_class(self):
        return self._tissue_class

    def set_class(self, x):
        self._tissue_class = x

    def get_bigwigs(self):
        return self._bigwigs

    def set_bigwigs(self, x):
        self._bigwigs = x

    def set_is_scalar(self, x):
        self._is_scalar = x

    def is_scalar(self):
        return self._is_scalar

    def __repr__(self):
        name = os.path.basename(self.get_file())  # stringify object
        return name + ' ; ' + self.get_tissue() + ' ; ' + self.get_class()
