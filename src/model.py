'''
Represents empirical classes that are used throughout application runtime.
Such classes allow for interaction between data-structures its underlying
BED contents.
'''


class BEDFile():
    def __init__(self):
        self._bedfile = None  # input filename
        self._fasta = None  # corresponding FASTA sequences
        self._data = None  # parsed BED contents
        self._tissue_name = None  # tissue BED file references
        self._tissue_class = None  # magnitude of tissue-specificity
        self._bigwigs = []  # BED graph files useful in expression analysis
        self._vector = []  # vector if BED file is vectorized (optional)

    def get_filename(self):
        return self._bedfile

    def set_filename(self, f):
        self._bedfile = f

    def get_fasta(self):
        return self._fasta

    def set_fasta(self, f):
        self._fasta = f

    def get_data(self):
        return self._data

    def set_data(self, x):
        self._data = x

    def get_tissue_name(self):
        return self._tissue_name

    def set_tissue_name(self, x):
        self._tissue_name = x

    def get_tissue_class(self):
        return self._tissue_class

    def set_tissue_class(self, x):
        self._tissue_class = x

    def get_bigwigs(self):
        return self._bigwigs

    def set_bigwigs(self, x):
        self._bigwigs = x

    def get_vector(self):
        return self._vector

    def set_vector(self, x):
        self._vector = x
