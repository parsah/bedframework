'''
Represents empirical classes that are used throughout application runtime.
Such classes allow for interaction between data-structures its underlying
BED contents.
'''


class BEDFile():
    def __init__(self):
        self._f = None  # input filename
        self._data = None  # parsed BED contents
        self._tissue_name = None  # tissue BED file references
        self._tissue_class = None  # magnitude of tissue-specificity
        self._vector = []  # BED vector; present if BED contains base-values
        self._bedgraphs = []  # BED graph files useful in expression analysis

    def get_filename(self):
        return self._f

    def set_filename(self, f):
        self._f = f

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

    def get_vector(self):
        return self._vector

    def set_vector(self, x):
        self._vector = x

    def get_bedgraphs(self):
        return self._bedgraphs

    def set_bedgraphs(self, x):
        self._bedgraphs = x

    def is_vectorized(self):
        if len(self.get_vector()):
            return True
        return False
