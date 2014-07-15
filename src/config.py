'''
System-wide variables that help application runtime.
'''


GENOMIC_BW = None  # references location of an optional BigWig file.
GENOMIC_BED = None  # BED file referencing all genes in this genome.
IS_SCALAR = False  # whether the BED files are all scalar or vectorized.
