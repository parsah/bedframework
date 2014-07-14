'''
Given a BED file and a corresponding set of all genes, this module helps
map each BED entry to its most proximal gene. The bulk of such logic is made
possible using the BEDtools package.
'''

import argparse
from src import parser


def closest():
    pass

if __name__ == '__main__':
    argsparser = argparse.ArgumentParser()
    argsparser.add_argument('-in', metavar='XML', required=True,
                        help='XML configuration file [req]')
