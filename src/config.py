'''
System-wide variables that help application runtime.
'''

import os
import argparse
from xml.etree import ElementTree
from xml.dom import minidom

TISSUE_SPEC = 'Tissue-Specific'
UBIQUITOUS = 'Ubiquitous'


def to_element(parent, tag, text=None, attrib={}):
    elem = ElementTree.SubElement(parent, tag, attrib)
    if not text:
        elem.text = ''
    else:
        elem.text = text
    elem.attrib = attrib
    return elem


def generate(args):
    '''
    Generates a skeleton configuration XML output file.
    @param args: dictionary of command-line arguments.
    '''
    beds = [i for i in os.listdir(args['dir']) if i.endswith('bed')]
    assert len(beds) > 0  # at least 1x BED file must be present
    root = ElementTree.Element('configuration')
    bedfiles = ElementTree.SubElement(root, 'bedfiles')  # add children node
    for fname in beds:  # create elements per BED file
        bed_elem = to_element(parent=bedfiles, tag='bed')
        to_element(parent=bed_elem, tag='file', text=args['dir'] + '/' + fname)
        to_element(parent=bed_elem, tag='class', text=' ')
        to_element(parent=bed_elem, tag='tissue', text=' ')
        to_element(parent=bed_elem, tag='fasta', text=' ')
        bwfiles = to_element(parent=bed_elem, tag='bigwigfiles')
        to_element(parent=bwfiles, tag='bw', text=' ')
        to_element(parent=bwfiles, tag='bw', text=' ')
        to_element(parent=bed_elem, tag='is_scalar', text=str(args['scalar']))

    # parse XML object and prettify string, output to standard-output
    parsed = minidom.parseString(ElementTree.tostring(root))
    print(parsed.toprettyxml())


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-dir', metavar='DIR', required=True,
                        help='Folder containing BED files [req]')
    parser.add_argument('--scalar', action='store_true',
                        help='BED files are scalar [False]')
    args = vars(parser.parse_args())
    generate(args)
