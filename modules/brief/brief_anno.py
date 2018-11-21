#!/usr/bin/env python
# -*- coding=utf-8 -*-
"""
add header and extract columns
"""
import re
import os
import sys
import argparse

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(os.path.dirname(BASE_DIR)))

from utils import mkdir_if_not_exists, safe_open


HEADER_DIR = os.path.join(BASE_DIR, 'header')

header_map = {
    'snpindel': [
        'Priority:FILTER', 'GeneName:Gene', 'GeneDetail:Gencode+1', 'cytoBand',
        'genomicSuperDups:avsnp', 'FORMAT:Ori_ALT+1', '1000g_Chinese:LRT',
        'CADD:MCAP+1', 'ClinVar', 'OMIM', 'HGMD_Disease_ID:'
    ],
    'acmg': ['ACMG', 'Heredity', 'Priority:avsnp', 'ClinVar:Heredity'],
    'sv_cnv': ['Priority:cpgIslandExt', 'genomicSuperDups', 'Repeat', 'encodeK562+1:StringentLib', 'dgvMerged', 'gwasCatalog', 'encodeGm12878:encodeK562+1', 'StringentLib:'],
    'noncoding': [
        'CHROM:FILTER', 'GeneName:cpgIslandExt', 'cytoBand', 'genomicSuperDups:avsnp',
        'FORMAT:Ori_ALT+1', '1000g_Chinese:NovoDb_WGS+1',
        'CADD:gerp++gt2+1', 'ClinVar', 'OMIM', 'HGMD_Disease_ID:'
    ],
}


def get_indices(headerlist, columns):

    # print headerlist, columns
    for column in columns:
        if ':' in column:
            start, stop = column.split(':')

            start_n = re.findall(r'^(.+?)\+(\d+)$', start)
            stop_n = re.findall(r'^(.+?)\+(\d+)$', stop)

            if start_n:
                start_idx = headerlist.index(start_n[0][0]) + int(start_n[0][1])
            else:
                start_idx = headerlist.index(start) if start else None

            if stop_n:
                stop_idx = headerlist.index(stop_n[0][0]) + int(stop_n[0][1])
            else:
                stop_idx = headerlist.index(stop) if stop else None

            yield slice(start_idx, stop_idx)
        else:
            yield headerlist.index(column)

def get_new_line(linelist, indices):

    new_linelist = []
    for indice in indices:
        if isinstance(indice, slice):
            new_linelist += linelist[indice]
        else:
            new_linelist.append(linelist[indice])

    return '\t'.join(new_linelist) + '\n'


def get_added_header(header_file):

    with safe_open(header_file) as h:
        for line in h:
            yield line


def main():

    infile = args['infile']
    filetype = args['type']
    add_header = args['add_header']

    outdir = args['outdir']

    # print args;exit()

    if len(infile) == 1:
        infile_list = re.split(r'\s+|;|:|,', infile[0])
    else:
        infile_list = infile

    for infile in infile_list:

        outfile = infile.replace('.xls', '.brief.xls')
        if outdir:
            outfile = os.path.join(outdir, os.path.basename(outfile))

        with safe_open(infile) as f, safe_open(outfile, 'w') as out:

            if add_header:
                header_file = os.path.join(BASE_DIR, 'header/{}.header'.format(filetype))
                print 'add header:', header_file
                for line in get_added_header(header_file):
                    # print line
                    out.write(line)

            for line in f:
                linelist = line.rstrip('\n').split('\t')
                # if all(h in linelist for h in ['CHROM', 'POS']):
                if linelist[0] in ('Priority', 'Chr', 'CHROM'):
                    indices = list(get_indices(linelist, header_map[filetype]))
                    new_header = get_new_line(linelist, indices)
                    # print new_header
                    out.write(new_header)
                    continue
                new_line = get_new_line(linelist, indices)
                # print new_line
                out.write(new_line)

        print 'write brief file:', outfile


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument(
        '-i', '--infile', help='The input file', required=True, nargs='*')
    parser.add_argument('-o', '--outfile', help='The output file', nargs='*')
    parser.add_argument('-O', '--outdir', help='The output directory')

    parser.add_argument(
        '-t',
        '--type',
        help='The type of infile',
        choices=header_map,
        required=True)

    parser.add_argument('-header', '--add-header', help='Add readme header', action='store_true')

    args = vars(parser.parse_args())

    main()
