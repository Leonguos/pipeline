#!/usr/bin/env python
# -*- coding=utf-8 -*-
"""
    Report for Disease Pipeline
"""
import argparse

__version__ = '4.7'
__author__ = 'suqingdong'


def get_args(config):

    analysis_mapping = '\n'.join(
        '{:-4}  {}'.format(code, config.ANALYSIS_CODE[code][0])
        for code in config.ANALYSIS_CODE)

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description=__doc__,
        epilog='contact: {0} <{0}@novogene.com>'.format(__author__))

    parser.add_argument(
        '-v', '--version', action='version', version=__version__, help="Show the version")

    parser.add_argument(
        '--analy_array',
        metavar='Str',
        help='The anaylsis code list seperated by comma\n{}'.format(analysis_mapping))

    parser.add_argument(
        '--project', help="the file with info, project name.", required=True)

    parser.add_argument('--projdir', help="the path of project", required=True)

    parser.add_argument(
        '--seqsty',
        help="sequencing strategy",
        choices=['WES_ag', 'WES_illu', 'WGS', 'TS'],
        default='WES_ag')

    parser.add_argument(
        '--ref', help="reference", choices=['b37', 'hg19', 'hg38'], default='b37')
    parser.add_argument(
        '--sample', help=" the information of the sample(sample_info)")
    parser.add_argument(
        '--date',
        help="The job name for this analysis")
    parser.add_argument(
        '--TR', help="A bed file for sequence region\n")
    parser.add_argument(
        '--repsty',
        help="primary analysis, qc report or mapping report",
        choices=['primary', 'qc', 'mapping', 'advance'],
        default='primary')
    parser.add_argument('--pre', help="The prefix of qc_list", required=True)
    parser.add_argument(
        '--english',
        help="Use english report template",
        choices=['N', 'Y'],
        default='N')
    parser.add_argument(
        '--WES_xten',
        help="Use WES by xten report template",
        choices=['N', 'Y'],
        default='Y')
    parser.add_argument(
        '--sf',
        help="Do single sample's filter or not",
        choices=['N', 'Y'],
        default='N')
    parser.add_argument(
        '--NGC', help="Show N in GC", choices=['N', 'Y'], default='Y')
    parser.add_argument(
        '--odir', help="the path of output files", required=True)
    parser.add_argument(
        '--mail', help="mail", default='suqingdong@novogene.com')
    parser.add_argument(
        '--HPAmode', help="HPA", choices=['brief', 'complex'], default='brief')
    parser.add_argument(
        '--datastat', help="in sample database?", choices=['N', 'Y'], default='Y')
    parser.add_argument(
        '--Dtype',
        help="the type of contract, Monogenic(M), Complex(C) or Unknown(N).",
        choices=['M', 'C', 'N'],
        default='N')
    parser.add_argument(
        '--network', help="Use phenolyzer, default=NONENULL", default='NONENULL')
    parser.add_argument(
        '--Funcenrichement', help="Funcenrichement report", choices=['N', 'Y'], default='N')

    parser.add_argument(
        '--pdf',
        help='generate pdf report or not[default=%(default)s]',
        choices=['Y', 'N'],
        default='Y')

    parser.add_argument(
        '--disease',
        help='generate disease report or not[default=%(default)s]',
        choices=['Y', 'N'],
        default='N')

    args = vars(parser.parse_args())

    return args
