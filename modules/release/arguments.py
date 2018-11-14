#!/usr/bin/env python
# -*- coding=utf-8 -*-
"""
    Data Release for Disease Pipeline
"""
import argparse

__version__ = '4.7'
__author__ = 'suqingdong'


def get_args(config):
    """
    argument parser for data_release
    """

    analysis_mapping = '\n'.join(
        '{:-4}  {}'.format(code, config.ANALYSIS_CODE[code][0])
        for code in config.ANALYSIS_CODE)
    # print analysis_mapping

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description=__doc__,
        epilog='contact: {}@novogene.com'.format(__author__))

    parser.add_argument(
        '-v', '--version', action='version', version=__version__, help="Show the version")

    parser.add_argument(
        '--analy_array',
        metavar='Str',
        help='the anaylsis code list seperated by comma\n{}'.format(analysis_mapping))

    parser.add_argument(
        '--samp_info', help="the sample_info file", required=True)

    parser.add_argument(
        '--samp_info_done', help="the sample_info_done file")

    parser.add_argument(
        '--qc_list', help="the qc_list file", required=True)

    parser.add_argument(
        '--pn', help="the file with info,project name.", required=True)

    parser.add_argument(
        '--analydir', help="the path of project", required=True)

    parser.add_argument(
        '--ER',
        help="generate english report",
        choices=['N', 'Y'],
        default='N')

    parser.add_argument(
        '--MT',
        help="Generate MT analysis result.",
        choices=['N', 'Y'],
        default='N')

    parser.add_argument(
        '--sf',
        help="Standard filter for every sample",
        choices=['N', 'Y'],
        default='N')

    parser.add_argument(
        '--newjob',
        help="the job name",
        required=True)
        
    parser.add_argument(
        '--odir', help="data release directory")

    args = vars(parser.parse_args())

    return args

def get_args_tar(config):
    """
    argument parser for tar_release
    """

    analysis_mapping = '\n'.join(
        '{:-4}  {}'.format(code, config.ANALYSIS_CODE[code][0])
        for code in config.ANALYSIS_CODE)
    # print analysis_mapping

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description=__doc__,
        epilog='contact: {}@novogene.com'.format(__author__))

    parser.add_argument(
        '-v', '--version', action='version', version=__version__, help="Show the version")

    parser.add_argument(
        '--analydir', help="the path of project", required=True)

    parser.add_argument(
        '--analy_array',
        help='the anaylsis code list seperated by comma\n{}'.format(analysis_mapping))

    parser.add_argument(
        '--newjob',
        help="the job name",
        required=True)

    parser.add_argument(
        '--fenqi', help="the fenqi number", required=True)

    parser.add_argument(
        '--pn', help="the pn.txt file", required=True)

    parser.add_argument(
        '--odir', help="the output directory")

    parser.add_argument(
        '--mail', help="the email of the xinxi")

    parser.add_argument(
        '--yymail', help="the email of the yunying")

    parser.add_argument(
        '--release',
        help="specifies the data to be released, default will judge from analy_array")

    parser.add_argument(
        '--fore-tar',
        action='store_true',
        help="ignore the data size and force packaging")

    parser.add_argument(
        '--no-record',
        action='store_true',
        help="do not store the release informations to database")

    args = vars(parser.parse_args())

    return args
