#!/usr/bin/env python
# -*- coding=utf-8 -*-
import os
import sys
import glob
import datetime
import gzip
import json
import commands
import textwrap
import chardet
from collections import defaultdict

import django
from django.template import Context, loader
from django.conf import settings

from arguments import get_args_tar

RESULT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(os.path.dirname(RESULT_DIR))
sys.path.append(ROOT_DIR)

from config import config
import utils


class TarRelease(object):

    def __init__(self, **args):
        self.args = args
        self.projectid = self.parse_pn()[0]
        self.ymd = utils.get_now_time('%Y%m%d')
        self.__dict__.update(**args)

    def start(self):

        print '> start tar release for {projectid} / {newjob} at {ymd}'.format(**self.__dict__)


    def parse_pn(self):

        with utils.safe_open(self.args['pn']) as f:
            projectid, projname = f.read().strip().split(None, 1)

        return projectid, projname

    @staticmethod
    @property
    def dir_map():

        dir_map = {
            "RawData": "ReleaseResult/Data/RawData",
            "QC": "ReleaseResult/Data/CleanData",
            "Mapping": "ReleaseResult/Data/BamData",
            "Variation": "ReleaseResult/PrimaryAnalysis",
            "Advance": "ReleaseResult/FinalResult",
            "Readme": "ReleaseResult/Readme"
        }


def main():

    args = get_args_tar(config)

    tr = TarRelease(**args)

    tr.start()


if __name__ == "__main__":

    main()
