#!/usr/bin/env python2
# -*- coding=utf-8 -*-
import os
import re
import sys
import time
import json
import datetime
import getpass
import commands
import argparse
# import signal
from pprint import pprint

from fuzzywuzzy import process

from subprocess import Popen, PIPE
from commands import getoutput
from commands import getstatusoutput
from ConfigParser import ConfigParser
from collections import OrderedDict

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(BASE_DIR)

from config import config

from QC import QC
from Mapping import Mapping
from Mutation import Mutation, SV, CNV
from Advance import FilterDB, FilterModel, Denovo, IntegrateResult, Other, Association, HLA, ROH, Linkage
from Report import Report
from Result import Result

from arguments import get_args

import utils

#  # test something
# print config.CONFIG.get('common', 'root_dir')
# print dict(config.CONFIG.items('genome_b37'))
# print 'color_text is used to ', utils.color_text('add color for part of sentence')
# utils.print_color('print_color is used to add color for a full sentence', fg='green', bg='aaa', style='bright')
# exit()

__version__ = '4.7'
__author__ = 'suqingdong'


class MainPipeline(object):

    def __init__(self, args):

        self.args = args
        self.analydir = os.path.abspath(args['pwd'])
        self.samp_list = os.path.abspath(args['samp_list'])
        self.samp_info = os.path.abspath(args['samp_info'])
        self.samp_info_done = args['samp_info_done']
        self.pn = os.path.abspath(args['pn'])
        self.seqstrag = args['seqstrag']
        self.analy_array = [
            float(each) if '.' in each else int(each)
            for each in args['analy_array'].split(',')
        ]
        self.startpoint = args['startpoint']
        self.refgenome = args['ref']
        self.newjob = args['newjob']
        self.rawTR = args['TR']
        self.callTR = args['callTR']
        # self.no_check = args['no_check']
        self.queues = args['queues'].split(',')

        self.qc_list = os.path.abspath('qc_list_{}'.format(args['qcsuffix']))

        self.username = getpass.getuser()
        self.email = args['email'] or '{}@novogene.com'.format(self.username)

        self.ROOT_DIR = config.CONFIG.get('common', 'root_dir')

        self.args.update(dict(config.CONFIG.items('genome_' + self.refgenome)))

        self.args.update(self.__dict__)
        # print self.args

        self.jobs = []
        self.orders = []

    def start(self):

        # temp
        advance_dirs = {
            'Merged_vcf': '{analydir}/Advance/{newjob}/Merged_vcf',
            'ACMG': '{analydir}/Advance/{newjob}/ACMG',
            'FilterSV': '{analydir}/Advance/{newjob}/FilterSV',
            'FilterCNV': '{analydir}/Advance/{newjob}/FilterCNV',
            'Noncoding': '{analydir}/Advance/{newjob}/Noncoding',
            'ModelF': '{analydir}/Advance/{newjob}/ModelF',
            'Share': '{analydir}/Advance/{newjob}/Share',
            'Denovo': '{analydir}/Advance/{newjob}/Denovo',
            'Linkage': '{analydir}/Advance/{newjob}/Linkage',
            'ROH': '{analydir}/Advance/{newjob}/ROH',
            'Network': '{analydir}/Advance/{newjob}/Network',
            'Pathway': '{analydir}/Advance/{newjob}/Pathway',
            'PPI': '{analydir}/Advance/{newjob}/PPI',
            'HLA': '{analydir}/Advance/{newjob}/HLA',
            'SiteAS': '{analydir}/Advance/{newjob}/SiteAS',
            'GeneAS': '{analydir}/Advance/{newjob}/GeneAS',
            'IntegrateResult': '{analydir}/Advance/{newjob}/IntegrateResult',
            'Disease': '{analydir}/Advance/{newjob}/Disease',
            'BriefResults': '{analydir}/Advance/{newjob}/BriefResults',
        }

        for k, v in advance_dirs.iteritems():
            self.args.update({
                k: v.format(**self.args)
            })

        # print self.args['SiteAS']
        # exit()

        # print self.analy_array
        print 'hello, {}'.format(self.username)

        # Require rawdata or not
        qc_status = utils.get_status('qc', self.startpoint, config.ANALYSIS_POINTS)
        mapping_status = utils.get_status('bwa_mem', self.startpoint, config.ANALYSIS_POINTS)

        print 'qc status:', qc_status
        print 'mapping status:', mapping_status

        ANALY_DICT = utils.get_analysis_dict(self.analy_array, config.ANALYSIS_CODE)
        self.args.update({'ANALY_DICT': ANALY_DICT})
        # print ANALY_DICT.keys();exit()

        softwares = utils.get_softwares(self.analy_array, self.args['ANALY_DICT'], self.args, self.seqstrag)
        # pprint(softwares);exit()
        self.args.update({'softwares': softwares})


        # check inputs
        self.queues = utils.check_queues(self.queues, self.username)
        self.args.update({'queues': self.queues})

        # use sentieon specific queues if needed
        if 'sentieon' in softwares.values():
            print 'add sentieon_queues'
            sentieon_queues = self.queues
            if config.CONFIG.has_option('resource', 'sentieon_queues'):
                sentieon_queues = config.CONFIG.get(
                    'resource', 'sentieon_queues').split(',')
                sentieon_queues = utils.check_queues(sentieon_queues,
                                                     self.username)
                if not sentieon_queues:
                    sentieon_queues = self.queues
            self.args.update({'sentieon_queues': sentieon_queues})

        # print self.args['sentieon_queues'];exit()
        # print sentieon_queues;exit()

        utils.check_analy_array(self.seqstrag, self.analy_array, config.ANALYSIS_CODE)
        utils.check_files(self.pn, self.samp_info, self.samp_list)
        newTR = utils.check_target_region(config.CONFIG, self.seqstrag, self.refgenome, self.rawTR)
        self.args.update({'TR': newTR})

        print 'analysis items:'
        for analysis_code in self.analy_array:
            print utils.color_text('{:4}  {}'.format(
                analysis_code, config.ANALYSIS_CODE[analysis_code][0]), 'yellow')

        # Analysis start point
        if self.startpoint:
            if self.startpoint in config.ANALYSIS_POINTS:
                print 'start point: {}'.format(
                    utils.color_text(self.startpoint))
            else:
                print '[error] invalid startpoint: {}'.format(
                    utils.color_text(self.startpoint))

                print 'maybe you want to choose: {}'.format(
                    utils.color_text(
                        process.extractOne(self.startpoint,
                                           config.ANALYSIS_POINTS.keys())[0], 'cyan'))

                print 'available startpoints are as follows:\n  {}'.format(
                    '  '.join(config.ANALYSIS_POINTS.keys()))
                exit(1)



        is_advance = max(self.analy_array) > 6.1
        project = utils.Project(self.analydir, self.samp_info, self.samp_info_done, self.samp_list,
                                self.qc_list, qc_status, mapping_status, is_advance)

        # Extract sample_info
        print 'extract sample informations...'

        fenqi, tissue, disease_name, sample_infos, sample_infos_all, sample_done = project.get_sample_infos(
            self.samp_list, self.samp_info, self.samp_info_done,
            is_advance)

        database = '{}/project/DisGeNet.json'.format(config.CONFIG.get('software', 'soft_dir'))
        disease_ids = utils.get_disease_id(disease_name, database)
        self.args.update({
            'disease_name': disease_name,
            'disease_ids': disease_ids,
            })

        sample_infos_waiting = {sampleid: infos for sampleid, infos in sample_infos.iteritems() if sampleid  not in sample_done}
        self.args.update({'sample_infos_waiting': sample_infos_waiting})
        # print sample_infos_waiting
        # exit()

        # print 'fenqi:', fenqi
        # print 'tissue:', tissue
        # exit()

        sample_lists = project.get_sample_lists
        # print sample_lists
        # print sample_infos.keys()
        # print sample_infos_all.keys()
        # for sample in sample_infos:
        #     print sample, sample_infos[sample]['familyid']
        # exit()

        if mapping_status == 'waiting':
            sample_lists = project.update_qc_list()

        print '  report number: {}'.format(utils.color_text(fenqi))
        if disease_name:
            print '  disease name: {}'.format(utils.color_text(disease_name))
            print '  disease id: {}'.format(utils.color_text(disease_ids))
        if tissue:
            print '  tissue: {}'.format(utils.color_text(tissue))
        print '  samples ({}): {}'.format(
            len(sample_infos), utils.color_text(sample_infos.keys()))

        if sample_done:
            print '  samples done({}): {}'.format(
                len(sample_done), utils.color_text(sample_done))

        # Update qc_list and extract sample_list
        # print 'update qc_list...'
        # print json.dumps(sample_lists, indent=2)

        # set memory according seqstrag
        print 'set analysis memory...'
        if self.seqstrag == 'WGS':
            print 'upate memory for WGS...'
            for analysis, memory in config.ANALYSIS_MEM_WGS.items():
                if analysis in config.ANALYSIS_POINTS:
                    config.ANALYSIS_POINTS[analysis][0] = memory
        # exit()

        # ===========================================================
        # ===========================================================
        print '>>> pipeline start...'

        mutation_soft, sv_soft, cnv_soft, denovo_soft = [softwares[each] for each in ('mutation', 'sv', 'cnv', 'denovo')]

        print '  mutation_soft:{}, sv_soft:{}, cnv_soft:{}, denovo_soft:{}'.format(
            mutation_soft, sv_soft, cnv_soft, denovo_soft)

        # QC
        if ANALY_DICT['quality_control'] and qc_status == 'waiting':
            utils.print_color('> QC', 'white')
            QC(self.args, self.jobs, self.orders, sample_lists, config).start()

        # Mapping
        if ANALY_DICT['mapping']:
            utils.print_color('> Mapping', 'white')
            Mapping(self.args, self.jobs, self.orders, sample_lists, sample_infos, config, qc_status, mapping_status).start()

        # Mutation
        if ANALY_DICT['snpindel_call']:
            utils.print_color('> Mutation', 'white')
            Mutation(self.args, self.jobs, self.orders, sample_lists, sample_infos, config).start()

        # SV
        if ANALY_DICT['sv_call']:
            utils.print_color('> SV', 'white')
            SV(self.args, self.jobs, self.orders, sample_infos, config).start()

        # CNV
        if ANALY_DICT['cnv_call']:
            utils.print_color('> CNV', 'white')
            CNV(self.args, self.jobs, self.orders, sample_infos, config).start()

        # FilterDB
        if ANALY_DICT['filter']:
            utils.print_color('> FilterDB', 'white')
            FilterDB(self.args, self.jobs, self.orders, mutation_soft, sv_soft,
                     cnv_soft, sample_infos, config, disease_name, tissue,
                     ANALY_DICT).start()

        # ModelF
        if ANALY_DICT['filter_model']:
            utils.print_color('> Model', 'white')
            FilterModel(self.args, self.jobs, self.orders, mutation_soft,
                        sv_soft, cnv_soft, sample_infos, config).start()

        # Denovo
        if ANALY_DICT['denovo']:
            utils.print_color('> Denovo', 'white')
            Denovo(self.args, self.jobs, self.orders, mutation_soft, sv_soft,
                   cnv_soft, denovo_soft, sample_infos, config,
                   ANALY_DICT).start()

        # Linkage
        if ANALY_DICT['linkage']:
            utils.print_color('> Linkage', 'white')
            Linkage(self.args, self.jobs, self.orders, mutation_soft, sv_soft,
                    cnv_soft, denovo_soft, sample_infos_all, config,
                    ANALY_DICT).start()

        # IntegrateResult
        if any(ANALY_DICT[analysis] for analysis in ['filter', 'filter_model', 'denovo', 'phenolyzer']):
            utils.print_color('> IntegrateResult', 'white')
            IntegrateResult(self.args, self.jobs, self.orders, config).start()

        # ROH
        if ANALY_DICT['roh']:
            utils.print_color('> ROH', 'white')
            ROH(self.args, self.jobs, self.orders, sample_infos, mutation_soft, config).start()

        # OTHER
        other = Other(self.args, self.jobs, self.orders, config, disease_name)

        # IBD
        if any(ANALY_DICT[each] for each in ['filter_model', 'linkage', 'denovo']) and len(sample_infos_waiting) > 1:
            utils.print_color('> IBD', 'white')
            other.ibd()

        # Network
        if ANALY_DICT['phenolyzer']:
            utils.print_color('> Phenolyzer', 'white')
            other.phenolyzer()

        # Pathway
        if ANALY_DICT['pathway']:
            utils.print_color('> Pathway', 'white')
            other.pathway()

        # PPI
        if ANALY_DICT['ppi']:
            utils.print_color('> PPI', 'white')
            other.ppi()

        # SiteAS
        if ANALY_DICT['site_association']:
            utils.print_color('> SiteAS', 'white')
            Association(self.args, self.jobs, self.orders, config).site_association()

        # GeneAS
        if ANALY_DICT['gene_association']:
            utils.print_color('> GeneAS', 'white')
            Association(self.args, self.jobs, self.orders, config).gene_association()

        # HLA
        if ANALY_DICT['hla']:
            utils.print_color('> HLA', 'white')
            HLA(self.args, self.jobs, self.orders, sample_lists, sample_infos, config, qc_status).start()

        # result and report
        utils.print_color('> Result', 'white')
        Result(self.args, self.jobs, self.orders, config).start()

        utils.print_color('> Report', 'white')
        Report(self.args, self.jobs, self.orders, config).start()

        # job summary
        print 'lenght of jobs waiting/total: {}/{}'.format(
            len([job for job in self.jobs if job.get('status') == 'waiting']), len(self.jobs))

        utils.write_job(self.analydir, self.newjob, self.jobs, self.orders)

        print '{:-^80}'.format(' all done ')


def main():

    args = get_args(config, utils, __version__)

    pipeline = MainPipeline(args)

    pipeline.start()


if __name__ == '__main__':

    main()

# the end
