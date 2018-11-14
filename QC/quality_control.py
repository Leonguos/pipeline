#!/usr/bin/env python
# -*- coding=utf-8 -*-
import os

import utils


class QC(object):

    def __init__(self, args, jobs, orders, sample_lists, config):

        self.args = args
        self.jobs = jobs
        self.orders = orders

        self.sample_lists = sample_lists

        self.analydir = args.get('analydir')
        self.queues = args.get('queues')

        self.ANALYSIS_POINTS = config.ANALYSIS_POINTS

    def start(self):

        print '>  qc start...'

        # print self.args.get('analy_array')
        if self.args['ANALY_DICT']['quality_control_no_clean']  and len(self.args.get('analy_array')) > 1:
            print '[error] analysis code 1.1 applies only to the self-built library project!'
            exit(1)
        elif self.args['ANALY_DICT']['quality_control_rm_clean']:
            print '>  clean data will be removed after mapping_check'
        else:
            print '>  clean data will be keeped'

        for sampleID, items in self.sample_lists.items():
            if sampleID not in self.args['sample_infos_waiting']:
                continue
            for lane in items.get('lanes'):
                self.link_raw(sampleID, lane)
                self.md5_raw(sampleID, lane)
                self.qc(sampleID, lane)
                if self.args['ANALY_DICT']['quality_control'] or self.args['ANALY_DICT']['quality_control_keep_clean']:
                    self.gzip_md5_clean(sampleID, lane)

            self.qc_check(sampleID)

    def link_raw(self, sampleID, lane):

        # print '  link rawdata...'
        cmd = '''
            set -eo pipefail
            mkdir -p {analydir}/RawData/{sampleID}
            cd {analydir}/RawData/{sampleID}
            ln -sf {path}/{libID}/{libID}_L{lane}_1.fq.gz {sampleID}_{novoid}_{flowcell}_L{lane}_1.fq.gz
            ln -sf {path}/{libID}/{libID}_L{lane}_1.adapter.list.gz {sampleID}_{novoid}_{flowcell}_L{lane}_1.adapter.list.gz
            ln -sf {path}/{libID}/{libID}_L{lane}_1.adap.stat {sampleID}_{novoid}_{flowcell}_L{lane}_1.adap.stat
            ln -sf {path}/{libID}/{libID}_L{lane}_2.fq.gz {sampleID}_{novoid}_{flowcell}_L{lane}_2.fq.gz
            ln -sf {path}/{libID}/{libID}_L{lane}_2.adapter.list.gz {sampleID}_{novoid}_{flowcell}_L{lane}_2.adapter.list.gz
            ln -sf {path}/{libID}/{libID}_L{lane}_2.adap.stat {sampleID}_{novoid}_{flowcell}_L{lane}_2.adap.stat
        '''.format(
            analydir=self.analydir, sampleID=sampleID, **lane)

        os.system(cmd)

    def md5_raw(self, sampleID, lane):

        # print '  md5sum rawdata...'
        # write shell
        cmd = '''
            set -eo pipefail
            echo MD5 rawdata for sample {sampleID} start: `date "+%F %T"`

            cd {analydir}/RawData/{sampleID}

            md5sum {sampleID}_{novoid}_{flowcell}_L{lane}_1.fq.gz |unix2dos > {sampleID}_{novoid}_{flowcell}_L{lane}_1.fq.gz.MD5.txt
            md5sum {sampleID}_{novoid}_{flowcell}_L{lane}_2.fq.gz |unix2dos > {sampleID}_{novoid}_{flowcell}_L{lane}_2.fq.gz.MD5.txt

            echo MD5 rawdata for sample {sampleID} done: `date "+%F %T"`
        '''.format(
            analydir=self.analydir, sampleID=sampleID, **lane)

        shell_path = '{analydir}/RawData/{sampleID}/md5_raw_{sampleID}_{novoid}_{flowcell}_L{lane}.sh'.format(
            analydir=self.analydir, sampleID=sampleID, **lane)

        utils.write_shell(shell_path, cmd)

        # add job
        now_point = 'md5_raw'
        job_name = 'md5_raw_{sampleID}_{novoid}_{flowcell}_L{lane}'.format(sampleID=sampleID, **lane)
        utils.add_job(self.jobs, now_point, self.args['startpoint'],
                      self.ANALYSIS_POINTS, job_name, shell_path, self.queues)

    def qc(self, sampleID, lane):

        # print '  quality control...'
        # write shell
        if self.args['softwares']['qc'] == 'fastp':
            print '>   qc with fastp ...'
            cmd = '''
                set -eo pipefail

                echo QC for sample {sampleID} start: `date "+%F %T"`
            
                cd {analydir}/QC/{sampleID}

                fastp \\
                    -i ../../RawData/{sampleID}/{sampleID}_{novoid}_{flowcell}_L{lane}_1.fq.gz \\
                    -I ../../RawData/{sampleID}/{sampleID}_{novoid}_{flowcell}_L{lane}_2.fq.gz \\
                    -o {sampleID}_{novoid}_{flowcell}_L{lane}_1.clean.fq.gz \\
                    -O {sampleID}_{novoid}_{flowcell}_L{lane}_2.clean.fq.gz \\
                    -j {sampleID}_{novoid}_{flowcell}_L{lane}.fastp.json \\
                    -h {sampleID}_{novoid}_{flowcell}_L{lane}.fastp.html \\
                    -q 5 -u 50 -n 15 \\
                    -w 4

                python {moduledir}/QC/fastp_convert.py \\
                    --json {sampleID}_{novoid}_{flowcell}_L{lane}.fastp.json \\
                    --identity {sampleID}_{novoid}_{flowcell}_L{lane}

                echo QC for sample {sampleID} done: `date "+%F %T"`
            '''
        else:
            no_clean = '-n ' if self.args['ANALY_DICT']['quality_control_no_clean'] else ''

            cmd = '''
                set -eo pipefail
                echo QC for sample {sampleID} start: `date "+%F %T"`

                cd {analydir}/QC/{sampleID}

                raw2clean_QC \\
                    -i ../../RawData/{sampleID}/{sampleID}_{novoid}_{flowcell}_L{lane}_1.fq.gz,../../RawData/{sampleID}/{sampleID}_{novoid}_{flowcell}_L{lane}_2.fq.gz \\
                    -a ../../RawData/{sampleID}/{sampleID}_{novoid}_{flowcell}_L{lane}_1.adapter.list.gz,../../RawData/{sampleID}/{sampleID}_{novoid}_{flowcell}_L{lane}_2.adapter.list.gz \\
                    {no_clean}-o .

                # Rscript {sampleID}_{novoid}_{flowcell}_L{lane}.QCplot.R

                echo QC for sample {sampleID} done: `date "+%F %T"`
            '''

        cmd = cmd.format(**dict(self.args, **dict(locals(), **lane)))

        shell_path = '{analydir}/QC/{sampleID}/qc_{sampleID}_{novoid}_{flowcell}_L{lane}.sh'.format(
            analydir=self.analydir, sampleID=sampleID, **lane)
    
        # add job
        now_point = 'qc'
        job_name = 'qc_{sampleID}_{novoid}_{flowcell}_L{lane}'.format(sampleID=sampleID, **lane)
        status = utils.add_job(self.jobs, now_point, self.args['startpoint'],
                      self.ANALYSIS_POINTS, job_name, shell_path, self.queues)

        # add order
        after_jobs = ['qc_check_{sampleID}'.format(**locals())]
        utils.add_order(self.orders, job_name, after_jobs=after_jobs)

        if status == 'waiting':

            utils.write_shell(shell_path, cmd)

    def gzip_md5_clean(self, sampleID, lane):

        # print '  gzip and md5sum clean data...'
        # write shell
        cmd = '''
            set -eo pipefail
            echo Compress and md5sum clean for sample {sampleID} start: `date "+%F %T"`\n

            cd {analydir}/QC/{sampleID}

            fq1={sampleID}_{novoid}_{flowcell}_L{lane}_1.clean.fq
            fq2={sampleID}_{novoid}_{flowcell}_L{lane}_2.clean.fq

            for fq in $fq1 $fq2;do
                if [ -s $fq.gz -a ! -s $fq ];then
                    echo $fq has been compressed.
                else
                    pigz -p4 -f $fq
                    md5sum $fq.gz | unix2dos > $fq.gz.MD5.txt
                fi
            done

            echo Compress and md5sum clean for sample {sampleID} done: `date "+%F %T"`
        '''.format(
            analydir=self.analydir, sampleID=sampleID, **lane)

        shell_path = '{analydir}/QC/{sampleID}/gzip_md5_clean_{sampleID}_{novoid}_{flowcell}_L{lane}.sh'.format(
            analydir=self.analydir, sampleID=sampleID, **lane)

        utils.write_shell(shell_path, cmd)

        # add job
        now_point = 'gzip_md5_clean'
        job_name = 'gzip_md5_clean_{sampleID}_{novoid}_{flowcell}_L{lane}'.format(sampleID=sampleID, **lane)

        utils.add_job(self.jobs, now_point, self.args['startpoint'],
                      self.ANALYSIS_POINTS, job_name, shell_path, self.queues)

        # add order
        before_jobs = ['qc_{sampleID}_{novoid}_{flowcell}_L{lane}'.format(sampleID=sampleID, **lane)]
        after_jobs = ['data_release']
        utils.add_order(self.orders, job_name, before_jobs=before_jobs, after_jobs=after_jobs)

    def qc_check(self, sampleID):

        # print '  qc check...'
        # write shell
        cmd = '''
            set -eo pipefail
            echo qc check for sample {sampleID} start: `date "+%F %T"`

            python2 {moduledir}/QC/auto_check.py \\
                --qc_list {qc_list} \\
                --sampid {sampleID} \\
                --pwd {analydir} \\
                --check qc \\
                --jobname {newjob} \\
                --seqstrag {seqstrag} \\
                --email {email} \\
                --PE {PE} \\
                --q20 {Q20} \\
                --q30 {Q30} \\
                --error {error} \\
                --raw {rawdata}
            
            echo qc check for sample {sampleID} done: `date "+%F %T"`
        '''.format(
            sampleID=sampleID,
            **self.args)

        shell_path = '{analydir}/QC/{sampleID}/qc_check_{sampleID}.sh'.format(
            analydir=self.analydir, sampleID=sampleID)

        utils.write_shell(shell_path, cmd)

        # add job
        now_point = 'qc_check'
        job_name = 'qc_check_{sampleID}'.format(sampleID=sampleID)
        utils.add_job(self.jobs, now_point, self.args['startpoint'],
                      self.ANALYSIS_POINTS, job_name, shell_path, self.queues)

        # add order
        after_jobs = ['qc_report']
        utils.add_order(self.orders, job_name, after_jobs=after_jobs)

# the end
