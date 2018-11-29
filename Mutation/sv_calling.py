#!/usr/bin/env python
# -*- coding=utf-8 -*-
import os
import utils
import string


class SV(object):

    def __init__(self, args, jobs, orders, sample_infos, config):

        self.args = args
        self.jobs = jobs
        self.orders = orders

        self.sv_soft = args['softwares']['sv']
        self.analydir = args.get('analydir')
        self.queues = args.get('queues')

        self.sample_infos = sample_infos
        self.sample_infos_waiting = args['sample_infos_waiting']
        self.ANALYSIS_POINTS = config.ANALYSIS_POINTS

        self.config = config

        self.soft_dir = config.CONFIG.get('software', 'soft_dir')

        self.ref = args.get('ref')
        self.annovar = self.config.CONFIG.get('genome_' + self.ref, 'annovar')

        self.__dict__.update(dict(config.CONFIG.items('genome_' + self.ref)))
        self.__dict__.update(args)

    def start(self):

        print '>  sv with {} ...'.format(self.sv_soft)
        for sampleID, infos in self.sample_infos_waiting.iteritems():
            sex = infos['sex']
            chrom_list = utils.get_chrom_list(self.ref, sex, self.args['MT'])
            # print sampleID, chrom_list
            if self.sv_soft == 'crest':
                self.crest_call(sampleID, chrom_list)
                self.crest_txt2gff(sampleID, chrom_list)
            elif self.sv_soft == 'breakdancer':
                self.breakdancer_config(sampleID)
                self.breakdancer_call(sampleID)
            elif self.sv_soft == 'lumpy':
                self.lumpy_call(sampleID)
            self.annotate_gff(sampleID)

    def crest_call(self, sampleID, chrom_list):

        for chrom in chrom_list:
            # 1 extract soft clip
            cmd = '''
                set -eo pipefail
                echo sv call with crest for {sampleID} start: `date "+%F %T"`\n

                # 1 Extract Softclip
                echo extract softclip...
                perl {soft_dir}/CREST/CREST/extractSClip.pl \\
                    -o {analydir}/SV/{sampleID}/crest/bychr \\
                    -i {analydir}/Mapping/{sampleID}.{sampleID}/{sampleID}.final.bam \\
                    -ref_genome {reffasta} \\
                    -r {chrom} \\
                    -p {sampleID}

                # 2 CREST Calling
                echo crest calling...
                perl {soft_dir}/CREST/crest_sv_calling_pipe4.6.pl \\
                    -outDir {analydir}/SV/{sampleID}/crest/bychr \\
                    -tumorBam {analydir}/Mapping/{sampleID}.{sampleID}/{sampleID}.final.bam \\
                    --min_one_side_reads 4 \\
                    -sampleID {sampleID}.{chrom} \\
                    -regionList {refbed} \\
                    -cover {analydir}/SV/{sampleID}/crest/bychr/{sampleID}.{chrom}.cover \\
                    -ref {reffasta} \\
                    -bit {reffasta2bit}

                echo sv call with crest for {sampleID} done: `date "+%F %T"`
            '''.format(
                sampleID=sampleID, chrom=chrom, **self.__dict__)

            shell_path = '{analydir}/SV/{sampleID}/crest/crest_call_chr_{chrom}_{sampleID}.sh'.format(
                sampleID=sampleID, chrom=chrom.strip('chr'), **self.__dict__)

            utils.write_shell(shell_path, cmd)

            # add job
            now_point = 'crest_call'
            job_name = 'crest_call_{}_{}'.format(chrom, sampleID)
            utils.add_job(self.jobs, now_point, self.args['startpoint'],
                          self.ANALYSIS_POINTS, job_name, shell_path,
                          self.queues)

            # add order
            before_jobs = ['final_bam_{sampleID}'.format(sampleID=sampleID)]
            after_jobs = ['crest_txt2gff_{sampleID}'.format(sampleID=sampleID)]
            utils.add_order(self.orders, job_name, before_jobs=before_jobs, after_jobs=after_jobs)

    def crest_txt2gff(self, sampleID, chrom_list):

        auto_chrom_last = filter(lambda x: x.strip('chr').isdigit(), chrom_list)[-1]

        other_chrom = filter(lambda x: not x.strip('chr').isdigit(), chrom_list)

        # print auto_chrom_last
        # print other_chrom

        crest_results = '{{1..%s},%s}' % (
            auto_chrom_last.strip('chr'),
            ','.join(map(lambda x: x.strip('chr'), other_chrom)))

        if 'chr' in chrom_list[-1]:
            crest_results = 'chr' + crest_results

        # print crest_results

        cmd = '''
            set -eo pipefail
            echo convert sv results to gff for {sampleID} start: `date "+%F %T"`\n

            cd {analydir}/SV/{sampleID}/crest

            cat bychr/{sampleID}.{crest_results}.predSV.txt |
                grep -vw hs37d5 > {sampleID}.predSV.txt

            perl {moduledir}/Varition/SV/crest/crest_txt2gff.pl \\
                {sampleID}.predSV.txt > {sampleID}.crest.gff

            echo convert sv results to gff for {sampleID} done: `date "+%F %T"`
        '''.format(
            crest_results=crest_results, sampleID=sampleID, **self.__dict__)

        shell_path = '{analydir}/SV/{sampleID}/crest/crest_txt2gff_{sampleID}.sh'.format(
            sampleID=sampleID, **self.__dict__)

        utils.write_shell(shell_path, cmd)

        # add job
        now_point = 'crest_txt2gff'
        job_name = 'crest_txt2gff_{}'.format(sampleID)
        utils.add_job(self.jobs, now_point, self.args['startpoint'],
                      self.ANALYSIS_POINTS, job_name, shell_path, self.queues)

        # add order
        after_jobs = ['annotate_gff_{sampleID}'.format(sampleID=sampleID)]
        utils.add_order(self.orders, job_name, after_jobs=after_jobs)

    def breakdancer_config(self, sampleID):

        cmd = '''
            set -eo pipefail
            echo breakdancer config for {sampleID} start: `date "+%F %T"`\n

            cd {analydir}/SV/{sampleID}/breakdancer

            perl {soft_dir}/breakdancer/current/bam2cfg.pl \\
                -g -h -n 100000 \\
                {analydir}/Mapping/{sampleID}.{sampleID}/{sampleID}.final.bam \\
                > {sampleID}.breakdancer.cfg

            echo breakdancer config for {sampleID} done: `date "+%F %T"`
        '''.format(
            sampleID=sampleID, **self.__dict__)

        shell_path = '{analydir}/SV/{sampleID}/breakdancer/breakdancer_config_{sampleID}.sh'.format(
            sampleID=sampleID, **self.__dict__)

        utils.write_shell(shell_path, cmd)

        # add job
        now_point = 'breakdancer_config'
        job_name = 'breakdancer_config_{}'.format(sampleID)
        utils.add_job(self.jobs, now_point, self.args['startpoint'],
                      self.ANALYSIS_POINTS, job_name, shell_path, self.queues)

        # add order
        before_jobs = ['final_bam_{sampleID}'.format(sampleID=sampleID)]
        utils.add_order(self.orders, job_name, before_jobs=before_jobs)

    def breakdancer_call(self, sampleID):

        cmd = '''
            set -eo pipefail
            echo breakdancer call for {sampleID} start: `date "+%F %T"`\n

            cd {analydir}/SV/{sampleID}/breakdancer

            {soft_dir}/breakdancer/current/breakdancer-max \\
                -h \\
                -d {sampleID}.breakdancer.SV-supporting \\
                -g {sampleID}.breakdancer.bed \\
                {sampleID}.breakdancer.cfg |
                grep -vwE 'hs37d5|GL000220' \\
                > {sampleID}.breakdancer.txt

            perl {moduledir}/Varition/SV/breakdancer/breakdancer_filter.pl \\
                -g {sex} \\
                -n 6 \\
                -a {sampleID}.breakdancer.txt \\
                > {sampleID}.breakdancer.flt.txt

            perl {moduledir}/Varition/SV/breakdancer/breakdancer_txt2gff.pl \\
                {sampleID}.breakdancer.flt.txt \\
                > {sampleID}.breakdancer.gff

            echo breakdancer call for {sampleID} done: `date "+%F %T"`
        '''.format(
            sampleID=sampleID,
            sex=self.sample_infos[sampleID]['sex'],
            **self.__dict__)

        shell_path = '{analydir}/SV/{sampleID}/breakdancer/breakdancer_call_{sampleID}.sh'.format(
            sampleID=sampleID, **self.__dict__)

        utils.write_shell(shell_path, cmd)

        # add job
        now_point = 'breakdancer_call'
        job_name = 'breakdancer_call_{}'.format(sampleID)
        utils.add_job(self.jobs, now_point, self.args['startpoint'],
                      self.ANALYSIS_POINTS, job_name, shell_path, self.queues)

        # add order
        before_jobs = ['breakdancer_config_{sampleID}'.format(sampleID=sampleID)]
        after_jobs = ['annotate_gff_{sampleID}'.format(sampleID=sampleID)]
        utils.add_order(self.orders, job_name, before_jobs=before_jobs, after_jobs=after_jobs)

    def lumpy_call(self, sampleID):

        sex = self.sample_infos[sampleID]['sex']

        cmd = '''
            set -eo pipefail

            echo lumpy call for {sampleID} start: `date "+%F %T"`\n

            cd {analydir}/SV/{sampleID}/lumpy

            python {moduledir}/Varition/SV/lumpy/lumpy.py \\
                -b {analydir}/Mapping/{sampleID}.{sampleID}/{sampleID}.final.bam \\
                -r {ref} \\
                -o {sampleID}

            rm -f *bam*

            echo lumpy call for {sampleID} done: `date "+%F %T"`
        '''.format(
            **dict(self.__dict__, **locals()))

        shell_path = '{analydir}/SV/{sampleID}/lumpy/lumpy_call_{sampleID}.sh'.format(
            **dict(self.__dict__, **locals()))

        utils.write_shell(shell_path, cmd)

        # add job
        now_point = 'lumpy_call'
        job_name = 'lumpy_call_{}'.format(sampleID)
        utils.add_job(self.jobs, now_point, self.args['startpoint'],
                      self.ANALYSIS_POINTS, job_name, shell_path, self.queues)

        # add order
        before_jobs = ['final_bam_{sampleID}'.format(sampleID=sampleID)]
        after_jobs = ['annotate_gff_{sampleID}'.format(sampleID=sampleID)]
        utils.add_order(self.orders, job_name, before_jobs=before_jobs, after_jobs=after_jobs)

    def annotate_gff(self, sampleID):

        cmd = '''
            set -eo pipefail
            echo annotate gff for {sampleID}: `date "+%F %T"`\n

            cd {analydir}/SV/{sampleID}/{sv_soft}

            sh {annovar} \\
                -t SVType \\
                {sampleID}.{sv_soft}.gff \\
                {sampleID}.{sv_soft}

            python {moduledir}/Varition/SV/sv_cnv_stat.py \\
                -i {sampleID}.{sv_soft}.hg19_multianno.xls \\
                -s {sampleID} \\
                -soft {sv_soft}

            echo annotate gff for {sampleID} `date "+%F %T"`
            '''.format(
                sampleID=sampleID, **self.__dict__)

        shell_path = '{analydir}/SV/{sampleID}/{sv_soft}/annotate_gff_{sampleID}.sh'.format(
            sampleID=sampleID, **self.__dict__)

        utils.write_shell(shell_path, cmd)

        # add job
        now_point = 'annotate_gff'
        job_name = 'annotate_gff_{}'.format(sampleID)
        utils.add_job(self.jobs, now_point, self.args['startpoint'],
                      self.ANALYSIS_POINTS, job_name, shell_path, self.queues)

        after_jobs = ['data_release', 'primary_report']
        utils.add_order(self.orders, job_name, after_jobs=after_jobs)

# the end
