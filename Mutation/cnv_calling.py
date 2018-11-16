#!/usr/bin/env python
# -*- coding=utf-8 -*-
import os
import utils


class CNV(object):

    def __init__(self, args, jobs, orders, sample_infos, config):

        self.args = args
        self.jobs = jobs
        self.orders = orders

        self.mutation_soft = args['softwares']['mutation']
        self.sv_soft = args['softwares']['sv']
        self.cnv_soft = args['softwares']['cnv']
        self.analydir = args.get('analydir')
        self.queues = args.get('queues')

        self.sample_infos = sample_infos
        self.sample_infos_waiting = args['sample_infos_waiting']
        self.ANALYSIS_POINTS = config.ANALYSIS_POINTS

        self.config = config
        self.soft_dir = config.CONFIG.get('software', 'soft_dir')
        self.circos_dir = config.CONFIG.get('software', 'circos_dir')

        self.__dict__.update(args)

    def start(self):

        print '>  cnv with {} ...'.format(self.cnv_soft)

        sampleIDs = self.sample_infos.keys()
        if self.cnv_soft == 'conifer':
            print '>  cnv call conifer'
            self.conifer_call(sampleIDs)
        elif self.cnv_soft == 'freec':
            print '>  cnv call freec'
            if self.args['seqstrag'] != 'WGS':
                print '[warn] freec is suitable for WGS, not {}'.format(
                    self.args['seqstrag'])
            for sampleID in self.sample_infos_waiting:
                self.freec_call(sampleID)

        # 同时做snpindel, sv, cnv时绘制circos图
        if all([self.mutation_soft, self.sv_soft, self.cnv_soft]) and self.cnv_soft == 'freec':
            print '>  circos start...'
            for sampleID in self.sample_infos_waiting:
                self.circos_draw(sampleID)

    def conifer_call(self, sampleIDs):

        if 'V5' in self.args['TR']:
            probe = 'V5'
        elif 'V6' in self.args['TR']:
            probe = 'V6'
        else:
            print '[Error] Only agilent V5 or V6 can do CoNIFER analysis. '
            exit(1)

        # prepare data for conifer
        outfile = '{analydir}/SV/CoNIFER_{newjob}/sample_for_cnv_call'.format(
            **self.args)
        with utils.safe_open(outfile, 'w') as out:
            for sampleID in sampleIDs:
                bam = '{analydir}/Mapping/{sampleID}.{sampleID}/{sampleID}.final.bam'.format(
                    sampleID=sampleID, analydir=self.analydir)
                out.write('{}\t{}\n'.format(sampleID, bam))

        cmd = '''
            set -eo pipefail
            echo cnv call with conifer start: `date "+%F %T"`\n
            cd {analydir}/SV/CoNIFER_{newjob}

            python {moduledir}/Varition/CNV/CoNIFER/conifer_v0.2.2/conifer.pipe4.7.py \\
                --svd 10 \\
                --probe {probe} \\
                --ref {ref} \\
                --in sample_for_cnv_call \\
                --suffix {newjob} \\
                --out {analydir}/SV

            echo cnv call with conifer done: `date "+%F %T"`
            '''.format(
                probe=probe, **self.__dict__)

        shell_path = '{analydir}/SV/CoNIFER_{newjob}/conifer_call.sh'.format(
            **self.args)

        utils.write_shell(shell_path, cmd)

        # add job
        now_point = 'conifer_call'
        job_name = 'conifer_call'
        utils.add_job(self.jobs, now_point, self.args['startpoint'],
                      self.ANALYSIS_POINTS, job_name, shell_path, self.queues)

        # add order
        before_jobs = ['final_bam_{sampleID}'.format(sampleID=sampleID) for sampleID in sampleIDs]
        after_jobs = ['primary_report']
        utils.add_order(self.orders, job_name, before_jobs=before_jobs, after_jobs=after_jobs)

    def freec_call(self, sampleID):

        seqtype = self.args['seqstrag'].split('_')[0]
        target = ''
        if seqtype != 'WGS':
            target = '\\\n{:16}--target {} '.format(' ', self.args['TR'])
        sex = 'XX' if self.sample_infos[sampleID]['sex'] == 'F' else 'XY'

        cmd = '''
            set -eo pipefail
            echo cnv call with freec for {sampleID} start: `date "+%F %T"`

            cd {analydir}/SV/{sampleID}/freec

            python {moduledir}/Varition/CNV/freec/freec_calling.py \\
                --type {seqtype} {target}\\
                --format BAM \\
                --loh 0 \\
                --contamination 0 \\
                --samName {sampleID} \\
                --sample {analydir}/Mapping/{sampleID}.{sampleID}/{sampleID}.final.bam \\
                --sex {sex} \\
                --ref {ref} \\
                --o .

            echo cnv call with freec for {sampleID} done: `date "+%F %T"`
            '''.format(
                sampleID=sampleID, sex=sex, seqtype=seqtype, target=target, **self.__dict__)

        shell_path = '{analydir}/SV/{sampleID}/freec/freec_call_{sampleID}.sh'.format(
            sampleID=sampleID, **self.args)

        utils.write_shell(shell_path, cmd)

        # add job
        now_point = 'freec_call'
        job_name = 'freec_call_{}'.format(sampleID)
        utils.add_job(self.jobs, now_point, self.args['startpoint'],
                      self.ANALYSIS_POINTS, job_name, shell_path, self.queues)

        # add order
        before_jobs = ['final_bam_{sampleID}'.format(sampleID=sampleID)]
        after_jobs = ['primary_report']
        utils.add_order(self.orders, job_name, before_jobs=before_jobs, after_jobs=after_jobs)

    # ======================== Circos ===============================
    def circos_draw(self, sampleID):

        cmd = '''
            set -eo pipefail
            echo circos drawing for {sampleID} start: `date "+%F %T"`

            cd {analydir}/SV/{sampleID}/freec/Circos

            # ================= Circos Configuration ======================
            # 1 depth info
            python {ROOT_DIR}/modules/scripts/circos_get_info.py \\
                -i {analydir}/Mapping/{sampleID}.{sampleID}/{sampleID}.final.bam \\
                -vt bam

            # 2 snp info
            python {ROOT_DIR}/modules/scripts/circos_get_info.py \\
                -i {analydir}/Mutation/{sampleID}.{mutation_soft}/{sampleID}.{mutation_soft}.snp.vcf.gz \\
                -vt snp

            # 3 indel info
            python {ROOT_DIR}/modules/scripts/circos_get_info.py \\
                -i {analydir}/Mutation/{sampleID}.{mutation_soft}/{sampleID}.{mutation_soft}.indel.vcf.gz \\
                -vt indel

            # 4 sv info
            python {ROOT_DIR}/modules/scripts/circos_get_info.py \\
                -i {analydir}/SV/{sampleID}/{sv_soft}/{sampleID}*.hg19_multianno.xls \\
                -vt {sv_soft}

            # 5 cnv info
            python {ROOT_DIR}/modules/scripts/circos_get_info.py \\
                -i {analydir}/SV/{sampleID}/{cnv_soft}/{sampleID}*.hg19_multianno.xls \\
                -vt {cnv_soft}

            # circos configure all
            cp {circos_dir}/conf/*.conf .
            python {circos_dir}/scripts/config.py \\
                --id {sampleID} \\
                --inpath . \\
                --outpath . \\
                --ref {ref} \\
                --vt depth,snp,indel,{sv_soft},{cnv_soft} \\
                --in {analydir}/Mapping/{sampleID}.{sampleID}/{sampleID}.final.bam,\\
            {analydir}/Mutation/{sampleID}.{mutation_soft}/{sampleID}.{mutation_soft}.snp.vcf.gz,\\
            {analydir}/Mutation/{sampleID}.{mutation_soft}/{sampleID}.{mutation_soft}.indel.vcf.gz,\\
            {analydir}/SV/{sampleID}/{sv_soft}/{sampleID}.{sv_soft}.hg19_multianno.xls,\\
            {analydir}/SV/{sampleID}/{cnv_soft}/{sampleID}.{cnv_soft}.hg19_multianno.xls

            # ================= Circos Drawing ======================
            {circos_dir}/circos-0.67-4/bin/circos \\
                -conf circos_{sampleID}.conf \\
                --outputdir . \\
                --outputfile {sampleID}

            echo circos drawing for {sampleID} done: `date "+%F %T"`
            '''.format(
                sampleID=sampleID, **self.__dict__)

        shell_path = '{analydir}/SV/{sampleID}/freec/Circos/circos_draw_{sampleID}.sh'.format(
            sampleID=sampleID, **self.args)

        utils.write_shell(shell_path, cmd)

        # add job
        now_point = 'circos_draw'
        job_name = 'circos_draw_{}'.format(sampleID)
        utils.add_job(self.jobs, now_point, self.args['startpoint'],
                      self.ANALYSIS_POINTS, job_name, shell_path, self.queues)

        # add order
        if self.mutation_soft == 'gatk':
            mutation_filter = 'gatk_filter_{sampleID}'.format(
                sampleID=sampleID)
        elif self.mutation_soft == 'sentieon':
            mutation_filter = 'extract_annotation_snp_{sampleID} extract_annotation_indel_{sampleID}'.format(
                sampleID=sampleID)
        else:
            mutation_filter = 'bcftools_filter_{sampleID}'.format(sampleID=sampleID)

        before_jobs = 'final_bam_{sampleID} {mutation_filter} annotate_gff_{sampleID} {cnv_soft}_call_{sampleID}'.format(
            cnv_soft=self.cnv_soft, **locals()).split()
        after_jobs = ['primary_report']
        utils.add_order(self.orders, job_name, before_jobs=before_jobs, after_jobs=after_jobs)

# the end
