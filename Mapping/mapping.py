#!/usr/bin/env python
# -*- coding=utf-8 -*-
import utils


class Mapping(object):
    '''
    Mapping Steps:
        clean.fq -> sort.bam -> nodup.bam -> final.bam
    '''

    def __init__(self, args, jobs, orders, sample_lists, sample_infos, config, qc_status, mapping_status):

        self.args = args
        self.__dict__.update(args)

        self.jobs = jobs
        self.orders = orders

        self.analydir = args.get('analydir')

        self.queues = args.get('queues')
        self.sentieon_queues = args.get('sentieon_queues')

        self.sample_lists = sample_lists
        self.sample_infos = sample_infos

        self.qc_status = qc_status
        self.mapping_status = mapping_status

        self.ANALYSIS_POINTS = config.ANALYSIS_POINTS

        self.gatkDIR = config.CONFIG.get('software', 'gatkDIR')

        self.refgenome = args.get('ref')
        self.__dict__.update(
            dict(config.CONFIG.items('genome_' + self.refgenome))
        )
        self.__dict__.update(dict(config.CONFIG.items('software')))

        self.alignment_soft = args['softwares']['alignment']
        self.merge_soft = args['softwares']['merge']
        self.markdup_soft = args['softwares']['markdup']

        self.threads = {
            'bwa_mem': 6,
            'sentieon_bwa_mem': 12,
            'merge': 8,
            'markdup': 8,
        }

        self.rm_clean = 'true' if self.args['ANALY_DICT']['quality_control_rm_clean'] else 'false'

    def start(self):

        print '>  mapping start...'

        for sampleID, items in self.sample_lists.iteritems():
            if sampleID  not in self.args['sample_infos_waiting']:
                continue

            if self.mapping_status == 'done':
                self.final_bam(sampleID, sampleID)
                continue

            is_jiace = items.get('jiace')
            # print 'Jiace', is_jiace
            sort_bams = []
            for lane in items.get('lanes'):
                # print lane
                sort_bam = '{sampleID}_{novoid}_{flowcell}_L{lane}.sort.bam'.format(sampleID=sampleID, **lane)
                patientID = lane['patientID']
                if sort_bam not in sort_bams:
                    sort_bams.append(sort_bam)

                if self.alignment_soft == 'bwa':
                    self.bwa_mem(sampleID, lane)
                elif self.alignment_soft == 'sentieon':
                    self.sentieon_bwa_mem(sampleID, lane)

            # print sampleID, sort_bams
            # merge sortbam and markdup
            if self.merge_soft == 'sambamba':
                self.sambamba_merge(patientID, sampleID, sort_bams, is_jiace)
            elif self.merge_soft == 'picard':
                print 'picard merge to be added...'
                exit(1)

            if self.markdup_soft == 'sambamba':
                self.sambamba_markdup(patientID, sampleID)
            elif self.markdup_soft == 'sentieon':
                self.sentieon_markdup(patientID, sampleID)

            # # sentieon realign and recal
            # if self.args['softwares']['recal'] == 'sentieon':
            #     self.sentieon_recal(patientID, sampleID)
            # elif self.args['softwares']['recal'] == 'gatk':

            # # gatk realign and recal
            # if self.mapping_code == 2.2:
            #     self.gatk_realign(patientID, sampleID)
            #     self.gatk_recal(patientID, sampleID)

            # # gatk just recal
            # if self.mapping_code == 2.3:
            #     self.gatk_recal(patientID, sampleID)

            # final bam
            self.final_bam(patientID, sampleID)

            # align stat
            self.stat_depth(patientID, sampleID)
            self.stat_flag(patientID, sampleID)
            self.stat_uncover(patientID, sampleID)
            self.combine_stat(patientID, sampleID)

            self.mapping_check(patientID, sampleID)

        # self.mapping_report()

    def bwa_mem(self, sampleID, lane):

        # print '  bwa mem...'
        # write shell
        cmd = '''
            set -eo pipefail
            echo bwa mem and samtools sort for {sampleID} start: `date "+%F %T"`

            cd {analydir}/Mapping/{patientID}.{sampleID}

            fq1={analydir}/QC/{sampleID}/{sampleID}_{novoid}_{flowcell}_L{lane}_1.clean.fq
            fq2={analydir}/QC/{sampleID}/{sampleID}_{novoid}_{flowcell}_L{lane}_2.clean.fq
            if [ ! -f $fq1 ];then
                fq1=$fq1.gz
                fq2=$fq2.gz
            fi

            bwa mem \\
                -t {bwa_threads} -M \\
                -R "@RG\\tID:{sampleID}_{novoid}_{flowcell}_L{lane}\\tSM:{sampleID}\\tLB:{sampleID}\\tPU:{novoid}_{flowcell}_L{lane}\\tPL:illumina\\tCN:novogene" \\
                {reffasta} \\
                $fq1 $fq2 |
            samtools-1.6 view \\
                -@ {bwa_threads} -b -S -t \\
                {reffasta}.fai |
            samtools-1.6 sort \\
                -@ {bwa_threads} -m 2G \\
                -T {sampleID}_{novoid}_{flowcell}_L{lane}.tmp \\
                -o {sampleID}_{novoid}_{flowcell}_L{lane}.sort.bam

            echo bwa mem and samtools sort for {sampleID} done: `date "+%F %T"`
        '''.format(
            sampleID=sampleID,
            bwa_threads=self.threads['bwa_mem'],
            **dict(lane, **self.__dict__))

        shell_path = '{analydir}/Mapping/{patientID}.{sampleID}/bwa_mem_{sampleID}_{novoid}_{flowcell}_L{lane}.sh'.format(
            analydir=self.analydir, sampleID=sampleID, **lane)

        utils.write_shell(shell_path, cmd)

        # add job
        now_point = 'bwa_mem'
        job_name = 'bwa_mem_{sampleID}_{novoid}_{flowcell}_L{lane}'.format(sampleID=sampleID, **lane)
        utils.add_job(
            self.jobs,
            now_point,
            self.args['startpoint'],
            self.ANALYSIS_POINTS,
            job_name,
            shell_path,
            self.queues,
            threads=self.threads['bwa_mem'])

        # add order
        before_jobs = []
        if self.qc_status == 'waiting':
            before_jobs = ['qc_{sampleID}_{novoid}_{flowcell}_L{lane}'.format(sampleID=sampleID, **lane)]
        after_jobs = ['{merge_soft}_merge_{sampleID}'.format(sampleID=sampleID, **self.__dict__)]

        if self.qc_status == 'waiting' and (self.args['ANALY_DICT']['quality_control'] or self.args['ANALY_DICT']['quality_control_keep_clean']):
            after_jobs += [
                'gzip_md5_clean_{sampleID}_{novoid}_{flowcell}_L{lane}'.format(sampleID=sampleID, **lane)
            ]

        utils.add_order(self.orders, job_name, before_jobs=before_jobs, after_jobs=after_jobs)

    def sentieon_bwa_mem(self, sampleID, lane):

        # print '  sentieon bwa mem and sort...'
        # write shell
        cmd = '''
            set -eo pipefail
            echo sentieon bwa mem and sort for {sampleID} start: `date "+%F %T"`

            cd {analydir}/Mapping/{patientID}.{sampleID}

            fq1={analydir}/QC/{sampleID}/{sampleID}_{novoid}_{flowcell}_L{lane}_1.clean.fq
            fq2={analydir}/QC/{sampleID}/{sampleID}_{novoid}_{flowcell}_L{lane}_2.clean.fq
            if [ ! -f $fq1 ];then
                fq1=$fq1.gz
                fq2=$fq2.gz
            fi

            sentieon bwa mem \\
                -t {bwa_threads} -M \\
                -K 10000000 \\
                -R "@RG\\tID:{sampleID}_{novoid}_{flowcell}_L{lane}\\tSM:{sampleID}\\tLB:{sampleID}\\tPU:{novoid}_{flowcell}_L{lane}\\tPL:illumina\\tCN:novogene" \\
                {reffasta} \\
                $fq1 $fq2 |
            sentieon util sort \\
                -t {bwa_threads} --sam2bam -i - \\
                -r {reffasta} \\
                -o {sampleID}_{novoid}_{flowcell}_L{lane}.sort.bam

            echo sentieon bwa mem and sort for {sampleID} done: `date "+%F %T"`
        '''.format(
            sampleID=sampleID,
            bwa_threads=self.threads['sentieon_bwa_mem'],
            **dict(lane, **self.__dict__))

        shell_path = '{analydir}/Mapping/{patientID}.{sampleID}/sentieon_bwa_mem_{sampleID}_{novoid}_{flowcell}_L{lane}.sh'.format(
            analydir=self.analydir, sampleID=sampleID, **lane)

        utils.write_shell(shell_path, cmd)

        # add job
        now_point = 'sentieon_bwa_mem'
        job_name = 'sentieon_bwa_mem_{sampleID}_{novoid}_{flowcell}_L{lane}'.format(sampleID=sampleID, **lane)

        utils.add_job(
            self.jobs,
            now_point,
            self.args['startpoint'],
            self.ANALYSIS_POINTS,
            job_name,
            shell_path,
            self.sentieon_queues,
            threads=self.threads['sentieon_bwa_mem'])

        # add order
        before_jobs = []
        if self.qc_status == 'waiting':
            before_jobs = ['qc_{sampleID}_{novoid}_{flowcell}_L{lane}'.format(sampleID=sampleID, **lane)]
        after_jobs = ['{merge_soft}_merge_{sampleID}'.format(sampleID=sampleID, **self.__dict__)]

        if self.qc_status == 'waiting' and (self.args['ANALY_DICT']['quality_control'] or self.args['ANALY_DICT']['quality_control_keep_clean']):
            after_jobs += [
                'gzip_md5_clean_{sampleID}_{novoid}_{flowcell}_L{lane}'.format(
                    sampleID=sampleID, **lane)
            ]

        utils.add_order(self.orders, job_name, before_jobs=before_jobs, after_jobs=after_jobs)

    def sambamba_merge(self, patientID, sampleID, sort_bams, is_jiace=False):

        # print '  sambamba merge...'
        # write shell
        print sort_bams

        merge_threads = self.threads['merge']

        sort_bams_list = ' '.join(sort_bams)
        if is_jiace:
            cmd = '''
                set -eo pipefail
                echo sambamba merge for {sampleID} start: `date "+%F %T"`
                
                cd {analydir}/Mapping/{patientID}.{sampleID}

                mv {sampleID}.sort.bam {sampleID}.old.bam

                sambamba merge \\
                    -t {merge_threads} \\
                    {sampleID}.sort.bam \\
                    {sampleID}.old.bam \\
                    {sort_bams_list}

                rm -f {sort_bams_list}

                echo sambamba merge for {sampleID} done: `date "+%F %T"`
            '''
        else:
            if len(sort_bams) == 1:
                cmd = '''
                    set -eo pipefail
                    echo rename sortbam for {sampleID} start: `date "+%F %T"`

                    cd {analydir}/Mapping/{patientID}.{sampleID}

                    mv {sort_bams_list} {sampleID}.sort.bam
                    
                    sambamba index {sampleID}.sort.bam

                    echo rename sortbam for {sampleID} done: `date "+%F %T"`
                '''
            else:
                cmd = '''
                    set -eo pipefail
                    echo sambamba merge for {sampleID} start: `date "+%F %T"`

                    cd {analydir}/Mapping/{patientID}.{sampleID}

                    sambamba merge \\
                        -t {merge_threads} \\
                        {sampleID}.sort.bam \\
                        {sort_bams_list}

                    rm -f {sort_bams_list}

                    echo sambamba merge for {sampleID} done: `date "+%F %T"`
                '''

        cmd = cmd.format(**dict(self.__dict__, **locals()))

        shell_path = '{analydir}/Mapping/{patientID}.{sampleID}/sambamba_merge_{sampleID}.sh'.format(**dict(self.__dict__, **locals()))

        utils.write_shell(shell_path, cmd)

        # add job
        now_point = 'sambamba_merge'
        job_name = 'sambamba_merge_{sampleID}'.format(sampleID=sampleID)
        utils.add_job(
            self.jobs,
            now_point,
            self.args['startpoint'],
            self.ANALYSIS_POINTS,
            job_name,
            shell_path,
            self.queues,
            threads=self.threads['merge'])

        # add order
        before_jobs = []
        after_jobs = ['stat_depth_{sampleID}'.format(**locals())]
        utils.add_order(self.orders, job_name, before_jobs=before_jobs, after_jobs=after_jobs)


    def sambamba_markdup(self, patientID, sampleID):

        # print '  sambamba merge...'
        # write shell
        markdup_threads = self.threads['markdup']

        cmd = '''
            set -eo pipefail
            echo sambamba markdup for {sampleID} start: `date "+%F %T"`

            cd {analydir}/Mapping/{patientID}.{sampleID}

            sambamba markdup \\
                -t {markdup_threads} \\
                --overflow-list-size=10000000 \\
                --tmpdir=tmp \\
                {sampleID}.sort.bam \\
                {sampleID}.nodup.bam

            rm -rf tmp

            echo sambamba markdup for {sampleID} done: `date "+%F %T"`
        '''.format(**dict(self.__dict__, **locals()))

        shell_path = '{analydir}/Mapping/{patientID}.{sampleID}/sambamba_markdup_{sampleID}.sh'.format(**dict(self.__dict__, **locals()))

        utils.write_shell(shell_path, cmd)

        # add job
        now_point = 'sambamba_markdup'
        job_name = 'sambamba_markdup_{sampleID}'.format(sampleID=sampleID)

        utils.add_job(
            self.jobs,
            now_point,
            self.args['startpoint'],
            self.ANALYSIS_POINTS,
            job_name,
            shell_path,
            self.queues,
            threads=markdup_threads)

        # add order
        before_jobs = ['{merge_soft}_merge_{sampleID}'.format(sampleID=sampleID, **self.__dict__)]
        after_jobs = 'mapping_check_{sampleID} stat_flag_{sampleID}'.format(sampleID=sampleID).split()
        utils.add_order(self.orders, job_name, before_jobs=before_jobs, after_jobs=after_jobs)

    def sentieon_markdup(self, patientID, sampleID):

        # print '  sentieon markdup...'
        # write shell
        markdup_threads = self.threads['markdup']

        cmd = '''
            set -eo pipefail            
            echo sentieon markdup for {sampleID} start: `date "+%F %T"`

            cd {analydir}/Mapping/{patientID}.{sampleID}

            sentieon driver \\
                -t {markdup_threads} \\
                -i {sampleID}.sort.bam \\
                --algo LocusCollector \\
                --fun score_info \\
                {sampleID}.score.txt

            sentieon driver \\
                -t {markdup_threads} \\
                -i {sampleID}.sort.bam \\
                --algo Dedup \\
                --score_info {sampleID}.score.txt \\
                --metrics {sampleID}.dedup.metrics.txt \\
                {sampleID}.nodup.bam

            echo sentieon markdup for {sampleID} done: `date "+%F %T"`
        '''.format(**dict(self.__dict__, **locals()))

        shell_path = '{analydir}/Mapping/{patientID}.{sampleID}/sentieon_markdup_{sampleID}.sh'.format(
            **dict(self.__dict__, **locals()))

        utils.write_shell(shell_path, cmd)

        # add job
        now_point = 'sentieon_markdup'
        job_name = 'sentieon_markdup_{sampleID}'.format(sampleID=sampleID)
        utils.add_job(
            self.jobs,
            now_point,
            self.args['startpoint'],
            self.ANALYSIS_POINTS,
            job_name,
            shell_path,
            self.sentieon_queues,
            threads=markdup_threads)

        # add order
        before_jobs = ['{merge_soft}_merge_{sampleID}'.format(sampleID=sampleID, **self.__dict__)]
        after_jobs = 'mapping_check_{sampleID} stat_flag_{sampleID}'.format(sampleID=sampleID).split()
        utils.add_order(self.orders, job_name, before_jobs=before_jobs, after_jobs=after_jobs)

    def final_bam(self, patientID, sampleID):

        # print '  fianl bam...'
        # write shell

        cmd = '''
            set -eo pipefail
            echo final bam for {sampleID} start: `date "+%F %T"`
            
            cd {analydir}/Mapping/{patientID}.{sampleID}
            
            ln -sf {sampleID}.nodup.bam {sampleID}.final.bam
            ln -sf {sampleID}.nodup.bam.bai {sampleID}.final.bam.bai
            
            echo final bam for {sampleID} done: `date "+%F %T"`
        '''.format(**dict(self.__dict__, **locals()))

        shell_path = '{analydir}/Mapping/{patientID}.{sampleID}/final_bam_{sampleID}.sh'.format(**dict(self.__dict__, **locals()))

        utils.write_shell(shell_path, cmd)

        # add job
        now_point = 'final_bam'
        job_name = 'final_bam_{sampleID}'.format(sampleID=sampleID)
        utils.add_job(self.jobs, now_point, self.args['startpoint'],
                      self.ANALYSIS_POINTS, job_name, shell_path, self.queues)

    # ===== Alignment Stat =====
    def stat_depth(self, patientID, sampleID):

        # based on sort.bam
        # print '  stat depth...'
        # write shell
        cmd = '''
            set -eo pipefail
            echo stat depth for {sampleID} start: `date "+%F %T"`
            
            cd {analydir}/Alnstat/{sampleID}
            
            perl {moduledir}/Alnstat/depthStat_pipe4.6.pl \\
                -s {sampleID} \\
                -g {ref} \\'''

        if self.args['seqstrag'] != 'WGS':
            cmd += '''
                -r {TR} \\'''

        cmd += '''
                {analydir}/Mapping/{patientID}.{sampleID}/{sampleID}.sort.bam \\
                .
    
            echo stat depth for {sampleID} done: `date "+%F %T"`
        '''
        
        cmd = cmd.format(
            patientID=patientID, sampleID=sampleID, **self.args)

        shell_path = '{analydir}/Alnstat/{sampleID}/stat_depth_{sampleID}.sh'.format(
            analydir=self.analydir, sampleID=sampleID)

        utils.write_shell(shell_path, cmd)

        # add job
        now_point = 'stat_depth'
        job_name = 'stat_depth_{sampleID}'.format(sampleID=sampleID)
        utils.add_job(self.jobs, now_point, self.args['startpoint'],
                      self.ANALYSIS_POINTS, job_name, shell_path, self.queues)

    def stat_uncover(self, patientID, sampleID):

        # based on sort.bam
        # print '  stat uncover...'
        # write shell
        if self.args['seqstrag'] != 'WGS':
            cmd = '''
                set -eo pipefail
                echo stat uncover for {sampleID} start: `date "+%F %T"`

                cd {analydir}/Alnstat/{sampleID}

                samtools-1.6 depth \\
                    -aa -q 0 -Q 0 \\
                    -b {TR} \\
                    {analydir}/Mapping/{patientID}.{sampleID}/{sampleID}.sort.bam |
                awk -F'\\t' '$3==0' |
                grep -vwf target_region.00.depth > target_region.0.depth

                python {moduledir}/Alnstat/uncover_pos_chr_pipe4.6.py \\
                    target_region.0.depth \\
                    {sampleID} \\
                    {sampleID}.uncovered_region.annovar.result.xls

                rm -f target_region.0.depth

                echo stat uncover for {sampleID} done: `date "+%F %T"`
            '''
        else:
            cmd = '''
                set -eo pipefail
                echo stat uncover for {sampleID} start: `date "+%F %T"`
            
                rm -f *.depth *.bed *.pdf* *.png

                echo stat uncover for {sampleID} done: `date "+%F %T"`
            '''

        cmd = cmd.format(
            patientID=patientID,
            sampleID=sampleID,
            **self.args)

        shell_path = '{analydir}/Alnstat/{sampleID}/stat_uncover_{sampleID}.sh'.format(
            analydir=self.analydir, sampleID=sampleID)

        utils.write_shell(shell_path, cmd)

        # add job
        now_point = 'stat_uncover'
        job_name = 'stat_uncover_{sampleID}'.format(sampleID=sampleID)
        utils.add_job(self.jobs, now_point, self.args['startpoint'],
                      self.ANALYSIS_POINTS, job_name, shell_path, self.queues)

        # add order
        before_jobs = ['stat_depth_{sampleID}'.format(sampleID=sampleID)]
        after_jobs = []
        utils.add_order(self.orders, job_name, before_jobs=before_jobs, after_jobs=after_jobs)

    def stat_flag(self, patientID, sampleID):

        # based on nodup.bam
        # print '  stat flag...'
        # write shell
        cmd = '''
            set -eo pipefail
            echo stat flag for {sampleID} start: `date "+%F %T"`

            cd {analydir}/Alnstat/{sampleID}

            python {moduledir}/Alnstat/sam_flagstat_pipe4.6.py \\
                --bam {analydir}/Mapping/{patientID}.{sampleID}/{sampleID}.nodup.bam \\
                > {sampleID}.flagstat

            echo stat flag for {sampleID} done: `date "+%F %T"`
        '''.format(
            patientID=patientID, sampleID=sampleID, **self.args)

        shell_path = '{analydir}/Alnstat/{sampleID}/stat_flag_{sampleID}.sh'.format(
            analydir=self.analydir, sampleID=sampleID)

        utils.write_shell(shell_path, cmd)

        # add job
        now_point = 'stat_flag'
        job_name = 'stat_flag_{sampleID}'.format(sampleID=sampleID)
        utils.add_job(self.jobs, now_point, self.args['startpoint'],
                      self.ANALYSIS_POINTS, job_name, shell_path, self.queues)

    def combine_stat(self, patientID, sampleID):

        # print '  combine stat result...'
        # write shell
        cmd = '''
            set -eo pipefail
            echo combine stat result for {sampleID} start: `date "+%F %T"`\n
            cd {analydir}/Alnstat/{sampleID}\n
            python {moduledir}/Alnstat/combine_pipe4.6.py \\
                information.xlsx \\
                {sampleID}.flagstat \\
                {sampleID} \\
                {seqstrag} \\
                > {sampleID}_mapping_coverage.txt\n
            echo combine stat result for {sampleID} done: `date "+%F %T"`
        '''.format(
            patientID=patientID, sampleID=sampleID, **self.args)

        shell_path = '{analydir}/Alnstat/{sampleID}/combine_stat_{sampleID}.sh'.format(
            analydir=self.analydir, sampleID=sampleID)

        utils.write_shell(shell_path, cmd)

        # add job
        now_point = 'combine_stat'
        job_name = 'combine_stat_{sampleID}'.format(sampleID=sampleID)

        utils.add_job(self.jobs, now_point, self.args['startpoint'],
                      self.ANALYSIS_POINTS, job_name, shell_path, self.queues)
        # add order
        before_jobs = 'stat_depth_{sampleID} stat_flag_{sampleID}'.format(
            sampleID=sampleID).split()

        after_jobs = ['mapping_check_{sampleID}'.format(sampleID=sampleID), 'mapping_report']
        utils.add_order(self.orders, job_name, before_jobs=before_jobs, after_jobs=after_jobs)

    def mapping_check(self, patientID, sampleID):

        # print '  mapping check...'
        # write shell
        sex = self.sample_infos[sampleID]['sex']

        cmd = '''
            set -eo pipefail
            echo mapping check for sample {sampleID} start: `date "+%F %T"`

            python2 {moduledir}/QC/auto_check.py \\
                --qc_list {qc_list} \\
                --sampid {sampleID} \\
                --pwd {analydir} \\
                --check map \\
                --jobname {newjob} \\
                --seqstrag {seqstrag} \\
                --email {email} \\
                --PE {PE} \\
                --gender {sex} \\
                --dup {dup} \\
                --depth {depth}
            
            # remove clean data if mapping check passed
            if {rm_clean};then
                rm -f {analydir}/QC/{sampleID}/{sampleID}_*.clean.fq*
            fi

            echo mapping check for sample {sampleID} done: `date "+%F %T"`
        '''.format(**dict(self.__dict__, **locals()))

        shell_path = '{analydir}/Mapping/{patientID}.{sampleID}/mapping_check_{sampleID}.sh'.format(**dict(self.__dict__, **locals()))

        utils.write_shell(shell_path, cmd)

        # add job
        now_point = 'mapping_check'
        job_name = 'mapping_check_{sampleID}'.format(sampleID=sampleID)
        utils.add_job(self.jobs, now_point, self.args['startpoint'],
                      self.ANALYSIS_POINTS, job_name, shell_path, self.queues)

        # add order
        before_jobs = []
        after_jobs = ['final_bam_{sampleID}'.format(sampleID=sampleID), 'data_release']
        utils.add_order(self.orders, job_name, before_jobs=before_jobs, after_jobs=after_jobs)

# the end
