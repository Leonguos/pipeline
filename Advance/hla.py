#!/usr/bin/env python
# -*- coding=utf-8 -*-
import glob
import utils


class HLA(object):

    def __init__(self, args, jobs, orders, sample_lists, sample_infos, config, qc_status):

        self.args = args
        self.jobs = jobs
        self.orders = orders
        self.queues = args.get('queues')

        self.sample_lists = sample_lists
        self.sample_infos = sample_infos

        self.qc_status = qc_status
        
        self.ANALYSIS_POINTS = config.ANALYSIS_POINTS

        self.__dict__.update(args)
        self.refgenome = args.get('ref')
        self.__dict__.update(
            dict(config.CONFIG.items('genome_' + self.refgenome)))
        self.__dict__.update(dict(config.CONFIG.items('software')))

    def start(self):

        print '>  hla start ...'

        hla_beds = glob.glob('{athlates_db_dir}/bed/*.non*.bed'.format(**self.__dict__))
        hla_genes = map(lambda x: x.rsplit('.')[-2].strip('non-'), hla_beds)
        
        if self.args['hla_gene']:
            new_hla_genes = []
            for gene in self.args['hla_gene'].split(','):
                if gene not in hla_genes:
                    print '[error] your supply gene not in database: {gene}'.format(**locals())
                    exit(1)
                new_hla_genes.append(gene)
            hla_genes = new_hla_genes

        print '>  hla genelist: {}'.format(sorted(hla_genes))

        for sampleid, items in self.sample_lists.iteritems():

            sort_bams = []
            lanes = items.get('lanes')
            # print sampleid, lanes
            for lane in lanes:
                # print lane
                sort_bam = '{sampleid}_{novoid}_{flowcell}_L{lane}.sort.bam'.format(
                    sampleid=sampleid, **lane)
                # patientid = lane['patientID']
                if sort_bam not in sort_bams:
                    sort_bams.append(sort_bam)

                self.hla_bwa_mem(sampleid, lane)

            self.hla_sambamba_merge(sampleid, sort_bams)
            # self.hla_sambamba_markdup(sampleid)
            self.hla_picard_markdup(sampleid)

            for gene in hla_genes:
                self.hla_sort_by_name(sampleid, gene)
                self.hla_athlates_typing(sampleid, gene)

            self.hla_hlahd_typing(sampleid, lanes)

    def hla_bwa_mem(self, sampleid, lane):

        # print '>  hla bwa mem ...'
        cmd = '''
            set -eo pipefail
            echo hla bwa mem and samtools sort for {sampleid} start: `date "+%F %T"`

            cd {analydir}/Advance/{newjob}/HLA/ATHLATES_typing/{sampleid}

            fq1={analydir}/QC/{sampleid}/{sampleid}_{novoid}_{flowcell}_L{lane}_1.clean.fq
            fq2={analydir}/QC/{sampleid}/{sampleid}_{novoid}_{flowcell}_L{lane}_2.clean.fq
            if [ ! -f $fq1 ];then
                fq1=$fq1.gz
                fq2=$fq2.gz
            fi

            bwa mem \\
                -t 6 -M \\
                -R "@RG\\tID:{sampleid}_{novoid}_{flowcell}_L{lane}\\tSM:{sampleid}\\tLB:{sampleid}\\tPU:{novoid}_{flowcell}_L{lane}\\tPL:illumina\\tCN:novogene" \\
                {athlates_db_dir}/ref/hla_nclean.fasta \\
                $fq1 $fq2 |
            samtools-1.6 view \\
                -@ 5 -b -S -F 4 -t \\
                {athlates_db_dir}/ref/hla_nclean.fasta.fai |
            samtools-1.6 sort \\
                -@ 3 -m 2G \\
                -T {sampleid}_{novoid}_{flowcell}_L{lane}.tmp \\
                -o {sampleid}_{novoid}_{flowcell}_L{lane}.sort.bam

            echo hla bwa mem and samtools sort for {sampleid} done: `date "+%F %T"`
        '''.format(
            sampleid=sampleid, **dict(lane, **self.__dict__))

        shell_path = '{analydir}/Advance/{newjob}/HLA/ATHLATES_typing/{sampleid}/hla_bwa_mem_{sampleid}_{novoid}_{flowcell}_L{lane}.sh'.format(
            sampleid=sampleid, **dict(lane, **self.__dict__))

        utils.write_shell(shell_path, cmd)

        # add job
        now_point = 'hla_bwa_mem'
        job_name = 'hla_bwa_mem_{sampleid}_{novoid}_{flowcell}_L{lane}'.format(
            sampleid=sampleid, **lane)
        utils.add_job(self.jobs, now_point, self.args['startpoint'],
                      self.ANALYSIS_POINTS, job_name, shell_path, self.queues)
        # add order
        before_jobs = []
        if self.qc_status == 'waiting':
            before_jobs = ['qc_{sampleid}_{novoid}_{flowcell}_L{lane}'.format(sampleid=sampleid, **lane)]

        after_jobs = ['hla_sambamba_merge_{sampleid}'.format(sampleid=sampleid, **lane)]
        utils.add_order(self.orders, job_name, before_jobs=before_jobs, after_jobs=after_jobs)

    def hla_sambamba_merge(self, sampleid, sort_bams):

        # print '  sambamba merge...'
        # write shell
        if len(sort_bams) == 1:
            cmd = '''
                set -eo pipefail
                echo rename sortbam for {sampleid} start: `date "+%F %T"`

                cd {analydir}/Advance/{newjob}/HLA/ATHLATES_typing/{sampleid}

                mv {sort_bam} {sampleid}.sort.bam

                samtools-1.6 index -@ 4 {sampleid}.sort.bam

                echo rename sortbam for {sampleid} done: `date "+%F %T"`
            '''.format(
                sampleid=sampleid,
                sort_bam=sort_bams[0],
                **self.__dict__)
        else:
            cmd = '''
                set -eo pipefail
                echo hla sambamba merge for {sampleid} start: `date "+%F %T"`

                cd {analydir}/Advance/{newjob}/HLA/ATHLATES_typing/{sampleid}

                sambamba merge \\
                    -t 4 \\
                    {sampleid}.sort.bam \\
                    {sort_bams}

                echo hla sambamba merge for {sampleid} done: `date "+%F %T"`
            '''.format(
                sampleid=sampleid,
                sort_bams=' '.join(sort_bams),
                **self.__dict__)

        shell_path = '{analydir}/Advance/{newjob}/HLA/ATHLATES_typing/{sampleid}/hla_sambamba_merge_{sampleid}.sh'.format(
            sampleid=sampleid, **self.__dict__)

        utils.write_shell(shell_path, cmd)

        # add job
        now_point = 'hla_sambamba_merge'
        job_name = 'hla_sambamba_merge_{sampleid}'.format(sampleid=sampleid)
        utils.add_job(self.jobs, now_point, self.args['startpoint'],
                      self.ANALYSIS_POINTS, job_name, shell_path, self.queues)

    def hla_sambamba_markdup(self, sampleid):

        # print '  sambamba merge...'
        # write shell
        cmd = '''
            set -eo pipefail
            echo sambamba markdup for {sampleid} start: `date "+%F %T"`

            cd {analydir}/Advance/{newjob}/HLA/ATHLATES_typing/{sampleid}

            sambamba markdup \\
                -t 5 \\
                --overflow-list-size=10000000 \\
                --tmpdir=tmp \\
                {sampleid}.sort.bam \\
                {sampleid}.nodup.bam

            rm -rf tmp

            echo sambamba markdup for {sampleid} done: `date "+%F %T"`
        '''.format(
            sampleid=sampleid, **self.__dict__)

        shell_path = '{analydir}/Advance/{newjob}/HLA/ATHLATES_typing/{sampleid}/hla_sambamba_markdup_{sampleid}.sh'.format(
            sampleid=sampleid, **self.__dict__)

        utils.write_shell(shell_path, cmd)

        # add job
        now_point = 'hla_sambamba_markdup'
        job_name = 'hla_sambamba_markdup_{sampleid}'.format(sampleid=sampleid)

        utils.add_job(self.jobs, now_point, self.args['startpoint'],
                      self.ANALYSIS_POINTS, job_name, shell_path, self.queues)

        # add order
        before_jobs = ['hla_sambamba_merge_{sampleid}'.format(sampleid=sampleid)]
        after_jobs = []
        utils.add_order(self.orders, job_name, before_jobs=before_jobs, after_jobs=after_jobs)

    def hla_picard_markdup(self, sampleid):

        # print '  picard markdup ...'
        cmd = '''
            set -eo pipefail
            echo picard markdup for {sampleid} start: `date "+%F %T"`

            cd {analydir}/Advance/{newjob}/HLA/ATHLATES_typing/{sampleid}

            java1.8.0 -Xmx5g -jar {picard_jar} \\
                MarkDuplicates \\
                TMP_DIR=TMP \\
                INPUT={sampleid}.sort.bam \\
                OUTPUT={sampleid}.nodup.bam \\
                METRICS_FILE={sampleid}.nodup.metrics.txt \\
                CREATE_INDEX=true \\
                ASSUME_SORTED=true

            echo picard markdup for {sampleid} done: `date "+%F %T"`
        '''.format(
            sampleid=sampleid, **self.__dict__)

        shell_path = '{analydir}/Advance/{newjob}/HLA/ATHLATES_typing/{sampleid}/hla_picard_markdup_{sampleid}.sh'.format(
            sampleid=sampleid, **self.__dict__)

        utils.write_shell(shell_path, cmd)

        # add job
        now_point = 'hla_picard_markdup'
        job_name = 'hla_picard_markdup_{sampleid}'.format(sampleid=sampleid)

        utils.add_job(self.jobs, now_point, self.args['startpoint'],
                      self.ANALYSIS_POINTS, job_name, shell_path, self.queues)

        # add order
        before_jobs = ['hla_sambamba_merge_{sampleid}'.format(sampleid=sampleid)]
        after_jobs = []
        utils.add_order(self.orders, job_name, before_jobs=before_jobs, after_jobs=after_jobs)

    def hla_sort_by_name(self, sampleid, gene):

        # print '  sort by name ...'
        cmd = '''
            set -eo pipefail
            echo hla sort by name for {sampleid} {gene} start: `date "+%F %T"`

            cd {analydir}/Advance/{newjob}/HLA/ATHLATES_typing/{sampleid}/gene/{gene}

            mkdir -p TMP

            for gene in {gene} non-{gene};do

                samtools-1.6 view \\
                    -b -L {athlates_db_dir}/bed/hla.$gene.bed \\
                    -o {sampleid}.$gene.bam \\
                    -@ 4 \\
                    ../../{sampleid}.nodup.bam

                (
                    samtools-1.6 view -H {sampleid}.$gene.bam
                    samtools-1.6 view {sampleid}.$gene.bam | sort -k1,1 -k3,3 -T TMP
                ) | samtools-1.6 view -bS -o {sampleid}.$gene.sort.bam -@ 4 -

                rm -f {sampleid}.$gene.bam

            done

            rm -rf TMP

            rm -f ../../{sampleid}.nodup.bam

            echo hla sort by name for {sampleid} {gene} done: `date "+%F %T"`
        '''.format(
            **dict(self.__dict__, **locals()))

        shell_path = '{analydir}/Advance/{newjob}/HLA/ATHLATES_typing/{sampleid}/gene/{gene}/hla_sort_by_name_{sampleid}_{gene}.sh'.format(
            **dict(self.__dict__, **locals()))

        utils.write_shell(shell_path, cmd)

        # add job
        now_point = 'hla_sort_by_name'
        job_name = 'hla_sort_by_name_{sampleid}_{gene}'.format(**locals())

        utils.add_job(self.jobs, now_point, self.args['startpoint'],
                      self.ANALYSIS_POINTS, job_name, shell_path, self.queues)

        # add order
        before_jobs = ['hla_picard_markdup_{sampleid}'.format(**locals())]
        after_jobs = []
        utils.add_order(self.orders, job_name, before_jobs=before_jobs, after_jobs=after_jobs)

    def hla_athlates_typing(self, sampleid, gene):

        # print '  athlates typing ...'
        cmd = '''
            set -eo pipefail
            echo hla athlates typing for {sampleid} {gene} start: `date "+%F %T"`

            cd {analydir}/Advance/{newjob}/HLA/ATHLATES_typing/{sampleid}/gene/{gene}

            typing \\
                -hd 2 \\
                -msa {athlates_db_dir}/msa/{gene}_nuc.txt \\
                -bam {sampleid}.{gene}.sort.bam \\
                -exlbam {sampleid}.non-{gene}.sort.bam \\
                -o {sampleid}.{gene}

            rm -f *bam

            # link result
            mkdir -p {analydir}/Advance/{newjob}/HLA/ATHLATES_typing/result

            cd {analydir}/Advance/{newjob}/HLA/ATHLATES_typing/result

            ln -sf ../{sampleid}/gene/{gene}/{sampleid}.{gene}.typing.txt .

            echo hla athlates typing for {sampleid} {gene} done: `date "+%F %T"`
        '''.format(
            **dict(self.__dict__, **locals()))

        shell_path = '{analydir}/Advance/{newjob}/HLA/ATHLATES_typing/{sampleid}/gene/{gene}/hla_athlates_typing_{sampleid}_{gene}.sh'.format(
            **dict(self.__dict__, **locals()))

        utils.write_shell(shell_path, cmd)

        # add job
        now_point = 'hla_athlates_typing'
        job_name = 'hla_athlates_typing_{sampleid}_{gene}'.format(**locals())

        utils.add_job(self.jobs, now_point, self.args['startpoint'],
                      self.ANALYSIS_POINTS, job_name, shell_path, self.queues)

        # add order
        before_jobs = ['hla_sort_by_name_{sampleid}_{gene}'.format(**locals())]
        after_jobs = ['data_release']
        utils.add_order(self.orders, job_name, before_jobs=before_jobs, after_jobs=after_jobs)

    def hla_hlahd_typing(self, sampleid, lanes):

        read1_list, merge_read1, merge_read2, merge_read1_gz, merge_read2_gz = utils.get_merge_read(
            self.args['analydir'], sampleid, lanes)

        # print '  hlahd typing ...'
        cmd = '''
            set -eo pipefail
            echo hla hlahd typing for {sampleid} start: `date "+%F %T"`

            cd {analydir}/Advance/{newjob}/HLA/HLAHD_typing/{sampleid}

            if [ -f {read1_1} ];then
                {merge_read1}
                {merge_read2}
            else
                {merge_read1_gz}
                {merge_read2_gz}
            fi

            export PATH={hlahd_dir}/bin:$PATH

            hlahd.sh \\
                -t 2 \\
                -m 150 \\
                -f {hlahd_dir}/freq_data \\
                {sampleid}.R1.fastq {sampleid}.R2.fastq \\
                {hlahd_dir}/HLA_gene.split.txt \\
                {hlahd_dir}/dictionary \\
                {sampleid} \\
                .

            rm -f {sampleid}.R[12].fastq

            # link result
            mkdir -p {analydir}/Advance/{newjob}/HLA/HLAHD_typing/result

            cd {analydir}/Advance/{newjob}/HLA/HLAHD_typing/result

            ln -sf ../{sampleid}/{sampleid}/result/{sampleid}_final.result.txt .

            echo hla hlahd typing for {sampleid} done: `date "+%F %T"`
        '''.format(
            read1_1=read1_list[0], **dict(self.__dict__, **locals()))

        shell_path = '{analydir}/Advance/{newjob}/HLA/HLAHD_typing/{sampleid}/hla_hlahd_typing_{sampleid}.sh'.format(
            sampleid=sampleid, **self.__dict__)

        utils.write_shell(shell_path, cmd)

        # add job
        now_point = 'hla_hlahd_typing'
        job_name = 'hla_hlahd_typing_{sampleid}'.format(sampleid=sampleid)

        utils.add_job(self.jobs, now_point, self.args['startpoint'],
                      self.ANALYSIS_POINTS, job_name, shell_path, self.queues)

        # add order
        before_jobs = []
        if self.qc_status == 'waiting':
            for lane in lanes:
                before_jobs.append('qc_{sampleid}_{novoid}_{flowcell}_L{lane}'.format(sampleid=sampleid, **lane))
        after_jobs = ['data_release']
        utils.add_order(self.orders, job_name, before_jobs=before_jobs, after_jobs=after_jobs)

# the end
