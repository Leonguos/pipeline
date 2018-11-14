#!/usr/bin/env python
# -*- coding=utf-8 -*-
import os
import utils


class Mutation(object):

    def __init__(self, args, jobs, orders, sample_lists, sample_infos, config):

        self.args = args
        self.jobs = jobs
        self.orders = orders

        self.mutation_soft = args['softwares']['mutation']
        self.analydir = args.get('analydir')
        self.queues = args.get('queues')
        self.sentieon_queues = args.get('sentieon_queues')

        self.sample_lists = sample_lists
        self.sample_infos = sample_infos

        self.ANALYSIS_POINTS = config.ANALYSIS_POINTS

        self.config = config

        self.ref = args.get('ref')
        self.REF='hg19' if self.ref == 'b37' else self.ref

        self.__dict__.update(dict(config.CONFIG.items('genome_' + self.ref)))
        self.__dict__.update(dict(config.CONFIG.items('software')))

        self.call_bed = self.config.CONFIG.get('genome_' + self.ref, 'refbed')
        if self.args['callTR']:
            self.call_bed = self.args['TR']

        self.annovar = self.config.CONFIG.get('genome_' + self.ref, 'annovar')
        self.moduledir = self.args['moduledir']
        self.newjob = self.args['newjob']

        self.mpileup_max_idepth = 9999  # INDEL检测允许的最大深度

        self.threads = {
            'sentieon_realign': 8,
            'sentieon_recal': 8,
            'sentieon_hc_call': 8,
            'sentieon_joint_call': 8,
            'sentieon_vqsr': 8,

            'annotate_merge': 4,
        }


    def start(self):

        print '>  mutation with {} ...'.format(self.mutation_soft)
        # ====== mutation calling ======

        # GATK Calling
        if self.mutation_soft == 'gatk':
            for sampleid, infos in self.args['sample_infos_waiting'].iteritems():
                self.gatk_recal(sampleid)
                self.gatk_hc_call(sampleid)
            self.gatk_consolidate()
            self.gatk_joint_call()
            self.gatk_vqsr()

        elif self.mutation_soft == 'sentieon':
            for sampleid, infos in self.args['sample_infos_waiting'].iteritems():
                sex = infos['sex']
                chrom_list = utils.get_chrom_list(self.ref, sex)
                self.sentieon_realign(sampleid)
                self.sentieon_recal(sampleid)
                self.sentieon_hc_call(sampleid, sex, chrom_list)
            self.sentieon_joint_call()
            self.sentieon_vqsr()

        for sampleid, infos in self.args['sample_infos_waiting'].iteritems():
            sex = infos['sex']
            chrom_list = utils.get_chrom_list(self.ref, sex)
            # print sampleid, sex, self.ref
            # print chrom_list
            if self.mutation_soft == 'samtools':
                self.samtools_mpileup(sampleid, chrom_list)
                self.bcftools_concat(sampleid, chrom_list)
                self.bcftools_filter(sampleid, sex)

            elif self.mutation_soft == 'mtoolbox':
                # print 'mutation with mtoolbox...'
                lanes = self.sample_lists[sampleid]['lanes']
                # print sampleid, lanes
                self.mtoolbox_call(sampleid, lanes)

        # merge
        if self.mutation_soft == 'samtools':
            self.bcftools_merge(self.sample_infos)

        # 注释和拆分
        if self.mutation_soft not in ('mtoolbox', ):
            for mtype in ('snp', 'indel'):
                self.annotate_merged_vcf(mtype)
                for sampleid, infos in self.sample_infos.iteritems():
                    sex = infos['sex']
                    self.extract_annotation_vcf(sampleid, mtype, sex)

        # 统计
        self.variation_summary()

# =============================================== samtools start =========================================================
    def samtools_mpileup(self, sampleid, chrom_list, bychr=True):

        # mutation calling by chrom
        if bychr:

            for chrom in chrom_list:
                cmd = '''
                    set -eo pipefail
                    echo samtools mpileup for {sampleid} start: `date "+%F %T"`\n
                    cd {analydir}/Mutation/{sampleid}.samtools\n
                    samtools-1.6 mpileup \\
                        -l {call_bed} \\
                        -r {chrom} \\
                        -q 1 -t DP,AD -C 50 -m 2 -F 0.002 -L {mpileup_max_idepth} -ugf \\
                        {reffasta} \\
                        {analydir}/Mapping/{sampleid}.{sampleid}/{sampleid}.final.bam |
                    bcftools-1.6 call \\
                        -g 5,10,30,50,100 \\
                        -m |
                    bcftools-1.6 view \\
                        -e \'type!="SNP" && type!="INDEL" && format/DP<5\' \\
                        -Ov -o {sampleid}.samtools.{chrom}.var.raw.vcf\n
                    echo samtools mpileup for {sampleid} done: `date "+%F %T"`
                '''.format(
                    sampleid=sampleid,
                    chrom=chrom,
                    **self.__dict__)

                shell_path = '{analydir}/Mutation/{sampleid}.samtools/samtools_mpileup_{chrom}_{sampleid}.sh'.format(
                    analydir=self.analydir,
                    sampleid=sampleid,
                    chrom=chrom
                )

                utils.write_shell(shell_path, cmd)

                # add job
                now_point = 'samtools_call'
                job_name = 'samtools_mpileup_{chrom}_{sampleid}'.format(chrom=chrom, sampleid=sampleid)
                utils.add_job(self.jobs, now_point, self.args['startpoint'],
                              self.ANALYSIS_POINTS, job_name, shell_path,
                              self.queues)

                # add order
                before_jobs = ['final_bam_{sampleid}'.format(sampleid=sampleid)]
                after_jobs = ['bcftools_concat_{sampleid}'.format(sampleid=sampleid)]
                utils.add_order(self.orders, job_name, before_jobs=before_jobs, after_jobs=after_jobs)

    def bcftools_concat(self, sampleid, chrom_list, bychr=True):

        # mutation calling by chrom
        if bychr:

            cmd = '''
                set -eo pipefail
                echo bcftools concat for {sampleid} start: `date "+%F %T"`
                
                cd {analydir}/Mutation/{sampleid}.samtools
                
                bcftools-1.6 concat \\
                    -Oz -o {sampleid}.samtools.var.raw.vcf.gz \\
                    {sampleid}.samtools.*.var.raw.vcf
                    
                bcftools-1.6 index -tf {sampleid}.samtools.var.raw.vcf.gz
                
                rm -f {sampleid}.samtools.*.var.raw.vcf
                
                echo bcftools concat for {sampleid} done: `date "+%F %T"`
            '''.format(
                analydir=self.analydir,
                sampleid=sampleid)

            shell_path = '{analydir}/Mutation/{sampleid}.samtools/bcftools_concat_{sampleid}.sh'.format(
                analydir=self.analydir,
                sampleid=sampleid
            )

            utils.write_shell(shell_path, cmd)

            # add job
            now_point = 'bcftools_concat'
            job_name = 'bcftools_concat_{sampleid}'.format(sampleid=sampleid)
            utils.add_job(self.jobs, now_point, self.args['startpoint'],
                          self.ANALYSIS_POINTS, job_name, shell_path,
                          self.queues)

    def bcftools_filter(self, sampleid, sex):

        # print '  bcftools filter...'
        # get all mutation for male in pseudoaotosomal chrY region:60001-2649520 and 59034050-59363566, other region just keep hom mutation.???
        filter_XY = ''
        if self.args['ref'] in ('b37', 'hg19', 'hg38') and sex in ('M',):
            filter_XY = ' filterXY - |'

        cmd = '''
            set -eo pipefail
            echo bcftools filter for {sampleid} start: `date "+%F %T"`

            cd {analydir}/Mutation/{sampleid}.{mutation_soft}
            
            # filter QUAL DP MQ, keep GT=0/0; then normalize
            bcftools-1.6 filter \\
                -i \'(%QUAL>20 && format/DP>4 && MQ>30) || (GT="0/0")\' \\
                {sampleid}.{mutation_soft}.var.raw.vcf.gz |{filter_XY}
            bcftools-1.6 norm \\
                -m -both |
            bcftools-1.6 norm \\
                -f {reffasta} \\
                -m +both \\
                -Oz -o {sampleid}.{mutation_soft}.var.flt.vcf.gz
                
            bcftools-1.6 index -tf {sampleid}.{mutation_soft}.var.flt.vcf.gz

            # separate snp
            bcftools-1.6 filter \\
                -i \'%TYPE="snp"\' \\
                -Oz -o {sampleid}.{mutation_soft}.snp.vcf.gz \\
                {sampleid}.{mutation_soft}.var.flt.vcf.gz
                
            bcftools-1.6 index -tf {sampleid}.{mutation_soft}.snp.vcf.gz
                          
            # separate indel
            bcftools-1.6 filter \\
                -i \'%TYPE="indel"\' \\
                -Oz -o {sampleid}.{mutation_soft}.indel.vcf.gz \\
                {sampleid}.{mutation_soft}.var.flt.vcf.gz
                
            bcftools-1.6 index -tf {sampleid}.{mutation_soft}.indel.vcf.gz
            
            echo bcftools filter for {sampleid} done: `date "+%F %T"`
        '''.format(**dict(self.__dict__, **locals()))

        shell_path = '{analydir}/Mutation/{sampleid}.{mutation_soft}/bcftools_filter_{sampleid}.sh'.format(
            **dict(self.__dict__, **locals()))

        utils.write_shell(shell_path, cmd)

        # add job
        now_point = 'bcftools_filter'
        job_name = 'bcftools_filter_{sampleid}'.format(sampleid=sampleid)
        utils.add_job(self.jobs, now_point, self.args['startpoint'],
                      self.ANALYSIS_POINTS, job_name, shell_path, self.queues)

        # add order
        before_jobs = ['bcftools_concat_{sampleid}'.format(sampleid=sampleid)]

        utils.add_order(self.orders, job_name, before_jobs=before_jobs)

    def bcftools_merge(self, sampleids):

        # print 'merge vcfs...'
        # preprae vcf list for merge
        cmd = '''
                set -eo pipefail
                echo merge gvcf start: `date "+%F %T"`
            
                cd {analydir}/Advance/{newjob}/Merged_vcf

                mkdir -p VCF
        '''

        if len(sampleids) == 1:
            sampleid = sampleids.keys()[0]
            cmd += '''
                ln -sf {analydir}/Mutation/{sampleid}.{mutation_soft}/{sampleid}.{mutation_soft}.var.flt.vcf.gz VCF/all.merged.vcf.gz
            '''
        else:
            vcf_list = '{analydir}/Advance/{newjob}/Merged_vcf/vcf.list'.format(**dict(self.__dict__, **locals()))
            with utils.safe_open(vcf_list, 'w') as out:
                for sampleid in sampleids:
                    vcf = '{analydir}/Mutation/{sampleid}.{mutation_soft}/{sampleid}.{mutation_soft}.var.flt.vcf.gz'.format(**dict(self.__dict__, **locals()))
                    out.write(vcf + '\n')

            cmd += '''
                bcftools-1.6 merge \\
                    -g {reffasta} \\
                    -l vcf.list \\
                    -Oz -o VCF/all.merged.vcf.gz
            '''

        cmd += '''
                bcftools-1.6 index -tf VCF/all.merged.vcf.gz
                    
                bcftools-1.6 view -i '%TYPE=="snp"' -Oz -o VCF/snp.merged.vcf.gz VCF/all.merged.vcf.gz
                
                bcftools-1.6 index -tf VCF/snp.merged.vcf.gz

                bcftools-1.6 view -i '%TYPE=="indel"' -Oz -o VCF/indel.merged.vcf.gz VCF/all.merged.vcf.gz

                bcftools-1.6 index -tf VCF/indel.merged.vcf.gz

                rm -f VCF/all.merged.vcf.gz

                echo merge gvcf done: `date "+%F %T"`
            '''

        cmd = cmd.format(**dict(self.__dict__, **locals()))

        shell_path = '{analydir}/Advance/{newjob}/Merged_vcf/bcftools_merge_vcf.sh'.format(**dict(self.__dict__, **locals()))

        utils.write_shell(shell_path, cmd)

        # # add job
        now_point = 'bcftools_merge'
        job_name = 'bcftools_merge'
        utils.add_job(self.jobs, now_point, self.args['startpoint'],
                      self.ANALYSIS_POINTS, job_name, shell_path, self.queues)

        # add order
        after_jobs = ['annotate_merged_snp', 'annotate_merged_indel']
        utils.add_order(self.orders, job_name, after_jobs=after_jobs)

        for sampleid in self.args['sample_infos_waiting']:
            before_jobs = ['bcftools_filter_{sampleid}'.format(sampleid=sampleid)]
            utils.add_order(self.orders, job_name, before_jobs=before_jobs)
# =============================================== samtools end =========================================================


# =============================================== gatk start =========================================================
    def gatk_recal(self, sampleid):

        # print '> bam recal with gatk...'
        cmd = '''
            set -eo pipefail
            echo gatk recal for {sampleid} start: `date "+%F %T"`

            cd {analydir}/Mutation/{sampleid}.gatk

            source {gatk4_dir}/source_this

            # step1
            gatk --java-options -Xmx10g BaseRecalibrator \\
                -R {reffasta} \\
                -I {analydir}/Mapping/{sampleid}.{sampleid}/{sampleid}.final.bam \\
                --known-sites {dbsnp} \\
                --known-sites {mills_indels} \\
                --TMP_DIR tmp \\
                -O {sampleid}.recal.table

            # step2
            gatk --java-options -Xmx10g ApplyBQSR \\
                -R {reffasta} \\
                -I {analydir}/Mapping/{sampleid}.{sampleid}/{sampleid}.final.bam \\
                -bqsr {sampleid}.recal.table \\
                --TMP_DIR tmp \\
                -O {sampleid}.recal.bam

            echo gatk recal for {sampleid} done: `date "+%F %T"`
        '''.format(
            sampleid=sampleid, **self.__dict__)

        shell_path = '{analydir}/Mutation/{sampleid}.gatk/gatk_recal_{sampleid}.sh'.format(
            analydir=self.analydir, sampleid=sampleid)

        utils.write_shell(shell_path, cmd)

        # add job
        now_point = 'gatk_recal'
        job_name = 'gatk_recal_{sampleid}'.format(sampleid=sampleid)
        utils.add_job(self.jobs, now_point, self.args['startpoint'],
                      self.ANALYSIS_POINTS, job_name, shell_path, self.queues)

        # add order
        before_jobs = ['final_bam_{sampleid}'.format(sampleid=sampleid)]
        utils.add_order(self.orders, job_name, before_jobs=before_jobs)

    def gatk_hc_call(self, sampleid):

        # print '> mutation call with gatk...'
        chroms = ','.join([str(i) for i in range(1, 23)] + ['X', 'Y'])

        cmd = '''
            set -eo pipefail
            echo gatk call variation for {sampleid} start: `date "+%F %T"`

            cd {analydir}/Mutation/{sampleid}.gatk

            source {gatk4_dir}/source_this

            gatk --java-options -Xmx10g HaplotypeCaller \\
                -R {reffasta} \\
                -I {sampleid}.recal.bam \\
                -O {sampleid}.gatk.g.vcf.gz \\
                -stand-call-conf 30 \\
                --native-pair-hmm-threads 8 \\
                -ERC GVCF

            # gvcf to vcf
            gatk --java-options -Xmx4g GenotypeGVCFs \\
                -R {reffasta} \\
                -V {sampleid}.gatk.g.vcf.gz \\
                -O {sampleid}.gatk.vcf.gz

            # separate
            bcftools-1.6 filter \\
                -i '%TYPE="SNP" && FORMAT/DP>4 && %QUAL>20 && MQ>30' \\
                -r {chroms} \\
                -Oz -o {sampleid}.gatk.snp.vcf.gz \\
                {sampleid}.gatk.vcf.gz

            bcftools-1.6 filter \\
                -i '%TYPE="INDEL" && FORMAT/DP>4 && %QUAL>20 && MQ>30' \\
                -r {chroms} \\
                -Oz -o {sampleid}.gatk.indel.vcf.gz \\
                {sampleid}.gatk.vcf.gz

            # norm
            bcftools-1.6 norm \\
                -f {reffasta} \\
                -m -both \\
                -Oz -o {sampleid}.gatk.snp_sn.vcf.gz \\
                {sampleid}.gatk.snp.vcf.gz

            bcftools-1.6 norm \\
                -f {reffasta} \\
                -m -both \\
                -Oz -o {sampleid}.gatk.indel_sn.vcf.gz \\
                {sampleid}.gatk.indel.vcf.gz

            # rm -f {sampleid}.recal.bam {sampleid}.recal.table

            echo gatk call variation for {sampleid} done: `date "+%F %T"`
        '''.format(**dict(self.__dict__, **locals()))

        shell_path = '{analydir}/Mutation/{sampleid}.gatk/gatk_hc_call_{sampleid}.sh'.format(
            **dict(self.__dict__, **locals()))
        utils.write_shell(shell_path, cmd)

        # add job
        now_point = 'gatk_hc_call'
        job_name = 'gatk_hc_call_{sampleid}'.format(sampleid=sampleid)
        utils.add_job(self.jobs, now_point, self.args['startpoint'],
                      self.ANALYSIS_POINTS, job_name, shell_path, self.queues)

        # add order
        before_jobs = ['gatk_recal_{sampleid}'.format(sampleid=sampleid)]
        after_jobs = ['gatk_consolidate']
        utils.add_order(self.orders, job_name, before_jobs=before_jobs, after_jobs=after_jobs)


    def gatk_consolidate(self):

        # print '> gatk consolidate gvcfs...'

        # with GenomicsDBImport
        # to be added ...

        # with CombineGVCFs
        input_gvcfs = ' \\\n                '.join(
            ['-V {analydir}/Mutation/{sampleid}.gatk/{sampleid}.gatk.g.vcf.gz'.format(
                **dict(self.__dict__, **locals())) for sampleid in self.sample_infos])

        cmd = '''
            set -eo pipefail
            echo gatk consolidate gvcfs start: `date "+%F %T"`

            cd {analydir}/Advance/{newjob}/Merged_vcf

            source {gatk4_dir}/source_this

            mkdir -p VCF

            gatk --java-options -Xmx10g CombineGVCFs \\
                -R {reffasta} \\
                -O VCF/all.merged.g.vcf.gz \\
                {input_gvcfs}

            echo gatk consolidate gvcfs done: `date "+%F %T"`
        '''.format(**dict(self.__dict__, **locals()))

        shell_path = '{analydir}/Advance/{newjob}/Merged_vcf/gatk_consolidate.sh'.format(**self.__dict__)

        utils.write_shell(shell_path, cmd)

        # add job
        now_point = 'gatk_consolidate'
        job_name = 'gatk_consolidate'
        utils.add_job(self.jobs, now_point, self.args['startpoint'],
                      self.ANALYSIS_POINTS, job_name, shell_path, self.queues)

        # add order

    def gatk_joint_call(self):

        # print '> gatk joint call...'
        cmd = '''
            set -eo pipefail
            echo gatk joint call start: `date "+%F %T"`

            cd {analydir}/Advance/{newjob}/Merged_vcf

            source {gatk4_dir}/source_this

            gatk --java-options -Xmx4g GenotypeGVCFs \\
                -R {reffasta} \\
                -V VCF/all.merged.g.vcf.gz \\
                -O VCF/all.merged.vcf.gz

            echo gatk joint call done: `date "+%F %T"`
        '''.format(**dict(self.__dict__, **locals()))

        shell_path = '{analydir}/Advance/{newjob}/Merged_vcf/gatk_joint_call.sh'.format(
            **dict(self.__dict__, **locals()))

        utils.write_shell(shell_path, cmd)

        # add job
        now_point = 'gatk_joint_call'
        job_name = 'gatk_joint_call'
        utils.add_job(self.jobs, now_point, self.args['startpoint'],
                      self.ANALYSIS_POINTS, job_name, shell_path, self.queues)

        # add order
        before_jobs = ['gatk_consolidate']
        utils.add_order(self.orders, job_name, before_jobs=before_jobs)

    def gatk_vqsr(self):

        # print '> gatk vqsr...'
        chroms = ','.join([str(i) for i in range(1, 23)] + ['X', 'Y'])

        cmd = '''
            set -eo pipefail
            echo gatk vqsr start: `date "+%F %T"`

            cd {analydir}/Advance/{newjob}/Merged_vcf

            source {gatk4_dir}/source_this

            # 1 SNP VQSR
            gatk --java-options -Xmx4g VariantRecalibrator \\
                -R {reffasta} \\
                -V VCF/all.merged.vcf.gz \\
                -mode SNP \\
                --resource hapmap,known=false,training=true,truth=true,prior=15.0:{hapmap} \\
                --resource omni,known=false,training=true,truth=false,prior=12.0:{1000g_omni} \\
                --resource 1000G,known=false,training=true,truth=false,prior=10.0:{1000g_phase1} \\
                --resource dbsnp,known=true,training=false,truth=false,prior=2.0:{dbsnp} \\
                -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \\
                -O VCF/all.merged.snp.recal \\
                --tranches-file VCF/all.merged.snp.tranches \\
                --rscript-file VCF/all.merged.snp.plots.R

            # 2 SNP Apply VQSR
            gatk --java-options -Xmx4g ApplyVQSR \\
                -R {reffasta} \\
                -V VCF/all.merged.vcf.gz \\
                -O VCF/all.merged.snp.recal.vcf.gz \\
                -mode SNP \\
                --recal-file VCF/all.merged.snp.recal \\
                --tranches-file VCF/all.merged.snp.tranches \\
                --truth-sensitivity-filter-level 99.0

            # 3 INDEL VQSR
            gatk --java-options -Xmx4g VariantRecalibrator \\
                -R {reffasta} \\
                -V VCF/all.merged.snp.recal.vcf.gz \\
                -mode INDEL \\
                --resource mills,known=false,training=true,truth=true,prior=12.0:{mills_indels} \\
                --resource dbsnp,known=true,training=false,truth=false,prior=2.0:{dbsnp} \\
                -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \\
                -O VCF/all.merged.snp.indel.recal \\
                --tranches-file VCF/all.merged.snp.indel.tranches \\
                --rscript-file VCF/all.merged.snp.indel.plots.R

            # 4 INDEL Apply VQSR
            gatk --java-options -Xmx4g ApplyVQSR \\
                -R {reffasta} \\
                -V VCF/all.merged.snp.vcf.gz \\
                -O VCF/all.merged.snp.indel.recal.vcf.gz \\
                -mode INDEL \\
                --recal-file VCF/all.merged.snp.indel.recal \\
                --tranches-file VCF/all.merged.snp.indel.tranches \\
                --truth-sensitivity-filter-level 99.0

            # Filter and Separate
            bcftools-1.6 filter \\
                -i '%FILTER=="PASS" && %TYPE="SNP"' \\
                -r {chroms} \\
                -Oz -o VCF/snp.merged.vcf.gz \\
                VCF/all.merged.snp.indel.recal.vcf.gz 

            bcftools-1.6 filter \\
                -i '%FILTER=="PASS" && %TYPE="INDEL"' \\
                -r {chroms} \\
                -Oz -o VCF/indel.merged.vcf.gz \\
                VCF/all.merged.snp.indel.recal.vcf.gz  

            echo gatk joint call done: `date "+%F %T"`
        '''.format(**dict(self.__dict__, **locals()))

        shell_path = '{analydir}/Advance/{newjob}/Merged_vcf/gatk_vqsr.sh'.format(
            **dict(self.__dict__, **locals()))

        utils.write_shell(shell_path, cmd)

        # add job
        now_point = 'gatk_vqsr'
        job_name = 'gatk_vqsr'
        utils.add_job(self.jobs, now_point, self.args['startpoint'],
                      self.ANALYSIS_POINTS, job_name, shell_path, self.queues)

        # add order
        before_jobs = ['gatk_joint_call']
        after_jobs = ['annotate_merged_snp', 'annotate_merged_indel']
        utils.add_order(self.orders, job_name, before_jobs=before_jobs, after_jobs=after_jobs)

    # def gatk_filter(self, sampleid):

    #     # print '> gatk filter...'
    #     cmd = '''
    #         set -eo pipefail
    #         echo gatk filter for {sampleid} start: `date "+%F %T"`\n
    #         cd {analydir}/Mutation/{sampleid}.gatk

    #         # ========================== Filter SNP ==========================
    #         # 1 Extract SNP Variants
    #         java1.8.0 -Xmx10g -jar {gatkdir}/GenomeAnalysisTK.jar \\
    #             --no_cmdline_in_header \\
    #             -T SelectVariants \\
    #             -R {reffasta} \\
    #             -selectType SNP \\
    #             -V {sampleid}.gatk.hc.vcf.gz \\
    #             -o {sampleid}.gatk.raw.snp.vcf

    #         # 2 Apply SNP Filter
    #         java1.8.0 -Xmx10g -jar {gatkdir}/GenomeAnalysisTK.jar \\
    #             --no_cmdline_in_header \\
    #             -T VariantFiltration \\
    #             -R {reffasta} \\
    #             --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \\
    #             --filterName "my_snp_filter" \\
    #             -V {sampleid}.gatk.raw.snp.vcf \\
    #             -o {sampleid}.gatk.filter.snp.vcf

    #         # 3 Filter Result
    #         bcftools-1.6 filter \\
    #             -i '%FILTER=="PASS" && DP>4' \\
    #             -Oz -o {sampleid}.gatk.snp.vcf.gz \\
    #             {sampleid}.gatk.filter.snp.vcf

    #         # 4 Norm
    #         bcftools-1.6 norm \\
    #             -f {reffasta} \\
    #             -m -both \\
    #             -Oz -o {sampleid}.gatk.snp_sn.vcf.gz \\
    #             {sampleid}.gatk.snp.vcf.gz

    #         # ========================== Filter INDEL ==========================
    #         # 1 Extract INDEL Variants
    #         java1.8.0 -Xmx10g -jar {gatkdir}/GenomeAnalysisTK.jar \\
    #             --no_cmdline_in_header \\
    #             -T SelectVariants \\
    #             -R {reffasta} \\
    #             -selectType INDEL \\
    #             -V {sampleid}.gatk.hc.vcf.gz \\
    #             -o {sampleid}.gatk.raw.indel.vcf

    #         # 2 Apply INDEL Filter
    #         java1.8.0 -Xmx10g -jar {gatkdir}/GenomeAnalysisTK.jar \\
    #             --no_cmdline_in_header \\
    #             -T VariantFiltration \\
    #             -R {reffasta} \\
    #             --filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \\
    #             --filterName "my_indel_filter" \\
    #             -V {sampleid}.gatk.raw.indel.vcf \\
    #             -o {sampleid}.gatk.filter.indel.vcf

    #         # 3 Filter Result
    #         bcftools-1.6 filter \\
    #             -i '%FILTER=="PASS" && DP>4' \\
    #             -Oz -o {sampleid}.gatk.indel.vcf.gz \\
    #             {sampleid}.gatk.filter.indel.vcf

    #         # 4 Norm
    #         bcftools-1.6 norm \\
    #             -f {reffasta} \\
    #             -m -both \\
    #             -Oz -o {sampleid}.gatk.indel_sn.vcf.gz \\
    #             {sampleid}.gatk.indel.vcf.gz

    #         # Make index
    #         bcftools-1.6 index -tf {sampleid}.gatk.snp.vcf.gz
    #         bcftools-1.6 index -tf {sampleid}.gatk.indel.vcf.gz

    #         # Remove temp file
    #         rm -f *.gatk.hc.g.vcf* *gatk.raw.*.vcf* *gatk.filter.*.vcf*\n
    #         echo gatk filter for {sampleid} done: `date "+%F %T"`
    #     '''.format(
    #         sampleid=sampleid, **self.__dict__)

    #     shell_path = '{analydir}/Mutation/{sampleid}.gatk/gatk_filter_{sampleid}.sh'.format(
    #         analydir=self.analydir, sampleid=sampleid)

    #     utils.write_shell(shell_path, cmd)

    #     # add job
    #     now_point = 'gatk_filter'
    #     job_name = 'gatk_filter_{sampleid}'.format(sampleid=sampleid)
    #     utils.add_job(self.jobs, now_point, self.args['startpoint'],
    #                   self.ANALYSIS_POINTS, job_name, shell_path, self.queues)

    #     # add order
    #     before_jobs = ['gatk_gvcf2vcf_{sampleid}'.format(sampleid=sampleid)]
    #     utils.add_order(self.orders, job_name, before_jobs=before_jobs)

# =============================================== gatk end =========================================================


# =============================================== sentieon start =========================================================
# ========================================== just for human b37 now ======================================================
    def sentieon_realign(self, sampleid):

        # print '  sentieon realign...'
        # write shell
        realign_threads = self.threads['sentieon_realign']

        cmd = '''
            set -eo pipefail
            echo sentieon realign for {sampleid} start: `date "+%F %T"`

            cd {analydir}/Mutation/{sampleid}.sentieon

            sentieon driver \\
                -t {realign_threads} \\
                -i {analydir}/Mapping/{sampleid}.{sampleid}/{sampleid}.final.bam \\
                -r {reffasta} \\
                --algo Realigner \\
                -k {mills_indels} \\
                {sampleid}.realn.bam

            echo sentieon realign for {sampleid} done: `date "+%F %T"`
        '''.format(**dict(self.__dict__, **locals()))

        shell_path = '{analydir}/Mutation/{sampleid}.sentieon/sentieon_realign_{sampleid}.sh'.format(
            **dict(self.__dict__, **locals()))

        utils.write_shell(shell_path, cmd)

        # add job
        now_point = 'sentieon_realign'
        job_name = 'sentieon_realign_{sampleid}'.format(sampleid=sampleid)
        utils.add_job(
            self.jobs,
            now_point,
            self.args['startpoint'],
            self.ANALYSIS_POINTS,
            job_name,
            shell_path,
            self.sentieon_queues,
            threads=realign_threads)
        # add order
        before_jobs = ['final_bam_{sampleid}'.format(sampleid=sampleid)]
        utils.add_order(self.orders, job_name, before_jobs=before_jobs)

    def sentieon_recal(self, sampleid):

        # print '> sentieon recal...'
        recal_threads = self.threads['sentieon_recal']

        cmd = '''
            set -eo pipefail
            echo "sentieon recal(BQSR) for {sampleid} start: `date "+%F %T"`"
            
            cd {analydir}/Mutation/{sampleid}.sentieon
            
            sentieon driver \\
                -t {recal_threads} \\
                -i {sampleid}.realn.bam \\
                -r {reffasta} \\
                --algo QualCal \\
                -k {mills_indels} \\
                -k {dbsnp} \\
                {sampleid}.recal.table
                
            echo "sentieon recal(BQSR) for {sampleid} done: `date "+%F %T"`"
        '''.format(**dict(self.__dict__, **locals()))

        shell_path = '{analydir}/Mutation/{sampleid}.sentieon/sentieon_recal_{sampleid}.sh'.format(
            **dict(self.__dict__, **locals()))

        utils.write_shell(shell_path, cmd)

        # add job
        now_point = 'sentieon_recal'
        job_name = 'sentieon_recal_{sampleid}'.format(sampleid=sampleid)
        utils.add_job(
            self.jobs,
            now_point,
            self.args['startpoint'],
            self.ANALYSIS_POINTS,
            job_name,
            shell_path,
            self.sentieon_queues,
            threads=recal_threads)

        # add order
        before_jobs = ['sentieon_realign_{sampleid}'.format(sampleid=sampleid)]
        utils.add_order(self.orders, job_name, before_jobs=before_jobs)

    def sentieon_hc_call(self, sampleid, sex, chrom_list):

        # print '> mutation call with sentieon...'
        hc_call_threads = self.threads['sentieon_hc_call']

        chroms = ','.join(chrom_list)

        cmd = '''
            set -eo pipefail
            echo sentieon call variation for {sampleid} start: `date "+%F %T"`
            
            cd {analydir}/Mutation/{sampleid}.sentieon
            
            sentieon driver \\
                -t {hc_call_threads} \\
                -r {reffasta} \\
                -i {sampleid}.realn.bam \\
                -q {sampleid}.recal.table \\
                --algo Haplotyper \\
                -d {dbsnp} \\
                --emit_conf=10  \\
                --call_conf=10 \\
                --emit_mode gvcf \\
                {sampleid}.sentieon.hc.g.vcf.gz
        '''

        if sex in ('F', 'Female'):
            cmd += '''
            bcftools-1.6 view -t ^Y \\
                -Oz -o {sampleid}.sentieon.hc.g.vcf.gz.tmp \\
                {sampleid}.sentieon.hc.g.vcf.gz

            mv {sampleid}.sentieon.hc.g.vcf.gz.tmp {sampleid}.sentieon.hc.g.vcf.gz

            bcftools-1.6 index -tf {sampleid}.sentieon.hc.g.vcf.gz
            '''

        if self.args['ref'] in ('b37', 'hg19', 'hg38') and sex in ('M', ):
            cmd += '''
            filterXY {sampleid}.sentieon.hc.g.vcf.gz -o {sampleid}.sentieon.hc.g.vcf

            bgzip -f {sampleid}.sentieon.hc.g.vcf

            bcftools-1.6 index -tf {sampleid}.sentieon.hc.g.vcf.gz
            '''

        cmd += '''
            rm -f {sampleid}.recal.bam {sampleid}.recal.table

            echo sentieon call variation for {sampleid} done: `date "+%F %T"`
        '''

        cmd = cmd.format(**dict(self.__dict__, **locals()))

        shell_path = '{analydir}/Mutation/{sampleid}.sentieon/sentieon_hc_call_{sampleid}.sh'.format(
            **dict(self.__dict__, **locals()))

        utils.write_shell(shell_path, cmd)

        # add job
        now_point = 'sentieon_hc_call'
        job_name = 'sentieon_hc_call_{sampleid}'.format(sampleid=sampleid)
        utils.add_job(
            self.jobs,
            now_point,
            self.args['startpoint'],
            self.ANALYSIS_POINTS,
            job_name,
            shell_path,
            self.sentieon_queues,
            threads=hc_call_threads)

        # add order
        before_jobs = ['sentieon_recal_{sampleid}'.format(sampleid=sampleid)]
        after_jobs = ['sentieon_joint_call']
        utils.add_order(self.orders, job_name, before_jobs=before_jobs, after_jobs=after_jobs)

    def sentieon_joint_call(self):

        # print '> sentieon gvcftyper...'
        joint_threads = self.threads['sentieon_joint_call']

        input_gvcfs = ' \\\n                '.join(
            ['-v {analydir}/Mutation/{sampleid}.sentieon/{sampleid}.sentieon.hc.g.vcf.gz'.format(
                **dict(self.__dict__, **locals())) for sampleid in self.sample_infos])

        cmd = '''
            set -eo pipefail
            echo sentieon joint call start: `date "+%F %T"`
            
            cd {analydir}/Advance/{newjob}/Merged_vcf

            mkdir -p VCF
            
            sentieon driver \\
                -t {joint_threads} \\
                -r {reffasta} \\
                --algo GVCFtyper \\
                -d {dbsnp} \\
                {input_gvcfs} \\
                VCF/all.merged.vcf.gz
                
            echo sentieon joint call done: `date "+%F %T"`
        '''.format(**dict(self.__dict__, **locals()))

        shell_path = '{analydir}/Advance/{newjob}/Merged_vcf/sentieon_joint_call.sh'.format(
            **self.__dict__)

        utils.write_shell(shell_path, cmd)

        # add job
        now_point = 'sentieon_joint_call'
        job_name = 'sentieon_joint_call'
        utils.add_job(
            self.jobs,
            now_point,
            self.args['startpoint'],
            self.ANALYSIS_POINTS,
            job_name,
            shell_path,
            self.sentieon_queues,
            threads=joint_threads)

    def sentieon_vqsr(self):

        # print sampleid, sex
        vqsr_threads = self.threads['sentieon_vqsr']

        chroms = ','.join([str(i) for i in range(1, 23)] + ['X', 'Y'])

        # print '> sentieon gvcftyper...'
        cmd = '''
            set -eo pipefail
            echo sentieon vqsr start: `date "+%F %T"`

            cd {analydir}/Advance/{newjob}/Merged_vcf

            # Step1: SNP VQSR
            sentieon driver \\
                -t {vqsr_threads} \\
                -r {reffasta} \\
                --algo VarCal \\
                --var_type SNP \\
                --resource {1000g_phase1} \\
                --resource_param 1000G,known=false,training=true,truth=false,prior=10.0 \\
                --resource {1000g_omni} \\
                --resource_param omni,known=false,training=true,truth=true,prior=12.0 \\
                --resource {dbsnp} \\
                --resource_param dbsnp,known=true,training=false,truth=false,prior=2.0 \\
                --resource {hapmap} \\
                --resource_param hapmap,known=false,training=true,truth=true,prior=15.0 \\
                --annotation QD --annotation MQ --annotation MQRankSum \\
                --annotation ReadPosRankSum --annotation FS \\
                --tranches_file VCF/all.merged.snp.tranches \\
                -v VCF/all.merged.vcf.gz \\
                VCF/all.merged.snp.recal 

            # Step2: Apply SNP VQSR
            sentieon driver \\
                -t {vqsr_threads} \\
                -r {reffasta} \\
                --algo ApplyVarCal \\
                --var_type SNP \\
                --recal VCF/all.merged.snp.recal \\
                --tranches_file VCF/all.merged.snp.tranches \\
                -v VCF/all.merged.vcf.gz \\
                VCF/all.merged.snp.recal.vcf.gz
                
            # Step3: INDEL VQSR
            sentieon driver \\
                -t {vqsr_threads} \\
                -r {reffasta} \\
                --algo VarCal \\
                --var_type INDEL \\
                --resource {dbsnp} \\
                --resource_param dbsnp,known=true,training=false,truth=false,prior=2.0 \\
                --resource {mills_indels} \\
                --resource_param Mills,known=false,training=true,truth=true,prior=12.0 \\
                --annotation QD --annotation MQ \\
                --annotation ReadPosRankSum --annotation FS \\
                --tranches_file VCF/all.merged.snp.indel.tranches \\
                -v VCF/all.merged.snp.recal.vcf.gz \\
                VCF/all.merged.snp.indel.recal
                
            # Step4: Apply INDEL VQSR
            sentieon driver \\
                -t {vqsr_threads} \\
                -r {reffasta} \\
                --algo ApplyVarCal \\
                --var_type INDEL \\
                --recal VCF/all.merged.snp.indel.recal \\
                --tranches_file VCF/all.merged.snp.indel.tranches \\
                -v VCF/all.merged.snp.recal.vcf.gz \\
                VCF/all.merged.snp.indel.recal.vcf.gz
                
            # Filter and Separate
            bcftools-1.6 filter \\
                -i '%FILTER=="PASS" && %TYPE="SNP" && (%QUAL>20 && format/DP>4 && MQ>30)' \\
                -r {chroms} \\
                -Oz -o VCF/snp.merged.vcf.gz \\
                VCF/all.merged.snp.indel.recal.vcf.gz 

            bcftools-1.6 index -tf VCF/snp.merged.vcf.gz

            bcftools-1.6 filter \\
                -i '%FILTER=="PASS" && %TYPE="INDEL" && (%QUAL>20 && format/DP>4 && MQ>30)' \\
                -r {chroms} \\
                -Oz -o VCF/indel.merged.vcf.gz \\
                VCF/all.merged.snp.indel.recal.vcf.gz 

            bcftools-1.6 index -tf VCF/indel.merged.vcf.gz
        '''

        if self.args['MT']:
            cmd += '''
            bcftools-1.6 filter \\
                -i '%FILTER=="PASS" && %TYPE="SNP" && (%QUAL>20 && format/DP>4 && MQ>30)' \\
                -r MT \\
                -Oz -o VCF/snp.MT.merged.vcf.gz \\
                VCF/all.merged.snp.indel.recal.vcf.gz

            bcftools-1.6 filter \\
                -i '%FILTER=="PASS" && %TYPE="INDEL" && (%QUAL>20 && format/DP>4 && MQ>30)' \\
                -r MT \\
                -Oz -o VCF/indel.MT.merged.vcf.gz \\
                VCF/all.merged.snp.indel.recal.vcf.gz
            '''

        cmd += '''        
            echo sentieon vqsr done: `date "+%F %T"`
        '''

        cmd = cmd.format(**dict(self.__dict__, **locals()))

        shell_path = '{analydir}/Advance/{newjob}/Merged_vcf/sentieon_vqsr.sh'.format(
            **dict(self.__dict__, **locals()))

        utils.write_shell(shell_path, cmd)

        # add job
        now_point = 'sentieon_vqsr'
        job_name = 'sentieon_vqsr'
        utils.add_job(
            self.jobs,
            now_point,
            self.args['startpoint'],
            self.ANALYSIS_POINTS,
            job_name,
            shell_path,
            self.sentieon_queues,
            threads=vqsr_threads)

        # add order
        before_jobs = ['sentieon_joint_call']
        after_jobs = ['annotate_merged_snp', 'annotate_merged_indel']
        utils.add_order(self.orders, job_name, before_jobs=before_jobs, after_jobs=after_jobs)

# =============================================== sentieon end =========================================================

# =============================================== mtoolbox start =========================================================
    def mtoolbox_call(self, sampleid, lanes):

        read1_list, merge_read1, merge_read2, merge_read1_gz, merge_read2_gz = utils.get_merge_read(
            self.args['analydir'], sampleid, lanes)

        config_file = '{analydir}/Mutation/{sampleid}.mtoolbox/{sampleid}.mt.config'.format(**dict(self.__dict__, **locals()))

        config_txt = '''
            mtdb_fasta=chrM.fa
            hg19_fasta=hg19RCRS.fa
            mtdb=chrM
            humandb=hg19RCRS
            input_path={analydir}/Mutation/{sampleid}.mtoolbox
            list={sampleid}.reads.list
            input_type=fastq
            ref=RCRS
            UseMarkDuplicates=true
            UseIndelRealigner=true
            MitoExtraction=true
        '''.format(**dict(self.__dict__, **locals()))

        utils.write_shell(config_file, config_txt)

        list_file = '{analydir}/Mutation/{sampleid}.mtoolbox/{sampleid}.reads.list'.format(
            **dict(self.__dict__, **locals()))

        list_txt = configtxt = '''
            {sampleid}.R1.fastq
            {sampleid}.R2.fastq
        '''.format(sampleid=sampleid)

        utils.write_shell(list_file, list_txt)

        cmd = '''
            set -eo pipefail
            echo mtoolbox call for {sampleid} start: `date "+%F %T"`

            cd {analydir}/Mutation/{sampleid}.mtoolbox

            if [ -f {read1_1} ];then
                {merge_read1}
                {merge_read2}
            else
                {merge_read1_gz}
                {merge_read2_gz}
            fi

            # calling 
            {mtoolbox_dir}/MToolBox.sh \\
                -i {sampleid}.mt.config \\
                -m "-t 20" \\
                -a "-t 10 -z 0.6"

            # annotation
            python {moduledir}/Varition/MT/anno2xls.py \\
                --anno_dir OUT_{sampleid} \\
                --outpre {sampleid} \\
                --mtref {mtoolbox_dir}/data/chrRCRS.fa

            echo mtoolbox call for {sampleid} done: `date "+%F %T"`
        '''.format(
            read1_1=read1_list[0], **dict(self.__dict__, **locals()))

        shell_path = '{analydir}/Mutation/{sampleid}.mtoolbox/mtoolbox_call_{sampleid}.sh'.format(
            sampleid=sampleid, **self.__dict__)

        utils.write_shell(shell_path, cmd)

        # add job
        now_point = 'mtoolbox_call'
        job_name = 'mtoolbox_call_{sampleid}'.format(sampleid=sampleid)
        utils.add_job(self.jobs, now_point, self.args['startpoint'],
                      self.ANALYSIS_POINTS, job_name, shell_path, self.queues)

        # add order
        before_jobs = []
        for lane in lanes:
            before_jobs.append('qc_{sampleid}_{novoid}_{flowcell}_L{lane}'.format(sampleid=sampleid, **lane))
        after_jobs = ['variation_summary']
        utils.add_order(self.orders, job_name, before_jobs=before_jobs, after_jobs=after_jobs)

# =============================================== mtoolbox end =========================================================

# =============================================== annotation start =========================================================
# def annotate_vcf(self, sampleid, mut_type):

#     # print self.annovar
#     # print self.args['moduledir']
#     # print '  annotate vcf...'
#     cmd = '''
#         set -eo pipefail
#         echo annotate {mut_type} vcf for {sampleid} start: `date "+%F %T"`

#         cd {analydir}/Mutation/{sampleid}.{mutation_soft}

#         # ===== 1 databases annotation with annovar =====
#         sh {annovar} \\
#             -b {MUT_TYPE} \\
#             -p 5 \\
#             -r {REF} \\
#             {sampleid}.{mutation_soft}.{mut_type}_sn.vcf.gz \\
#             {sampleid}

#         # ===== 2 add PubmedID and reformat HGMD =====
#         python {moduledir}/AddOMIM_HGMD/AddHGMD_OMIM_Priority_pipe4.6.py \\
#             {sampleid}.{mutation_soft}.{mut_type}_sn.annovar.{REF}_multianno.xls \\
#             {sampleid}.{mutation_soft}.{mut_type}.annovar.{REF}_multianno_mid.xls \\
#             {REF}

#         # ===== 3 add HPA annotation =====
#         python {moduledir}/HPA.v15/annotatExpression_for_multiannofile.py \\
#             -i {sampleid}.{mutation_soft}.{mut_type}.annovar.{REF}_multianno_mid.xls \\
#             -o {sampleid}.{mutation_soft}.{mut_type}.annovar.{REF}_multianno.xls

#         gzip -f {sampleid}.{mutation_soft}.{mut_type}.annovar.{REF}_multianno.xls

#         rm -f *{mutation_soft}.{mut_type}*multianno.xls.bak *{mutation_soft}.{mut_type}*multianno_mid.xls

#         rm -f {sampleid}.{mutation_soft}.{mut_type}_sn.vcf.gz

#         echo annotate {mut_type} vcf for {sampleid} done: `date "+%F %T"`
#     '''.format(
#         sampleid=sampleid,
#         mut_type=mut_type,
#         MUT_TYPE=mut_type.upper(),
#         **self.__dict__)

#     shell_path = '{analydir}/Mutation/{sampleid}.{mutation_soft}/annotate_{mut_type}_{mutation_soft}_{sampleid}.sh'.format(
#         sampleid=sampleid, mut_type=mut_type, **self.__dict__)

#     utils.write_shell(shell_path, cmd)

#     # add job
#     now_point = 'annotate_vcf'
#     job_name = 'annotate_{mut_type}_{sampleid}'.format(sampleid=sampleid, mut_type=mut_type)
#     utils.add_job(self.jobs, now_point, self.args['startpoint'],
#                   self.ANALYSIS_POINTS, job_name, shell_path, self.queues)

#     # add order
#     filter_soft = 'gatk' if self.mutation_soft == 'gatk' else 'bcftools'
#     before_jobs = ['{filter_soft}_filter_{sampleid}'.format(**locals())]
#     after_jobs = ['variation_summary']
#     utils.add_order(self.orders, job_name, before_jobs=before_jobs, after_jobs=after_jobs)

    def annotate_merged_vcf(self, mut_type):

        # print self.annovar
        # print self.args['moduledir']
        # print '  annotate merged vcf...'
        annotate_merge_threads = self.threads['annotate_merge']

        sample_list = ','.join(self.sample_infos)
        MUT_TYPE = mut_type.upper()

        cmd = '''
            set -eo pipefail
            echo annotate merged {mut_type} vcf start: `date "+%F %T"`

            cd {analydir}/Advance/{newjob}/Merged_vcf
    
            # ===== 1 databases annotation with annovar =====
            sh {annovar} \\
                -b {MUT_TYPE} \\
                -p {annotate_merge_threads} \\
                -r {ref} \\
                VCF/{mut_type}.merged.vcf.gz \\
                {sample_list}

            # ===== 2 add PubmedID and reformat HGMD =====
            python {moduledir}/AddOMIM_HGMD/AddHGMD_OMIM_Priority_pipe4.6.py \\
                VCF/{mut_type}.merged_sn.annovar.{REF}_multianno.xls \\
                VCF/{mut_type}.merged.annovar.{REF}_multianno_mid.xls \\
                {REF}

            # ===== 3 add HPA annotation =====
            python {moduledir}/HPA.v15/annotatExpression_for_multiannofile.py \\
                -i VCF/{mut_type}.merged.annovar.{REF}_multianno_mid.xls \\
                -o VCF/{mut_type}.merged.annovar.{REF}_multianno.xls

            gzip -f VCF/{mut_type}.merged.annovar.{REF}_multianno.xls

            rm -f VCF/{mut_type}*_sn*gz VCF/{mut_type}*multianno.xls.bak VCF/{mut_type}*multianno_mid.xls
        '''

        cmd += '''
            echo annotate merged {mut_type} vcf done: `date "+%F %T"`
        '''

        cmd = cmd.format(**dict(self.__dict__, **locals()))

        shell_path = '{analydir}/Advance/{newjob}/Merged_vcf/annotate_merged_{mut_type}_vcf.sh'.format(**dict(self.__dict__, **locals()))

        utils.write_shell(shell_path, cmd)

        # add job
        now_point = 'annotate_merged_vcf'
        job_name = 'annotate_merged_{mut_type}'.format(mut_type=mut_type)
        utils.add_job(
            self.jobs,
            now_point,
            self.args['startpoint'],
            self.ANALYSIS_POINTS,
            job_name,
            shell_path,
            self.queues,
            threads=annotate_merge_threads)

    def extract_annotation_vcf(self, sampleid, mut_type, sex):

        # print '  extract annotation from merged...'
        chroms = ','.join(map(str, range(1, 23)) + ['X'])
        if sex.lower() in ('m', 'male'):
            chroms += ',Y'

        cmd = '''
            set -eo pipefail
            echo extract annotation for {sampleid} start: `date "+%F %T"`

            cd {analydir}/Mutation/{sampleid}.{mutation_soft}

            if [ '{mutation_soft}' == 'gatk' -o '{mutation_soft}' == 'sentieon' ];then
                bcftools-1.6 view \\
                    -s {sampleid} \\
                    -r {chroms} \\
                    {analydir}/Advance/{newjob}/Merged_vcf/VCF/{mut_type}.merged.vcf.gz |
                awk '$1 ~ /^#/ || $NF !~ /\.\/\./' > {sampleid}.{mutation_soft}.{mut_type}.vcf
                bgzip -f {sampleid}.{mutation_soft}.{mut_type}.vcf
            fi

            # norm
            bcftools-1.6 norm \\
                -f {reffasta} \\
                -m -both \\
                -Oz -o {sampleid}.{mutation_soft}.{mut_type}_sn.vcf.gz \\
                {sampleid}.{mutation_soft}.{mut_type}.vcf.gz

            python {moduledir}/Varition/Ann/GetAnnoFromMerge_pipe4.6.py \\
                -M {analydir}/Advance/{newjob}/Merged_vcf/VCF/{mut_type}.merged.annovar.{REF}_multianno.xls.gz \\
                -nv {sampleid}.{mutation_soft}.{mut_type}_sn.vcf.gz \\
                -rv {sampleid}.{mutation_soft}.{mut_type}.vcf.gz \\
                -S {sampleid} \\
                -T {mut_type} \\
                -ref {REF} \\
                -O {sampleid}.{mutation_soft}.{mut_type}.annovar.{REF}_multianno.xls

            gzip -f {sampleid}.{mutation_soft}.{mut_type}.annovar.{REF}_multianno.xls

            # rm -f {sampleid}.{mutation_soft}.{mut_type}_sn.vcf.gz
        '''

        if self.args['MT']:

            if self.args['ref'] == 'b37':
                buildver = 'GRCh37_MT'
            elif self.args['ref'] == 'hg19':
                buildver = 'hg19_MT'
            else:
                print 'only hg19 or b37 can do MT annotation!'
                exit(1)

            cmd += '''
            if [ '{mutation_soft}' == 'gatk' -o '{mutation_soft}' == 'sentieon' ];then
                bcftools-1.6 view \\
                    -s {sampleid} \\
                    {analydir}/Advance/{newjob}/Merged_vcf/VCF/{mut_type}.MT.merged.vcf.gz |
                awk '$1 ~ /^#/ || $NF !~ /\.\/\./' > {sampleid}.{mutation_soft}.{mut_type}.MT.vcf
                bgzip -f {sampleid}.{mutation_soft}.{mut_type}.MT.vcf
            fi

            perl {annovar_dir}/table_annovar.pl \\
                {sampleid}.{mutation_soft}.{mut_type}.MT.vcf.gz \\
                {humandb_dir} \\
                --buildver {buildver} \\
                -otherinfo \\
                -nastring . \\
                -protocol ensGene \\
                --operation g \\
                --vcfinput \\
                --outfile {sampleid}.{mutation_soft}.{mut_type}.MT.annovar

            perl {annovar_dir}/reformat_annovar.pl \\
                {sampleid}.{mutation_soft}.{mut_type}.MT.annovar.{buildver}_multianno.txt \\
                -v vcf4 \\
                -id {sampleid} \\
                > {sampleid}.{mutation_soft}.{mut_type}.MT.annovar.{buildver}_multianno.xls

            gzip -f {sampleid}.{mutation_soft}.{mut_type}.MT.annovar.{buildver}_multianno.xls

            rm -f {sampleid}.{mutation_soft}.{mut_type}.MT.annovar.ensGene* {sampleid}.{mutation_soft}.{mut_type}.MT.annovar.avinput
            rm -f {sampleid}.{mutation_soft}.{mut_type}.MT.annovar.*.txt {sampleid}.{mutation_soft}.{mut_type}.MT.annovar.*.vcf
            '''

        cmd += '''
            echo extract annotation for {sampleid} done: `date "+%F %T"`
        '''

        cmd = cmd.format(**dict(self.__dict__, **locals()))

        shell_path = '{analydir}/Mutation/{sampleid}.{mutation_soft}/extract_annotation_{mut_type}_{sampleid}.sh'.format(**dict(self.__dict__, **locals()))

        utils.write_shell(shell_path, cmd)

        # add job
        now_point = 'extract_annotation'
        job_name = 'extract_annotation_{mut_type}_{sampleid}'.format(
            sampleid=sampleid, mut_type=mut_type)
        utils.add_job(self.jobs, now_point, self.args['startpoint'],
                      self.ANALYSIS_POINTS, job_name, shell_path, self.queues)

        # add order
        before_jobs = ['annotate_merged_{mut_type}'.format(**locals())]
        after_jobs = ['variation_summary']
        utils.add_order(self.orders, job_name, before_jobs=before_jobs, after_jobs=after_jobs)

# =============================================== annotation end =========================================================

    def variation_summary(self):

        print '>  summary'

        countlist = '{analydir}/Advance/{newjob}/Summary/countlist'.format(
            **self.args)
        with utils.safe_open(countlist, 'w') as out:
            for sampleid in self.sample_infos:
                line = '{analydir}/Mutation/{sampleid}.{mutation_soft}/{sampleid}.{mutation_soft}.snp.annovar.hg19_multianno.xls.gz\n'
                if self.mutation_soft == 'mtoolbox':
                    line = '{analydir}/Mutation/{sampleid}.{mutation_soft}/{sampleid}.{mutation_soft}.snp.multianno.xls\n'
                line = line.format(
                    sampleid=sampleid,
                    mutation_soft=self.mutation_soft,
                    **self.args)
                out.write(line)

        vcf_counts = 'Varition/VCF/vcf_count/vcf_counts_v2.py'
        if self.mutation_soft == 'mtoolbox':
            vcf_counts = 'Varition/MT/vcf_counts_v2.py'

        cmd = '''
            set -eo pipefail
            echo variation summary start: `date "+%F %T"`

            cd {analydir}/Advance/{newjob}/Summary

            python {moduledir}/{vcf_counts} \\
                --list countlist

            echo variation summary done: `date "+%F %T"`
        '''.format(
            vcf_counts=vcf_counts, **self.args)

        shell_path = '{analydir}/Advance/{newjob}/Summary/variation_summary.sh'.format(
            **self.args)

        utils.write_shell(shell_path, cmd)

        # add job
        now_point = job_name = 'variation_summary'
        utils.add_job(self.jobs, now_point, self.args['startpoint'],
                      self.ANALYSIS_POINTS, job_name, shell_path, self.queues)

        after_jobs = ['data_release', 'primary_report']
        utils.add_order(self.orders, job_name, after_jobs=after_jobs)

# the end
