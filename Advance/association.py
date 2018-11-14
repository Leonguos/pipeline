#!/usr/bin/env python
# -*- coding=utf-8 -*-
import utils


class Association(object):

    def __init__(self, args, jobs, orders, config):

        self.args = args
        self.jobs = jobs
        self.orders = orders
        self.queues = args.get('queues')
        self.ANALYSIS_POINTS = config.ANALYSIS_POINTS

        self.__dict__.update(**args)

        if ('V6' in args['TR']) or args['seqstrag'] == 'WGS':
            self.__dict__.update(
                dict(config.CONFIG.items('burden_v6')))
        elif 'V5' in args['TR']:
            self.__dict__.update(
                dict(config.CONFIG.items('burden_v5')))
        else:
            print '[error] burden analysis only support: WES V5, WES V6 or WGS'
            exit(1)

    def site_association(self):

        print '>  site association ...'
        # write shell
        cmd = '''
            set -eo pipefail
            echo site association start: `date "+%F %T"`

            cd {analydir}/Advance/{newjob}/SiteAS

            python {moduledir}/Association/site_AS/v3/SiteAS.py \\
                --pwd . \\
                --infile ../Merged_vcf/VCF/snp.merged.annovar.hg19_multianno.xls.gz \\
                --type allele \\
                --db NOVO \\
                --out AS
            
            sh AS_siteAS_Allele_NOVO.sh

            echo site association done: `date "+%F %T"`
        '''.format(**self.args)

        shell_path = '{analydir}/Advance/{newjob}/SiteAS/site_association.sh'.format(
            **self.args)

        utils.write_shell(shell_path, cmd)

        # add job
        now_point = job_name = 'site_association'
        utils.add_job(self.jobs, now_point, self.args['startpoint'],
                      self.ANALYSIS_POINTS, job_name, shell_path, self.queues)

        # add order
        before_jobs = ['annotate_merged_snp', 'annotate_merged_indel']
        after_jobs = ['data_release']
        utils.add_order(self.orders, job_name, before_jobs=before_jobs, after_jobs=after_jobs)

    def gene_as_filter(self):

        print '>  filter for gene association ...'
        # write shell
        cmd = '''
            set -eo pipefail
            echo filter for gene association start: `date "+%F %T"`

            cd {analydir}/Advance/{newjob}/GeneAS

            ## 1 Filter
            # filter more strict
            for mtype in snp indel;do
                python {moduledir}/Association/Burden/VariantFilter.py \\
                    --in ../Merged_vcf/VCF/$mtype.merged.annovar.hg19_multianno.xls.gz \\
                    --freqli '1000g_ALL 0.01;GnomAD_ALL_AF 0.001;GnomAD_EAS_AF 0.001;NovoDb_WES 0.01' \\
                    --func \\
                    --sys \\
                    --loss \\
                    --sp 2 \\
                    --dam P \\
                    --gerp \\
                    --repeat \\
                    --out Filter/$mtype
            done

            suffix=func.sys.1000g_ALL.GnomAD_ALL_AF.GnomAD_EAS_AF.NovoDb_WES.repeat.xls

            # remove XY for snp and snp.indel
            awk '$2!~/X/ && $2!~/Y/' Filter/snp.$suffix |
                gzip - > snp.filter.noXY.xls.gz
            
            awk '(NR==1 || FNR!=1) && $2!~/X$/ && $2!~/Y$/' Filter/snp.$suffix Filter/indel.$suffix |
                gzip - > snp.indel.filter.noXY.xls.gz

            ## 2 Extract Bed
            python {moduledir}/Association/Burden/ExtractAgilentBed.py \\
                snp.filter.noXY.xls.gz \\
                {pad_100_bed} \\
                snp.filter.noXY.V6Pad100.xls
            gzip -f snp.filter.noXY.V6Pad100.xls

            python {moduledir}/Association/Burden/ExtractAgilentBed.py \\
                snp.indel.filter.noXY.xls.gz \\
                {pad_100_bed} \\
                snp.indel.filter.noXY.V6Pad100.xls
            gzip -f snp.indel.filter.noXY.V6Pad100.xls

            echo filter for gene association done: `date "+%F %T"`
        '''.format(**self.__dict__)



        shell_path = '{analydir}/Advance/{newjob}/GeneAS/gene_as_filter.sh'.format(
            **self.args)

        utils.write_shell(shell_path, cmd)

        # add job
        now_point = job_name = 'gene_as_filter'
        utils.add_job(self.jobs, now_point, self.args['startpoint'],
                      self.ANALYSIS_POINTS, job_name, shell_path, self.queues)

        # add order
        before_jobs = ['annotate_merged_snp', 'annotate_merged_indel']
        after_jobs = []
        utils.add_order(self.orders, job_name, before_jobs=before_jobs, after_jobs=after_jobs)

    def gene_association(self):

        print '>  gene association ...'
        self.gene_as_filter()
        # write shell
        cmd = '''
            set -eo pipefail
            echo gene association start: `date "+%F %T"`

            cd {analydir}/Advance/{newjob}/GeneAS

            for mtype in snp snp.indel;do
                python {moduledir}/Association/Burden/GetBurdenFre.py \\
                    -case $mtype.filter.noXY.V6Pad100.xls.gz \\
                    -control {moduledir}/{pad_100_stat} \\
                    -cr 0.95 \\
                    -nr 0.6 \\
                    -cc N \\
                    -Num 2827 \\
                    -out $mtype.burden.stat.xls
                
                rows=`wc -l $mtype.burden.stat.xls | cut -d' ' -f1`
                if [ $rows -eq 1 ];then
                    echo "[error] no data in $mtype.burden.stat.xls"
                    exit 1
                fi

                Rscript {moduledir}/Association/Burden/GeneFisherPlot.R \\
                    --infile $mtype.burden.stat.xls \\
                    --outpre $mtype.burden

                paste $mtype.burden.fisher.xls $mtype.burden.stat.samstat.xls > $mtype.burden.result.xls
            done

            echo gene association done: `date "+%F %T"`
        '''.format(**self.__dict__)



        shell_path = '{analydir}/Advance/{newjob}/GeneAS/gene_association.sh'.format(
            **self.args)

        utils.write_shell(shell_path, cmd)

        # add job
        now_point = job_name = 'gene_association'
        utils.add_job(self.jobs, now_point, self.args['startpoint'],
                      self.ANALYSIS_POINTS, job_name, shell_path, self.queues)

        # add order
        before_jobs = ['gene_as_filter']
        after_jobs = ['data_release']
        utils.add_order(self.orders, job_name, before_jobs=before_jobs, after_jobs=after_jobs)

# the end
