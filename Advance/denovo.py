#!/usr/bin/env python
# -*- coding=utf-8 -*-
import os
import string
import utils


class Denovo(object):

    def __init__(self, args, jobs, orders, mutation_soft, sv_soft, cnv_soft, denovo_soft,
                 sample_infos, config, ANALY_DICT):

        self.args = args
        self.jobs = jobs
        self.orders = orders

        self.mutation_soft = mutation_soft

        self.sv_cnv_soft = {
            'sv': sv_soft,
            'cnv': cnv_soft,
        }

        self.denovo_soft = denovo_soft

        self.queues = args.get('queues')

        self.sample_infos = sample_infos

        self.ANALY_DICT = ANALY_DICT

        self.ANALYSIS_POINTS = config.ANALYSIS_POINTS

        self.config = config
        self.__dict__.update(self.args)
        self.ref = args.get('ref')
        self.REF = 'hg19' if self.ref == 'b37' else self.ref
        self.__dict__.update(dict(config.CONFIG.items('genome_' + self.ref)))
        self.__dict__.update(dict(config.CONFIG.items('software')))

        self.denovo_dir_map = {
            'samtools': 'DenovoSam',
            'denovogear': 'DenovoGear',
            'triodenovo': 'DenovoTrio',
            'sv': 'DenovoSV',
            'cnv': 'DenovoCNV'
        }

    def start(self):

        print '>  denovo start...'

        # 当sample_info中有指定data列为2或3的，则只做标有2或3的家系
        # 否则默认sample_info中符合的家系都做
        denovo_infos = dict(
            filter(lambda (k, v): v['data'] in ('2', '3'),
                   self.sample_infos.items()))
        denovo_samples = map(lambda (k, v): v['sampleid'],
                             denovo_infos.items())
        # print denovo_infos
        # print denovo_samples

        # print self.sample_infos;exit()

        peds_raw = utils.sampleinfo2ped(self.sample_infos)

        peds = []
        for ped in peds_raw:
            # print ped
            sampleid = ped['sampleid']

            if denovo_samples and sampleid not in denovo_samples:
                continue

            # 判断父母是否有测序数据
            if (ped['pa'] not in self.sample_infos) or (ped['ma'] not in self.sample_infos):
                continue

            phenotype = ped['phenotype']
            pa_phenotype = ped['pa_context'].get('phenotype')
            ma_phenotype = ped['ma_context'].get('phenotype')
            # print ped['familyid'], ped['sampleid'], phenotype, pa_phenotype, ma_phenotype
            # 孩子和父母均有测序数据，且孩子患病，父母正常
            if (phenotype == '2') and (pa_phenotype == ma_phenotype == '1'):
                familyid_sampleid = '{familyid}_{sampleid}'.format(**ped)
                ped['familyid'] = familyid_sampleid
                ped['pa_context']['familyid'] = familyid_sampleid
                ped['ma_context']['familyid'] = familyid_sampleid
                peds.append(ped)
        if not peds:
            print '[error] no sample can do denovo analysis according to your sample_info, please check!'
            exit(1)

        # SNP/INDEL
        if any([soft in self.denovo_soft for soft in ('samtools', 'triodenovo')]):
            for context in peds:
                # print context
                self.samtools_call_trio(context)
                self.denovo_call(context)
                self.denovo_annotate(context, 'snp')
                self.denovo_annotate(context, 'indel')
                self.denovo_intersect(context)
        self.denovo_rate(peds)

        # SV/CNV
        if self.ANALY_DICT['denovo_sv']:
            for context in peds:
                self.denovo_sv_cnv(context, 'sv')

        if self.ANALY_DICT['denovo_cnv']:
            for context in peds:
                self.denovo_sv_cnv(context, 'cnv')

    def write_ped(self, pedfile, context):

        with utils.safe_open(pedfile, 'w') as ped:
            line = '{familyid}\t{sampleid}\t{pa}\t{ma}\t{sex}\t{phenotype}\n'.format(**context)
            line += '{familyid}\t{sampleid}\t{pa}\t{ma}\t{sex}\t{phenotype}\n'.format(**context['pa_context'])
            line += '{familyid}\t{sampleid}\t{pa}\t{ma}\t{sex}\t{phenotype}\n'.format(**context['ma_context'])
            ped.write(line)

    def samtools_call_trio(self, context):

        print '>   samtools call trio for {familyid}'.format(**context)

        chrom_list = utils.get_chrom_list(
            self.ref, self.sample_infos[context['sampleid']]['sex'], self.args['MT'])

        pedfile = '{analydir}/Advance/{newjob}/Denovo/CallTrio/{familyid}/{familyid}.ped'.format(
            **dict(self.__dict__, **context))

        self.write_ped(pedfile, context)

        for chrom in chrom_list:
            cmd = '''
                set -eo pipefail
                echo samtools call trio for {familyid} start: `date "+%F %T"`

                cd {analydir}/Advance/{newjob}/Denovo/CallTrio/{familyid}/bychr

                # call
                samtools-1.6 mpileup \\
                    -r {chrom} \\
                    -q 1 -C 50 -t DP,AD -m 2 -F 0.002 \\
                    -ugf {reffasta} \\
                    {analydir}/Mapping/{sampleid}.{sampleid}/{sampleid}.final.bam \\
                    {analydir}/Mapping/{pa}.{pa}/{pa}.final.bam \\
                    {analydir}/Mapping/{ma}.{ma}/{ma}.final.bam |
                bcftools-1.6 call  \\
                    -C trio \\
                    -S ../{familyid}.ped \\
                    -v -m \\
                    -o {familyid}.raw.trio.{chrom}.vcf

                echo samtools call trio for {familyid} done: `date "+%F %T"`
            '''.format(chrom=chrom, **dict(self.__dict__, **context))

            shell_path = '{analydir}/Advance/{newjob}/Denovo/CallTrio/{familyid}/bychr/samtools_call_trio_{familyid}_{chrom}.sh'.format(
                chrom=chrom, **dict(self.__dict__, **context))

            utils.write_shell(shell_path, cmd)

            # add job
            now_point = 'samtools_call_trio'
            job_name = 'samtools_call_trio_{familyid}_{chrom}'.format(chrom=chrom, **context)
            utils.add_job(self.jobs, now_point, self.args['startpoint'],
                          self.ANALYSIS_POINTS, job_name, shell_path,
                          self.queues)

            # add order
            before_jobs = ['final_bam_{sampleid}'.format(**context)]
            after_jobs = ['bcftools_concat_trio_{familyid}'.format(**context)]
            utils.add_order(self.orders, job_name, before_jobs=before_jobs, after_jobs=after_jobs)

        # concat and filter
        vcf_list = 'bychr/%s.raw.trio.{{1..22}%s}.vcf' % (
            context['familyid'], ',' + ','.join(chrom_list[22:]))
        # print vcf_list

        cmd = '''
            set -eo pipefail
            echo bcftools concat for {familyid} start: `date "+%F %T"`
            
            cd {analydir}/Advance/{newjob}/Denovo/CallTrio/{familyid}

            bcftools-1.6 concat \\
                -Oz -o {familyid}.raw.trio.vcf.gz \\
                {vcf_list}

            # rm -f bychr/*.vcf

            echo bcftools concat for {familyid} done: `date "+%F %T"`
        '''.format(
            vcf_list=vcf_list, **dict(self.__dict__, **context))

        shell_path = '{analydir}/Advance/{newjob}/Denovo/CallTrio/{familyid}/bcftools_concat_trio_{familyid}.sh'.format(
            **dict(self.__dict__, **context))

        utils.write_shell(shell_path, cmd)

        # add job
        now_point = 'bcftools_concat_trio'
        job_name = 'bcftools_concat_trio_{familyid}'.format(**context)
        utils.add_job(self.jobs, now_point, self.args['startpoint'],
                      self.ANALYSIS_POINTS, job_name, shell_path, self.queues)

    def denovo_call(self, context):

        for soft in self.denovo_soft:

            print '>     denovo call {} for {}'.format(soft, context['familyid'])

            denovo_dir = self.denovo_dir_map[soft]

            if soft == 'samtools':
                cmd = string.Template('''\
                    set -eo pipefail
                    echo denovo call samtools for ${familyid} start: `date "+%F %T"`

                    cd ${analydir}/Advance/${newjob}/Denovo/DenovoSam/${familyid}/

                    bcftools-1.6 filter \\
                        -i '%QUAL>15 && DP>12' \\
                        ../../CallTrio/${familyid}/${familyid}.raw.trio.vcf.gz |
                    awk '$1~/^#/ || ($10!~/0\/0/ && $10!~/\.\/\./ && $11~/0\/0/ && $12~/0\/0/)' \\
                        > ${familyid}.denovo.vcf

                    bcftools-1.6 filter \\
                        -i '%TYPE=="snp"' \\
                        -Oz -o ${familyid}.denovo.snp.vcf.gz \\
                        ${familyid}.denovo.vcf 

                    bcftools-1.6 index -tf ${familyid}.denovo.snp.vcf.gz 

                    bcftools-1.6 filter \\
                        -i '%TYPE=="indel"' \\
                        -Oz -o ${familyid}.denovo.indel.vcf.gz \\
                        ${familyid}.denovo.vcf 

                    bcftools-1.6 index -tf ${familyid}.denovo.indel.vcf.gz 

                    echo denovo call samtools for ${familyid} done: `date "+%F %T"`
                ''')

            # 注意triodenovo输出的顺序为：Father, Mother, Child
            if soft == 'triodenovo':
                cmd = string.Template('''\
                    set -eo pipefail
                    echo denovo call triodenovo for ${familyid} start: `date "+%F %T"`

                    cd ${analydir}/Advance/${newjob}/Denovo/DenovoTrio/${familyid}

                    triodenovo \\
                        --ped ../../CallTrio/${familyid}/${familyid}.ped \\
                        --in_vcf ../../CallTrio/${familyid}/${familyid}.raw.trio.vcf.gz \\
                        --minDepth 10 \\
                        --out_vcf ${familyid}.triodenovo.raw.vcf
                    
                    # change sample order
                    awk -v OFS='\\t' '$1~/^##/{print} $1!~/^##/{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$12,$10,$11}' \\
                        ${familyid}.triodenovo.raw.vcf \\
                        > ${familyid}.triodenovo.mid.vcf

                    grep -E '^#|INDEL' ${familyid}.triodenovo.mid.vcf > ${familyid}.triodenovo.mid.indel.vcf
                    grep -vE 'INDEL' ${familyid}.triodenovo.mid.vcf > ${familyid}.triodenovo.mid.snp.vcf

                    ${moduledir}/Varition/DNM/change_triodenovo_GT.py \\
                        ${familyid}.triodenovo.mid.indel.vcf \\
                        ${familyid}.denovo.indel.vcf \\
                        indel

                    bgzip -f ${familyid}.denovo.indel.vcf

                    ${moduledir}/Varition/DNM/change_triodenovo_GT.py \\
                        ${familyid}.triodenovo.mid.snp.vcf \\
                        ${familyid}.denovo.snp.vcf \\
                        snp

                    bgzip -f ${familyid}.denovo.snp.vcf

                    echo denovo call triodenovo for ${familyid} done: `date "+%F %T"`
                ''')

            cmd = cmd.safe_substitute(**dict(self.__dict__, **context))

            shell_path = '{analydir}/Advance/{newjob}/Denovo/{denovo_dir}/{familyid}/denovo_call_{soft}_{familyid}.sh'.format(
                soft=soft,
                denovo_dir=denovo_dir,
                **dict(self.__dict__, **context))

            utils.write_shell(shell_path, cmd)

            # add job
            now_point = 'denovo_call'
            job_name = 'denovo_call_{soft}_{familyid}'.format(soft=soft, **context)

            utils.add_job(self.jobs, now_point, self.args['startpoint'],
                          self.ANALYSIS_POINTS, job_name, shell_path,
                          self.queues)

            # add order
            before_jobs = ['bcftools_concat_trio_{familyid}'.format(**context)]
            utils.add_order(self.orders, job_name, before_jobs=before_jobs)

    def denovo_annotate(self, context, mtype):

        for soft in self.denovo_soft:

            denovo_dir = self.denovo_dir_map[soft]

            print '>     denovo {} annotate {} for {}'.format(soft, mtype, context['familyid'])

            cmd = '''
                set -eo pipefail
                echo denovo annotate {soft} {mtype} vcf for {familyid} start: `date "+%F %T"`

                cd {analydir}/Advance/{newjob}/Denovo/{denovo_dir}/{familyid}

                # norm
                bcftools-1.6 norm \\
                    -m -both \\
                    -f {reffasta} \\
                    -o {familyid}.denovo.{mtype}_sn.vcf \\
                    {familyid}.denovo.{mtype}.vcf.gz

                # 1 Annotate
                sh {annovar} \\
                    -b {MTYPE} \\
                    -r {REF} \\
                    {familyid}.denovo.{mtype}_sn.vcf \\
                    {sampleid},{pa},{ma}

                # 2 Add PubmedID and reformat HGMD
                python {moduledir}/AddOMIM_HGMD/AddHGMD_OMIM_Priority_pipe4.6.py \\
                    {familyid}.denovo.{mtype}_sn.annovar.{REF}_multianno.xls \\
                    {familyid}.denovo.{mtype}.annovar.{REF}_multianno_mid.xls \\
                    {REF}

                # 3 Add HPA annotation
                python {moduledir}/HPA.v15/annotatExpression_for_multiannofile.py \\
                    -i {familyid}.denovo.{mtype}.annovar.{REF}_multianno_mid.xls \\
                    -o {familyid}.denovo.{mtype}.annovar.{REF}_multianno.xls

                # gzip -f {familyid}.denovo.{mtype}.annovar.{REF}_multianno.xls

                # 4 Filter
                perl {moduledir}/Varition/Filter/filter_annovar_hg19_multianno_xls_pipe4.6.pl \\
                    -pipeline Y \\
                    -input {familyid}.denovo.{mtype}.annovar.{REF}_multianno.xls \\
                    -output {familyid}.denovo.{mtype}

                rm -f *{mtype}*multianno.xls.bak *.{mtype}*multianno_mid.xls
                rm -f {familyid}.denovo.{mtype}_sn.vcf.gz *_sn.annovar.{REF}_multianno.xls

                echo denovo annotate {soft} {mtype} vcf for {familyid} done: `date "+%F %T"`
            '''.format(
                soft=soft,
                denovo_dir=denovo_dir,
                mtype=mtype,
                MTYPE=mtype.upper(),
                **dict(self.__dict__, **context))

            shell_path = '{analydir}/Advance/{newjob}/Denovo/{denovo_dir}/{familyid}/denovo_annotate_{soft}_{mtype}_{familyid}.sh'.format(
                soft=soft,
                mtype=mtype,
                denovo_dir=denovo_dir,
                **dict(self.__dict__, **context))

            utils.write_shell(shell_path, cmd)

            # add job
            now_point = 'denovo_annotate'
            job_name = 'denovo_annotate_{soft}_{mtype}_{familyid}'.format(soft=soft, mtype=mtype, **context)

            utils.add_job(self.jobs, now_point, self.args['startpoint'],
                          self.ANALYSIS_POINTS, job_name, shell_path,
                          self.queues)

            # add order
            before_jobs = ['denovo_call_{soft}_{familyid}'.format(soft=soft, **context)]
            after_jobs = ['denovo_intersect_{familyid}'.format(**context)]
            utils.add_order(self.orders, job_name, before_jobs=before_jobs, after_jobs=after_jobs)

    def denovo_intersect(self, context):

        print '>   intersect result for {familyid}'.format(**context)

        snps = []
        for soft in self.denovo_soft:
            snp = '../../{denovo_dir}/{familyid}/{familyid}.denovo.snp.annovar.{REF}_multianno.xls'.format(
                denovo_dir=self.denovo_dir_map[soft], **dict(self.__dict__, **context))
            snps.append(snp)

        snps = ','.join(snps)
        indels = snps.replace('snp', 'indel')

        cmd = '''
            set -eo pipefail

            echo denovo intersect for {familyid} start: `date "+%F %T"`

            cd {analydir}/Advance/{newjob}/Denovo/Intersect/{familyid}

            # Intersect
            python {moduledir}/Varition/Filter/IntegrateFile_pipe4.6.py \\
                -inputs {snps} \\
                -type intersect \\
                -out SNP/{familyid}.denovo.snp.intersect.xls

            python {moduledir}/Varition/Filter/IntegrateFile_pipe4.6.py \\
                -inputs {indels} \\
                -type intersect \\
                -out INDEL/{familyid}.denovo.indel.intersect.xls

            for mtype in snp indel;do
                MTYPE=`echo $mtype | tr a-z A-Z`

                # Filter
                perl {moduledir}/Varition/Filter/filter_annovar_hg19_multianno_xls_pipe4.6.pl \\
                    -pipeline Y \\
                    -input $MTYPE/{familyid}.denovo.$mtype.intersect.xls \\
                    -output $MTYPE/{familyid}.denovo.$mtype.intersect

                # Brief Result
                echo generate brief results for $MTYPE

                python {ROOT_DIR}/modules/brief/brief_anno.py \\
                    -i $MTYPE/{familyid}.denovo.$mtype.intersect.freq.func.syn.deleterious.xls \\
                    -O {BriefResults}/Denovo/{familyid} \\
                    -t snpindel

                python {ROOT_DIR}/modules/brief/text2excel.py \\
                    {BriefResults}/Denovo/{familyid}/{familyid}.denovo.$mtype.intersect.freq.func.syn.deleterious.xlsx \\
                    {ROOT_DIR}/modules/brief/readme/denovo.readme.xls \\
                    $MTYPE/{familyid}.denovo.$mtype.intersect.filter.stat.xls \\
                    {BriefResults}/Denovo/{familyid}/{familyid}.denovo.$mtype.intersect.freq.func.syn.deleterious.brief.xls
            done

            echo denovo intersect for {familyid} done: `date "+%F %T"`
        '''.format(
            snps=snps, indels=indels, **dict(self.__dict__, **context))

        shell_path = '{analydir}/Advance/{newjob}/Denovo/Intersect/{familyid}/denovo_intersect_{familyid}.sh'.format(
            **dict(self.__dict__, **context))

        utils.write_shell(shell_path, cmd)

        # add job
        now_point = 'denovo_intersect'
        job_name = 'denovo_intersect_{familyid}'.format(**context)
        utils.add_job(self.jobs, now_point, self.args['startpoint'],
                      self.ANALYSIS_POINTS, job_name, shell_path, self.queues)

        # add order
        before_jobs = []
        after_jobs = ['denovo_rate', 'integrate_result']
        utils.add_order(
            self.orders,
            job_name,
            before_jobs=before_jobs,
            after_jobs=after_jobs)

    def denovo_rate(self, peds):

        print '>   denovo rate calculate'

        familyids = [ped['familyid'] for ped in peds]
        # print familyids

        for mtype in ('snp', 'indel'):
            MTYPE=mtype.upper()
            configfile = '{analydir}/Advance/{newjob}/Denovo/DenovoRate/denovo_stat_{mtype}.conf'.format(
                mtype=mtype, **self.__dict__)
            # print configfile
            with utils.safe_open(configfile, 'w') as conf:
                for familyid in familyids:
                    sampleid = familyid.split('_', 1)[1]
                    linelist = [
                        '{familyid}',
                        '{analydir}/Alnstat/{sampleid}/{sampleid}.sample_cumulative_coverage_counts',
                        '../Intersect/{familyid}/{MTYPE}/{familyid}.denovo.{mtype}.intersect.xls'
                    ]
                    line = '\t'.join(linelist).format(**dict(self.__dict__, **locals())) + '\n'
                    conf.write(line)

        cmd = '''
            set -eo pipefail

            echo denovo rate calculate start: `date "+%F %T"`

            cd {analydir}/Advance/{newjob}/Denovo/DenovoRate/

            python {moduledir}/Varition/DNM/DNM_Rate.py \\
                -config denovo_stat_snp.conf \\
                -out denovo_rate_snp.xls

            python {moduledir}/Varition/DNM/DNM_Rate.py \\
                -config denovo_stat_indel.conf \\
                -out denovo_rate_indel.xls
            
            echo denovo rate calculate done: `date "+%F %T"`
        '''

        cmd = cmd.format(**self.__dict__)

        shell_path = '{analydir}/Advance/{newjob}/Denovo/DenovoRate/denovo_rate.sh'.format(**self.__dict__)

        utils.write_shell(shell_path, cmd)

        # add job
        now_point = 'denovo_rate'
        job_name = 'denovo_rate'

        utils.add_job(self.jobs, now_point, self.args['startpoint'],
                      self.ANALYSIS_POINTS, job_name, shell_path, self.queues)

        # add order
        before_jobs = []
        after_jobs = ['data_release']
        utils.add_order(
            self.orders,
            job_name,
            before_jobs=before_jobs,
            after_jobs=after_jobs)

    def denovo_sv_cnv(self, context, svtype):

        print '>    denovo {svtype} for {familyid}'.format(svtype=svtype, **context)

        denovo_dir = self.denovo_dir_map[svtype]

        pedfile = '{analydir}/Advance/{newjob}/Denovo/{denovo_dir}/{familyid}/{familyid}.ped'.format(
            denovo_dir=denovo_dir, **dict(self.__dict__, **context))

        self.write_ped(pedfile, context)

        cmd = '''
            set -eo pipefail

            echo denovo {svtype} start: `date "+%F %T"`

            cd {analydir}/Advance/{newjob}/Denovo/{denovo_dir}/{familyid}

            python {moduledir}/Varition/DNM/denovo_sv_cnv.py \\
                --proj {analydir} \\
                --ped {familyid}.ped \\
                --outdir {analydir}/Advance/{newjob}/Denovo \\
                --soft {soft}

            # Brief Result
            echo generate brief results

            python {ROOT_DIR}/modules/brief/text2excel.py \\
                {BriefResults}/Denovo/{familyid}/{familyid}.denovo{SVTYPE}.hg19_multianno.xlsx \\
                {ROOT_DIR}/modules/brief/readme/denovo_sv_cnv.readme.xls \\
                {familyid}.denovo{SVTYPE}.hg19_multianno.xls

            echo denovo {svtype} done: `date "+%F %T"`
        '''

        cmd = cmd.format(
            svtype=svtype,
            SVTYPE=svtype.upper(),
            soft=self.sv_cnv_soft[svtype],
            denovo_dir=denovo_dir,
            **dict(self.__dict__, **context))

        shell_path = '{analydir}/Advance/{newjob}/Denovo/{denovo_dir}/{familyid}/denovo_{svtype}_{familyid}.sh'.format(
            svtype=svtype,
            denovo_dir=denovo_dir,
            **dict(self.__dict__, **context))

        utils.write_shell(shell_path, cmd)

        # add job
        now_point = 'denovo_{svtype}'.format(svtype=svtype)
        job_name = 'denovo_{svtype}_{familyid}'.format(svtype=svtype, **context)

        utils.add_job(self.jobs, now_point, self.args['startpoint'],
                      self.ANALYSIS_POINTS, job_name, shell_path, self.queues)

        # add order
        before_jobs = []
        after_jobs = ['data_release']
        utils.add_order(
            self.orders,
            job_name,
            before_jobs=before_jobs,
            after_jobs=after_jobs)

# the end
