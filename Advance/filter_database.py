#!/usr/bin/env python
# -*- coding=utf-8 -*-
from string import Template
from collections import defaultdict
import utils

class FilterDB(object):

    def __init__(self, args, jobs, orders, mutation_soft, sv_soft, cnv_soft,
                 sample_infos, config, disease_name, tissue, ANALY_DICT):

        self.args = args
        self.__dict__.update(self.args)
        
        self.jobs = jobs
        self.orders = orders

        self.mutation_soft = mutation_soft
        self.sv_soft = sv_soft
        self.cnv_soft = cnv_soft

        self.queues = args.get('queues')

        self.disease_name = disease_name

        self.tissue = tissue

        self.sample_infos = sample_infos

        self.ANALY_DICT = ANALY_DICT

        self.familyids = defaultdict(set)
        for sampleid, items in sample_infos.iteritems():
            self.familyids[items['familyid']].add(sampleid)

        self.ANALYSIS_POINTS = config.ANALYSIS_POINTS

        self.config = config

        self.ref = args.get('ref')
        self.__dict__.update(dict(config.CONFIG.items('genome_' + self.ref)))
        self.__dict__.update(dict(config.CONFIG.items('software')))

    def start(self):

        print '>  filter_database start...'
        # ====== filter_db calling ======
        if self.ANALY_DICT['filter_db'] or self.ANALY_DICT['filter_acmg']:
            self.filter_snpindel()

        if self.ANALY_DICT['filter_acmg']:
            self.filter_acmg()
        if self.ANALY_DICT['filter_sv']:
            self.filter_sv()
        if self.ANALY_DICT['filter_cnv']:
            self.filter_cnv()
        if self.ANALY_DICT['filter_noncoding']:
            self.filter_noncoding()

    def filter_snpindel(self):

        print '>   filter_snpindel'
        cmd = '''
            set -eo pipefail
            echo filter database start: `date "+%F %T"`

            cd {Merged_vcf}

            mkdir -p Filter/SNP Filter/INDEL

            # Filter SNP
            perl {moduledir}/Varition/Filter/filter_annovar_hg19_multianno_xls_pipe4.6.pl \\
                -i VCF/snp.merged.annovar.hg19_multianno.xls.gz \\
                -pipeline Y \\
                -output Filter/SNP/snp.merged

            # Filter INDEL
            perl {moduledir}/Varition/Filter/filter_annovar_hg19_multianno_xls_pipe4.6.pl \\
                -i VCF/indel.merged.annovar.hg19_multianno.xls.gz \\
                -pipeline Y \\
                -output Filter/INDEL/indel.merged

            # Merge SNP INDEL
            head -1 Filter/SNP/snp.merged.freq.func.syn.deleterious.xls > Filter/snp.indel.merged.freq.func.syn.deleterious.xls
            tail -n +2 Filter/SNP/snp.merged.freq.func.syn.deleterious.xls > Filter/snp.tmp
            tail -n +2 Filter/INDEL/indel.merged.freq.func.syn.deleterious.xls > Filter/indel.tmp
            sort -k2,3V -m Filter/snp.tmp Filter/indel.tmp >> Filter/snp.indel.merged.freq.func.syn.deleterious.xls
            rm -f Filter/snp.tmp Filter/indel.tmp

            # Brief Result
            echo generate brief results

            for mtype in snp indel;do
                MTYPE=`echo $mtype | tr a-z A-Z`

                python {ROOT_DIR}/modules/brief/brief_anno.py \\
                    -i Filter/$MTYPE/$mtype.merged.freq.func.syn.deleterious.xls \\
                    -O {BriefResults}/FilterDB \\
                    -t snpindel

                python {ROOT_DIR}/modules/brief/text2excel.py \\
                    {BriefResults}/FilterDB/$mtype.merged.freq.func.syn.deleteriou.xlsx \\
                    {ROOT_DIR}/modules/brief/readme/filterdb.readme.xls \\
                    Filter/$MTYPE/$mtype.merged.filter.stat.xls \\
                    {BriefResults}/FilterDB/$mtype.merged.freq.func.syn.deleterious.brief.xls
            done

            echo filter database done: `date "+%F %T"`
        '''.format(**self.__dict__)

        shell_path = '{analydir}/Advance/{newjob}/Merged_vcf/filter_database.sh'.format(
            **self.args)

        utils.write_shell(shell_path, cmd)

        # add job
        now_point = 'filter_snpindel'
        job_name = 'filter_snpindel'
        utils.add_job(self.jobs, now_point, self.args['startpoint'],
                      self.ANALYSIS_POINTS, job_name, shell_path, self.queues)

        # add order
        before_jobs = ['annotate_merged_snp', 'annotate_merged_indel']
        after_jobs = ['data_release']
        utils.add_order(self.orders, job_name, before_jobs=before_jobs, after_jobs=after_jobs)

    def filter_acmg(self):

        print '>   filter_acmg'

        for familyid in self.familyids:
            # print familyid
            cmd = Template('''
                set -eo pipefail
                echo acmg class start: `date "+%F %T"`

                cd ${ACMG}/${familyid}

                python ${moduledir}/Varition/Filter/acmg.class.py \\
                    -in ${Merged_vcf}/VCF/snp.merged.annovar.hg19_multianno.xls.gz,${Merged_vcf}/VCF/indel.merged.annovar.hg19_multianno.xls.gz \\
                    -ed F \\
                    -fid ${familyid} \\
                    -samp_info ${samp_info} \\
                    -type snp.indel \\
                    -prefix snp.indel \\
                    -o .
            
                # Brief Result
                echo generate brief results

                python ${ROOT_DIR}/modules/brief/brief_anno.py \\
                    -i ${familyid}/${familyid}.snp.indel.{Pathogenic,LikelyPathogenic,VUS,LikelyBenign,Benign}.xls \\
                    -O ${BriefResults}/ACMG/${familyid} \\
                    -t acmg

                python ${ROOT_DIR}/modules/brief/text2excel.py \\
                    ${BriefResults}/ACMG/${familyid}/${familyid}.snp.indel.Pathogenic.LikelyPathogenic.VUS.LikelyBenign.Benign.xlsx \\
                    ${ROOT_DIR}/modules/brief/readme/acmg.readme.xls \\
                    ${familyid}/${familyid}.snp.indel.stat.xls \\
                    ${BriefResults}/ACMG/${familyid}/${familyid}.snp.indel.{Pathogenic,LikelyPathogenic,VUS,LikelyBenign,Benign}.brief.xls

                rm -rf Merge_snp_indel Extract_effective ModelFilter

                echo acmg class done: `date "+%F %T"`
            ''')

            cmd = cmd.safe_substitute(familyid=familyid, **self.__dict__)

            shell_path = '{analydir}/Advance/{newjob}/ACMG/{familyid}/filter_acmg_{familyid}.sh'.format(
                familyid=familyid, **self.args)

            utils.write_shell(shell_path, cmd)

            # add job
            now_point = 'filter_acmg'
            job_name = 'filter_acmg_{}'.format(familyid)
            utils.add_job(self.jobs, now_point, self.args['startpoint'],
                          self.ANALYSIS_POINTS, job_name, shell_path,
                          self.queues)

            # add order
            before_jobs = ['annotate_merged_snp', 'annotate_merged_indel']
            after_jobs = ['data_release', 'integrate_result']
            utils.add_order(self.orders, job_name, before_jobs=before_jobs, after_jobs=after_jobs)

    def filter_sv(self):

        print '>   filter_sv'
        for sampleid in self.sample_infos:
            cmd = '''
                set -eo pipefail
                echo filter sv for {sampleid} start: `date "+%F %T"`

                cd {FilterSV}/{sampleid}

                python {moduledir}/Varition/Filter/filter_sv_cnv.py \\
                    --proj {analydir} \\
                    --overlap 0.7 \\
                    --sampleID {sampleid} \\
                    --outdir {FilterSV} \\
                    --soft {sv_soft} \\
                    --lib StringentLib,InclusiveLib,DGV.GoldStandard.July2015,DGV,CNVD

                # Brief Result
                echo generate brief results

                python {ROOT_DIR}/modules/brief/brief_anno.py \\
                    -i {sampleid}.LikelyDeleterious.SV.xls \\
                    -O {BriefResults}/FilterSV \\
                    -t sv_cnv

                python {ROOT_DIR}/modules/brief/text2excel.py \\
                    {BriefResults}/FilterSV/{sampleid}.LikelyDeleterious.SV.xlsx \\
                    {ROOT_DIR}/modules/brief/readme/filter_sv_cnv.readme.xls \\
                    {BriefResults}/FilterSV/{sampleid}.LikelyDeleterious.SV.brief.xls

                echo filter sv for {sampleid} done: `date "+%F %T"`
            '''.format(
                sampleid=sampleid,
                **self.__dict__)

            shell_path = '{FilterSV}/{sampleid}/filter_sv_{sampleid}.sh'.format(
                sampleid=sampleid, **self.args)

            utils.write_shell(shell_path, cmd)

            # add job
            now_point = 'filter_sv'
            job_name = 'filter_sv_{}'.format(sampleid)
            utils.add_job(self.jobs, now_point, self.args['startpoint'],
                          self.ANALYSIS_POINTS, job_name, shell_path,
                          self.queues)

            # add order
            if self.sv_soft == 'crest':
                sv_last = 'crest_txt2gff_{sampleid}'
            elif self.sv_soft in ('breakdancer', 'lumpy'):
                sv_last = '{sv_soft}_call_{sampleid}'

            before_jobs = [sv_last.format(**dict(locals(), **self.__dict__))]
            after_jobs = ['data_release']
            utils.add_order(self.orders, job_name, before_jobs=before_jobs, after_jobs=after_jobs)

    def filter_cnv(self):

        print '>   filter_cnv'
        for sampleid in self.sample_infos:
            cmd = '''
                set -eo pipefail
                echo filter cnv for {sampleid} start: `date "+%F %T"`

                cd {FilterCNV}/{sampleid}

                python {moduledir}/Varition/Filter/filter_sv_cnv.py \\
                    --proj {analydir} \\
                    --overlap 0.7 \\
                    --sampleID {sampleid} \\
                    --outdir {FilterCNV} \\
                    --soft {cnv_soft} \\
                    --lib StringentLib,InclusiveLib,DGV.GoldStandard.July2015,DGV,CNVD

                # Brief Result
                echo generate brief results

                python {ROOT_DIR}/modules/brief/brief_anno.py \\
                    -i {sampleid}.LikelyDeleterious.CNV.xls \\
                    -O {BriefResults}/FilterCNV \\
                    -t sv_cnv

                python {ROOT_DIR}/modules/brief/text2excel.py \\
                    {BriefResults}/FilterCNV/{sampleid}.LikelyDeleterious.CNV.xlsx \\
                    {ROOT_DIR}/modules/brief/readme/filter_sv_cnv.readme.xls \\
                    {BriefResults}/FilterCNV/{sampleid}.LikelyDeleterious.CNV.brief.xls

                echo filter cnv for {sampleid} done: `date "+%F %T"`
            '''.format(
                sampleid=sampleid,
                **self.__dict__)

            shell_path = '{FilterCNV}/{sampleid}/filter_cnv_{sampleid}.sh'.format(
                sampleid=sampleid, **self.args)

            utils.write_shell(shell_path, cmd)

            # add job
            now_point = 'filter_cnv'
            job_name = 'filter_cnv_{}'.format(sampleid)
            utils.add_job(self.jobs, now_point, self.args['startpoint'],
                          self.ANALYSIS_POINTS, job_name, shell_path,
                          self.queues)

            # add order
            if self.cnv_soft == 'freec':
                cnv_last = 'freec_call_{}'.format(sampleid)
            elif self.cnv_soft == 'conifer':
                cnv_last = 'conifer_call'

            before_jobs = [cnv_last]
            after_jobs = ['data_release']
            utils.add_order(self.orders, job_name, before_jobs=before_jobs, after_jobs=after_jobs)

    def filter_noncoding(self):

        print '>   filter noncoding'

        if not self.disease_name:
            print '[error] noncoding needs disease name in your sample_info'
            exit(1)

        for familyid in self.familyids:
            # print familyid
            samples = ','.join(self.familyids[familyid])
            pedfile = '{Noncoding}/{familyid}/ped.xls'.format(familyid=familyid, **self.__dict__)
            with utils.safe_open(pedfile, 'w') as out:
                for sample in samples.split(','):
                    line = 'N\t{sample}\t0\t0\t0\t0\n'.format(sample=sample)
                    out.write(line)

            tissue = ''
            model = ''
            model_list = []
            if self.tissue:
                tissue += " --tissue '{}' ".format(self.tissue)
            if self.ANALY_DICT['model_dominant']:
                model_list += ['D']
            if self.ANALY_DICT['model_recessive']:
                model_list += ['R', 'C']
            if model_list:
                model = " --model '{}' ".format(';'.join(model_list))

            cmd = '''
                set -eo pipefail
                echo filter noncoding for {familyid} start: `date "+%F %T"`

                cd {Noncoding}/{familyid}

                for mtype in snp indel;do

                    bcftools-1.6 norm \\
                        -m -both \\
                        -f {reffasta} \\
                        {Merged_vcf}/VCF/$mtype.merged.vcf.gz |
                    bcftools-1.6 view \\
                        -x -s {samples} \\
                        -Oz -o {familyid}.$mtype.merged_sn.vcf.gz

                    bcftools-1.6 norm \\
                        -m +both \\
                        -Oz -o {familyid}.$mtype.merged.vcf.gz \\
                        {familyid}.$mtype.merged_sn.vcf.gz

                    python {moduledir}/Varition/Ann/GetAnnoFromMerge_v4.6.0.2_multiple.py \\
                        -M {Merged_vcf}/VCF/$mtype.merged.annovar.hg19_multianno.xls.gz \\
                        -nv {familyid}.$mtype.merged_sn.vcf.gz \\
                        -rv {familyid}.$mtype.merged.vcf.gz \\
                        -S {samples} \\
                        -T $mtype \\
                        -ref {ref} \\
                        -O {familyid}.$mtype.merged.annovar.hg19_multianno.xls

                    gzip -f {familyid}.$mtype.merged.annovar.hg19_multianno.xls
                    rm -f {familyid}.$mtype.merged_sn.vcf.gz
                done

                python {moduledir}/Nocoding_Genomiser/noncoding_genomiser_GTEx_Epigenome.Filter.py \\
                    --vcf '{familyid}.snp.merged.vcf.gz;{familyid}.indel.merged.vcf.gz' \\
                    --anno '{familyid}.snp.merged.annovar.hg19_multianno.xls.gz;{familyid}.indel.merged.annovar.hg19_multianno.xls.gz' \\
                    --ped ped.xls \\
                    --samp_info {samp_info} \\
                    --diseases '{disease_name}' {tissue}{model} \\
                    --Fid {familyid} \\
                    --type 'snp;indel' \\
                    --freq 'GnomAD_EAS_AF 0.01' \\
                    --out .
                
                # Brief Result
                echo generate brief results

                for mtype in snp indel;do
                    python {ROOT_DIR}/modules/brief/brief_anno.py \\
                        -i {familyid}.$mtype.anno.noncoding.GnomAD_EAS_AF.conserv.epigenome.xls.gz \\
                        -O {BriefResults}/Noncoding/{familyid} \\
                        -t noncoding

                    python {ROOT_DIR}/modules/brief/text2excel.py \\
                        {BriefResults}/Noncoding/{familyid}/{familyid}.$mtype.anno.noncoding.GnomAD_EAS_AF.conserv.epigenome.xlsx \\
                        {ROOT_DIR}/modules/brief/readme/noncoding.readme.xls \\
                        {familyid}.$mtype.anno.noncoding.filter.stat.xls.gz \\
                        {BriefResults}/Noncoding/{familyid}/{familyid}.$mtype.anno.noncoding.GnomAD_EAS_AF.conserv.epigenome.brief.xls.gz
                done

                echo filter noncoding for {familyid} done: `date "+%F %T"`
            '''.format(**dict(self.__dict__, **locals()))

            shell_path = '{Noncoding}/{familyid}/filter_noncoding_{familyid}.sh'.format(
                familyid=familyid, **self.args)

            utils.write_shell(shell_path, cmd)

            # add job
            now_point = 'filter_noncoding'
            job_name = 'filter_noncoding_{}'.format(familyid)
            utils.add_job(self.jobs, now_point, self.args['startpoint'],
                          self.ANALYSIS_POINTS, job_name, shell_path,
                          self.queues)

            # add order
            before_jobs = ['annotate_merged_snp', 'annotate_merged_indel']
            after_jobs = ['data_release']
            utils.add_order(self.orders, job_name, before_jobs=before_jobs, after_jobs=after_jobs)

# the end
