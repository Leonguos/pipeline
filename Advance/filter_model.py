#!/usr/bin/env python
# -*- coding=utf-8 -*-
import utils


class FilterModel(object):

    def __init__(self, args, jobs, orders, mutation_soft, sv_soft, cnv_soft, sample_infos, config):

        self.args = args
        self.jobs = jobs
        self.orders = orders

        self.mutation_soft = mutation_soft
        self.sv_soft = sv_soft
        self.cnv_soft = cnv_soft

        self.queues = args.get('queues')

        self.sample_infos = sample_infos

        self.patients = [
            sample for sample in sample_infos
            if sample_infos[sample]['phenotype'] == 'P'
        ]
        self.normals = [
            sample for sample in sample_infos
            if sample_infos[sample]['phenotype'] == 'N'
        ]

        self.familyids = set(
            each.get('familyid') for each in sample_infos.itervalues())

        self.ANALYSIS_POINTS = config.ANALYSIS_POINTS

        self.config = config
        self.__dict__.update(self.args)

        self.ref = args.get('ref')
        self.__dict__.update(dict(config.CONFIG.items('genome_' + self.ref)))
        self.__dict__.update(dict(config.CONFIG.items('software')))

    def start(self):

        print '>  filter_database start...'
        # ====== filter_db calling ======
        if 7.1 in self.args['analy_array']:
            self.filter_model('dominant')
        if 7.2 in self.args['analy_array']:
            self.filter_model('recessive')
            self.filter_model('compound_heterozygous')
        if 7.3 in self.args['analy_array']:
            self.share_compare()

    def filter_model(self, model):

        print '>   filter_model {}'.format(model)

        for familyid in self.familyids:
            if model in ('dominant', 'recessive'):
                cmd = '''
                    set -eo pipefail
                    echo filter mendel{model} for {familyid} start: `date "+%F %T"`

                    cd {ModelF}/{familyid}

                    for mtype in snp indel;do

                        MTYPE=`echo $mtype | tr a-z A-Z`

                        python {moduledir}/ModelFilter/filter_model.py \\
                            -M {MODEL} \\
                            -Fid {familyid} \\
                            -S {samp_info} \\
                            -I {Merged_vcf}/Filter/$MTYPE/$mtype.merged.freq.func.syn.deleterious.xls \\
                            -T $mtype \\
                            -O {ModelF}

                        # Brief Result
                        echo generate brief results

                        python {ROOT_DIR}/modules/brief/brief_anno.py \\
                            -i {familyid}.$mtype.{model}.xls \\
                            -O {BriefResults}/ModelF/{familyid} \\
                            -t snpindel

                        python {ROOT_DIR}/modules/brief/text2excel.py \\
                            {BriefResults}/ModelF/{familyid}/{familyid}.$mtype.{model}.xlsx \\
                            {ROOT_DIR}/modules/brief/readme/model.{model}.readme.xls \\
                            {familyid}.$mtype.{model}.xls.CandidateGene.xls \\
                            {familyid}.$mtype.{model}.xls
                            
                    done

                    echo filter mendel{model} for {familyid} done: `date "+%F %T"`
                '''
            else:
                # 符合杂合是以基因为单位，可以是SNP和INDEL的组合
                cmd = '''
                    set -eo pipefail
                    echo filter mendel{model} for {familyid} start: `date "+%F %T"`

                    cd {ModelF}/{familyid}

                    # SNP and INDEL
                    python {moduledir}/ModelFilter/filter_model.py \\
                        -M {MODEL} \\
                        -Fid {familyid} \\
                        -S {samp_info} \\
                        -I {Merged_vcf}/Filter/snp.indel.merged.freq.func.syn.deleterious.xls \\
                        -T snp.indel \\
                        -O {ModelF}

                        # Brief Result
                        echo generate brief results

                        python {ROOT_DIR}/modules/brief/brief_anno.py \\
                            -i {familyid}.snp.indel.{model}.xls \\
                            -O {BriefResults}/ModelF/{familyid} \\
                            -t snpindel

                        python {ROOT_DIR}/modules/brief/text2excel.py \\
                            {BriefResults}/ModelF/{familyid}/{familyid}.snp.indel.{model}.xlsx \\
                            {ROOT_DIR}/modules/brief/readme/model.{model}.readme.xls \\
                            {familyid}.snp.indel.{model}.xls.CandidateGene.xls \\
                            {BriefResults}/ModelF/{familyid}/{familyid}.snp.indel.{model}.brief.xls
                            
                    echo filter mendel{model} for {familyid} done: `date "+%F %T"`
                '''

            cmd = cmd.format(
                model=model,
                MODEL=model[0].upper(),
                familyid=familyid,
                **self.args)

            shell_path = '{ModelF}/{familyid}/filter_model_{model}_{familyid}.sh'.format(
                model=model, familyid=familyid, **self.args)

            utils.write_shell(shell_path, cmd)

            # add job
            now_point = 'filter_model'
            job_name = 'filter_model_{model}_{familyid}'.format(model=model, familyid=familyid)
            utils.add_job(self.jobs, now_point, self.args['startpoint'],
                          self.ANALYSIS_POINTS, job_name, shell_path,
                          self.queues)

            # add order
            before_jobs = ['filter_snpindel']
            after_jobs = ['integrate_result', 'data_release']
            utils.add_order(
                self.orders,
                job_name,
                before_jobs=before_jobs,
                after_jobs=after_jobs)

    def share_compare(self):

        print '>   share_compare'
        if not self.patients:
            print '[Error] No patients in your sample_info, please check!'
            exit(1)

        controls = ''
        if self.normals:
            controls = '\\\n                -cL {} -cf 0.9 '.format(','.join(self.normals))

        cmd = '''
            set -eo pipefail
            echo share compare start: `date "+%F %T"`\n

            cd {analydir}/Advance/{newjob}/Share

            python {moduledir}/Varition/Filter/CompareCaseControl.py \\
                -in {Merged_vcf}/Filter/snp.indel.merged.freq.func.syn.deleterious.xls \\
                -t gene \\
                -sL {patients} -sf 0 {controls}\\
                -o snp.indel.filter.ShareGene.PatientShare.NotInNormal.xls

            python {moduledir}/Varition/Filter/MutationLand_pipe4.6.py \\
                -I snp.indel.filter.ShareGene.PatientShare.NotInNormal.xls \\
                -S {patients} \\
                -F Y

            # Brief Result
            echo generate brief results

            python {ROOT_DIR}/modules/brief/brief_anno.py \\
                -i snp.indel.filter.ShareGene.PatientShare.NotInNormal.xls \\
                -O {BriefResults}/Share \\
                -t snpindel

            python {ROOT_DIR}/modules/brief/text2excel.py \\
                {BriefResults}/Share/snp.indel.filter.ShareGene.PatientShare.NotInNormal.xlsx \\
                {ROOT_DIR}/modules/brief/readme/share.readme.xls \\
                snp.indel.filter.ShareGene.PatientShare.NotInNormal.stat.xls \\
                {BriefResults}/Share/snp.indel.filter.ShareGene.PatientShare.NotInNormal.brief.xls

            gene=`head -1 {BriefResults}/Share/snp.indel.filter.ShareGene.PatientShare.NotInNormal.brief.xls | tr '\\t' '\\n' | grep -n GeneName | cut -d : -f 1`
            tail -n +2 {BriefResults}/Share/snp.indel.filter.ShareGene.PatientShare.NotInNormal.brief.xls |
              cut -f $gene | tr , '\\n' > {BriefResults}/Share/snp.filter.ShareGene.PatientShare.brief.xls

            echo share compare done: `date "+%F %T"`
        '''.format(
            patients=','.join(self.patients), controls=controls, **self.args)

        shell_path = '{analydir}/Advance/{newjob}/Share/share_compare.sh'.format(**self.args)

        utils.write_shell(shell_path, cmd)

        # add job
        now_point = 'share_compare'
        job_name = 'share_compare'
        utils.add_job(self.jobs, now_point, self.args['startpoint'],
                      self.ANALYSIS_POINTS, job_name, shell_path, self.queues)

        # add order
        before_jobs = ['filter_snpindel']
        after_jobs = ['integrate_result', 'data_release']
        utils.add_order(
            self.orders,
            job_name,
            before_jobs=before_jobs,
            after_jobs=after_jobs)

# the end
