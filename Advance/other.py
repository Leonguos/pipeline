#!/usr/bin/env python
# -*- coding=utf-8 -*-
from string import Template
import utils


class Other(object):

    def __init__(self, args, jobs, orders, config, disease):

        self.__dict__.update(args)
        self.args = args
        self.jobs = jobs
        self.orders = orders
        self.queues = args.get('queues')
        self.ANALYSIS_POINTS = config.ANALYSIS_POINTS

        self.__dict__.update(dict(config.CONFIG.items('software')))

    def phenolyzer(self):

        print '>  phenolyzer ...'
        # write shell
        if not self.args['disease_name']:
            print '[error] phenolyzer needs disease name in your sample_info'
            exit(1)

        cmd = '''
            set -eo pipefail
            echo phenolyzer start: `date "+%F %T"`

            cd {analydir}/Advance/{newjob}/Network

            # Phenolyzer
            python {moduledir}/Phenolyzer/phenolyzer-0.1.5/phenolyzer_pipe4.7.py \\
                --dir {analydir} \\
                --disease "{disease_name}" \\
                --genelist {analydir}/Advance/{newjob}/IntegrateResult/total.candidate.gene.xls \\
                --job {newjob}

            # DisGeNet
            python {moduledir}/DisGeNet/disgenet.py \\
                --id '{disease_ids}' \\
                --glist {analydir}/Advance/{newjob}/IntegrateResult/total.candidate.gene.xls \\
                --out_dir .

            # Brief Result
            echo generate brief results

            python {ROOT_DIR}/modules/brief/text2excel.py \\
                {BriefResults}/Network/phenolyzer.xlsx \\
                {ROOT_DIR}/modules/brief/readme/phenolyzer.readme.xls \\
                AllGene_list.xls \\
                CandidateGene_list.xls \\
                CandidateGene_score.xls

            python {ROOT_DIR}/modules/brief/text2excel.py \\
                {BriefResults}/Network/disgenet.xlsx \\
                {ROOT_DIR}/modules/brief/readme/disgenet.readme.xls \\
                DisGeNet_shared_gene.xls

            echo phenolyzer done: `date "+%F %T"`
        '''.format(**self.__dict__)

        shell_path = '{analydir}/Advance/{newjob}/Network/phenolyzer.sh'.format(
            **self.__dict__)

        utils.write_shell(shell_path, cmd)

        # add job
        now_point = job_name = 'phenolyzer'
        utils.add_job(self.jobs, now_point, self.args['startpoint'],
                      self.ANALYSIS_POINTS, job_name, shell_path, self.queues)

        # add order
        before_jobs = ['integrate_result']
        after_jobs = ['data_release']
        utils.add_order(self.orders, job_name, before_jobs=before_jobs, after_jobs=after_jobs)

    def ibd(self):

        if self.args['seqstrag'] == 'WGS':
            region = '-r 1'
        else:
            region = '-R {TR}'.format(**self.__dict__)

        print '>   IBD'
        cmd = '''
            set -eo pipefail

            echo IBD start: `date "+%F %T"`

            cd ${Merged_vcf}/IBD

            # extract region
            bcftools-1.6 view \\
                ${region} \\
                ../VCF/snp.merged.vcf.gz |
            awk '$5!~/*/' > snp.merged.bed.vcf

            # extract sample_info
            awk -F '\\t' -v OFS='\\t' '$1!~/^#/{print $2, $2, $1, $2}' \\
                ${samp_info} \\
                > sample.ped

            # plink 
            plink --vcf snp.merged.bed.vcf --double-id --update-ids sample.ped --make-bed -out plink

            plink --bfile plink --genome -out all
            plink --bfile plink --genome  --rel-check -out family

            # result
            awk -v OFS='\\t' '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' all.genome > all.IBD.xls
            awk -v OFS='\\t' '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' family.genome > family.IBD.xls

            ln -sf ${moduledir}/IBD/readme.txt IBD.readme.txt

            rm -f plink.* *.{log,nosex,genome} *.vcf *.ped

            echo IBD done: `date "+%F %T"`
        '''
        
        cmd = Template(cmd).safe_substitute(**dict(self.__dict__, **locals()))

        shell_path = '{Merged_vcf}/IBD/IBD.sh'.format(**self.__dict__)

        utils.write_shell(shell_path, cmd)

        # add job
        now_point = job_name = 'ibd'
        utils.add_job(self.jobs, now_point, self.args['startpoint'],
                      self.ANALYSIS_POINTS, job_name, shell_path, self.queues)

        # add order
        before_jobs = ['annotate_merged_snp']
        after_jobs = ['data_release']
        utils.add_order(self.orders, job_name, before_jobs=before_jobs, after_jobs=after_jobs)

    def pathway(self):

        print '>  pathway ...'
        # write shell
        cmd = '''
            set -eo pipefail
            echo pathway start: `date "+%F %T"`

            cd {analydir}/Advance/{newjob}/Pathway

            python {moduledir}/Enrich_R/bin/stat_phyper_table.py \\
                -i {analydir}/Advance/{newjob}/IntegrateResult/total.candidate.gene.xls \\
                -o Pathway \\
                -f 1

            # extract gene
            awk -F '\\t' 'NR>1 && $5<0.05 {print8}' Pathway/Pathway_kegg.xls |
                tr '/' '\\n' |
                sort -u > KEGG.xls

            # extract path
            awk -F '\\t' 'NR>1 && $5<0.05 {print1}' Pathway/Pathway_kegg.xls |
                sort -u >> KEGG.xls

            python {moduledir}/KEGG/kegg_svg.py KEGG.xls

            # Brief Result
            echo generate brief results

            python {ROOT_DIR}/modules/brief/text2excel.py \\
                {BriefResults}/Pathway/pathway.xlsx \\
                {ROOT_DIR}/modules/brief/readme/pathway.readme.xls \\
                Pathway/Pathway_go_MF.xls \\
                Pathway/Pathway_go_BP.xls \\
                Pathway/Pathway_go_CC.xls \\
                Pathway/Pathway_kegg.xls

            echo pathway done: `date "+%F %T"`
        '''.format(
            print1='{print $1}',
            print8='{print $8}',
            **self.args)

        shell_path = '{analydir}/Advance/{newjob}/Pathway/pathway_enrichment.sh'.format(
            **self.args)

        utils.write_shell(shell_path, cmd)

        # add job
        now_point = job_name = 'pathway'
        utils.add_job(self.jobs, now_point, self.args['startpoint'],
                      self.ANALYSIS_POINTS, job_name, shell_path, self.queues)

        # add order
        before_jobs = ['integrate_result']
        after_jobs = ['data_release']
        utils.add_order(self.orders, job_name, before_jobs=before_jobs, after_jobs=after_jobs)

    def ppi(self):

        print '>  ppi ...'
        # write shell
        cmd = '''
            set -eo pipefail
            echo ppi start: `date "+%F %T"`

            cd {PPI}

            echo 9606 > PPI_genes.xls

            cat {analydir}/Advance/{newjob}/IntegrateResult/total.candidate.gene.xls |
                tr ',' '\\n' |
                tr '\\n' '\\t' >> PPI_genes.xls

            echo -e "\\nall\\n20" >> PPI_genes.xls

            java -Xmx6G -jar {genemania_jar} QueryRunner \\
                --data {genemania_data} \\
                --out flat \\
                --results . \\
                PPI_genes.xls

            python {moduledir}/PPI/SplitPPI_Result.py .

            # Brief Result
            echo generate brief results

            python {ROOT_DIR}/modules/brief/text2excel.py \\
                {BriefResults}/PPI/PPI.xlsx \\
                {ROOT_DIR}/modules/brief/readme/ppi.readme.xls \\
                Gene_interactions.xls \\
                Gene_description.xls \\
                Networks.description.xls

            echo ppi done: `date "+%F %T"`
        '''.format(**self.__dict__)

        shell_path = '{PPI}/ppi.sh'.format(
            **self.args)

        utils.write_shell(shell_path, cmd)

        # add job
        now_point = job_name = 'ppi'
        utils.add_job(self.jobs, now_point, self.args['startpoint'],
                      self.ANALYSIS_POINTS, job_name, shell_path, self.queues)

        # add order
        before_jobs = ['integrate_result']
        after_jobs = ['data_release']
        utils.add_order(self.orders, job_name, before_jobs=before_jobs, after_jobs=after_jobs)

# the end
