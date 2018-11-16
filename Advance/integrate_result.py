#!/usr/bin/env python
# -*- coding=utf-8 -*-
import glob
import utils
from string import Template


class IntegrateResult(object):

    def __init__(self, args, jobs, orders, config):

        self.args = args
        self.jobs = jobs
        self.orders = orders
        self.queues = args.get('queues')
        self.ANALYSIS_POINTS = config.ANALYSIS_POINTS
        
        self.disease_db = config.CONFIG.get('database', 'disease_db')

    def start(self):

        self.integrate()

        if self.args['disease'] == 'Y':
            find_disease_id = False
            for disease_id in self.args['disease_ids'].split(';'):

                if glob.glob('{disease_db}/Disease_BackGround/{disease_id}__*txt'.format(**dict(self.__dict__, **locals()))):
                    find_disease_id = True
                    utils.print_color('> Disease Analysis', 'white')
                    self.disease(disease_id)
                    break
            if not find_disease_id:
                print '[warn] disease analysis will not to do, cause the disease "{}" not in database yet'.format(self.args['disease_name'])


    def integrate(self):

        print '>  integrate result ...'
        # write shell
        cmd = Template('''
            set -eo pipefail
            echo integrate result start: `date "+%F %T"`

            cd ${IntegrateResult}

            python ${ROOT_DIR}/modules/integrate/integrate.py \\
                --analydir ${analydir} \\
                --samp-info ${samp_info} \\
                --newjob ${newjob} \\
                --moduledir ${moduledir} \\
                --confidence N \\
                --analy-array ${array} \\
                --out .

            grep -v '\\.' total.candidate.gene.xls | tr  ',' '\\n' | sort -u > temp
            mv temp total.candidate.gene.xls

            awk -v OFS='\\t' 'NR==1 || $1!="GeneName" {print $1, $2, $3}' Integrate.xls  | sort -u > CandidateGene.xls
            python ${ROOT_DIR}/modules/brief/text2excel.py \\
                ${BriefResults}/IntegrateResult/candidate_gene.xlsx \\
                ${ROOT_DIR}/modules/brief/readme/candidate_gene.readme.xls \\
                CandidateGene.xls

            rm -f CandidateGene.xls

            echo integrate result done: `date "+%F %T"`
        ''')
        
        cmd = cmd.safe_substitute(
            array=','.join(map(str, self.args['analy_array'])), **self.args)

        shell_path = '{IntegrateResult}/integrate_result.sh'.format(
            **self.args)

        utils.write_shell(shell_path, cmd)

        # add job
        now_point = job_name = 'integrate_result'
        utils.add_job(self.jobs, now_point, self.args['startpoint'],
                      self.ANALYSIS_POINTS, job_name, shell_path, self.queues)

    def disease(self, disease_id):

        print '>  disease analysis ...'
        # write shell
        array = ','.join(map(str, self.args['analy_array']))
        cmd = '''
            set -eo pipefail
            echo disease analysis start: `date "+%F %T"`

            cd {Disease}

            python {moduledir}/Disease/disease.py \\
                -i {IntegrateResult}/Integrate.candidate.xls \\
                -o Integrate.disease.xls \\
                -id {disease_id} \\
                -enc utf8

            python {ROOT_DIR}/modules/brief/text2excel.py \\
                {BriefResults}/Disease/Integrate.disease.xlsx \\
                Integrate.disease.xls

            echo disease analysis done: `date "+%F %T"`
        '''.format(**dict(self.args, **locals()))

        shell_path = '{Disease}/disease_analysis.sh'.format(
            **self.args)

        utils.write_shell(shell_path, cmd)

        # add job
        now_point = job_name = 'disease_analysis'
        utils.add_job(self.jobs, now_point, self.args['startpoint'],
                      self.ANALYSIS_POINTS, job_name, shell_path, self.queues)

        # add order
        before_jobs = ['integrate_result']
        after_jobs = ['data_release']
        utils.add_order(self.orders, job_name, before_jobs=before_jobs, after_jobs=after_jobs)

# the end
