#!/usr/bin/env python
# -*- coding=utf-8 -*-
import utils


class Report(object):

    def __init__(self, args, jobs, orders, config):

        self.args = args
        self.jobs = jobs
        self.orders = orders
        self.queues = args.get('queues')
        self.ANALYSIS_POINTS = config.ANALYSIS_POINTS
        

    def start(self):

        qc_array = filter(lambda x: x < 2, self.args['analy_array'])
        mapping_array = filter(lambda x: 2 <= x < 3, self.args['analy_array'])
        primary_array = filter(lambda x: 3 <= x < 6, self.args['analy_array'])
        advance_array = filter(lambda x: x >= 6, self.args['analy_array'])

        print '>  qc report start...'
        self.report('qc', qc_array)

        if mapping_array:
            print '>  mapping report start...'
            self.report('mapping', qc_array + mapping_array)

        if primary_array:
            print '>  primary report start...'
            self.report('primary', qc_array + mapping_array + primary_array)

        if advance_array:
            print '>  advance report start...'
            self.report('advance', self.args['analy_array'])


    def report(self, report_type, array):
        """
        report_type: qc, mapping, primary, advance
        """

        # print '  qc report...'
        # write shell
        cmd = '''
            set -eo pipefail
            echo {report_type} report start: `date "+%F %T"`

            python {ROOT_DIR}/modules/report/report.py \\
                --repsty {report_type} \\
                --projdir {analydir} \\
                --project {pn} \\
                --sample {samp_info} \\
                --seqsty {seqstrag} \\
                --date {newjob} \\
                --TR {TR} \\
                --mail {email} \\
                --pre {qcsuffix} \\
                --ref {ref} \\
                --Dtype {disease_type} \\
                --datastat {datastat} \\
                --WES_xten {WES_xten} \\
                --disease {disease} \\
                --pdf {pdf} \\
                --odir {analydir}/Report/{newjob}/{outdir} \\
                --analy_array {array}

            echo {report_type} report done: `date "+%F %T"`
        '''.format(
            report_type=report_type,
            array=','.join(map(str, array)),
            outdir='QC' if report_type == 'qc' else report_type.title(),
            **self.args)

        shell_path = '{analydir}/Report/{newjob}/{report_type}_report.sh'.format(
            report_type=report_type, **self.args)

        utils.write_shell(shell_path, cmd)

        # add job
        now_point = job_name = '{}_report'.format(report_type)
        utils.add_job(self.jobs, now_point, self.args['startpoint'],
                      self.ANALYSIS_POINTS, job_name, shell_path, self.queues)

        # add order
        if report_type == 'advance':
            before_jobs = ['data_release']
            after_jobs = []
            utils.add_order(
                self.orders,
                job_name,
                before_jobs=before_jobs,
                after_jobs=after_jobs)
        else:
            before_jobs = []
            after_jobs = ['data_release']
            utils.add_order(
                self.orders,
                job_name,
                before_jobs=before_jobs,
                after_jobs=after_jobs)
            
