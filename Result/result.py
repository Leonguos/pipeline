#!/usr/bin/env python
# -*- coding=utf-8 -*-
import utils


class Result(object):

    def __init__(self, args, jobs, orders, config):

        self.args = args
        self.jobs = jobs
        self.orders = orders
        self.queues = args.get('queues')
        self.ANALYSIS_POINTS = config.ANALYSIS_POINTS
        

        self.__dict__.update(args)

    def start(self):

        self.release()

    def release(self):

        # print '  data release ...'
        array=','.join(map(str, self.args['analy_array']))

        sample_info_done = ''
        if self.args['samp_info_done']:
            sample_info_done = '\\\n{:16}--samp_info_done {samp_info_done}'.format('', **self.args)

        cmd = '''
            set -eo pipefail
            echo data release start: `date "+%F %T"`

            python2 {ROOT_DIR}/modules/release/data_release.py \\
               --analydir {analydir} \\
               --qc_list {qc_list} \\
               --samp_info {samp_info} {sample_info_done} \\
               --pn {pn} \\
               --odir {analydir}/Result/{newjob} \\
               --newjob {newjob} \\
               --analy_array {array}

            echo data release done: `date "+%F %T"`
        '''.format(
            **dict(self.__dict__, **locals()))

        shell_path = '{analydir}/Result/{newjob}/release.{newjob}.sh'.format(**self.args)

        utils.write_shell(shell_path, cmd)

        # add job
        now_point = job_name = 'data_release'
        utils.add_job(self.jobs, now_point, self.args['startpoint'],
                      self.ANALYSIS_POINTS, job_name, shell_path, self.queues)

        # tar release
        cmd = '''
            set -eo pipefail
            echo tar release start: `date "+%F %T"`

            python2 {ROOT_DIR}/modules/release/tar_release.old.py \\
                --projdir {analydir} \\
                --analy_array {array} \\
                --odir {analydir}/Result/{newjob} \\
                --date {newjob} \\
                --pre {qcsuffix} \\
                --pn {pn} \\
                --mail {email} \\
                --yymail {yymail}

            echo tar release done: `date "+%F %T"`
        '''.format(
            **dict(self.__dict__, **locals()))

        shell_path = '{analydir}/Result/{newjob}/tar_release.{newjob}.sh'.format(**self.args)

        utils.write_shell(shell_path, cmd)
# the end
