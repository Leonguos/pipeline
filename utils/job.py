#!/usr/bin/env python
# -*- coding=utf-8 -*-
"""
Job associate tools
"""
import os
import textwrap
import utils


def make_job(name, status, memory, shell_path, queues, threads=None, moniter=False):
    """
    docstring here
        :param queues: a list of queue, eg. ['disease1.q', 'disease2.q']
    """

    # print queues

    qs = ' '.join(['-q ' + q for q in queues])
    sched = '-V -cwd -S /bin/bash {qs}'
    # multi threading for
    if threads and str(threads).isdigit():
        cpu = int(threads * 0.8)
        sched += ' -l p={cpu}'
    sched = sched.format(**locals())

    cmd = 'sh '
    # =============== 用于监控资源消耗情况，结果存于stderr中 ============================
    # %S: Total number of CPU-seconds that the process spent in kernel mode.
    # %U: Total number of CPU-seconds that the process spent in user mode.
    # %E: Elapsed real time (in [hours:]minutes:seconds).
    # %P: Percentage of the CPU that this job got, computed as (%U + %S) / %E.
    # %M: Maximum resident set size of the process during its lifetime, in Kbytes.
    # ================================================================================
    if moniter:
        cmd = r"/usr/bin/time '-f\nElapsed\t%Es\nCPU\t%P\nMemory\t%Mkb\n\nExit\t%x' sh "

    cmd += shell_path

    new_job = {
        'name': name,
        'status': status,
        'sched': sched,
        'memory': memory,
        'cmd': cmd
    }

    return new_job


def write_shell(shell_path, cmd):

    shell_dirname = os.path.dirname(shell_path)
    create_dir(shell_dirname)

    cmd = textwrap.dedent(cmd).lstrip('\n')

    with open(shell_path, 'w') as out:
        out.write(cmd)

    os.chmod(shell_path, 0744)


def create_dir(path):

    if not os.path.exists(path):
        os.makedirs(path)


def write_job(analydir, jobname, jobs, orders):

    job_dir = os.path.join(analydir, 'job')
    log_dir = os.path.join(analydir, 'log', jobname)
    create_dir(job_dir)
    create_dir(log_dir)
    log_path = os.path.join(job_dir, jobname)
    with open(log_path, 'w') as log:
        # write jobs
        for job in jobs:
            job_part = '''\
                job_begin
                  name {name}
                  memory {memory}
                  status {status}
                  sched_options {sched}
                  cmd_begin
                    {cmd}
                  cmd_end
                job_end\n
            '''.format(**job)

            log.write(textwrap.dedent(job_part))
        # write orders
        for order in orders:
            log.write(order + '\n')
        # write log_dir
        log.write('\nlog_dir %s\n' % log_dir)

    print 'successfully generated job, and you can start with: %s' % utils.color_text('sjm job/' + jobname, 'green', style='bright')


def get_status(now_point, start_point, ANALYSIS_POINTS):

    if not start_point:
        return 'waiting'

    now_point = ANALYSIS_POINTS[now_point][1]

    # for point in re.split(r',|;', start_point):

    start_point = ANALYSIS_POINTS[start_point][1]

    if now_point[0] < start_point[0]:  # 第1级不同，直接比较
        return 'done'
    elif now_point[0] > start_point[0]:
        return 'waiting'
    elif now_point[1] != start_point[1]:  # 第1级相同，第2级不同，返回'waiting'?
        return 'done'
    elif now_point[2] < start_point[2]:  # 第1,2级相同，比较第3级
        return 'done'
    else:
        return 'waiting'


def add_job(jobs,
            now_point,
            startpoint,
            ANALYSIS_POINTS,
            job_name,
            shell_path,
            queues,
            threads=None):

    status = get_status(now_point, startpoint, ANALYSIS_POINTS)
    memory = ANALYSIS_POINTS[now_point][0]

    job = make_job(job_name, status, memory, shell_path, queues, threads)

    jobs.append(job)

    return status


def add_order(orders, job_name, before_jobs=(), after_jobs=()):

    for before_job in before_jobs:
        orders.append('order {job_name} after {before_job}'.format(**locals()))

    for after_job in after_jobs:
        orders.append('order {after_job} after {job_name}'.format(**locals()))
