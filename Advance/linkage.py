#!/usr/bin/env python
# -*- coding=utf-8 -*-
from collections import defaultdict
import string
import utils


class Linkage(object):

    def __init__(self, args, jobs, orders, mutation_soft, sv_soft, cnv_soft,
                 denovo_soft, sample_infos_all, config, ANALY_DICT):

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

        self.sample_infos_all = sample_infos_all

        self.ANALY_DICT = ANALY_DICT

        self.ANALYSIS_POINTS = config.ANALYSIS_POINTS

        self.config = config
        self.__dict__.update(self.args)
        self.ref = args.get('ref')
        self.REF = 'hg19' if self.ref == 'b37' else self.ref
        self.__dict__.update(dict(config.CONFIG.items('genome_' + self.ref)))
        self.__dict__.update(dict(config.CONFIG.items('software')))

    def start(self):

        print '>  linkage start...'

        # 当sample_info中有指定data列为2或3的，则只做标有2或3的家系
        # 否则默认sample_info中符合的家系都做
        linkage_infos = dict(
            filter(lambda (k, v): v['data'] in ('1', '3'),
                   self.sample_infos_all.items()))
        linkage_families = set(map(lambda (k, v): v['familyid'], linkage_infos.items()))

        peds_raw = utils.sampleinfo2ped(self.sample_infos_all)

        linkage_peds = defaultdict(list)
        for ped in peds_raw:
            familyid = ped['familyid']
            # print familyid, ped

            if linkage_families and familyid not in linkage_families:
                continue

            linkage_peds[familyid].append(ped)

        if not linkage_peds:
            print '[error] no sample can do linkage analysis according to your sample_info, please check!'
            exit(1)

        # print linkage_peds
        # exit()
        for familyid, ped in linkage_peds.items():
            # print familyid, ped
            samples = [each['sampleid'] for each in ped]

            # 至少两代
            flag = False
            for sample in samples:
                if (self.sample_infos_all[sample]['pa'] in samples) or (
                    (self.sample_infos_all[sample]['ma'] in samples)):
                    flag = True
            if not flag:
                continue

            # print familyid, samples
            pedfile = '{analydir}/Advance/{newjob}/Linkage/{familyid}/{familyid}.ped'.format(
                **dict(self.__dict__, **locals()))
            wsfile = '{analydir}/Advance/{newjob}/Linkage/{familyid}/{familyid}.ws'.format(
                **dict(self.__dict__, **locals()))
            samples_with_data = self.write_ped_ws(pedfile, wsfile, ped)
            # print samples_with_data
            self.samtools_call_hapmap(familyid, samples_with_data)
            self.linkdatagen(familyid)
            self.merlin(familyid)

    def write_ped_ws(self, pedfile, wsfile, pedlist):

        samples_with_data = []

        with utils.safe_open(pedfile, 'w') as pf, utils.safe_open(wsfile, 'w') as wf:
            ws_count = []
            n = 0
            for ped in pedlist:
                sampleid = ped['sampleid']
                if self.sample_infos_all[sampleid]['data'] != '0':
                    samples_with_data.append(sampleid)
                    n += 1
                    ws_count.append(n)
                else:
                    ws_count.append(0)
                ped_text = '{familyid}\t{sampleid}\t{pa}\t{ma}\t{sex}\t{phenotype}\n'.format(
                    **ped)
                pf.write(ped_text)
            # print ws_count
            ws_text = ' '.join(map(str, ws_count)) + '\n'
            wf.write(ws_text)

        return samples_with_data

    def samtools_call_hapmap(self, familyid, samples_with_data):

        vcf_list = '{analydir}/Advance/{newjob}/Linkage/{familyid}/vcf_{familyid}.list'.format(
            **dict(self.__dict__, **locals()))
        with utils.safe_open(vcf_list, 'w') as out:
            for sampleid in samples_with_data:
                out.write('{}.vcf\n'.format(sampleid))

        for sampleid in samples_with_data:
            print '>    samtools call hapmap for', sampleid

            cmd = '''
                set -eo pipefail
                echo samtools call hapmap for {sampleid} start: `date "+%F %T"`

                cd {analydir}/Advance/{newjob}/Linkage/{familyid}

                samtoolsv0.1.19 mpileup \\
                    -d 10000 -C 50 -D -S -m 2 -F 0.02 -q 13 -Q 13 \\
                    -gf {reffasta} \\
                    -l {moduledir}/Linkage/annotHapMap2L.txt \\
                    {analydir}/Mapping/{sampleid}.{sampleid}/{sampleid}.final.bam |
                bcftools_lh view \\
                    -cg -t 0.5 \\
                    -> {sampleid}.vcf

                echo samtools call hapmap for {sampleid} done: `date "+%F %T"`
            '''.format(**dict(self.__dict__, **locals()))

            shell_path = '{analydir}/Advance/{newjob}/Linkage/{familyid}/samtools_call_hapmap_{sampleid}.sh'.format(
                **dict(self.__dict__, **locals()))

            utils.write_shell(shell_path, cmd)

            # add job
            now_point = 'samtools_call_hapmap'
            job_name = 'samtools_call_hapmap_{sampleid}'.format(
                **locals())
            utils.add_job(self.jobs, now_point, self.args['startpoint'],
                          self.ANALYSIS_POINTS, job_name, shell_path,
                          self.queues)

            # add order
            before_jobs = ['final_bam_{sampleid}'.format(**locals())]
            after_jobs = ['linkdatagen_{familyid}'.format(**locals())]
            utils.add_order(
                self.orders,
                job_name,
                before_jobs=before_jobs,
                after_jobs=after_jobs)

    def linkdatagen(self, familyid):

        print '>   linkdatagen for', familyid

        cmd = '''
            set -eo pipefail
            echo linkdatagen for {familyid} start: `date "+%F %T"`

            cd {analydir}/Advance/{newjob}/Linkage/{familyid}

            # vcf2linkdatagen: VCF  => BRLMM genotype file
            perl {moduledir}/Linkage/vcf2linkdatagen.pl \\
                -annotfile {moduledir}/Linkage/annotHapMap2.txt \\
                -pop CHB -mindepth 10 -missingness 0 \\
                -idlist vcf_{familyid}.list \\
                > {familyid}.brlmm
            
            # linkdatagen: generate datasets for linkage anallysis by MERLIN
            perl {moduledir}/Linkage/linkdatagen.pl \\
                -data m \\
                -pedfile {familyid}.ped \\
                -whichSamplesFile {familyid}.ws \\
                -callFile {familyid}.brlmm \\
                -annotFile {moduledir}/Linkage/annotHapMap2.txt \\
                -pop CHB -binsize 0.3 -prog me \\
                -outputDir {familyid}_HapMap2 \\
                > {familyid}_HapMap2.out

            echo linkdatagen for {familyid} done: `date "+%F %T"`
        '''.format(**dict(self.__dict__, **locals()))

        shell_path = '{analydir}/Advance/{newjob}/Linkage/{familyid}/linkdatagen_{familyid}.sh'.format(
            **dict(self.__dict__, **locals()))

        utils.write_shell(shell_path, cmd)

        # add job
        now_point = 'linkdatagen'
        job_name = 'linkdatagen_{familyid}'.format(
            **locals())
        utils.add_job(self.jobs, now_point, self.args['startpoint'],
                      self.ANALYSIS_POINTS, job_name, shell_path, self.queues)

    def merlin(self, familyid):

        print '>   merlin for', familyid

        cmd = '''
            set -eo pipefail
            echo merlin for {familyid} start: `date "+%F %T"`

            cd {analydir}/Advance/{newjob}/Linkage/{familyid}/{familyid}_HapMap2_merlin/genome

            # ========================================
            # === 1 merlin non-parametric analysis ===
            # ========================================
            merlin \\
                -d merlin_autosome_{familyid}.dat \\
                -p merlin_genome_{familyid}.ped \\
                -f merlin_genome_{familyid}.freq \\
                -m merlin_genome_{familyid}.map \\
                --prefix merlin_autosome_{familyid} \\
                --error

            pedwipe \\
                -d merlin_autosome_{familyid}.dat \\
                -p merlin_genome_{familyid}.ped \\
                -e merlin_autosome_{familyid}.err
            mv wiped.dat merlin_autosome_wiped_{familyid}.dat
            mv wiped.ped merlin_autosome_wiped_{familyid}.ped

            minx \\
                -d merlin_X_{familyid}.dat \\
                -p merlin_genome_{familyid}.ped \\
                -f merlin_genome_{familyid}.freq \\
                -m merlin_genome_{familyid}.map \\
                --prefix merlin_X_{familyid} \\
                --error

            pedwipe \\
                -d merlin_X_{familyid}.dat \\
                -p merlin_genome_{familyid}.ped \\
                -e merlin_X_{familyid}.err

            mv wiped.dat merlin_X_wiped_{familyid}.dat
            mv wiped.ped merlin_X_wiped_{familyid}.ped

            merlin \\
                -d merlin_autosome_wiped_{familyid}.dat \\
                -p merlin_autosome_wiped_{familyid}.ped \\
                --Bits 2 -f merlin_genome_{familyid}.freq -m merlin_genome_{familyid}.map \\
                --smallswap \\
                --megabytes:9999 \\
                --founders \\
                --markerNames \\
                --pairs \\
                --exp \\
                --pdf \\
                --tabulate \\
                --prefix merlin_autosome_{familyid} \\
                > merlin_autosome_{familyid}.out

            minx \\
                -d merlin_X_wiped_{familyid}.dat \\
                -p merlin_X_wiped_{familyid}.ped \\
                -f merlin_genome_{familyid}.freq \\
                -m merlin_genome_{familyid}.map \\
                --smallswap \\
                --megabytes:9999 \\
                --founders \\
                --markerNames \\
                --pairs \\
                --exp \\
                --pdf \\
                --tabulate \\
                --prefix merlin_X_{familyid} \\
                > merlin_X_{familyid}.out

            # ========================================
            # ===            2 merlin2R            ===
            # ========================================
            perl {moduledir}/Linkage/merlin2r.pl \\
                -id {familyid} \\
                -i . \\
                -o {analydir}/Advance/{newjob}/Linkage/{familyid}

            cd {analydir}/Advance/{newjob}/Linkage/{familyid}
            Rscript merlin_{familyid}_plot.R


            # ========================================
            # ===         3 merlin2excel           ===
            # ========================================
            perl {moduledir}//Linkage/merlin2Excel.pl \\
                -id {familyid} \\
                -lod 0.4 \\
                -i merlin_R_{familyid}-nonparametric.tbl \\
                -o .

            # rm -rf {familyid}_HapMap2_tables {familyid}_HapMap2_merlin

            # Brief Result
            echo generate brief results

            python {ROOT_DIR}/modules/brief/text2excel.py \\
                {BriefResults}/Linkage/{familyid}/{familyid}.LinkageAnalysis.xlsx \\
                {ROOT_DIR}/modules/brief/readme/linkage.readme.xls \\
                LinkageAnalysis_{familyid}.xls

            echo merlin for {familyid} done: `date "+%F %T"`
        '''.format(**dict(self.__dict__, **locals()))

        shell_path = '{analydir}/Advance/{newjob}/Linkage/{familyid}/merlin_{familyid}.sh'.format(
            **dict(self.__dict__, **locals()))

        utils.write_shell(shell_path, cmd)

        # add job
        now_point = 'merlin'
        job_name = 'merlin_{familyid}'.format(
            **locals())
        utils.add_job(self.jobs, now_point, self.args['startpoint'],
                      self.ANALYSIS_POINTS, job_name, shell_path, self.queues)

        # add order
        before_jobs = ['linkdatagen_{familyid}'.format(**locals())]
        after_jobs = ['data_release']
        utils.add_order(
            self.orders,
            job_name,
            before_jobs=before_jobs,
            after_jobs=after_jobs)

# the end
