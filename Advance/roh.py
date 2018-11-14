#!/usr/bin/env python
# -*- coding=utf-8 -*-
from collections import defaultdict
import utils


class ROH(object):

    def __init__(self, args, jobs, orders, sample_infos, mutation_soft, config):

        self.args = args
        self.jobs = jobs
        self.orders = orders

        self.sample_infos = sample_infos

        self.queues = args.get('queues')
        self.ANALYSIS_POINTS = config.ANALYSIS_POINTS

        self.familyids = utils.get_familyids(sample_infos)

        self.mutation_soft = mutation_soft

        self.__dict__.update(**self.args)
        self.ref = args.get('ref')
        self.__dict__.update(dict(config.CONFIG.items('genome_' + self.ref)))
        self.__dict__.update(dict(config.CONFIG.items('software')))

    def start(self):

        if self.args['softwares']['roh'] == 'h3m2':
            self.roh_h3m2()
        elif self.args['softwares']['roh'] == 'plink':
            self.roh_plink()

    
    def roh_h3m2(self):

        print '>  roh with H3M2 ...'
        # 筛选患者有，正常人没有的区间，忽略患病未知的样本
        # 每个家系至少有一个患者才做
        for familyid, sampleids in self.familyids.iteritems():
            normals = []
            patients = []
            for sampleid in sampleids:
                phenotype = self.sample_infos[sampleid]['phenotype']
                if phenotype.lower() in ('n', 'normal'):
                    normals.append(sampleid)
                elif phenotype.lower() in ('p', 'patient'):
                    patients.append(sampleid)
                else:
                    # skip the unknown samples
                    continue
            if not patients:
                print '[warn] no patient in family "{familyid}", skip the family'.format(**locals())
                continue

            # 1 对每个正常人和患者进行ROH分析
            for sampleid in normals + patients:
                cmd = '''
                    set -eo pipefail
                    echo roh for {sampleid} start: `date "+%F %T"`

                    cd {analydir}/Advance/{newjob}/ROH/{familyid}

                    # vcf to bed
                    bcftools-1.6 view -H \\
                        {analydir}/Mutation/{sampleid}.{mutation_soft}/{sampleid}.{mutation_soft}.snp.vcf.gz |
                    cut -f 1,2 > {sampleid}.snp.bed

                    if [ -d H3M2/{sampleid} ];then
                        rm -rf H3M2/{sampleid}
                    fi

                    # parsing
                    {h3m2tool_dir}/H3M2BamParsing.sh \\
                        {h3m2tool_dir} \\
                        {analydir}/Mapping/{sampleid}.{sampleid} \\
                        {sampleid}.final.bam \\
                        . \\
                        H3M2 {sampleid} \\
                        {reffasta} \\
                        {sampleid}.snp.bed

                    # analyze
                    {h3m2tool_dir}/H3M2Analyze.sh \\
                        {h3m2tool_dir} \\
                        . \\
                        H3M2 {sampleid} \\
                        {sampleid}.snp.bed \\
                        100000 0.1 0.1 5

                    # result
                    cut -f 1-3 H3M2/{sampleid}/Results/{sampleid}_1e+05_0.1_0.1_HomozygosityTableCall.bed \\
                        > {sampleid}_roh.txt

                    bedtools intersect \\
                        -a {sampleid}_roh.txt \\
                        -b {refgene_cds_bed} \\
                        -wo > {sampleid}_roh_anno

                    perl {moduledir}/ROH/ROH_anno_v0.2.pl \\
                        {sampleid}_roh_anno \\
                        {analydir}/Mutation/{sampleid}.{mutation_soft}/{sampleid}.{mutation_soft}.snp.annovar.hg19_multianno.xls.gz

                    echo roh for {sampleid} done: `date "+%F %T"`
                '''.format(**dict(self.__dict__, **locals()))

                shell_path = '{analydir}/Advance/{newjob}/ROH/{familyid}/roh_{familyid}_{sampleid}.sh'.format(
                    **dict(self.__dict__, **locals()))
                utils.write_shell(shell_path, cmd)

                # add job
                now_point = 'roh'
                job_name = 'roh_{familyid}_{sampleid}'.format(**locals())
                utils.add_job(self.jobs, now_point, self.args['startpoint'],
                              self.ANALYSIS_POINTS, job_name, shell_path, self.queues)

                # add order
                if self.mutation_soft in ('samtools', 'sentieon'):
                    before_jobs = ['bcftools_filter_{sampleid}'.format(**locals())]
                elif self.mutation_soft in ('gatk', ):
                    before_jobs = ['gatk_filter_{sampleid}'.format(**locals())]
                elif self.mutation_soft in ('mtoolbox', ):
                    before_jobs = ['mtoolbox_call_{sampleid}'.format(**locals())]

                after_jobs = [
                    'data_release', 'roh_share_{familyid}'.format(**locals())
                ]

                utils.add_order(
                    self.orders,
                    job_name,
                    before_jobs=before_jobs,
                    after_jobs=after_jobs)

            # 2 筛选患者共有，正常人没有的区间
            normals_roh = ' '.join(
                '{}_roh.txt'.format(sampleid) for sampleid in normals)
            if not normals:
                normals_all = '# no normals'
            elif len(normals) == 1:
                normals_all = 'cp -f {normals_roh} {familyid}.normal.roh.xls'
            else:
                normals_all = 'cat {normals_roh} > {familyid}.normal.roh.xls'
            normals_all = normals_all.format(**locals())


            patients_roh = ['{}_roh.txt'.format(sampleid) for sampleid in patients]
            if len(patients) == 1:
                patient_roh = patients_roh[0]
                patients_share = ['cp -f {patient_roh} {familyid}.patient.roh.xls'.format(**locals())]
            elif len(patients) == 2:
                patient1_roh, patient2_roh = patients_roh
                patients_share = ['bedtools intersect -a {patient1_roh} -b {patient2_roh} -f 1E-9 -r > {familyid}.patient.roh.xls'.format(**locals())]
            elif len(patients) >= 3:
                patient1, patient2 = patients[:2]
                patients_share = ['bedtools intersect -a {patient1}_roh.txt -b {patient2}_roh.txt -f 1E-9 -r > share.2.roh.xls'.format(**locals())]

                n = 2
                for patient in patients[2:]:
                    n_1 = n
                    n += 1
                    patients_share.append(
                        'bedtools intersect -a share.{n_1}.roh.xls -b {patient}_roh.txt -f 1E-9 -r > share.{n}.roh.xls'.
                        format(**locals()))
                patients_share.append(
                    'cp -f share.{n}.roh.xls {familyid}.patient.roh.xls'.format(**locals()))
            patients_share = '; '.join(patients_share)

            if not normals:
                compare = 'cp -f {familyid}.patient.roh.xls {familyid}_Family_roh.txt'
            else:
                compare = 'bedtools intersect -a {familyid}.patient.roh.xls -b {familyid}.normal.roh.xls -v > {familyid}_Family_roh.txt'
            compare_family = compare.format(**locals())

            cmd = '''
                set -eo pipefail
                echo roh share for {familyid} start: `date "+%F %T"`

                cd {analydir}/Advance/{newjob}/ROH/{familyid}

                # patients share
                {patients_share}

                # normals all
                {normals_all}
                
                # compare
                {compare_family}

                bedtools intersect \\
                    -a {familyid}_Family_roh.txt \\
                    -b {refgene_cds_bed} \\
                    -wo \\
                    > {familyid}_Family_roh_anno

                perl /NJPROJ2/DISEASE/share/Disease/ROH/ROH_anno_v0.2.pl \\
                    {familyid}_Family_roh_anno \\
                    {analydir}/Mutation/{sampleid}.{mutation_soft}/{sampleid}.{mutation_soft}.snp.annovar.hg19_multianno.xls.gz

                # rm -rf H3M2 *.bed *_anno *.txt share.*.xls

                # Brief Result
                echo generate brief results

                python {ROOT_DIR}/modules/brief/text2excel.py \\
                    {BriefResults}/ROH/{familyid}/{familyid}.roh.xlsx \\
                    {ROOT_DIR}/modules/brief/readme/roh.readme.xls \\
                    {familyid}_Family_roh_anno.xls \\
                    {familyid}_Family_roh_anno_snp.xls

                echo roh share for {familyid} done: `date "+%F %T"`
            '''.format(**dict(self.__dict__, **locals()))

            shell_path = '{analydir}/Advance/{newjob}/ROH/{familyid}/roh_share_{familyid}.sh'.format(
                **dict(self.__dict__, **locals()))
            utils.write_shell(shell_path, cmd)

            # add job
            now_point = 'roh_share'
            job_name = 'roh_share_{familyid}'.format(**locals())
            utils.add_job(self.jobs, now_point, self.args['startpoint'],
                          self.ANALYSIS_POINTS, job_name, shell_path,
                          self.queues)

            # add order
            if len(self.sample_infos) == 1:
                pass
            else:
                before_jobs = []
                for sampleid in normals + patients:
                    before_jobs += 'extract_annotation_snp_{sampleid} extract_annotation_indel_{sampleid}'.format(**locals()).split()
            after_jobs = ['data_release']

            utils.add_order(
                self.orders,
                job_name,
                before_jobs=before_jobs,
                after_jobs=after_jobs)

    def roh_plink(self):

        print '>  roh with plink ...'
        print 'roh with plink not supported yet...'
        exit(1)

# the end
