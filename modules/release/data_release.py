#!/usr/bin/env python
# -*- coding=utf-8 -*-
import os
import sys
import glob
import time
import gzip
import json
import commands
import chardet
from collections import defaultdict

import django
from django.template import Context, loader
from django.conf import settings

from arguments import get_args

RESULT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(os.path.dirname(RESULT_DIR))
sys.path.append(ROOT_DIR)

from config import config
import utils


# print utils.get_software_version()
# exit()

reload(sys)
sys.setdefaultencoding('utf-8')


class DataRelease(object):

    def __init__(self, args):
        self.args = args
        self.analy_list = args['analy_array'].split(',')
        self.qc_list = args['qc_list']
        # self.samp_info = args['samp_info']
        # self.samp_info_done = args['samp_info_done']
        self.outdir = args['odir'] or os.path.join(args['analydir'], 'Result', args['newjob'])

        if args['ER'] == 'Y':
            self.readme_dir = os.path.join(RESULT_DIR, 'readme', 'en')
        else:
            self.readme_dir = os.path.join(RESULT_DIR, 'readme', 'chs')

        # A set of tuples (sampleid, patientid)
        self.ANALY_DICT = utils.get_analysis_dict(self.analy_list, config.ANALYSIS_CODE)
        self.qc_lists = utils.Project.get_qc_lists(self.qc_list)
        self.softwares =  utils.get_softwares(self.analy_list, self.ANALY_DICT)

        self.__dict__.update(**args)

        # 记录高级分析的数目
        self.final_result_counter = 0

        # 传数据到模板
        self.context = defaultdict(dict)
        self.context.update(self.ANALY_DICT)

        self.context['advance_counter'] = {}

        self.context['software'] = dict(utils.get_software_version(), **self.softwares)

        # print self.context['software']
        # exit()

        # 旧版本(1.6?)和新版本有所区别：settings.configure, template.render
        self.django_old = False
        if django.VERSION < (1, 9):
            self.django_old = True

        # exit()

    def start(self):

        # # tree of release directory
        # {job}/ReleaseResult
        # ├── Data
        # │   ├── BamData
        # │   └── RawData
        # ├── FinalResult
        # │   ├── 1.FilterDB-...
        # │   ├── 2. ...
        # ├── PrimaryAnalysis
        # │   ├── FilterAnalysis
        # │   └── SampleVariation
        # └── Readme

        dir_counter = iter('一 二 三 四 五 六 七 八 九 十'.split())

        dir_map = {
            'CandidateGene': 'CandidateGene-候选基因列表',
            'FilterDB': 'FilterDB-突变位点筛选',
            'ACMG': 'ACMG-突变位点有害性分类',
            'FilterSV_CNV': 'FilterSV_CNV-结构变异有害性分析',
            'Noncoding': 'Noncoding-非编码区突变位点筛选',
            'ModelF': 'ModelF-基于家系样本筛选',
            'Share': 'Share-共有突变基因筛选',
            'Denovo': 'Denovo-新生突变筛选',
            'Linkage': 'Linkage-连锁分析',
            'ROH': 'ROH-纯合子区域分析',
            'Network': 'Network-候选基因相关性排序',
            'Pathway': 'Pathway-候选基因富集分析',
            'PPI': 'PPI-蛋白相互作用分析',
            'SiteAS': 'SiteAS-基于位点关联分析',
            'GeneAS': 'GeneAS-基于基因的关联分析',

            'BriefResults': '{analydir}/Advance/{newjob}/BriefResults'.format(**self.__dict__),
            'BriefResults': '{analydir}/Advance/{newjob}/BriefResults'.format(**self.__dict__),

            'Data': '{outdir}/ReleaseResult/Data'.format(**self.__dict__),
            'PrimaryAnalysis': '{outdir}/ReleaseResult/PrimaryAnalysis'.format(**self.__dict__),
            'FinalResult': '{outdir}/ReleaseResult/FinalResult'.format(**self.__dict__),
            'Readme': '{outdir}/ReleaseResult/Readme/'.format(**self.__dict__),

        }
        self.__dict__.update(**dir_map)

        # print self.ANALY_DICT.keys()
        # exit()
        utils.print_color('release {} samples:\n{}'.format(len(self.qc_lists), self.qc_lists.keys()))

        # RawData
        if True:
            self.release_fastq('raw')
            self.context['raw_data'] = dir_counter.next()

        # QC
        if self.ANALY_DICT['quality_control_keep_clean']:
            self.release_fastq('clean')
            self.context['clean_data'] = dir_counter.next()

        # Mapping
        if self.ANALY_DICT['mapping']:
            self.release_mapping()
            self.context['bam_data'] = dir_counter.next()

        # SNP/INDEL
        if self.ANALY_DICT['snpindel_call']:
            self.release_snp_indel('snp')
            self.release_snp_indel('indel')
            if 'primary_result' not in self.context:
                self.context['primary_result'] = dir_counter.next()
            self.context['snpindel_result'] = True

        # SV
        if self.ANALY_DICT['sv_call']:
            self.release_sv_cnv('sv')
            if 'primary_result' not in self.context:
                self.context['primary_result'] = dir_counter.next()
            self.context['sv_result'] = True

        # CNV
        if self.ANALY_DICT['cnv_call']:
            self.release_sv_cnv('cnv')
            if 'primary_result' not in self.context:
                self.context['primary_result'] = dir_counter.next()
            self.context['cnv_result'] = True

        # Circos
        if self.ANALY_DICT['cnv_call_freec'] and self.ANALY_DICT['sv_call']:
            self.release_circos()
            self.context['circos_result'] = True

        # Adcance

        if any((self.ANALY_DICT['filter_acmg'], self.ANALY_DICT['filter_model'], self.ANALY_DICT['share_compare'], self.ANALY_DICT['denovo'])):
            self.release_candidate_gene()
            self.context['candidate_gene'] = {'name': '{CandidateGene}'.format(**self.__dict__)}
            self.context['candidate_gene'].update({'counter': self.final_result_counter})

        # FilterDB
        if self.ANALY_DICT['filter_db'] or self.ANALY_DICT['filter_acmg']:
            self.final_result_counter += 1
            self.release_filter_db('snp')
            self.release_filter_db('indel')
            if 'filter_analysis' not in self.context:
                self.context['filter_analysis'] = dir_counter.next()
            self.context['filterdb'] = {'name': '{FilterDB}'.format(**self.__dict__)}
            self.context['filterdb'].update({'counter': self.final_result_counter})

        # ACMG
        if self.ANALY_DICT['filter_acmg']:
            self.final_result_counter += 1
            self.release_acmg()
            self.context['filter_acmg'] = {'name': '{ACMG}'.format(**self.__dict__)}
            self.context['filter_acmg'].update({'counter': self.final_result_counter})
            
        # FilterSV/FilterCNV
        if self.ANALY_DICT['filter_sv'] or self.ANALY_DICT['filter_cnv']:
            self.final_result_counter += 1
            if 'filter_analysis' not in self.context:
                self.context['filter_analysis'] = dir_counter.next()
            self.context['filter_sv_cnv'] = {'name': '{FilterSV_CNV}'.format(**self.__dict__)}
            self.context['filter_sv_cnv'].update({'counter': self.final_result_counter})
            if self.ANALY_DICT['filter_sv']:
                self.release_filter_sv_cnv('sv')
            if self.ANALY_DICT['filter_cnv']:
                self.release_filter_sv_cnv('cnv')

        # Noncoding
        if self.ANALY_DICT['filter_noncoding']:
            self.final_result_counter += 1
            self.release_filter_noncoding()
            self.context['filter_noncoding'] = {'name': '{Noncoding}'.format(**self.__dict__)}
            self.context['filter_noncoding'].update({'counter': self.final_result_counter})

        # ModelF
        if self.ANALY_DICT['filter_model'] and not self.ANALY_DICT['share_compare']:
            self.final_result_counter += 1
            self.release_filter_model()
            self.context['filter_model'] = {'name': '{ModelF}'.format(**self.__dict__)}
            self.context['filter_model'].update({'counter': self.final_result_counter})

        # Share
        if self.ANALY_DICT['share_compare']:
            self.final_result_counter += 1
            self.release_share_compare()
            self.context['share_compare'] = {'name': '{Share}'.format(**self.__dict__)}
            self.context['share_compare'].update({'counter': self.final_result_counter})

        # Denovo
        if self.ANALY_DICT['denovo']:
            self.final_result_counter += 1
            self.context['denovo'] = {'name': '{Denovo}'.format(**self.__dict__)}
            self.context['denovo'].update({'counter': self.final_result_counter})
            if any((self.ANALY_DICT['denovo_samtools'], self.ANALY_DICT['denovo_triodenovo'], self.ANALY_DICT['denovo_denovogear'])):
                self.release_denovo()

            if self.ANALY_DICT['denovo_sv']:
                self.release_denovo_sv_cnv('sv')

            if self.ANALY_DICT['denovo_cnv']:
                self.release_denovo_sv_cnv('cnv')

        # Linkage
        if self.ANALY_DICT['linkage']:
            self.final_result_counter += 1
            self.release_linkage()
            self.context['linkage'] = {'name': '{Linkage}'.format(**self.__dict__)}
            self.context['linkage'].update({'counter': self.final_result_counter})

        # ROH
        if self.ANALY_DICT['roh']:
            self.final_result_counter += 1
            self.release_roh()
            self.context['roh'] = {'name': '{ROH}'.format(**self.__dict__)}
            self.context['roh'].update({'counter': self.final_result_counter})

        # Network
        if self.ANALY_DICT['phenolyzer']:
            self.final_result_counter += 1
            self.release_network()
            self.context['network'] = {'name': '{Network}'.format(**self.__dict__)}
            self.context['network'].update({'counter': self.final_result_counter})

        # Pathway
        if self.ANALY_DICT['pathway']:
            self.final_result_counter += 1
            self.release_pathway()
            self.context['pathway'] = {'name': '{Pathway}'.format(**self.__dict__)}
            self.context['pathway'].update({'counter': self.final_result_counter})

        # PPI
        if self.ANALY_DICT['ppi']:
            self.final_result_counter += 1
            self.release_ppi()
            self.context['ppi'] = {'name': '{PPI}'.format(**self.__dict__)}
            self.context['ppi'].update({'counter': self.final_result_counter})

        # SiteAS
        if self.ANALY_DICT['site_association']:
            self.final_result_counter += 1
            # self.release_site_as()
            self.context['site_as'] = {'name': '{SiteAS}'.format(**self.__dict__)}
            self.context['site_as'].update({'counter': self.final_result_counter})

        # GeneAS
        if self.ANALY_DICT['gene_association']:
            self.final_result_counter += 1
            # self.release_gene_as()
            self.context['gene_as'] = {'name': '{GeneAS}'.format(**self.__dict__)}
            self.context['gene_as'].update({'counter': self.final_result_counter})

        self.context['final_result'] = dir_counter.next()
        self.context['appendix'] = dir_counter.next()

        # Readme
        self.make_readme()

    # ======================== Readme ===========================
    def django_configure(self):

        print 'configure django ...'
        if self.django_old:

            settings.configure(
                DEBUG=True,
                TEMPLATE_DEBUG=True,
                TEMPLATE_DIRS=(
                    os.path.join(RESULT_DIR, 'templates'),))
        else:
            settings.configure(
                DEBUG=True,
                TEMPLATE_DEBUG=True,
                TEMPLATES=[
                    {
                        'BACKEND': 'django.template.backends.django.DjangoTemplates',
                        'DIRS': [
                            RESULT_DIR,
                            os.path.join(RESULT_DIR, 'templates')
                        ]
                    }
                ]
            )
            django.setup()

    def make_readme(self):
        
        self.django_configure()

        title = open(self.args['pn']).read().strip()
        encoding = chardet.detect(title)['encoding']
        if encoding != 'utf8':
            title = title.decode(encoding)
        self.context['title'] = title

        # self.context['software'] = self.softwares

        src = os.path.join(RESULT_DIR, 'src')
        dest = '{Readme}'.format(**self.__dict__)
        self.link_data(src, dest)

        max_code = max(map(float, self.analy_list))
        if max_code < 2:
            report_type = 'qc'
        elif 2 <= max_code < 3:
            report_type = 'mapping'
        elif 3 <= max_code < 6.2:
            report_type = 'primary'
        elif max_code >= 6.2:
            report_type = 'advance'
        # print report_type
        self.context['report_type'] = report_type

        print json.dumps(self.context, ensure_ascii=False, indent=2)
        # print os.path.join(RESULT_DIR, 'templates')
        # template = loader.get_template('test.html')
        # template = loader.get_template('readme_template_chs.html')
        template = loader.get_template('index.html')
        if self.django_old:
            html = template.render(Context(self.context))
        else:
            html = template.render(self.context)
        # print html
        dest_html = os.path.join(dest, 'index.html')
        with utils.safe_open(dest_html, 'w') as out:
            out.write(html)


    def release_fastq(self, fq_type):

        print '> release {} ...'.format(fq_type)

        if fq_type == 'raw':
            data_dir = 'RawData'
        elif fq_type == 'clean':
            data_dir = 'CleanData'
        else:
            exit('error fq_type')

        for sample in self.qc_lists:
            md5_list = []

            dest = '{Data}/{data_dir}/{sample}/'.format(
                **dict(self.__dict__, **locals()))
            dest_md5 = '{Data}/{data_dir}/{sample}/MD5.txt'.format(
                **dict(self.__dict__, **locals()))

            for lane in self.qc_lists[sample]['lanes']:
                for read in (1, 2):
                    fastq = '{analydir}/QC/{sample}/{sample}_{novoid}_{flowcell_lane}_{read}.clean.fq.gz'.format(
                        sample=sample,
                        read=read,
                        analydir=self.analydir,
                        **lane)

                    if fq_type == 'raw':
                        fastq = fastq.replace('clean.fq.gz', 'fq.gz').replace('QC', 'RawData')

                    self.link_data(fastq, dest)
                    fastq_md5 = fastq + '.MD5.txt'

                    if utils.check_file_exists(fastq_md5):
                        md5_list.append(fastq_md5)
            if md5_list:
                self.cat_md5(md5_list, dest_md5)

    def release_mapping(self):

        print '> release mapping ...'

        for sample in self.qc_lists:
            dest = '{Data}/BamData/{sample}/'.format(
                **dict(self.__dict__, **locals()))
            bam = '{analydir}/Mapping/{sample}.{sample}/{sample}.final.bam'.format(
                **dict(self.__dict__, **locals()))
            bam_bai = bam + '.bai'
            self.link_data(bam, dest)
            self.link_data(bam_bai, dest)

        # link readme

    def release_snp_indel(self, mtype):

        print '> release {} ...'.format(mtype)

        if mtype not in ('snp', 'indel'):
            exit('error mtype')

        data_dir = mtype.upper()

        soft = self.softwares['mutation']

        for sample in self.qc_lists:
            dest = '{PrimaryAnalysis}/SampleVariation/{sample}/{data_dir}/'.format(
                **dict(self.__dict__, **locals()))
            vcf = '{analydir}/Mutation/{sample}.{soft}/{sample}.{soft}.{mtype}.vcf.gz'.format(
                **dict(self.__dict__, **locals()))
            anno = '{analydir}/Mutation/{sample}.{soft}/{sample}.{soft}.{mtype}.annovar.*_multianno.xls.gz'.format(
                **dict(self.__dict__, **locals()))
            self.link_data(vcf, dest)
            self.link_data(anno, dest)

    def release_sv_cnv(self, svtype):

        print '> release {} ...'.format(svtype)

        if svtype not in ('sv', 'cnv'):
            exit('error svtype')

        data_dir = svtype.upper()

        soft = self.softwares[svtype]

        for sample in self.qc_lists:
            dest = '{PrimaryAnalysis}/SampleVariation/{sample}/{data_dir}/'.format(
                **dict(self.__dict__, **locals()))
            gff = '{analydir}/SV/{sample}/{soft}/{sample}.{soft}.gff'.format(
                **dict(self.__dict__, **locals()))
            anno = '{analydir}/SV/{sample}/{soft}/{sample}.{soft}.*_multianno.xls'.format(
                **dict(self.__dict__, **locals()))
            gene_info = '{analydir}/SV/{sample}/{soft}/{sample}.{soft}.geneInfo.xls'.format(
                **dict(self.__dict__, **locals()))

            self.link_data(gff, dest)
            self.link_data(gene_info, dest)
            self.link_data(anno, dest)

    def release_circos(self):

        print '> release circos ...'

        soft = self.softwares['sv']

        for sample in self.qc_lists:
            dest = '{PrimaryAnalysis}/SampleVariation/{sample}/Circos/'.format(
                **dict(self.__dict__, **locals()))
            png = '{analydir}/SV/{sample}/freec/Circos/{sample}.png'.format(
                **dict(self.__dict__, **locals()))

            self.link_data(png, dest)

    # ========================================
    # =             Advance                  =
    # ========================================
    def release_candidate_gene(self):

        print '> release candidate gene ...'

        dest = '{FinalResult}/{final_result_counter}.{CandidateGene}/'.format(**dict(self.__dict__, **locals()))
        gene = '{BriefResults}/IntegrateResult/candidate_gene.xlsx'.format(**dict(self.__dict__, **locals()))

        self.link_data(gene, dest)


    def release_filter_db(self, mtype):

        print '> release filter database for {} ...'.format(mtype)

        if mtype not in ('snp', 'indel'):
            exit('error mtype')

        data_dir = mtype.upper()

        # Filter Result Normal
        dest = '{PrimaryAnalysis}/FilterAnalysis/FilterDB/{data_dir}/'.format(
            **dict(self.__dict__, **locals()))
        xls_files = '{analydir}/Advance/{newjob}/Merged_vcf/Filter/{data_dir}/{mtype}.merged.*.xls'.format(
            **dict(self.__dict__, **locals()))
        self.link_data(xls_files, dest)

        # Merged VCF
        dest = '{FinalResult}/{final_result_counter}.{FilterDB}/VCF/'.format(
            **dict(self.__dict__, **locals()))
        vcf = '{analydir}/Advance/{newjob}/Merged_vcf/VCF/{mtype}.merged.vcf.gz'.format(
            **dict(self.__dict__, **locals()))
        anno = '{analydir}/Advance/{newjob}/Merged_vcf/VCF/{mtype}.merged.annovar.*_multianno.xls.gz'.format(
            **dict(self.__dict__, **locals()))
        self.link_data(vcf, dest)
        self.link_data(anno, dest)

        # Filter Result Excel
        dest = '{FinalResult}/{final_result_counter}.{FilterDB}/Filter/'.format(
            **dict(self.__dict__, **locals()))
        brief_xlsx = '{BriefResults}/FilterDB/{mtype}.*.xlsx'.format(
            **dict(self.__dict__, **locals()))
        self.link_data(brief_xlsx, dest)

    def release_acmg(self):

        print '> release acmg ...'

        for familyid in self.listdir('{BriefResults}/ACMG'):
            dest = '{FinalResult}/{final_result_counter}.{ACMG}/{familyid}/'.format(
                **dict(self.__dict__, **locals()))

            brief_xlsx = '{BriefResults}/ACMG/{familyid}/{familyid}.*.xlsx'.format(
                **dict(self.__dict__, **locals()))
            self.link_data(brief_xlsx, dest)

    def release_filter_sv_cnv(self, svtype):

        print '> release filter {} ...'.format(svtype)

        SVTYPE = svtype.upper()
        soft = self.softwares[svtype]

        for sample in self.qc_lists:

            # for primary
            dest = '{PrimaryAnalysis}/FilterAnalysis//Filter{SVTYPE}/{sample}/'.format(
                **dict(self.__dict__, **locals()))

            # priority_xls = '{analydir}/Advance/{newjob}/Filter{SVTYPE}/{sample}/{sample}.{soft}*.priority.xls'.format(
            #     **dict(self.__dict__, **locals()))
            # self.link_data(priority_xls, dest)

            deleterious_xls = '{analydir}/Advance/{newjob}/Filter{SVTYPE}/{sample}/{sample}.LikelyDeleterious.{SVTYPE}.xls'.format(
                **dict(self.__dict__, **locals()))

            self.link_data(deleterious_xls, dest)

            # for advance
            dest = '{FinalResult}/{final_result_counter}.{FilterSV_CNV}/{sample}/'.format(
                **dict(self.__dict__, **locals()))

            brief_xlsx = '{BriefResults}/Filter{SVTYPE}/{sample}.LikelyDeleterious.{SVTYPE}.xlsx'.format(
                **dict(self.__dict__, **locals()))
            self.link_data(brief_xlsx, dest)

    def release_filter_noncoding(self):
    
        print '> release filter noncoding ...'

        for familyid in self.listdir('{BriefResults}/Noncoding'):
            dest = '{FinalResult}/{final_result_counter}.{Noncoding}/{familyid}/'.format(
                **dict(self.__dict__, **locals()))

            brief_xlsx = '{BriefResults}/Noncoding/{familyid}/{familyid}.*.xlsx'.format(
                **dict(self.__dict__, **locals()))
            self.link_data(brief_xlsx, dest)

    def release_filter_model(self):

        print '> release filter model ...'

        for familyid in self.listdir('{BriefResults}/ModelF'):
            dest = '{FinalResult}/{final_result_counter}.{ModelF}/{familyid}/'.format(
                **dict(self.__dict__, **locals()))

            brief_xlsx = '{BriefResults}/ModelF/{familyid}/{familyid}*.xlsx'.format(
                **dict(self.__dict__, **locals()))
            self.link_data(brief_xlsx, dest)

    def release_share_compare(self):

        print '> release share compare ...'

        dest = '{FinalResult}/{final_result_counter}.{Share}/'.format(
            **dict(self.__dict__, **locals()))

        brief_xlsx = '{BriefResults}/Share/snp.indel.*.xlsx'.format(
            **dict(self.__dict__, **locals()))
        self.link_data(brief_xlsx, dest)

    def release_denovo(self):
    
        print '> release denovo snpindel ...'

        for familyid in self.listdir('{BriefResults}/Denovo'):

            dest = '{FinalResult}/{final_result_counter}.{Denovo}/{familyid}/'.format(
                **dict(self.__dict__, **locals()))

            snp_xlsx = '{BriefResults}/Denovo/{familyid}/*snp*.xlsx'.format(
                **dict(self.__dict__, **locals()))
            indel_xlsx = '{BriefResults}/Denovo/{familyid}/*indel*.xlsx'.format(
                **dict(self.__dict__, **locals()))
            self.link_data(snp_xlsx, dest)
            self.link_data(indel_xlsx, dest)

    def release_denovo_sv_cnv(self, svtype):
    
        print '> release denovo {} ...'.format(svtype)

        if svtype == 'sv':
            denovo_type = 'denovoSV'
        elif svtype == 'cnv':
            denovo_type = 'denovoCNV'

        for familyid in self.listdir('{BriefResults}/Denovo'):
            dest = '{FinalResult}/{final_result_counter}.{Denovo}/{familyid}/'.format(
                **dict(self.__dict__, **locals()))

            brief_xlsx = '{BriefResults}/Denovo/{familyid}/{familyid}.*{denovo_type}*.xlsx'.format(
                **dict(self.__dict__, **locals()))
            self.link_data(brief_xlsx, dest)

    def release_linkage(self):

        print '> release linkage ...'

        for familyid in self.listdir('{BriefResults}/Linkage'):
            dest = '{FinalResult}/{final_result_counter}.{Linkage}/{familyid}/'.format(
                **dict(self.__dict__, **locals()))

            brief_xlsx = '{BriefResults}/Linkage/{familyid}/{familyid}.*.xlsx'.format(
                **dict(self.__dict__, **locals()))
            self.link_data(brief_xlsx, dest)

            png = '{analydir}/Advance/{newjob}/Linkage/{familyid}/*png'.format(
                **dict(self.__dict__, **locals()))
            self.link_data(png, dest)

    def release_roh(self):

        print '> release roh ...'

        for familyid in self.listdir('{BriefResults}/ROH'):
            dest = '{FinalResult}/{final_result_counter}.{ROH}/{familyid}/'.format(
                **dict(self.__dict__, **locals()))

            brief_xlsx = '{BriefResults}/ROH/{familyid}/{familyid}.roh.xlsx'.format(
                **dict(self.__dict__, **locals()))
            self.link_data(brief_xlsx, dest)

    def release_network(self):
    
        print '> release network ...'

        dest = '{FinalResult}/{final_result_counter}.{Network}/'.format(
            **dict(self.__dict__, **locals()))

        brief_xlsx = '{BriefResults}/Network/phenolyzer.xlsx {BriefResults}/Network/disgenet.xlsx'.format(
            **dict(self.__dict__, **locals())).split()

        self.link_data(brief_xlsx, dest)

    def release_pathway(self):

        print '> release pathway ...'

        dest = '{FinalResult}/{final_result_counter}.{Pathway}/'.format(
            **dict(self.__dict__, **locals()))

        brief_xlsx = '{BriefResults}/Pathway/pathway.xlsx'.format(
            **dict(self.__dict__, **locals()))

        self.link_data(brief_xlsx, dest)

    def release_ppi(self):

        print '> release ppi ...'

        dest = '{FinalResult}/{final_result_counter}.{PPI}/'.format(
            **dict(self.__dict__, **locals()))

        brief_xlsx = '{BriefResults}/PPI/PPI.xlsx'.format(
            **dict(self.__dict__, **locals()))

        self.link_data(brief_xlsx, dest)

    def listdir(self, path):

        return os.listdir(path.format(**self.__dict__))

    @staticmethod
    def cat_md5(md5_list, dest_md5):

        md5s = ' '.join(md5_list)
        cmd = 'cat {md5s} | unix2dos > {dest_md5}'.format(**locals())
        # print cmd
        os.system(cmd)

    @staticmethod
    def link_data(source, dest):

        if isinstance(source, list):
            sources = source
        elif '*' in source:
            sources = glob.glob(source)
        else:
            sources = [source]
        # print sources

        if not sources:
            print '[warn] no data from: {}'.format(source)

        for s in sources:
            if not os.path.exists(s):
                print '[warn] file not exists: {}'.format(s)

            else:
                dest_dir = os.path.dirname(dest)
                utils.mkdir_if_not_exists(dest_dir)

                cmd = 'ln -sf {} {}'.format(s, dest)

                os.system(cmd)


def main():

    args = get_args(config)
    DataRelease(args).start()


if __name__ == "__main__":

    main()
