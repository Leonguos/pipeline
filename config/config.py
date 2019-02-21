#!/usr/bin/env python
# -*- coding=utf-8 -*-
import os
from ConfigParser import ConfigParser
from collections import OrderedDict


# Read config.ini to CONFIG object
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
config_default = os.path.join(BASE_DIR, 'config_default.ini')
CONFIG = ConfigParser()
CONFIG.read(config_default)

# CONFIG.get('software', 'database_dir')

# Override the default configration
config_nanjing = os.path.join(BASE_DIR, 'config_nanjing.ini')
config_tianjin = os.path.join(BASE_DIR, 'config_tianjin.ini')

if BASE_DIR.startswith('/NJPROJ') and os.path.exists(config_nanjing):
    print 'use configration: %s' % config_nanjing
    CONFIG.read(config_nanjing)
elif BASE_DIR.startswith('/ifs/TJPROJ') and os.path.exists(config_tianjin):
    print 'use configration: %s' % config_tianjin
    CONFIG.read(config_tianjin)
else:
    print 'use configration: %s' % config_default


# Set root dir
CONFIG.add_section('common')
CONFIG.set('common', 'root_dir', os.path.dirname(BASE_DIR))
# print CONFIG.get('common', 'root_dir')

# ANALYSIS_CODE[code] = 'name': [dependencies]
ANALYSIS_CODE = OrderedDict()

ANALYSIS_CODE[1]   = 'quality_control', []
ANALYSIS_CODE[1.1] = 'quality_control_no_clean', []
ANALYSIS_CODE[1.2] = 'quality_control_rm_clean', []
ANALYSIS_CODE[1.3] = 'quality_control_keep_clean', []

ANALYSIS_CODE[2]   = 'mapping', [1]
ANALYSIS_CODE[2.1] = 'mapping_with_default', [1]
ANALYSIS_CODE[2.2] = 'mapping_with_sentieon', [1]
# ANALYSIS_CODE[2.1] = 'mapping_no_gatk', [1]
# ANALYSIS_CODE[2.2] = 'mapping_with_gatk', [1]
# ANALYSIS_CODE[2.3] = 'mapping_with_gatk_no_realign', [1]
# ANALYSIS_CODE[2.4] = 'mapping_with_sentieon', [1]

ANALYSIS_CODE[3]   = 'snpindel_call', [1, 2]
ANALYSIS_CODE[3.1] = 'snpindel_call_samtools', [1, 2]
ANALYSIS_CODE[3.2] = 'snpindel_call_gatk', [1, 2]
ANALYSIS_CODE[3.3] = 'snpindel_call_sentieon', [1, 2]
ANALYSIS_CODE[3.4] = 'snpindel_call_mtoolbox', [1]  # for MT
# ANALYSIS_CODE[3.2] = 'snpindel_call_samtools_multi', [1, 2]
# ANALYSIS_CODE[3.3] = 'snpindel_call_gatk', [1, 2]
# ANALYSIS_CODE[3.4] = 'snpindel_call_sentieon', [1, 2]
# ANALYSIS_CODE[3.5] = 'snpindel_call_mtoolbox', [1]  # for MT

ANALYSIS_CODE[4]   = 'sv_call', [1, 2]
ANALYSIS_CODE[4.1] = 'sv_call_crest', [1, 2]
ANALYSIS_CODE[4.2] = 'sv_call_breakdancer', [1, 2]
ANALYSIS_CODE[4.3] = 'sv_call_breakmer', [1, 2]
ANALYSIS_CODE[4.4] = 'sv_call_lumpy', [1, 2]

ANALYSIS_CODE[5]   = 'cnv_call', [1, 2]
ANALYSIS_CODE[5.1] = 'cnv_call_freec', [1, 2]
ANALYSIS_CODE[5.2] = 'cnv_call_cnvnator', [1, 2]
ANALYSIS_CODE[5.3] = 'cnv_call_conifer', [1, 2]

ANALYSIS_CODE[6]   = 'filter', [1, 2, 3]
ANALYSIS_CODE[6.1] = 'filter_db', [1, 2, 3]
ANALYSIS_CODE[6.2] = 'filter_acmg', [1, 2, 3]
ANALYSIS_CODE[6.3] = 'filter_sv', [1, 2, 4]
ANALYSIS_CODE[6.4] = 'filter_cnv', [1, 2, 5]
ANALYSIS_CODE[6.5] = 'filter_noncoding', [1, 2, 3]
# ANALYSIS_CODE[6.6] = 'filter_drug', [1, 2, 3, 6]

ANALYSIS_CODE[7]   = 'filter_model', [1, 2, 3, 6]
ANALYSIS_CODE[7.1] = 'model_dominant', [1, 2, 3, 6]
ANALYSIS_CODE[7.2] = 'model_recessive', [1, 2, 3, 6]
ANALYSIS_CODE[7.3] = 'share_compare', [1, 2, 3, 6]

ANALYSIS_CODE[8]   = 'denovo', [1, 2]
ANALYSIS_CODE[8.1] = 'denovo_samtools', [1, 2]
ANALYSIS_CODE[8.2] = 'denovo_denovogear', [1, 2]
ANALYSIS_CODE[8.3] = 'denovo_triodenovo', [1, 2]
ANALYSIS_CODE[8.4] = 'denovo_sv', [1, 2, 4]
ANALYSIS_CODE[8.5] = 'denovo_cnv', [1, 2, 5]

ANALYSIS_CODE[9]   = 'linkage', [1, 2]
ANALYSIS_CODE[9.1] = 'merlinkage', [1, 2]

ANALYSIS_CODE[10]   = 'other', [1, 2]
ANALYSIS_CODE[10.1] = 'roh', [1, 2, 3]
ANALYSIS_CODE[10.2] = 'phenolyzer', [1, 2, 3, 6]
ANALYSIS_CODE[10.3] = 'pathway', [1, 2, 3, 6]
ANALYSIS_CODE[10.4] = 'ppi', [1, 2, 3, 6]
ANALYSIS_CODE[10.5] = 'site_association', [1, 2, 3, 6]
ANALYSIS_CODE[10.6] = 'gene_association', [1, 2, 3, 6]
ANALYSIS_CODE[10.7] = 'hla', [1]

# Analysis Points
# {
#   'analysis_point': [memory, (s1, s2, s3)]
# }
# 指定内存和分析顺序，用于startpoint判断
ANALYSIS_POINTS = {
    # QC
    'md5_raw':   ['1G', (1, 0, 0)],
    'qc':        ['1G', (1, 0, 0)],
    'qc_check':  ['1G', (1, 0, 1)],
    'qc_report': ['1G', (1, 0, 1)],

    # Mapping
    'mapping':          ['22G', (2, 0, 0)],           # *
    'bwa_mem':          ['22G', (2, 0, 0)],
    'sentieon_bwa_mem': ['24G', (2, 0, 0)],
    'gzip_md5_clean':   ['1G', (2, 0, 1)],  
    # 'samtools_sort':  ['8G', (2, 0, 1)],
    'sambamba_merge':   ['4G', (2, 0, 2)],
    'sambamba_markdup': ['2G', (2, 0, 3)],
    'sentieon_markdup': ['4G', (2, 0, 3)],

    # picard merge and markdup to be added
    # 'sentieon_realign': ['4G', (2, 0, 4)],
    # 'sentieon_recal': ['4G', (2, 0, 5)],
    # 'gatk_realign': ['15G', (2, 0, 4)],
    # 'gatk_recal': ['15G', (2, 0, 5)],

    'final_bam': ['1G', (2, 0, 6)],
    'finalbam':  ['1G', (2, 0, 6)],              # *

    # Alnstat
    'stat_depth':     ['1G', (2, 0, 2)],
    'stat_flag':      ['1G', (2, 0, 2)],
    'stat_uncover':   ['1G', (2, 0, 3)],
    'combine_stat':   ['1G', (2, 0, 3)],
    'mapping_check':  ['1G', (2, 0, 4)],
    'mapping_report': ['1G', (2, 0, 4)],

    # Mutation
    # call samtools
    'mutation':        ['3G', (3, 0, 0)],            # *
    'samtools_call':   ['3G', (3, 0, 0)],
    'bcftools_concat': ['1G', (3, 0, 1)],
    'bcftools_filter': ['1G', (3, 0, 2)],

    # call sentieon
    'sentieon_realign':    ['1G', (3, 0, 0)],
    'sentieon_recal':      ['1G', (3, 0, 1)],
    'sentieon_hc_call':    ['4G', (3, 0, 2)],
    'sentieon_joint_call': ['1G', (3, 0, 3)],
    'sentieon_vqsr':       ['1G', (3, 0, 4)],

    # call gatk
    # 'gatk_recal':       ['16G', (3, 0, 0)],
    'gatk_hc_call':     ['15G', (3, 0, 0)],
    'gatk_concat':     ['15G', (3, 0, 1)],
    'gatk_consolidate': ['15G', (3, 0, 2)],
    # 'gatk_joint_call':  ['4G', (3, 0, 3)],
    'gatk_vqsr':        ['15G', (3, 0, 3)],

    # annotation or merge
    'annotate_vcf':        ['3G', (3, 0, 3)],
    'bcftools_merge':      ['1G', (3, 0, 3)],
    'annotate_merged_vcf_snp': ['8G', (3, 0, 4)],
    'annotate_merged_vcf_indel': ['6G', (3, 0, 4)],
    'extract_annotation':  ['1G', (3, 0, 5)],
    'variation_summary':   ['1G', (3, 0, 6)],

    # SV
    'sv':                 ['5G', (3, 2, 0)],                   # *
    'crest_call':         ['5G', (3, 2, 0)],
    'crest_txt2gff':      ['1G', (3, 2, 1)],
    'breakdancer_config': ['1G', (3, 2, 0)],
    'breakdancer_call':   ['2G', (3, 2, 1)],
    'annotate_gff':       ['2G', (3, 2, 5)],
    'lumpy_call':         ['3G', (3, 2, 0)],

    # CNV
    'cnv':          ['5G', (3, 3, 0)],                   # *
    'conifer_call': ['5G', (3, 3, 0)],
    'freec_call':   ['5G', (3, 3, 0)],

    # Circos
    'circos_draw':    ['5G', (4, 0, 0)],
    'primary_report': ['1G', (4, 0, 0)],

    # FilterDB
    'filterdb':         ['1G', (5, 0, 0)],               # *
    'filter_snpindel':  ['1G', (5, 0, 0)],
    'filter_acmg':      ['1G', (5, 1, 0)],
    'filter_noncoding': ['1G', (5, 1, 0)],
    'filter_model':     ['1G', (5, 2, 0)],
    'share_compare':    ['1G', (5, 3, 0)],

    'integrate_result': ['1G', (7, 0, 0)],
    'disease_analysis': ['1G', (7, 0, 1)],

    # Denovo
    'denovo':               ['4G', (5, 4, 0)],                 # *
    'samtools_call_trio':   ['4G', (5, 4, 0)],
    'bcftools_concat_trio': ['1G', (5, 4, 1)],
    'denovo_call':          ['2G', (5, 4, 2)],
    'denovo_annotate':      ['3G', (5, 4, 3)],
    'denovo_intersect':     ['1G', (5, 4, 4)],
    'denovo_rate':          ['1G', (5, 4, 5)],
    'denovo_sv':            ['2G', (5, 4, 5)],
    'denovo_cnv':           ['2G', (5, 4, 5)],
    'filter_sv':            ['1G', (3, 2, 6)],
    'filter_cnv':           ['1G', (3, 3, 6)],

    # Linkage
    'linkage':              ['1G', (3, 4, 0)],                 # *
    'samtools_call_hapmap': ['1G', (3, 4, 0)],
    'linkdatagen':          ['3G', (3, 4, 1)],
    'merlin':               ['4G', (3, 4, 2)],

    # Other
    'ibd': ['1G', (5, 5, 0)],

    'phenolyzer': ['1G', (6, 0, 1)],
    'pathway':    ['1G', (6, 0, 1)],
    'ppi':        ['1G', (6, 0, 1)],

    # AS
    'site_association': ['1G', (3, 0, 5)],
    'gene_as_filter':   ['1G', (5, 0, 1)],
    'gene_association': ['1G', (5, 0, 2)],

    # Release and Report
    'data_release':   ['1G', (7, 0, 3)],
    'advance_report': ['1G', (7, 0, 4)],

    # HLA
    'hla':                  ['6G', (6, 5, 0)],      # *
    'hla_bwa_mem':          ['6G', (6, 5, 0)],
    'hla_sambamba_merge':   ['6G', (6, 5, 1)],
    'hla_sambamba_markdup': ['6G', (6, 5, 2)],
    'hla_picard_markdup':   ['6G', (6, 5, 2)],
    'hla_sort_by_name':     ['6G', (6, 5, 3)],
    'hla_athlates_typing':  ['6G', (6, 5, 4)],
    
    'hla_hlahd_typing':     ['2G', (6, 6, 0)],

    # ROH
    'roh':       ['4G', (5, 4, 0)],
    'roh_share': ['1G', (5, 4, 1)],

    # MT
    'mtoolbox_call': ['6G', (2, 3, 0)],
}

ANALYSIS_MEM_WGS = {
    'qc': '18G',
    'mapping': '20G',
    'bwa_mem': '20G',
    'sentieon_bwa_mem': '20G',
    'gatk_recal': '15G',
    'gatk_hc_call': '15G',
    'gatk_vqsr': '15G',
    'annotate_merged_vcf_snp': '15G',
    'annotate_merged_vcf_indel': '10G',
}

# the end
