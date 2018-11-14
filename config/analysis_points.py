#!/usr/bin/env python
# -*- coding=utf-8 -*-

# 第1位表示优先级，数值大的在后
# 第2位表示分支，分支相同比较先后，分支不同不处理

ANALYSIS_POINTS = {
    # QC
    'md5_raw': ['1G', (0, 0, 0, 0)],

    'qc': ['1G', 0, 0, 0, 0)],

    'mapping': ['16G', (1, 0, 0, 0)],           # *
    'qc_check': ['1G', (1, 1, 0, 0)],           #
    'qc_report': ['1G', (1, 2, 0, 0)],          # 
    'hla': ['6G', (1, 3, 0, 0)],                # *

    # Mapping

    # HLA

    'mutation': ['3G', (1, 0, 2, 0)],           # *
    'filterdb': ['3G', (1, 0, 2, 1)],           # *

    'release': ['3G', (9, 0, 0, 0)],             # *
    'report': ['3G', (9, 0, 0, 0)],              # *


    'bwa_mem': ['16G', (2, 0, 0)],
    'sentieon_bwa_mem': ['10G', (2, 0, 0)],
    'gzip_md5_clean': ['1G', (2, 0, 1)],  
    'samtools_sort': ['8G', (2, 0, 1)],
    'sambamba_merge': ['4G', (2, 0, 2)],
    'sambamba_markdup': ['2G', (2, 0, 3)],
    # picard merge and markdup to be added
    'sentieon_markdup': ['4G', (2, 0, 3)],
    'sentieon_realign': ['4G', (2, 0, 4)],
    'sentieon_recal': ['4G', (2, 0, 5)],
    'gatk_realign': ['15G', (2, 0, 4)],
    'gatk_recal': ['15G', (2, 0, 5)],
    'final_bam': ['1G', (2, 0, 6)],
    'finalbam': ['1G', (2, 0, 6)],

}