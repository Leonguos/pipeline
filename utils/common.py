#!/usr/bin/env python
# -*- coding=utf-8 -*-
"""
Some common tools for Disease pipeline
"""
import os
import sys
import gzip
import commands
from datetime import datetime
from collections import defaultdict
from ConfigParser import ConfigParser

try:
    import colorama
except ImportError:
    sys.path.append('/ifs/TJPROJ3/DISEASE/share/Software/python/site-packages')
    import colorama

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

colorama.init()


fore_map = {
    'red': colorama.Fore.RED,
    'green': colorama.Fore.GREEN,
    'blue': colorama.Fore.BLUE,
    'cyan': colorama.Fore.CYAN,
    'yellow': colorama.Fore.YELLOW,
    'white': colorama.Fore.WHITE,
    'black': colorama.Fore.BLACK,
}

back_map = {
    'red': colorama.Back.RED,
    'green': colorama.Back.GREEN,
    'blue': colorama.Back.BLUE,
    'cyan': colorama.Back.CYAN,
    'yellow': colorama.Back.YELLOW,
    'white': colorama.Back.WHITE,
    'black': colorama.Back.BLACK,
}

style_map = {
    'bright': colorama.Style.BRIGHT,
    'dim': colorama.Style.DIM,
    'normal': colorama.Style.NORMAL,
}

software_avail = {
    'qc': ['raw2clean', 'fastp'],
    'alignment': ['bwa', 'sentieon'],
    'sort': ['samtools', 'sambamba', 'sentieon'],
    'merge': ['sambamba', 'picard', 'sentieon'],
    'markdup': ['sambamba', 'picard', 'sentieon'],
    'recal': ['gatk', 'sentieon'],
    'mutation': ['samtools', 'gatk', 'sentieon'],
    'sv': ['lumpy', 'crest', 'breakdancer', 'breakmer'],
    'cnv': ['freec', 'cnvnator', 'conifer'],
    'roh': ['h3m2', 'plink'],
}


def color_text(text, fg='red', bg='', style=''):

    text = str(text)

    fore = fore_map.get(fg) or fore_map.get('red')
    back = back_map.get(bg)
    style = style_map.get(style)

    fore_reset = colorama.Fore.RESET
    back_reset = colorama.Back.RESET
    style_reset = colorama.Style.RESET_ALL

    text_new = '{fore}{text}{fore_reset}'
    if back:
        text_new = '{back}' + text_new + '{back_reset}'
    if style:
        text_new = '{style}' + text_new + '{style_reset}'

    return text_new.format(**locals())


def print_color(text, fg='red', bg='', style=''):

    print color_text(text, fg, bg, style)


def safe_open(filename, mode='r'):

    try:
        if mode == 'w':
            dirname = os.path.dirname(filename)
            if dirname and not os.path.exists(dirname):
                os.makedirs(dirname)
        if filename.endswith('.gz'):
            mode += 'b'
            return gzip.open(filename, mode)

        return open(filename, mode)
    except Exception as e:
        print '[error] fail to open file:', filename
        exit(1)


def get_now_time(fmt='%Y%m%d%H%M%S'):

    return datetime.now().strftime(fmt)


def get_index(name, title, name2=''):

    ntitle = [t.lower() for t in title]

    for n in (name.lower(), name2.lower()):
        # for each in ntitle:
        #     if n in each:
        #         return ntitle.index(each)
        if n in ntitle:
            return ntitle.index(n)


def get_value(linelist, idx, default=''):

    if idx is not None and idx < len(linelist):
        return linelist[idx]
    else:
        return default


def get_flowcell_lib(path, lib, index, lane):

    lane = lane[-1]
    fq1 = os.path.join(path, lib, '{}_L{}_1.fq.gz'.format(lib, lane))

    if not os.path.exists(fq1):
        if '-' not in lib:
            lib += '-' + index
        new_fq1 = os.path.join(path, lib, '{}_L{}_1.fq.gz'.format(lib, lane))
        if not os.path.exists(new_fq1):

            return '', lib

    with safe_open(fq1) as f:
        flowcell = f.readline().split(':')[2]

    return flowcell, lib


def get_chrom_list(refgenome, sex, MT=False):

    if refgenome.lower() in ('b37', ):
        chrom_list = [str(i) for i in range(1, 23)] + ['X']
        chrom_list += ['Y'] if sex.lower() in ('u', 'm', 'male') else []
        chrom_list += ['MT'] if MT else []
    elif refgenome.lower() in ('hg19', 'hg38'):
        chrom_list = ['chr' + str(i) for i in range(1, 23)] + ['chrX']
        chrom_list += ['Y'] if sex.lower() in ('u', 'm', 'male') else ['chrY']
        chrom_list += ['chrM'] if MT else []
    elif refgenome.lower() in ('mm9', 'mm10'):
        chrom_list = ['chr' + str(i) for i in range(1, 20)] + ['chrX']
        chrom_list += ['Y'] if sex.lower() in ('u', 'm', 'male') else ['chrY']
        chrom_list += ['chrM'] if MT else []
    else:
        print '[error] invalid refgenome: {}'.format(refgenome)
        exit(1)

    return chrom_list

def check_queues(queues, username):
    '''
    Check whether the queues are available for `qsub`, remove the unavailable queues
    '''
    print 'check queues...'
    cmd = 'qselect -U {} | cut -d @ -f 1 | sort -u'.format(username)
    unavailable_queues = commands.getoutput(cmd).split('\n')
    new_queues = []
    bad_queues = []
    for queue in queues:
        if queue not in unavailable_queues:
            bad_queues.append(queue)
        else:
            new_queues.append(queue)
    if not new_queues:
        print '[error] no available queues for qsub from {}'.format(queues)
        exit(1)
    if bad_queues:
        print '  unavailable queues: {}'.format(bad_queues)
    print '  used queues: {}'.format(new_queues)

    return new_queues


def check_analy_array(seqstrag, analy_array, ANALYSIS_CODE):

    print 'check analy_array...'
    # analysis_basic = set([int(str(each).split('.')[0]) for each in analy_array])
    analysis_basic = set([int(each) for each in analy_array])

    for analysis in analy_array:

        if analysis not in ANALYSIS_CODE:
            print '[error]: invalid analysis number: {}'.format(analysis)
            exit(1)
        for each in ANALYSIS_CODE[analysis][1]:
            if each not in analysis_basic:
                print '[error]: {} depends on {}'.format(analysis, ANALYSIS_CODE[analysis][1])
                exit(1)

    wgs_only = [4.1, 4.2, 4.3, 5.1, 5.2, 6.3, 8.5]
    if seqstrag in ('WES_ag', 'WES_illu', 'TS') and any(str(each) in analy_array for each in wgs_only):
        print '[error] WES/TS can not include {}'.format(wgs_only)
        exit(1)

    # if 3.4 in analy_array and 2.4 not in analy_array:
    #     print '[error] mutation with sentieon needs mapping with code 2.4'
    #     exit(1)


def check_files(pn, samp_info, samp_list):

    print 'Check files...'
    for each in (pn, samp_info, samp_list):
        if not os.path.isfile(each):
            print '[error] file not exists: "{}"'.format(os.path.basename(each))
            exit(1)


def check_target_region(CONFIG, seqstrag, refgenome, rawTR):

    print 'check target region...'

    if seqstrag == 'WES_illu':

        if refgenome not in ('hg19', 'b37'):
            print '  [error] WES_illu can only choose genome reference from [hg19, b37]'
            exit(1)

        default_TR = CONFIG.get('bed_wes_illu', refgenome)
        if rawTR is None:
            print '  [warn] no TR was supplied for reference {}\n  default TR was used: {}'.format(refgenome, default_TR)
            TR = default_TR
        elif rawTR not in [each[1] for each in CONFIG.items('bed_wes_illu')]:
            print '  [warn] wrong TR for reference {}: {}\n  default TR was used: {}'.format(refgenome, TR, default_TR)
            TR = default_TR
        else:
            TR = rawTR

    elif seqstrag == 'WES_ag':

        # print CONFIG.items('bed_wes_ag')

        if refgenome not in ('b37', 'hg19', 'hg38', 'mm9', 'mm10'):
            print '[error] WES_ag can only choose genome reference from [b37, hg19, hg38, mm9, mm10]'
            exit(1)

        if refgenome in ('b37', 'hg19'):

            default_TR = CONFIG.get('bed_wes_ag', refgenome + '_V6')
            if rawTR in ('V4', 'V5', 'V5_UTR', 'V6'):
                TR_version = rawTR
                TR = CONFIG.get('bed_wes_ag', refgenome + '_' + TR_version)
                print '  used {} TR for {}: {}'.format(TR_version, refgenome, TR)
            elif rawTR is None:
                TR = default_TR
                print '  [warn] no TR was supplied for reference {}\n  default V6 TR was used: {}'.format(refgenome, default_TR)
            elif rawTR not in [each[1] for each in CONFIG.items('bed_wes_ag')]:
                print '  [warn] wrong TR for reference {}: {}\n  default V6 TR was used: {}'.format(refgenome, rawTR, default_TR)
                TR = default_TR
            else:
                TR = rawTR

        elif refgenome == 'hg38':

            default_TR = CONFIG.get('bed_wes_ag', refgenome + '_V6')
            if rawTR in ('V5', 'V6'):
                TR_version = rawTR
                TR = CONFIG.get('bed_wes_ag', refgenome + '_' + TR_version)
                print '  used {} TR for {}: {}'.format(TR_version, refgenome, TR)
            elif rawTR is None:
                rawTR = default_TR
                print '  [warn] no TR was supplied for reference {}\n  default V6 TR was used: {}'.format(refgenome, default_TR)
            elif rawTR not in [each[1] for each in CONFIG.items('bed_wes_ag')]:
                print '  [warn] wrong TR for reference {}: {}\n  default V6 TR was used: {}'.format(refgenome, rawTR, default_TR)
                TR = default_TR
            else:
                TR = rawTR

        elif refgenome in ('mm9', 'mm10'):

            default_TR = CONFIG.get('bed_wes_ag', refgenome)
            if rawTR is None:
                print '  [warn] no TR was supplied for reference {}\n  default TR was used: {}'.format(refgenome, default_TR)
                TR = default_TR
            elif rawTR not in [each[1] for each in CONFIG.items('bed_wes_illu')]:
                print '  [warn] wrong TR for reference {}: {}\n  default TR was used: {}'.format(refgenome, rawTR, default_TR)
                TR = default_TR
            else:
                TR = rawTR

    elif seqstrag == 'WGS':

        if refgenome not in ('hg19', 'b37', 'hg38', 'mm9', 'mm10'):
            print '  [error] WGS can only choose genome reference from [hg19, b37, hg38]'
            exit(1)

        default_TR = CONFIG.get('bed_wgs', refgenome)
        if rawTR is None:
            print '  [info] WGS use default TR for {}: {}'.format(refgenome, default_TR)
            TR = default_TR
        elif rawTR not in [each[1] for each in CONFIG.items('bed_wes_illu')]:
            print '  [warn] wrong TR for reference {}: {}\n  default TR was used: {}'.format(refgenome, rawTR, default_TR)
            TR = default_TR
        else:
            TR = rawTR

    return TR


def get_softwares(analy_array, ANALY_DICT, args=None, seqstrag=None):

    # analy_array = map(str, analy_array)

    qc_soft = 'raw2clean'
    mutation_soft = sv_soft = cnv_soft = ''
    denovo_soft = []
    roh_soft = ''

    alignment_soft = sort_soft = merge_soft = markdup_soft = recal_soft = ''
    if ANALY_DICT['mapping']:
        alignment_soft = 'bwa'
        sort_soft = 'samtools'
        merge_soft = 'sambamba'
        markdup_soft = 'sambamba'
    if ANALY_DICT['mapping_with_sentieon']:
        alignment_soft = 'sentieon'
        sort_soft = 'sentieon'
        merge_soft = 'sambamba'
        markdup_soft = 'sentieon'

    # default: samtools
    if ANALY_DICT['snpindel_call']:
        mutation_soft = 'samtools'
    if ANALY_DICT['snpindel_call_gatk']:
        mutation_soft = 'gatk'
        recal_soft = 'gatk'
    if ANALY_DICT['snpindel_call_sentieon']:
        mutation_soft = 'sentieon'
        recal_soft = 'sentieon'
    if ANALY_DICT['snpindel_call_mtoolbox']:
        mutation_soft = 'mtoolbox'

    # default: lumpy
    if ANALY_DICT['sv_call']:
        sv_soft = 'lumpy'
    if ANALY_DICT['sv_call_breakdancer']:
        sv_soft = 'breakdancer'
    if ANALY_DICT['sv_call_breakmer']:
        sv_soft = 'breakmer'
    if ANALY_DICT['sv_call_crest']:
        sv_soft = 'crest'

    # default: WGS - freec, WES - conifer
    if ANALY_DICT['cnv_call']:
        if seqstrag in ('WGS', ):
            cnv_soft = 'freec'
        elif seqstrag in ('WES_ag', 'WES_illu'):
            cnv_soft = 'conifer'
    if ANALY_DICT['cnv_call_cnvnator']:
        cnv_soft = 'cnvnator'
    if ANALY_DICT['cnv_call_conifer']:
        cnv_soft = 'conifer'
    if ANALY_DICT['cnv_call_freec']:
        cnv_soft = 'freec'

    # default: samtools + triodenovo
    if ANALY_DICT['denovo_samtools']:
        denovo_soft.append('samtools')
    if ANALY_DICT['denovo_denovogear']:
        denovo_soft.append('denovogear')
    if ANALY_DICT['denovo_triodenovo']:
        denovo_soft.append('triodenovo')

    if (not denovo_soft) and ANALY_DICT['denovo']:
        denovo_soft += ['samtools', 'triodenovo']

    if ANALY_DICT['roh']:
        roh_soft = 'h3m2'

    softwares = {
        'qc': qc_soft,
        'alignment': alignment_soft,
        'sort': sort_soft,
        'merge': merge_soft,
        'markdup': markdup_soft,
        'recal': recal_soft,
        'mutation': mutation_soft,
        'sv': sv_soft,
        'cnv': cnv_soft,
        'denovo': denovo_soft,
        'roh': roh_soft
    }

    # eg.  alignment=sentieon;merge=picard;denovo=samtools,denovogear
    if args and args.get('software'):
        try:
            print 'specify software:', args['software']
            for k, soft in dict(each.split('=') for each in args['software'].split(';')).items():
                if ',' in soft and all(s in software_avail for s in soft.split(',')):
                    softwares.update({k: soft.split(',')})
                elif k in software_avail and soft in software_avail[k]:
                    softwares.update({k: soft})
        except Exception as e:
            print 'invalid software format, ignore'

    print softwares
    return softwares


def get_software_version(ini=None):

    default_ini = os.path.join(BASE_DIR, '../config', 'software_version.ini')

    ini = ini or default_ini

    conf = ConfigParser()

    with safe_open(ini) as f:
        conf.readfp(f)

    software_version = {section: dict(conf.items(section)) for section in conf.sections()}

    return software_version


def check_sample_info(sample_info):

    pass


def check_sample_list(sample_list):

    pass


def info2ped(info):

    sex_map = {
        'm': '1',
        'f': '2'
    }

    pheno_map = {
        'n': '1',
        'p': '2'
    }

    # print info

    sampleid = info.get('sampleid')
    familyid = info.get('familyid')
    pa = info.get('pa') or '0'
    ma = info.get('ma') or '0'
    sex = info.get('sex').lower()
    phenotype = info.get('phenotype').lower()

    sex = sex_map.get(sex, '0')
    phenotype = pheno_map.get(phenotype, '0')

    return familyid, pa, ma, sex, phenotype


def sampleinfo2ped(sample_infos):
    """
    Ped Format:
        1 FamilyID
        2 SampleID
        3 Pa
        4 Ma
        5 Sex             1=>male   2=>female  other=>unknown
        6 Patient/Normal  1=>normal 2=>patient other=>unknown
    """


    ped = []

    for sampleid, info in sample_infos.items():
        familyid, pa, ma, sex, phenotype = info2ped(info)

        # print sampleid, sample_infos[pa], sample_infos[ma]

        pa_context = {}
        if sample_infos.get(pa):
            pa_familyid, pa_pa, pa_ma, pa_sex, pa_phenotype = info2ped(sample_infos[pa])
            pa_context = {
                'familyid': familyid,
                'sampleid': pa,
                'pa': pa_pa,
                'ma': pa_ma,
                'sex': pa_sex,
                'phenotype': pa_phenotype
            }

        ma_context = {}
        if sample_infos.get(ma):
            ma_familyid, ma_pa, ma_ma, ma_sex, ma_phenotype = info2ped(sample_infos[ma])
            ma_context = {
                'familyid': familyid,
                'sampleid': ma,
                'pa': ma_pa,
                'ma': ma_ma,
                'sex': ma_sex,
                'phenotype': ma_phenotype
            }

        context = {
            'familyid': familyid,
            'sampleid': sampleid,
            'pa': pa,
            'ma': ma,
            'sex': sex,
            'phenotype': phenotype,
            'pa_context': pa_context,
            'ma_context': ma_context
        }

        ped.append(context)

    return ped

def get_analysis_dict(analy_array, ANALYSIS_CODE):

    ANALYSIS_DICT = {}

    for code in ANALYSIS_CODE:
        ANALYSIS_DICT[ANALYSIS_CODE[code][0]] = False

    if isinstance(analy_array, str):
        analy_array = analy_array.split(',')

    for code in analy_array:
        if isinstance(code, str):
            if '.' in code:
                code = float(code)
            else:
                code = int(code)
        ANALYSIS_DICT[ANALYSIS_CODE[code][0]] = True
        ANALYSIS_DICT[ANALYSIS_CODE[int(code)][0]] = True

    # print ANALYSIS_DICT
    return ANALYSIS_DICT


def get_merge_read(analydir, sampleid, lanes):

    read1_list = [
        '{analydir}/QC/{sampleid}/{sampleid}_{novoid}_{flowcell}_L{lane}_1.clean.fq'.
        format(sampleid=sampleid, analydir=analydir, **lane) for lane in lanes
    ]
    read1_list_gz = [read + '.gz' for read in read1_list]
    read2_list = [
        '{analydir}/QC/{sampleid}/{sampleid}_{novoid}_{flowcell}_L{lane}_2.clean.fq'.
        format(sampleid=sampleid, analydir=analydir, **lane) for lane in lanes
    ]
    read2_list_gz = [read + '.gz' for read in read2_list]

    if len(lanes) == 1:
        merge_read1 = 'ln -sf {reads} {sampleid}.R1.fastq'.format(
            sampleid=sampleid, reads=read1_list[0])
        merge_read1_gz = 'zcat {reads} > {sampleid}.R1.fastq'.format(
            sampleid=sampleid, reads=read1_list_gz[0])
        merge_read2 = 'ln -sf {reads} {sampleid}.R2.fastq'.format(
            sampleid=sampleid, reads=read2_list[0])
        merge_read2_gz = 'zcat {reads} > {sampleid}.R2.fastq'.format(
            sampleid=sampleid, reads=read2_list_gz[0])
    else:
        merge_read1 = 'cat {reads} > {sampleid}.R1.fastq'.format(
            sampleid=sampleid, reads=' '.join(read1_list))
        merge_read1_gz = 'zcat {reads} > {sampleid}.R1.fastq'.format(
            sampleid=sampleid, reads=' '.join(read1_list_gz))
        merge_read2 = 'cat {reads} > {sampleid}.R2.fastq'.format(
            sampleid=sampleid, reads=' '.join(read2_list))
        merge_read2_gz = 'zcat {reads} > {sampleid}.R2.fastq'.format(
            sampleid=sampleid, reads=' '.join(read2_list_gz))

    return read1_list, merge_read1, merge_read2, merge_read1_gz, merge_read2_gz


def get_familyids(sample_infos):

    familyids = defaultdict(list)

    for sampleid, infos in sample_infos.iteritems():
        familyid = infos['familyid']
        familyids[familyid].append(sampleid)

    return dict(familyids)


def mkdir_if_not_exists(path):

    if not os.path.exists(path):
        # print 'make a new dir:', path
        # os.makedirs(path, mode=0o777, exist_ok=False)
        os.makedirs(path)

def check_file_exists(infile, not_exists_ok=True):

    if not os.path.isfile(infile):
        if not_exists_ok:
            print '[warn] file not exists: {}'.format(infile)
        else:
            print '[error] file not exists: {}'.format(infile)
            exit()
        return False
    return True

# the end
