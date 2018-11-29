#!/usr/bin/env python
# -*- coding=utf-8 -*-
import os
import argparse


def get_args(config, utils, __version__):

    analysis_mapping = '\n'.join(
        '{}\t{}'.format(code, config.ANALYSIS_CODE[code][0])
        for code in config.ANALYSIS_CODE)
    # print analysis_mapping

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description="Human reseq pipeline for Disease",
        epilog='contact: suqingdong@novogene.com',
        version=__version__)

    parser.add_argument(
        '--pwd',
        metavar='Str',
        help='The path of analysis directory[defalut="."]',
        default=os.getcwd())
    parser.add_argument(
        '--ref',
        metavar='Str',
        help="The version of reference(choose from [%(choices)s], default=%(default)s)",
        choices=['b37', 'hg19', 'hg38', 'mm9', 'mm10'],
        default='b37')
    parser.add_argument(
        '--samp_list',
        metavar='File',
        help="The filename of sample_list\n"
             "The format is as follow:\n"
             "#Lane\tPatientID\tSampleID\tLibID\tNovoID\tIndex\tPath"
    )
    parser.add_argument(
        '--samp_info',
        metavar='File',
        help="The filename of sample_info\n"
        "The file has at least tab delimited fields(if do DNM or Merlin analysis,need Data(just for merlin),Pa,Ma additional):\n"
        "#FamilyID\tSampleID\tSEX\tNormal/Patient\tPN\tData\tPa\tMa\n"
        "format for content(tab delimited):./string,string,F/M/U,N/P/U,Project number,0/1,0/string,0/string\n'.'"
        "for no family info, for no data or no father/mother infomation.")
    parser.add_argument(
        '--samp_info_done',
        metavar='File',
        help="The sample in this file will not do primary analysis again")
    parser.add_argument(
        '--pn',
        metavar='File',
        help=
        'The file of Project name(pn) and contract number(cn) for Web report name[default=%(default)s]',
        default='pn.txt')
    parser.add_argument(
        '--seqstrag',
        metavar='Str',
        help=
        'Sequencing strategy (choose from [%(choices)s], default=%(default)s)',
        choices=['WES_ag', 'WES_illu', 'WGS', 'TS'],
        default='WES_ag')
    parser.add_argument(
        '--TR',
        metavar='File',
        help='The bed file for sequence region',
    )
    parser.add_argument(
        '--rmdup',
        metavar='Str',
        help=
        'Mark duplication reads in bam files (choice from [%(choices)s], default=%(default)s)',
        choices=['N', 'Y'],
        default='Y')
    parser.add_argument(
        '--analy_array',
        metavar='Str',
        help=
        'The anaylsis code list seperated by comma [default=%(default)s]\n{}'.
        format(analysis_mapping),
        default='1,2.1,3.1')
    parser.add_argument(
        '--rawdata',
        metavar='Float',
        type=float,
        help=
        'Required rawdata size for all sequence, measured by "G"[Default is 90G for WGS, 10G for WES, TS',
        default=0)
    parser.add_argument(
        '--depth',
        metavar='Float',
        type=float,
        help=
        'Required data depth for TR sequencing, measured by "x"[Default is 0x, 100x, 0x for WGS, WES, TS]',
        default=0)
    parser.add_argument(
        '--PE',
        metavar='Int',
        type=int,
        help='PE length of the sequence[default=%(default)s]',
        default=150)
    parser.add_argument(
        '--Q20',
        metavar='Float',
        type=float,
        help='Min tolerate Q20[default=%(default)s%%]',
        default=90)
    parser.add_argument(
        '--Q30',
        metavar='Float',
        type=float,
        help='Min tolerate Q30[default=%(default)s%%]',
        default=80)
    parser.add_argument(
        '--error',
        metavar='Float',
        type=float,
        help='Max tolerate error rate[default=%(default)s%%]',
        default=0.1)
    parser.add_argument(
        '--dup',
        metavar='Float',
        type=float,
        help='Max tolerate duplicate[default=%(default)s%%]',
        default=30)
    parser.add_argument(
        '--qcsuffix',
        metavar='Str',
        help='The suffix of qc_list[Required]')
    parser.add_argument(
        '--startpoint',
        metavar='Str',
        help="The start point to analysis if needed",
    )
    parser.add_argument(
        '--queues',
        metavar='Str',
        help='The specific queues for qsub [default=%(default)s]',
        default=config.CONFIG.get('resource', 'queues'))
    parser.add_argument(
        '--newjob',
        metavar='Str',
        help='The job file for SJM Job Manager [default=%(default)s]',
        default='{}_V{}.job'.format(utils.get_now_time('%Y%m%d'), __version__))
    parser.add_argument(
        '--callTR',
        help=
        'Whether only call the variation in the target region(TR) or not, default=%(default)s]',
        action='store_true',
        default=False)
    parser.add_argument(
        '--moduledir',
        metavar='Str',
        help='The directory where our modules are stored [default=%(default)s]',
        default=config.CONFIG.get('software', 'moduleDIR'))
    parser.add_argument(
        '--email', metavar='Str', help='The email for receiving report')
    parser.add_argument(
        '--yymail', metavar='Str', help='The email for yunying manager')
    # parser.add_argument(
    #     '--partQC',
    #     metavar='Float',
    #     type=float,
    #     help="Extract part reads for QC(GB)[default=%(default)s]",
    #     default=-1)
    # parser.add_argument(
    #     '--bwath',
    #     metavar='Str',
    #     help="The CPU number bwa-mem will use[default=%(default)s]",
    #     default='5')
    # parser.add_argument(
    #     '--reAnaly',
    #     metavar='Boolean',
    #     type=bool,
    #     help=
    #     "Whether some samples have replenished datas, calling results will refresh or not [True or False, default=%(default)s]",
    #     default=True)
    # parser.add_argument(
    #     '--jobstat',
    #     metavar='File',
    #     help="The job.statues file to run previously uncompleted tasks")
    # parser.add_argument(
    #     '--mergef',
    #     help=
    #     'Whether merge vcfs by family or not (not merge all sample in one vcf) [choice from [%(choices)s], default=%(default)s]',
    #     choices=['N', 'Y'],
    #     default='N')
    # parser.add_argument(
    #     '--ER',
    #     help='Whether generate English report or not',
    #     choices=['N', 'Y'],
    #     default='N')
    parser.add_argument(
        '--WES_xten',
        help=
        "Use WES by xten report template or not (choice from [N, Y], default=Y)",
        choices=['N', 'Y'],
        default='Y')
    parser.add_argument(
        '--MT',
        help="Call ",
        action='store_true')
    # parser.add_argument(
    #     '--sf',
    #     help=
    #     "Standard deleterious filter(6.1) for every sample (choice from [N, Y], default=N)",
    #     choices=['N', 'Y'],
    #     default='N')
    # parser.add_argument(
    #     '--ACMG_EearlyD',
    #     help='Eearly disease',
    #     choices=['T', 'F'],
    #     default='F')
    # parser.add_argument(
    #     '--ACMG_dmodel',
    #     help='Disease inheritance model in this family or sample',
    #     choices=['D', 'R', 'None'],
    #     default='None')
    # parser.add_argument(
    #     '--ACMG_CDS',
    #     help='Only class mutation in CDS (include splicing) or not,default=T',
    #     choices=['T', 'F'],
    #     default='T')
    # parser.add_argument(
    #     '--ACMG_type',
    #     help='snp/indel,requried',
    #     choices=['snp', 'indel', 'snp.indel'],
    #     default='snp.indel')
    # parser.add_argument(
    #     '--ACMG_dfreq', help='Disease freuency reported', default='None')
    parser.add_argument(
        '--disease_type',
        help="the type of contract,Monogenic(M) or Complex(C).",
        choices=['M', 'C', 'N'],
        default='N')
    # parser.add_argument(
    #     '--svd_value',
    #     help="the value of svd of CoNIFER parameter",
    #     default='10')
    # parser.add_argument(
    #     '--snpindel_to_advance',
    #     help='if analysis from snpindel,choices T',
    #     choices=['T', 'F'],
    #     default='F')
    # parser.add_argument(
    #     '--AStype',
    #     help="the type of Association_Analysis",
    #     choices=['allele', 'genotype'],
    #     default='allele')
    # parser.add_argument(
    #     '--ASDB',
    #     help=
    #     'choose the control database:NOVO or ExAC_ALL or ExAC_EAS or other,or you can specify the databse yourself(other),'
    #     'but should use with the DB_file',
    #     choices=['NOVO', 'ExAC_ALL', 'ExAC_EAS', 'other'],
    #     default='NOVO')
    # parser.add_argument(
    #     '--dbfile',
    #     help=
    #     'if control database is other(client provided),need to provide the allele count file for ASDB',
    #     default='None')
    # parser.add_argument(
    #     '--ASdb_N', help='the sample number of DB', default='60706')
    # parser.add_argument(
    #     '--no_check',
    #     help="if library is self-built-library,choices Y",
    #     choices=['N', 'Y'],
    #     default='N')
    # parser.add_argument(
    #     '--ploidy', help="多样本混样，需要采用gatk call，并设置此参数为样本数乘2", default='2')
    # parser.add_argument(
    #     '--bam_merge_and_rmdup_with',
    #     help="Choose the software for bam merge and remove dup",
    #     choices=['sambamba', 'picard'],
    #     default='sambamba')

    parser.add_argument(
        '--datastat',
        help="是否记录样本信息，默认记录；非正常项目、测试项目需要设置为N",
        choices=['N', 'Y'],
        default='Y')

    parser.add_argument(
        '--pdf',
        help='generate pdf report or not[default=%(default)s]',
        choices=['Y', 'N'],
        default='Y')

    parser.add_argument(
        '--disease',
        help='generate disease report or not[default=%(default)s]',
        choices=['Y', 'N'],
        default='N')

    parser.add_argument(
        '--confidence',
        help='mark confidence for integrate results or not[default=%(default)s]',
        choices=['Y', 'N'],
        default='N')

    parser.add_argument(
        '--hla-gene',
        help='the gene to do HLA typing for ATHLATES, default will do all genes')

    parser.add_argument(
        '-sps',
        '--show-startpoints',
        help='Show the available startpoints',
        action='store_true')

    parser.add_argument(
        '--software',
        help='speify software, like "roh=plink;", "aligment=sentieon;merge=picard"')

    args = vars(parser.parse_args())

    if all((args['samp_info'], args['samp_list'])):
        return args
    elif args['show_startpoints']:
        for k, v in sorted(
                config.ANALYSIS_POINTS.iteritems(), key=lambda (k, v): v[1]):
            print '{} -- {}'.format('.'.join(map(str, v[1])), k)
    else:
        parser.print_help()

        print '\033[1;31;42mYou need supply these two arguments at least: samp_info and samp_list\033[0m'
    exit(0)
