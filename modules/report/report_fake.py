#!/usr/bin/env python
# -*- coding=utf-8 -*-
import os
import re
import sys
import gzip
import time
import shutil
import string
import glob
import textwrap
import argparse
import numpy as np
from itertools import islice
import socket

from jinja2 import Environment, FileSystemLoader


REPORT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(os.path.dirname(REPORT_DIR))
sys.path.append(ROOT_DIR)

from arguments import get_args
from config import config
import utils

__version__ = '4.7'
__author__ = 'suqingdong'

disease_db = config.CONFIG.get('database', 'disease_db')

templates_dir = [
    os.path.join(REPORT_DIR, 'templates'),
    os.path.join(disease_db, 'Disease_BackGround')
]

env = Environment(loader=FileSystemLoader(templates_dir))


moduledir = config.CONFIG.get('software', 'moduledir')

# 用于解决中文字符写入问题
reload(sys)
sys.setdefaultencoding('utf-8')


report_host = '192.168.20.19:8086'
remote_dir = '/VM/Report'


args = get_args(config)

proj_n = args['project']
projdir = args['projdir']
mail = args['mail']
newjob = args['date']
TR = args['TR']
NGC = args['NGC']
ref = args['ref']
hpa_mode = args['HPAmode']
analy_array = args['analy_array'].split(',')
seq_ty = args['seqsty']
rep_ty = args['repsty']
odir = os.path.abspath(args['odir'])
english = args['english']
WES_xten = args['WES_xten']
sf = args['sf']
#pheno = args['network']
suffix = args['pre']
datastat = args['datastat']
Dtype = args['Dtype']

pdf_report = args['pdf']

# use for transport data to template
context = {}

ANALY_DICT = utils.get_analysis_dict(analy_array, config.ANALYSIS_CODE)
context.update(ANALY_DICT)

context['software'] = dict(utils.get_software_version(), **utils.get_softwares(analy_array, ANALY_DICT))
# print context
# exit()


ANALYSIS = []
for code in analy_array:
    if '.' in code:
        code = float(code)
    else:
        code = int(code)
    ANALYSIS.append(config.ANALYSIS_CODE[code][0])
    ANALYSIS.append(config.ANALYSIS_CODE[int(code)][0])
# print ANALYSIS
# exit()

## for out dir
if not os.path.exists(odir):
    os.mkdir(odir)

def txts2tab(infiles, start=0, mv=[]):
    i = 0
    tables = []
    for file in infiles:
        A = utils.safe_open(file, 'r')
        j = 0
        if i == 0:
            for n in A:
                if n.startswith('##'):
                    continue
                if len(mv) == 0:
                    tables.append(n.strip().split('\t'))
                else:
                    b = n.strip().split('\t')
                    for k in mv:
                        del b[k - 1]
                    tables.append(b)
        else:
            for n in A:
                tables[j] = tables[j] + n.strip().split('\t')[start:]
                j += 1
        A.close()
        i += 1
    return tables


def svtxts2tab(infiles):
    i = 0
    tables = []
    for file in infiles:
        A = utils.safe_open(file, 'r')
        j = 0
        t = []
        if i == 0:
            for n in A:
                if j == 0:
                    tables.append(n.strip().strip('#').split('\t'))
                    j += 1
                else:
                    if j == 1:
                        t.append(n.strip().split('\t'))
                        j += 1
                    else:
                        t.append(n.strip().split('\t')[1:])
        else:
            for n in A:
                j += 1
                if n.startswith('#'):
                    continue
                if j == 2:
                    t.append(n.strip().split('\t'))
                else:
                    t.append(n.strip().split('\t')[1:])
        tables.append(t)
        A.close()
        i += 1
    return tables


def getta(table, samlist):
    tableu = []
    h = 0
    for i in table:
        h += 1
        if h == 1:
            tableu.append(i)
        else:
            if i[0] in samlist:
                tableu.append(i)
    return tableu


def gettainfo(table, samlist):
    # print samlist
    # print table
    tableu = []
    totals = {}
    for i in table:
        if i[1] not in totals:
            totals[i[1]] = i
    for i in samlist:
        # print i
        if totals.get(i):
            tableu.append(totals[i])
        else:
            exit('%s not in sample_info' % i)
    return tableu

with open(proj_n) as f:
    pn = f.next().strip()
    try:
        context['title'] = pn.decode('utf8')
    except:
        context['title'] = pn.decode('gbk').encode('utf8')

context['datet'] = utils.get_now_time('%Y年%m月%d日')

if Dtype == 'M':
    context['typem'] = True
elif Dtype == 'C':
    context['typec'] = True

##
if rep_ty == 'qc':
    context['rep_qc'] = True
    context['samp_info'] = True
    context['analy_flow'] = True
elif rep_ty == 'mapping':
    context['rep_mapping'] = True
    context['samp_info'] = True
    context['analy_flow'] = True
    context['dep_res'] = True
elif rep_ty == 'primary':
    context['rep_primary'] = True
    context['samp_info'] = True
    context['analy_flow'] = True
    context['dep_res'] = True
    context['muts_res'] = True

elif rep_ty == 'advance':
    context['rep_advance'] = True
    context['samp_info'] = True
    context['analy_flow'] = True
    context['dep_res'] = True
    context['muts_res'] = True

if 'snpindel_call' in ANALYSIS:
    context['snp_indel_res'] = True

if 'sv_call' in ANALYSIS:
    context['sv_res'] = True

if 'cnv_call' in ANALYSIS:
    context['cnv_res'] = True

if 'filter_acmg' in ANALYSIS:
    ANALYSIS.append('filter_db')

if 'filter_db' in ANALYSIS:
    context['filter_res'] = True
    advanceDir = os.path.join(projdir, 'Advance', newjob)
    advanceBriDir = os.path.join(projdir, 'Advance', newjob, 'BriefResults')
    gene_dir = os.path.join(advanceDir, 'IntegrateResult', 'total.candidate.gene.xls')

if 'filter_model' in ANALYSIS:
    context['model_res'] = True

if 'denovo' in ANALYSIS:
    context['denovo_res'] = True

if 'linkage' in ANALYSIS:
    context['linkage_res'] = True

if 'filter_noncoding' in ANALYSIS:
    context['Noncoding'] = True

if 'roh' in ANALYSIS:
    context['roh_res'] = True

if 'pathway' in ANALYSIS:
    context['Pathway'] = True

context[seq_ty] = True

if 'denovo_sv' in ANALYSIS:
    context['DenovoSV'] = True

if 'denovo_cnv' in ANALYSIS:
    context['DenovoCNV'] = True

if 'ppi' in ANALYSIS:
    context['ppinter'] = True

if 'hpa' in ANALYSIS:
    context['HPA'] = True

if 'site_association' in ANALYSIS:
    context['asso_site'] = True

if 'filter_drug' in ANALYSIS:
    context['rep_person'] = True

infile = os.path.join(projdir, 'qc_list_' + suffix)

context['table_qc'] = []
context['table_qcb'] = []

A = utils.safe_open(infile, 'r')
samples = {}
library = {}
datasize = {}
rawsize = {}
effesize = {}
q20size = {}
q30size = {}
errsize = {}
gcsize = {}
lanen = {}
novoid = {}
qc_sta = {}
qcstb = {}

for n in A:
    if n.startswith('#'): continue
    l = n.strip().split()
    if n.strip().split()[5].endswith('.bam'): continue
    qc_stat = os.path.join(
        projdir, 'QC/' + l[2] + '/' + l[2] + '_' + l[4] + '_' + l[0] + '.stat')
    B = utils.safe_open(qc_stat, 'r')
    #x = [l[1],l[2],l[4],l[6]]
    if l[1] not in samples:

        samples[l[1]] = {}
        samples[l[1]][l[2]] = {}
        samples[l[1]][l[2]][l[4] + '_' + l[0]] = qc_stat

    elif l[2] not in samples[l[1]]:
        samples[l[1]][l[2]] = {}
        samples[l[1]][l[2]][l[4] + '_' + l[0]] = qc_stat
    else:
        if l[4] + '_' + l[0] in samples[l[1]][l[2]]: continue
        else:
            samples[l[1]][l[2]][l[4] + '_' + l[0]] = qc_stat
    items = [l[2], l[4], l[0]]
    lanename = l[2] + '_' + l[4] + '_' + l[0]
    qc_sta[lanename] = []
    B.next()
    qc_sta1 = B.next().strip().split('\t')[1]
    if l[2] not in novoid.keys():
        novoid[l[2]] = l[4]
    if l[2] not in rawsize.keys():
        lanen[l[2]] = 1
        rawsize[l[2]] = int(qc_sta1)
    else:
        lanen[l[2]] += 1
        rawsize[l[2]] += int(qc_sta1)
    qc_sta[lanename].append(qc_sta1)
    items.append(qc_sta1)
    x = B.next().strip()
    qc_sta2 = int(x.split('\t')[1])
    qc_sta[lanename].append(qc_sta2)
    qc_sta[lanename].append(qc_sta2)
    if l[2] not in datasize:
        datasize[l[2]] = int(x.split('\t')[1])
    else:
        datasize[l[2]] += int(x.split('\t')[1])
    qc_sta3 = x.split('\t')[2].split('(')[1][:-2]
    if l[2] not in effesize.keys():
        effesize[l[2]] = int(qc_sta1) * float(qc_sta3)
    else:
        effesize[l[2]] += int(qc_sta1) * float(qc_sta3)
    qc_sta[lanename].append(qc_sta3)
    items.append(qc_sta3)
    B.next(), B.next(), B.next(), B.next(), B.next(), B.next()
    Q20 = B.next().strip().split('\t')[2][:-1] + ';' + B.next().strip().split(
        '\t')[2][:-1]
    Q30 = B.next().strip().split('\t')[2][:-1] + ';' + B.next().strip().split(
        '\t')[2][:-1]
    GC = B.next().strip().split('\t')[2][:-1] + ';' + B.next().strip().split(
        '\t')[2][:-1]
    Err = B.next().strip().split('\t')[2][:-1] + ';' + B.next().strip().split(
        '\t')[2][:-1]
    Q20_use = round((float(Q20.split(';')[0]) + float(Q20.split(';')[1])) / 2,
                    2)
    Q30_use = round((float(Q30.split(';')[0]) + float(Q30.split(';')[1])) / 2,
                    2)
    GC_use = round((float(GC.split(';')[0]) + float(GC.split(';')[1])) / 2, 2)
    Err_use = round((float(Err.split(';')[0]) + float(Err.split(';')[1])) / 2,
                    2)
    items.append(Err_use)
    items.append(Q20_use)
    items.append(Q30_use)
    items.append(GC_use)
    if l[2] not in q20size.keys():
        q20size[l[2]] = int(qc_sta1) * Q20_use
    else:
        q20size[l[2]] += int(qc_sta1) * Q20_use
    if l[2] not in q30size.keys():
        q30size[l[2]] = int(qc_sta1) * Q30_use
    else:
        q30size[l[2]] += int(qc_sta1) * Q30_use
    if l[2] not in errsize.keys():
        errsize[l[2]] = int(qc_sta1) * Err_use
    else:
        errsize[l[2]] += int(qc_sta1) * Err_use
    if l[2] not in gcsize.keys():
        gcsize[l[2]] = int(qc_sta1) * GC_use
    else:
        gcsize[l[2]] += int(qc_sta1) * GC_use
    qc_sta[lanename].append(str(Err_use))
    qc_sta[lanename].append(str(Q20_use))
    qc_sta[lanename].append(str(Q30_use))
    qc_sta[lanename].append(str(GC_use))
    if l[2] not in library:

        library[l[2]] = [items]

    else:
        library[l[2]].append(items)
    B.close()
A.close()




samli = []
sam = args['sample']
A = utils.safe_open(sam, 'r')

for n in A:
    if n.lower().startswith('#familyid'):
        ns = n.lower().strip().split('\t')
        continue
    if n.startswith('#'):
        continue
    n = n.strip().split('\t')
    if 'data' in ns:
        datapos = ns.index('data')
        if len(n) >= (datapos + 1):
            if n[datapos] != "0":
                if n[1] not in samli:
                    samli.append(n[1])
        else:
            if n[1] not in samli:
                samli.append(n[1])
    else:
        if n[1] not in samli:
            samli.append(n[1])
    nsfamily = ns.index('#familyid')
    nssample = ns.index('sampleid')
    nssex = ns.index('sex')
    nsnormal = ns.index('normal/patient')
A.close()
if sam != 'Null':
    tusam = []
    samf = open(sam, 'r')
    nin = 0
    for i in samf:
        i = i.strip()
        if i.startswith('#'):
            nin += 1
    tsam = txts2tab([sam])[nin:]
    for i in tsam:
        if i not in tusam:
            tusam.append([i[nsfamily], i[nssample], i[nssex], i[nsnormal]])
    context['table_sam_info'] = gettainfo(tusam, samli)
    samf.close()
else:
    context['table_samp_info'] = False

nnn = 1
for i in context['table_sam_info']:
    i.insert(0, str(nnn))
    nnn += 1

bg = 'B1'
if sam != 'Null':
    diseases = ''
    genelist = ''
    with open(sam, 'r') as f:
        for i in f:
            i = i.strip()
            if i.lower().startswith('#disease'):
                diseases = i.split(':')[1]
            if i.lower().startswith('#tissue'):
                tissue = i.split(':')[1]
            elif i.lower().startswith('#gene'):
                genelist = i.split(':')[1]
            elif i.startswith('#B'):
                bg = i[1:]

if 'phenolyzer' in ANALYSIS and diseases != '' and rep_ty == "advance":
    context['net_work'] = True
else:
    context['net_work'] = False

try:
    pn = pn.decode('utf8')
except:
    pn = pn.decode('gbk').encode('utf8')

reportnum = pn.split()[0] + '-' + bg + '-5'
hetongname = pn.split()[1]
hetongnum = pn.split()[0]
hetongdir = os.path.join(projdir, 'Report', hetongnum)
context['reportnum'] = reportnum
context['hetongnum'] = hetongnum
context['hetongname'] = hetongname





lenTR = 0
#sqchoices=['WES_ag','WES_illu','WGS']
if seq_ty == 'TS':  ##for TS
    for i in utils.safe_open(TR, 'r'):
        i = i.strip().split()
        lenTR += int(i[2]) - int(i[1])
    print 'The base on the bed file is:\t', lenTR
elif seq_ty == 'WGS':  ##for WGS
    lenTR = 2894323861
elif seq_ty == 'WES_illu':  ##for illumina
    lenTR = 62085295
else:  ##for agilent
    if 'V6' in TR:
        lenTR = 60456963
        context['v6'] = True
    elif 'V5' in TR:  ##for v5
        lenTR = 50390601
        context['v5'] = True

for v in library:
    total = round(datasize[v] / 1000000000.0, 2)
    library[v][0][4:4] = [total]
    #if seq_ty=='WGS':
    #    library[v][0][5:5]=[round(total/2.894323861,2)]
    #if seq_ty=='WES_ag':
    #    library[v][0][5:5]=[round(total/0.050390601,2)]
    #if seq_ty=='WES_illu':
    #    library[v][0][5:5]=[round(total/0.062085295,2)]
    #if seq_ty=='TS':
    #    library[v][0][5:5]=[round(total/(lenTR*0.000000001),2)]
    library[v][0][5:5] = [round(total / (lenTR * 0.000000001), 2)]

for i in qc_sta:
    total = round(qc_sta[i][1] / 1000000000.0, 2)
    qc_sta[i][1] = str(total)
    #if seq_ty=='WGS':
    #    qc_sta[i][2]=str(round(total/2.894323861,2))
    #if seq_ty=='WES_ag':
    #    qc_sta[i][2]=str(round(total/0.050390601,2))
    #if seq_ty=='WES_illu':
    #    qc_sta[i][2]=str(round(total/0.062085295,2))
    #if seq_ty=='TS':
    #qc_sta[i][2]=str(round(total/(lenTR*0.000000001),2))
    qc_sta[i][2] = str(round(total / (lenTR * 0.000000001), 2))
nnn = 1
for i in samli:
    qcstb[i] = []
    qcstb[i].append(i)
    qcstb[i].append(novoid[i])
    qcstb[i].append(str(rawsize[i]))
    total = round(datasize[i] / 1000000000.0, 2)
    qcstb[i].append(str(total))
    totaldep = str(round(total / (lenTR * 0.000000001), 2))
    qcstb[i].append(str(totaldep))
    effe = round(effesize[i] / rawsize[i], 2)
    qcstb[i].append(str(effe))
    erro = round(errsize[i] / rawsize[i], 2)
    qcstb[i].append(str(erro))
    q20o = round(q20size[i] / rawsize[i], 2)
    qcstb[i].append(str(q20o))
    q30o = round(q30size[i] / rawsize[i], 2)
    qcstb[i].append(str(q30o))
    gco = round(gcsize[i] / rawsize[i], 2)
    qcstb[i].append(str(gco))
    qcstb[i].append(str(nnn))
    nnn += 1
for i in samli:
    context['table_qcb'].append(qcstb[i])
for k in samli:
    context['table_qc'].append(library[k])
#print context['table_qc']

####init
#A = utils.safe_open(infile,'r')
#samples = {}
#librarys = {}
#for n in A:
#    if n.startswith('#'):continue
#    l = n.strip().split()
#
#    ##
#    if l[1] not in samples:
#
#        samples[l[1]] = {}
#        samples[l[1]][l[2]] = {}
#        samples[l[1]][l[2]][l[4]+'_'+l[0]] = l[4]
#
#    elif l[2] not in samples[l[1]]:
#        samples[l[1]][l[2]] = {}
#        samples[l[1]][l[2]][l[4]+'_'+l[0]] = l[4]
#    else:
#        samples[l[1]][l[2]][l[4]+'_'+l[0]] = l[4]
#A.close()
#

## depth table
context['gaoliang'] = [
    'Mapped:',
    'Average_sequencing_depth_on_target:',
    'Average_sequencing_depth:',
    'Coverage:',
    'Coverage_at_least_10X:'
    'Coverage_of_target_region:',
    'Fraction_of_target_covered_with_at_least_10x:'
]
##add fig_class
depth_files = []
context['fig_gc'] = []
context['fig_error'] = []
context['fig_qm'] = []
context['fig_class'] = []

##copy report dirs
if english == 'Y':
    print 'English report not available yet'
else:
    assert not os.system('cp -r %s/src %s' % (REPORT_DIR, odir))

if NGC == 'N':
    context['NGC'] = '2GCz'
else:
    context['NGC'] = 'GCz'

## depth| yi gai hao
depth_files = []
for m in samli:
    depth_files.append(
        os.path.join(projdir, 'Alnstat', m, m + '_mapping_coverage.txt'))
if rep_ty != 'qc':
    context['table_depth'] = txts2tab(depth_files, 1)
    mmm = []
    nnn = 1
    for m in samli:
        mmm.append(str(nnn))
        nnn += 1
    mmm.insert(0, 'ReportID:')
    context['table_depth'].insert(0, mmm)
    context['table_dept'] = ((np.array(context['table_depth'])).T).tolist()
    depthavg = 0
    if seq_ty == 'TS' or seq_ty == 'WES_ag':
        context['table_dept'][0] = [
            "ReportID", "Sample", "Total", "Duplicate", "Mapped",
            "Properly mapped", "PE mapped", "SE mapped",
            "Initial_base_on_target", "Initial_base_near_target",
            "Total_effective_yield(Mb)", "Effective_yield_on_target(Mb)",
            "Fraction_of_effective_bases_on_target",
            "Fraction_of_effective_bases_on_near_target",
            "Average_sequencing_depth_on_target", "Base_covered_on_target",
            "Coverage_of_target_region", "100X", "50X", "20X", "10X", "4X"
        ]
        for i in context['table_depth'][14][1:]:
            depthavg += float(i)
        depthavg = depthavg * 5 / len(context['table_depth'][14][1:])
    elif seq_ty == 'WGS':
        context['table_dept'][0] = [
            "ReportID", "Sample", "Total", "Duplicate", "Mapped",
            "Properly mapped", "PE mapped", "SE mapped", "diff chr",
            "diff chr(mapQ≥5)", "Average_Depth", "Coverage", "4X", "10X", "20X"
        ]
        for i in context['table_depth'][10][1:]:
            depthavg += float(i)
        depthavg = depthavg * 5 / len(context['table_depth'][10][1:])
    elif seq_ty == 'WES_illu':
        context['table_dept'][0] = [""]
        for i in context['table_depth'][20][1:]:
            depthavg += float(i)
        depthavg = depthavg * 5 / len(context['table_depth'][20][1:])

for k in sorted(samples.keys()):
    for m in sorted(samples[k].keys()):
        for n in sorted(samples[k][m].keys()):
            n = m + '_' + n
            context['fig_gc'].append([
                '"src/pictures/GC/' + n + '.GC.png"',
                '"src/pictures/GC/' + n + '.GC.JPEG"'
            ])
            context['fig_error'].append([
                '"src/pictures/Error/' + n + '.Error.png"',
                '"src/pictures/Error/' + n + '.Error.JPEG"'
            ])
            context['fig_qm'].append([
                '"src/pictures/Quality/' + n + '.QM.png"',
                '"src/pictures/Quality/' + n + '.QM.JPEG"'
            ])
            context['fig_class'].append([
                '"src/pictures/Class/' + n + '.pie3d.png"',
                '"src/pictures/Class/' + n + '.pie3d.JPEG"'
            ])

###depth pics paint
if rep_ty != 'qc':
    # new file: depth.list
    depthList = utils.safe_open(odir + "/depth.list", 'w')
    cumuList = utils.safe_open(odir + "/cumu.list", 'w')
    covByChrList = utils.safe_open(odir + "/covByChr.list", 'w')

    for k in samples.keys():
        for l in samples[k].keys():
            depthPicDir = os.path.join(projdir, 'Alnstat', l)
            depthList.write(l + "\t" + depthPicDir + '/depth_frequency.xls\n')
            cumuList.write(l + "\t" + depthPicDir + '/cumu.xls\n')
            covByChrList.write(
                l + "\t" + depthPicDir + '/' + l + '.coverage.bychr.txt\n')

    depthList.close()
    cumuList.close()
    covByChrList.close()

    # run to draw depth pic
    assert not os.system('cd %s' % odir)
    assert not os.system(
        '%s/Mapping/Stat/collectDepthPlot.pl -i=0 -j=1 %s > %s' %
        (moduledir, odir + "/depth.list", odir + "/depth.out"))
    assert not os.system(
        '%s/Mapping/Stat/collectDepthPlot.pl -i=0 -j=1 %s > %s' %
        (moduledir, odir + "/cumu.list", odir + "/cumu.out"))
    assert not os.system('%s/Mapping/Stat/depth.plot.R %s %s %s' %
                         (moduledir, odir + "/depth.out", odir + "/cumu.out",
                          odir + '/src/pictures/depth/depth.png'))

    assert not os.system('collectCovByChr.sh %s %s' %
                         (odir + "/covByChr.list",
                          odir + '/src/pictures/depth/covByChr.png'))
    assert not os.system('rm %s %s %s %s %s' %
                         (odir + "/depth.list", odir + "/cumu.list",
                          odir + "/covByChr.list", odir + "/depth.out",
                          odir + "/cumu.out"))

if rep_ty == 'primary' or rep_ty == 'advance':
    k = samples.keys()[0]
    m = samples[k].keys()[0]
    advanceDir = os.path.join(projdir, 'Advance', newjob)
    advanceBriDir = os.path.join(projdir, 'Advance', newjob, 'BriefResults')
    SumDir = os.path.join(advanceDir, 'Summary')
    # only regonize GATK calling
    if 'snpindel_call_gatk' in ANALYSIS:
        if sf == 'Y':
            single_snp_filter = os.path.join(projdir, "Mutation", m + '.gatk',
                                             "FilterSNP",
                                             m + '.SnpFilter.filter.stat.xls')
            single_indel_filter = os.path.join(
                projdir, "Mutation", m + '.gatk', "FilterInDel",
                m + '.InDelFilter.filter.stat.xls')
            context['table_single_snp_filter'] = txts2tab([single_snp_filter])
            context['table_single_indel_filter'] = txts2tab(
                [single_indel_filter])
        snp_anno = os.path.join(projdir, "Mutation", m + '.gatk',
                                m + '.gatk.snp.annovar.hg19_multianno.xls.gz')
        indel_anno = os.path.join(
            projdir, "Mutation", m + '.gatk',
            m + '.gatk.indel.annovar.hg19_multianno.xls.gz')
        context['samtools'] = False

    elif ('snpindel_call_samtools' in ANALYSIS or 'snpindel_call_samtools_multi' in ANALYSIS):
        if sf == 'Y':
            single_snp_filter = os.path.join(projdir, "Mutation",
                                             m + '.samtools', "FilterSNP",
                                             m + '.SnpFilter.filter.stat.xls')
            single_indel_filter = os.path.join(
                projdir, "Mutation", m + '.samtools', "FilterInDel",
                m + '.InDelFilter.filter.stat.xls')
            context['table_single_snp_filter'] = txts2tab([single_snp_filter])
            context['table_single_indel_filter'] = txts2tab(
                [single_indel_filter])
        snp_anno = os.path.join(
            projdir, "Mutation", m + '.samtools',
            m + '.samtools.snp.annovar.hg19_multianno.xls.gz')
        indel_anno = os.path.join(
            projdir, "Mutation", m + '.samtools',
            m + '.samtools.indel.annovar.hg19_multianno.xls.gz')
        context['samtools'] = True

    elif 'snpindel_call_sentieon' in ANALYSIS:
        snp_anno = os.path.join(
            projdir, "Mutation", m + '.sentieon',
            m + '.sentieon.snp.annovar.hg19_multianno.xls.gz')
        indel_anno = os.path.join(
            projdir, "Mutation", m + '.sentieon',
            m + '.sentieon.indel.annovar.hg19_multianno.xls.gz')
        context['sentieon'] = True

    if 'snpindel_call' in ANALYSIS:
        context['table_anno'] = txts2tab([snp_anno])[:6]
        context['table_anno_indel'] = txts2tab([indel_anno])[:6]

        func_snp = txts2tab([os.path.join(SumDir, 'counts/snp_func.xls')], 0,
                            [12, 12])
        exfunc_snp = txts2tab(
            [os.path.join(SumDir, 'counts/snp_exonicfunc.xls')])
        ts_tv = txts2tab([os.path.join(SumDir, 'counts/ts_tv.xls')])
        trait_snp = txts2tab([os.path.join(SumDir, 'counts/snp_other.xls')])
        func_indel = txts2tab([os.path.join(SumDir, 'counts/indel_func.xls')],
                              0, [12, 12])
        exfunc_indel = txts2tab(
            [os.path.join(SumDir, 'counts/indel_exonicfunc.xls')])
        trait_indel = txts2tab(
            [os.path.join(SumDir, 'counts/indel_other.xls')])
        context['table_func_snp'] = getta(func_snp, samli)
        context['table_exfunc_snp'] = getta(exfunc_snp, samli)
        context['table_ts_tv'] = getta(ts_tv, samli)
        context['table_trait_snp'] = getta(trait_snp, samli)
        context['table_func_indel'] = getta(func_indel, samli)
        context['table_exfunc_indel'] = getta(exfunc_indel, samli)
        context['table_trait_indel'] = getta(trait_indel, samli)

    ##snp/indel figs
    images_dir = os.path.join(SumDir, 'pics')
    src_mut_images_dir = os.path.join(odir,
                                      'src/pictures/mutation_stat_images')
    context['fig_snp_func'] = []
    context['fig_snp_trait'] = []
    context['fig_indel_len'] = []
    context['fig_indel_func'] = []
    context['fig_indel_trait'] = []
    context['fig_CoNIFER'] = []
    context['fig_CNV'] = []

    if 'sv_call' in ANALYSIS:
        sv_stat, sv_anno = [], []
    if 'cnv_call' in ANALYSIS:
        cnv_stat, cnv_anno = [], []
        context['fig_cnv_circos'] = []

    sv_datatra = {}
    sv_tra = {}
    sv_datains = {}
    sv_ins = {}
    sv_datadel = {}
    sv_del = {}
    sv_datainv = {}
    sv_inv = {}
    sv_datains2 = {}
    sv_ins2 = {}
    sv_datainv2 = {}
    sv_inv2 = {}
    sv_datadel2 = {}
    sv_del2 = {}
    sv_dataitx2 = {}
    sv_itx2 = {}
    sv_datactx2 = {}
    sv_ctx2 = {}
    sv_datatdp2 = {}
    sv_tdp2 = {}
    cnv_datagai = {}
    cnv_gai = {}
    cnv_datalos = {}
    cnv_los = {}
    context['table_cnv_anno'] = []
    for k2 in sorted(samples.keys()):
        for m2 in sorted(samples[k2].keys()):
            sv_datatra[m2] = ':[{name:"Translocation",value:0},'
            sv_tra[
                m2] = ':[{name:"cds",value:0},{name:"splicing",value:0},{name:"utr5",value:0},{name:"utr3",value:0},{name:"intron",value:0},{name:"upstream",value:0},{name:"downstream",value:0},{name:"ncRNA",value:0},{name:"intergenic",value:0},{name:"unknown",value:0}],'
            sv_datadel[m2] = '{name:"Deletion",value:0},'
            sv_del[
                m2] = ':[{name:"cds",value:0},{name:"splicing",value:0},{name:"utr5",value:0},{name:"utr3",value:0},{name:"intron",value:0},{name:"upstream",value:0},{name:"downstream",value:0},{name:"ncRNA",value:0},{name:"intergenic",value:0},{name:"unknown",value:0}],'
            sv_datains[m2] = '{name:"Insertion",value:0},'
            sv_ins[
                m2] = ':[{name:"cds",value:0},{name:"splicing",value:0},{name:"utr5",value:0},{name:"utr3",value:0},{name:"intron",value:0},{name:"upstream",value:0},{name:"downstream",value:0},{name:"ncRNA",value:0},{name:"intergenic",value:0},{name:"unknown",value:0}],'
            sv_datainv[m2] = '{name:"Inversion",value:0}],'
            sv_inv[
                m2] = ':[{name:"cds",value:0},{name:"splicing",value:0},{name:"utr5",value:0},{name:"utr3",value:0},{name:"intron",value:0},{name:"upstream",value:0},{name:"downstream",value:0},{name:"ncRNA",value:0},{name:"intergenic",value:0},{name:"unknown",value:0}],'
            cnv_datagai[m2] = ':[{name:"gain",value:0},'
            cnv_gai[
                m2] = ':[{name:"cds",value:0},{name:"splicing",value:0},{name:"utr5",value:0},{name:"utr3",value:0},{name:"intron",value:0},{name:"upstream",value:0},{name:"downstream",value:0},{name:"ncRNA",value:0},{name:"intergenic",value:0},{name:"unknown",value:0}],'
            cnv_datalos[m2] = '{name:"loss",value:0}],'
            cnv_los[
                m2] = ':[{name:"cds",value:0},{name:"splicing",value:0},{name:"utr5",value:0},{name:"utr3",value:0},{name:"intron",value:0},{name:"upstream",value:0},{name:"downstream",value:0},{name:"ncRNA",value:0},{name:"intergenic",value:0},{name:"unknown",value:0}],'
            sv_datains2[m2] = ':[{name:"INS",value:0},'
            sv_ins2[
                m2] = ':[{name:"cds",value:0},{name:"splicing",value:0},{name:"utr5",value:0},{name:"utr3",value:0},{name:"intron",value:0},{name:"upstream",value:0},{name:"downstream",value:0},{name:"ncRNA",value:0},{name:"intergenic",value:0},{name:"unknown",value:0}],'
            sv_datadel2[m2] = '{name:"DEL",value:0},'
            sv_del2[
                m2] = ':[{name:"cds",value:0},{name:"splicing",value:0},{name:"utr5",value:0},{name:"utr3",value:0},{name:"intron",value:0},{name:"upstream",value:0},{name:"downstream",value:0},{name:"ncRNA",value:0},{name:"intergenic",value:0},{name:"unknown",value:0}],'
            sv_datainv2[m2] = '{name:"INV",value:0},'
            sv_inv2[
                m2] = ':[{name:"cds",value:0},{name:"splicing",value:0},{name:"utr5",value:0},{name:"utr3",value:0},{name:"intron",value:0},{name:"upstream",value:0},{name:"downstream",value:0},{name:"ncRNA",value:0},{name:"intergenic",value:0},{name:"unknown",value:0}],'
            sv_dataitx2[m2] = '{name:"ITX",value:0},'
            sv_itx2[
                m2] = ':[{name:"cds",value:0},{name:"splicing",value:0},{name:"utr5",value:0},{name:"utr3",value:0},{name:"intron",value:0},{name:"upstream",value:0},{name:"downstream",value:0},{name:"ncRNA",value:0},{name:"intergenic",value:0},{name:"unknown",value:0}],'
            sv_datactx2[m2] = '{name:"CTX",value:0}],'
            sv_ctx2[
                m2] = ':[{name:"cds",value:0},{name:"splicing",value:0},{name:"utr5",value:0},{name:"utr3",value:0},{name:"intron",value:0},{name:"upstream",value:0},{name:"downstream",value:0},{name:"ncRNA",value:0},{name:"intergenic",value:0},{name:"unknown",value:0}],'
            sv_datatdp2[m2] = '{name:"TDP",value:0}],'
            sv_tdp2[
                m2] = ':[{name:"cds",value:0},{name:"splicing",value:0},{name:"utr5",value:0},{name:"utr3",value:0},{name:"intron",value:0},{name:"upstream",value:0},{name:"downstream",value:0},{name:"ncRNA",value:0},{name:"intergenic",value:0},{name:"unknown",value:0}],'

            if 'sv_call_breakdancer' in ANALYSIS:
                sv_stat_file = os.path.join(
                    projdir, 'SV', m2, 'breakdancer',
                    m2 + '.breakdancer.flt.gff.ann.stat.txt')
                sv_stat.append(sv_stat_file)
                sv_anno_file = os.path.join(
                    projdir, 'SV', m2, 'breakdancer',
                    m2 + '.breakdancer.flt.hg19_multianno.xls')
                sv_anno.append(sv_anno_file)
                context['breakdancer'] = True
                context['crest'] = False
                context['BreaKmer'] = False
                context['table_sv'] = svtxts2tab(sv_stat)
                context['table_sv_anno'] = txts2tab([sv_anno[0]])[:6]
                #context['table_sv_anno_head']=['染色体','SV的起始位置','SV的终止位置','基因名称','SV所在的区域','SV相关的转录本',' 变异情况描述','Gencode注释','CpG岛名称','染色体区段','变异相关的microRNA和snoRNA','预测microRNA的靶点','变异保守性预测','转录因子结合位点保守分值','重复片段','dgv数据库','gwasCatalog数据库','重复序列','基因组中区域的功能注释','基因组中区域的功能注释','基因组中区域的功能注释','基因组中区域的功能注释','基因组中区域的功能注释','基因组中区域的功能注释','SV第一个断点处在锚定区域中比对到正链或负链的reads数量','SV另一断点处在锚定区域中比对到正链或负链的reads数量','SV的分值','SV的长度','支持SV的reads的数目','每个样本中支持该SV的reads数（多样本的情况）','易位类型','对于易位（translocation），另一个断点所在的染色体','对于易位（translocation），另一个断点所在染色体的坐标','SV的ID号','SV的类型']
                #context['table_filter_sv_head']=['染色体','SV的起始位置','SV的终止位置','基因名称','SV所在的区域','SV相关的转录本',' 变异情况描述','Gencode注释','CpG岛名称','染色体区段','变异相关的microRNA和snoRNA','预测microRNA的靶点','变异保守性预测','转录因子结合位点保守分值','重复片段','dgv数据库','gwasCatalog数据库','重复序列','基因组中区域的功能注释','基因组中区域的功能注释','基因组中区域的功能注释','基因组中区域的功能注释','基因组中区域的功能注释','基因组中区域的功能注释','SV第一个断点处在锚定区域中比对到正链或负链的reads数量','SV另一断点处在锚定区域中比对到正链或负链的reads数量','SV的分值','SV的长度','支持SV的reads的数目','每个样本中支持该SV的reads数（多样本的情况）','易位类型','对于易位（translocation），另一个断点所在的染色体','对于易位（translocation），另一个断点所在染色体的坐标','SV的ID号','SV的类型']
                with open(sv_stat_file) as f:
                    for k in islice(f, 1, None):
                        k = k.strip()
                        k = k.split('\t')
                        if k[1] == 'Inversion':
                            sv_datainv[
                                m2] = '{name:"Inversion",value:' + k[2] + '}],'
                            sv_inv[
                                m2] = ':[{name:"cds",value:' + k[3] + '},{name:"splicing",value:' + k[4] + '},{name:"utr5",value:' + k[5] + '},{name:"utr3",value:' + k[6] + '},{name:"intron",value:' + k[7] + '},{name:"upstream",value:' + k[8] + '},{name:"downstream",value:' + k[9] + '},{name:"ncRNA",value:' + k[10] + '},{name:"intergenic",value:' + k[11] + '},{name:"unknown",value:' + k[12] + '}],'
                        elif k[1] == 'Insertion':
                            sv_datains[
                                m2] = '{name:"Insertion",value:' + k[2] + '},'
                            sv_ins[
                                m2] = ':[{name:"cds",value:' + k[3] + '},{name:"splicing",value:' + k[4] + '},{name:"utr5",value:' + k[5] + '},{name:"utr3",value:' + k[6] + '},{name:"intron",value:' + k[7] + '},{name:"upstream",value:' + k[8] + '},{name:"downstream",value:' + k[9] + '},{name:"ncRNA",value:' + k[10] + '},{name:"intergenic",value:' + k[11] + '},{name:"unknown",value:' + k[12] + '}],'
                        elif k[1] == 'Translocation':
                            sv_datatra[
                                m2] = ':[{name:"Translocation",value:' + k[2] + '},'
                            sv_tra[
                                m2] = ':[{name:"cds",value:' + k[3] + '},{name:"splicing",value:' + k[4] + '},{name:"utr5",value:' + k[5] + '},{name:"utr3",value:' + k[6] + '},{name:"intron",value:' + k[7] + '},{name:"upstream",value:' + k[8] + '},{name:"downstream",value:' + k[9] + '},{name:"ncRNA",value:' + k[10] + '},{name:"intergenic",value:' + k[11] + '},{name:"unknown",value:' + k[12] + '}],'
                        elif k[1] == 'Deletion':
                            sv_datadel[
                                m2] = '{name:"Deletion",value:' + k[2] + '},'
                            sv_del[
                                m2] = ':[{name:"cds",value:' + k[3] + '},{name:"splicing",value:' + k[4] + '},{name:"utr5",value:' + k[5] + '},{name:"utr3",value:' + k[6] + '},{name:"intron",value:' + k[7] + '},{name:"upstream",value:' + k[8] + '},{name:"downstream",value:' + k[9] + '},{name:"ncRNA",value:' + k[10] + '},{name:"intergenic",value:' + k[11] + '},{name:"unknown",value:' + k[12] + '}],'

            if 'sv_call_crest' in ANALYSIS:
                sv_stat_file = os.path.join(projdir, 'SV', m2, 'crest',
                                            m2 + '.crest.gff.ann.stat.xls')
                sv_stat.append(sv_stat_file)
                sv_anno_file = os.path.join(projdir, 'SV', m2, 'crest',
                                            m2 + '.crest.hg19_multianno.xls')
                sv_anno.append(sv_anno_file)
                context['crest'] = True
                context['breakdancer'] = False
                context['BreaKmer'] = False
                context['table_sv'] = svtxts2tab(sv_stat)
                context['table_sv_anno'] = txts2tab([sv_anno[0]])[:6]
                #context['table_sv_anno_head']=['染色体','SV的起始位置','SV的终止位置','基因名称','SV所在的区域','SV相关的转录本',' 变异情况描述','Gencode注释','CpG岛名称','染色体区段','变异相关的microRNA和snoRNA','预测microRNA的靶点','变异保守性预测','转录因子结合位点保守分值','重复片段','dgv数据库','gwasCatalog数据库','重复序列','基因组中区域的功能注释','基因组中区域的功能注释','基因组中区域的功能注释','基因组中区域的功能注释','基因组中区域的功能注释','基因组中区域的功能注释','SV左断点位置及所在染色体','SV右断点位置及所在染色体','左断点所在链方向','右断点所在链方向','覆盖左断点的reads 数目','覆盖右断点的reads数目','左断点位置reads覆盖数目','右断点位置reads覆盖数目','在左断点位置被修剪掉的reads组装长度','在右断点位置被修剪掉的reads组装长度','左断点contig平均相似度百分比','左断点位置非unique reads覆盖数目占总覆盖数目的比例','右断点contig平均相似度百分比','右断点位置非unique reads覆盖数目占总覆盖数目的比例',' 覆盖左右断点的所有reads联合序列的相对起始位置','支持该SV类型的有效序列起始染色体','覆盖左右断点的所有reads联合序列绝对起始位置','覆盖左右断点的所有reads联合序列的相对终止位置','支持该SV类型的有效序列终止染色体','覆盖左右断点的所有reads联合序列绝对终止位置','覆盖左右断点的所有reads联合序列','对于易位，另一个断点所在的染色体','对于易位，另一个断点所在染色体的坐标','SV的ID号','SV的类型']
                #context['table_filter_sv_head']=['染色体','SV的起始位置','SV的终止位置','基因名称','SV所在的区域','SV相关的转录本',' 变异情况描述','Gencode注释','CpG岛名称','染色体区段','变异相关的microRNA和snoRNA','预测microRNA的靶点','变异保守性预测','转录因子结合位点保守分值','重复片段','dgv数据库','gwasCatalog数据库','重复序列','基因组中区域的功能注释','基因组中区域的功能注释','基因组中区域的功能注释','基因组中区域的功能注释','基因组中区域的功能注释','基因组中区域的功能注释','SV第一个断点处在锚定区域中比对到正链或负链的reads数量','SV另一断点处在锚定区域中比对到正链或负链的reads数量','SV的分值','SV的长度','支持SV的reads的数目','每个样本中支持该SV的reads数（多样本的情况）','易位类型','对于易位（translocation），另一个断点所在的染色体','对于易位（translocation），另一个断点所在染色体的坐标','SV的ID号','SV的类型']
                with open(sv_stat_file, 'r') as f:
                    for k in islice(f, 1, None):
                        k = k.strip()
                        k = k.split('\t')
                        if k[1] == 'INV':
                            sv_datainv2[
                                m2] = '{name:"INV",value:' + k[2] + '},'
                            sv_inv2[
                                m2] = ':[{name:"cds",value:' + k[3] + '},{name:"splicing",value:' + k[4] + '},{name:"utr5",value:' + k[5] + '},{name:"utr3",value:' + k[6] + '},{name:"intron",value:' + k[7] + '},{name:"upstream",value:' + k[8] + '},{name:"downstream",value:' + k[9] + '},{name:"ncRNA",value:' + k[10] + '},{name:"intergenic",value:' + k[11] + '},{name:"unknown",value:' + k[12] + '}],'
                        elif k[1] == 'ITX':
                            sv_dataitx2[
                                m2] = '{name:"ITX",value:' + k[2] + '},'
                            sv_itx2[
                                m2] = ':[{name:"cds",value:' + k[3] + '},{name:"splicing",value:' + k[4] + '},{name:"utr5",value:' + k[5] + '},{name:"utr3",value:' + k[6] + '},{name:"intron",value:' + k[7] + '},{name:"upstream",value:' + k[8] + '},{name:"downstream",value:' + k[9] + '},{name:"ncRNA",value:' + k[10] + '},{name:"intergenic",value:' + k[11] + '},{name:"unknown",value:' + k[12] + '}],'
                        elif k[1] == 'INS':
                            sv_datains2[
                                m2] = ':[{name:"INS",value:' + k[2] + '},'
                            sv_ins2[
                                m2] = ':[{name:"cds",value:' + k[3] + '},{name:"splicing",value:' + k[4] + '},{name:"utr5",value:' + k[5] + '},{name:"utr3",value:' + k[6] + '},{name:"intron",value:' + k[7] + '},{name:"upstream",value:' + k[8] + '},{name:"downstream",value:' + k[9] + '},{name:"ncRNA",value:' + k[10] + '},{name:"intergenic",value:' + k[11] + '},{name:"unknown",value:' + k[12] + '}],'
                        elif k[1] == 'DEL':
                            sv_datadel2[
                                m2] = '{name:"DEL",value:' + k[2] + '},'
                            sv_del2[
                                m2] = ':[{name:"cds",value:' + k[3] + '},{name:"splicing",value:' + k[4] + '},{name:"utr5",value:' + k[5] + '},{name:"utr3",value:' + k[6] + '},{name:"intron",value:' + k[7] + '},{name:"upstream",value:' + k[8] + '},{name:"downstream",value:' + k[9] + '},{name:"ncRNA",value:' + k[10] + '},{name:"intergenic",value:' + k[11] + '},{name:"unknown",value:' + k[12] + '}],'
                        elif k[1] == 'CTX':
                            sv_datactx2[
                                m2] = '{name:"CTX",value:' + k[2] + '}],'
                            sv_ctx2[
                                m2] = ':[{name:"cds",value:' + k[3] + '},{name:"splicing",value:' + k[4] + '},{name:"utr5",value:' + k[5] + '},{name:"utr3",value:' + k[6] + '},{name:"intron",value:' + k[7] + '},{name:"upstream",value:' + k[8] + '},{name:"downstream",value:' + k[9] + '},{name:"ncRNA",value:' + k[10] + '},{name:"intergenic",value:' + k[11] + '},{name:"unknown",value:' + k[12] + '}],'
            if 'sv_call_lumpy' in ANALYSIS:
                sv_stat_file = os.path.join(projdir, 'SV', m2, 'lumpy',
                                            m2 + '.lumpy.ann.stat.xls')
                sv_stat.append(sv_stat_file)
                sv_anno_file = os.path.join(projdir, 'SV', m2, 'lumpy',
                                            m2 + '.lumpy.hg19_multianno.xls')
                sv_anno.append(sv_anno_file)
                context['lumpy'] = True
                context['breakdancer'] = False
                context['BreaKmer'] = False
                context['table_sv'] = svtxts2tab(sv_stat)
                context['table_sv_anno'] = txts2tab([sv_anno[0]])[:6]
                #context['table_sv_anno_head']=['染色体','SV的起始位置','SV的终止位置','基因名称','SV所在的区域','SV相关的转录本',' 变异情况描述','Gencode注释','CpG岛名称','染色体区段','变异相关的microRNA和snoRNA','预测microRNA的靶点','变异保守性预测','转录因子结合位点保守分值','重复片段','dgv数据库','gwasCatalog数据库','重复序列','基因组中区域的功能注释','基因组中区域的功能注释','基因组中区域的功能注释','基因组中区域的功能注释','基因组中区域的功能注释','基因组中区域的功能注释','SV左断点位置及所在染色体','SV右断点位置及所在染色体','左断点所在链方向','右断点所在链方向','覆盖左断点的reads 数目','覆盖右断点的reads数目','左断点位置reads覆盖数目','右断点位置reads覆盖数目','在左断点位置被修剪掉的reads组装长度','在右断点位置被修剪掉的reads组装长度','左断点contig平均相似度百分比','左断点位置非unique reads覆盖数目占总覆盖数目的比例','右断点contig平均相似度百分比','右断点位置非unique reads覆盖数目占总覆盖数目的比例',' 覆盖左右断点的所有reads联合序列的相对起始位置','支持该SV类型的有效序列起始染色体','覆盖左右断点的所有reads联合序列绝对起始位置','覆盖左右断点的所有reads联合序列的相对终止位置','支持该SV类型的有效序列终止染色体','覆盖左右断点的所有reads联合序列绝对终止位置','覆盖左右断点的所有reads联合序列','对于易位，另一个断点所在的染色体','对于易位，另一个断点所在染色体的坐标','SV的ID号','SV的类型']
                #context['table_filter_sv_head']=['染色体','SV的起始位置','SV的终止位置','基因名称','SV所在的区域','SV相关的转录本',' 变异情况描述','Gencode注释','CpG岛名称','染色体区段','变异相关的microRNA和snoRNA','预测microRNA的靶点','变异保守性预测','转录因子结合位点保守分值','重复片段','dgv数据库','gwasCatalog数据库','重复序列','基因组中区域的功能注释','基因组中区域的功能注释','基因组中区域的功能注释','基因组中区域的功能注释','基因组中区域的功能注释','基因组中区域的功能注释','SV第一个断点处在锚定区域中比对到正链或负链的reads数量','SV另一断点处在锚定区域中比对到正链或负链的reads数量','SV的分值','SV的长度','支持SV的reads的数目','每个样本中支持该SV的reads数（多样本的情况）','易位类型','对于易位（translocation），另一个断点所在的染色体','对于易位（translocation），另一个断点所在染色体的坐标','SV的ID号','SV的类型']
                with open(sv_stat_file, 'r') as f:
                    for k in islice(f, 1, None):
                        k = k.strip()
                        k = k.split('\t')
                        if k[1] == 'INV':
                            sv_datainv2[
                                m2] = '{name:"INV",value:' + k[2] + '},'
                            sv_inv2[
                                m2] = ':[{name:"cds",value:' + k[3] + '},{name:"splicing",value:' + k[4] + '},{name:"utr5",value:' + k[5] + '},{name:"utr3",value:' + k[6] + '},{name:"intron",value:' + k[7] + '},{name:"upstream",value:' + k[8] + '},{name:"downstream",value:' + k[9] + '},{name:"ncRNA",value:' + k[10] + '},{name:"intergenic",value:' + k[11] + '},{name:"unknown",value:' + k[12] + '}],'
                        elif k[1] == 'ITX':
                            sv_dataitx2[
                                m2] = '{name:"ITX",value:' + k[2] + '},'
                            sv_itx2[
                                m2] = ':[{name:"cds",value:' + k[3] + '},{name:"splicing",value:' + k[4] + '},{name:"utr5",value:' + k[5] + '},{name:"utr3",value:' + k[6] + '},{name:"intron",value:' + k[7] + '},{name:"upstream",value:' + k[8] + '},{name:"downstream",value:' + k[9] + '},{name:"ncRNA",value:' + k[10] + '},{name:"intergenic",value:' + k[11] + '},{name:"unknown",value:' + k[12] + '}],'
                        elif k[1] == 'INS':
                            sv_datains2[
                                m2] = ':[{name:"INS",value:' + k[2] + '},'
                            sv_ins2[
                                m2] = ':[{name:"cds",value:' + k[3] + '},{name:"splicing",value:' + k[4] + '},{name:"utr5",value:' + k[5] + '},{name:"utr3",value:' + k[6] + '},{name:"intron",value:' + k[7] + '},{name:"upstream",value:' + k[8] + '},{name:"downstream",value:' + k[9] + '},{name:"ncRNA",value:' + k[10] + '},{name:"intergenic",value:' + k[11] + '},{name:"unknown",value:' + k[12] + '}],'
                        elif k[1] == 'DEL':
                            sv_datadel2[
                                m2] = '{name:"DEL",value:' + k[2] + '},'
                            sv_del2[
                                m2] = ':[{name:"cds",value:' + k[3] + '},{name:"splicing",value:' + k[4] + '},{name:"utr5",value:' + k[5] + '},{name:"utr3",value:' + k[6] + '},{name:"intron",value:' + k[7] + '},{name:"upstream",value:' + k[8] + '},{name:"downstream",value:' + k[9] + '},{name:"ncRNA",value:' + k[10] + '},{name:"intergenic",value:' + k[11] + '},{name:"unknown",value:' + k[12] + '}],'
                        elif k[1] == 'CTX':
                            sv_datactx2[
                                m2] = '{name:"CTX",value:' + k[2] + '}],'
                            sv_ctx2[
                                m2] = ':[{name:"cds",value:' + k[3] + '},{name:"splicing",value:' + k[4] + '},{name:"utr5",value:' + k[5] + '},{name:"utr3",value:' + k[6] + '},{name:"intron",value:' + k[7] + '},{name:"upstream",value:' + k[8] + '},{name:"downstream",value:' + k[9] + '},{name:"ncRNA",value:' + k[10] + '},{name:"intergenic",value:' + k[11] + '},{name:"unknown",value:' + k[12] + '}],'

            if 'sv_call_breakmer' in ANALYSIS:
                sv_stat_file = os.path.join(projdir, 'SV', m2, 'breakmer',
                                            m2 + '.breakmer.gff.ann.stat.xls')
                sv_stat.append(sv_stat_file)
                sv_anno_file = os.path.join(
                    projdir, 'SV', m2, 'breakmer',
                    m2 + '.breakmer.hg19_multianno.xls')
                sv_anno.append(sv_anno_file)
                context['BreaKmer'] = True
                context['breakdancer'] = False
                context['crest'] = False
                context['table_sv'] = svtxts2tab(sv_stat)
                context['table_sv_anno'] = txts2tab([sv_anno[0]])[:6]
                with open(sv_stat_file, 'r') as f:
                    for k in islice(f, 1, None):
                        k = k.strip()
                        k = k.split('\t')
                        if k[1] == 'INV':
                            sv_datainv2[
                                m2] = '{name:"INV",value:' + k[2] + '},'
                            sv_inv2[
                                m2] = ':[{name:"cds",value:' + k[3] + '},{name:"splicing",value:' + k[4] + '},{name:"utr5",value:' + k[5] + '},{name:"utr3",value:' + k[6] + '},{name:"intron",value:' + k[7] + '},{name:"upstream",value:' + k[8] + '},{name:"downstream",value:' + k[9] + '},{name:"ncRNA",value:' + k[10] + '},{name:"intergenic",value:' + k[11] + '},{name:"unknown",value:' + k[12] + '}],'
                        elif k[1] == 'ITX':
                            sv_dataitx2[
                                m2] = '{name:"ITX",value:' + k[2] + '},'
                            sv_itx2[
                                m2] = ':[{name:"cds",value:' + k[3] + '},{name:"splicing",value:' + k[4] + '},{name:"utr5",value:' + k[5] + '},{name:"utr3",value:' + k[6] + '},{name:"intron",value:' + k[7] + '},{name:"upstream",value:' + k[8] + '},{name:"downstream",value:' + k[9] + '},{name:"ncRNA",value:' + k[10] + '},{name:"intergenic",value:' + k[11] + '},{name:"unknown",value:' + k[12] + '}],'
                        elif k[1] == 'INS':
                            sv_datains2[
                                m2] = ':[{name:"INS",value:' + k[2] + '},'
                            sv_ins2[
                                m2] = ':[{name:"cds",value:' + k[3] + '},{name:"splicing",value:' + k[4] + '},{name:"utr5",value:' + k[5] + '},{name:"utr3",value:' + k[6] + '},{name:"intron",value:' + k[7] + '},{name:"upstream",value:' + k[8] + '},{name:"downstream",value:' + k[9] + '},{name:"ncRNA",value:' + k[10] + '},{name:"intergenic",value:' + k[11] + '},{name:"unknown",value:' + k[12] + '}],'
                        elif k[1] == 'DEL':
                            sv_datadel2[
                                m2] = '{name:"DEL",value:' + k[2] + '},'
                            sv_del2[
                                m2] = ':[{name:"cds",value:' + k[3] + '},{name:"splicing",value:' + k[4] + '},{name:"utr5",value:' + k[5] + '},{name:"utr3",value:' + k[6] + '},{name:"intron",value:' + k[7] + '},{name:"upstream",value:' + k[8] + '},{name:"downstream",value:' + k[9] + '},{name:"ncRNA",value:' + k[10] + '},{name:"intergenic",value:' + k[11] + '},{name:"unknown",value:' + k[12] + '}],'
                        elif k[1] == 'CTX':
                            sv_datactx2[
                                m2] = '{name:"CTX",value:' + k[2] + '},'
                            sv_ctx2[
                                m2] = ':[{name:"cds",value:' + k[3] + '},{name:"splicing",value:' + k[4] + '},{name:"utr5",value:' + k[5] + '},{name:"utr3",value:' + k[6] + '},{name:"intron",value:' + k[7] + '},{name:"upstream",value:' + k[8] + '},{name:"downstream",value:' + k[9] + '},{name:"ncRNA",value:' + k[10] + '},{name:"intergenic",value:' + k[11] + '},{name:"unknown",value:' + k[12] + '}],'
                        elif k[1] == 'TDP':
                            sv_datatdp2[
                                m2] = '{name:"TDP",value:' + k[2] + '}],'
                            sv_tdp2[
                                m2] = ':[{name:"cds",value:' + k[3] + '},{name:"splicing",value:' + k[4] + '},{name:"utr5",value:' + k[5] + '},{name:"utr3",value:' + k[6] + '},{name:"intron",value:' + k[7] + '},{name:"upstream",value:' + k[8] + '},{name:"downstream",value:' + k[9] + '},{name:"ncRNA",value:' + k[10] + '},{name:"intergenic",value:' + k[11] + '},{name:"unknown",value:' + k[12] + '}],'

            if ('cnv_call_cnvnator' in ANALYSIS) and ('cnv_call_freec' not in ANALYSIS):
                cnv_stat_file = os.path.join(projdir, 'SV', m2, 'cnvnator',
                                             m2 + '.cnvnator.gff.ann.stat.xls')
                cnv_stat.append(cnv_stat_file)
                cnv_anno_file = os.path.join(
                    projdir, 'SV', m2, 'cnvnator',
                    m2 + '.cnvnator.hg19_multianno.xls')
                cnv_anno.append(cnv_anno_file)
                cnv_stat_tmp = svtxts2tab(cnv_stat)
                #del cnv_stat_tmp[1][1][0]
                context['cnvnator'] = True
                context['table_cnv'] = cnv_stat_tmp
                context['table_cnv_anno'] = txts2tab([cnv_anno[0]])[:6]
                #context['table_cnv_anno_head']=['染色体','CNV起始位置','CNV终止位置','基因名称','CNV所在的区域','CNV相关的转录本',' 变异情况描述','Gencode注释','CpG岛名称','染色体区段','变异相关的microRNA和snoRNA','预测microRNA的靶点','变异保守性预测','转录因子结合位点保守分值','重复片段','dgv数据库','gwasCatalog数据库','重复序列','基因组中区域的功能注释','基因组中区域的功能注释','基因组中区域的功能注释','基因组中区域的功能注释','基因组中区域的功能注释','基因组中区域的功能注释','拷贝数个数','发生CNV的区域的大小','CNV的编号','CNV的类型']
                with open(cnv_stat_file, 'r') as f:
                    for k in islice(f, 1, None):
                        k = k.strip()
                        k = k.split('\t')
                        if k[1] == 'gain':
                            cnv_datagai[
                                m2] = ':[{name:"gain",value:' + k[2] + '},'
                            cnv_gai[
                                m2] = ':[{name:"cds",value:' + k[3] + '},{name:"splicing",value:' + k[4] + '},{name:"utr5",value:' + k[5] + '},{name:"utr3",value:' + k[6] + '},{name:"intron",value:' + k[7] + '},{name:"upstream",value:' + k[8] + '},{name:"downstream",value:' + k[9] + '},{name:"ncRNA",value:' + k[10] + '},{name:"intergenic",value:' + k[11] + '},{name:"unknown",value:' + k[12] + '}],'
                        elif k[1] == 'loss':
                            cnv_datalos[
                                m2] = '{name:"loss",value:' + k[2] + '}],'
                            cnv_los[
                                m2] = ':[{name:"cds",value:' + k[3] + '},{name:"splicing",value:' + k[4] + '},{name:"utr5",value:' + k[5] + '},{name:"utr3",value:' + k[6] + '},{name:"intron",value:' + k[7] + '},{name:"upstream",value:' + k[8] + '},{name:"downstream",value:' + k[9] + '},{name:"ncRNA",value:' + k[10] + '},{name:"intergenic",value:' + k[11] + '},{name:"unknown",value:' + k[12] + '}],'

            if 'cnv_call_freec' in ANALYSIS:
                cnv_stat_file = os.path.join(projdir, 'SV', m2, 'freec',
                                             m2 + '.freec.ann.stat.xls')
                cnv_stat.append(cnv_stat_file)
                cnv_anno_file = os.path.join(projdir, 'SV', m2, 'freec',
                                             m2 + '.freec.hg19_multianno.xls')
                cnv_anno.append(cnv_anno_file)
                cnv_stat_tmp = svtxts2tab(cnv_stat)
                #del cnv_stat_tmp[1][1][0]
                context['freec'] = True
                # assert not os.system(
                #     "python %s/Varition/CNV/freec/Chr_CNV_freec_pipe4.5.py --i %s --ref %s --sample_info %s"
                #     % (moduledir, cnv_anno_file, ref, sam))
                assert not os.system('convert -resize 800 %s %s' %
                                     (os.path.join(projdir, 'SV', m2, 'freec',
                                                   m2 + '.Chr_CNV.png'),
                                      os.path.join(odir, 'src/pictures/Circos',
                                                   m2 + '.Chr_CNV.jpg')))
                context['cnvnator'] = False
                context['table_cnv'] = cnv_stat_tmp
                context['table_cnv_anno'] = txts2tab([cnv_anno[0]])[:6]
                #context['table_cnv_anno_head']=['染色体','CNV起始位置','CNV终止位置','基因名称','CNV所在的区域','CNV相关的转录本',' 变异情况描述','Gencode注释','CpG岛名称','染色体区段','变异相关的microRNA和snoRNA','预测microRNA的靶点','变异保守性预测','转录因子结合位点保守分值','重复片段','dgv数据库','gwasCatalog数据库','重复序列','基因组中区域的功能注释','基因组中区域的功能注释','基因组中区域的功能注释','基因组中区域的功能注释','基因组中区域的功能注释','基因组中区域的功能注释','发生CNV的区域的大小','cnv区域内标准化的深度信号','落在该CNV区域的reads中，mapping quality为0的reads比例。','CNV的ID号','CNV的类型']
                with open(cnv_stat_file, 'r') as f:
                    for k in islice(f, 1, None):
                        k = k.strip()
                        k = k.split('\t')
                        if k[1] == 'gain':
                            cnv_datagai[
                                m2] = ':[{name:"gain",value:' + k[2] + '},'
                            cnv_gai[
                                m2] = ':[{name:"cds",value:' + k[3] + '},{name:"splicing",value:' + k[4] + '},{name:"utr5",value:' + k[5] + '},{name:"utr3",value:' + k[6] + '},{name:"intron",value:' + k[7] + '},{name:"upstream",value:' + k[8] + '},{name:"downstream",value:' + k[9] + '},{name:"ncRNA",value:' + k[10] + '},{name:"intergenic",value:' + k[11] + '},{name:"unknown",value:' + k[12] + '}],'
                        elif k[1] == 'loss':
                            cnv_datalos[
                                m2] = '{name:"loss",value:' + k[2] + '}],'
                            cnv_los[
                                m2] = ':[{name:"cds",value:' + k[3] + '},{name:"splicing",value:' + k[4] + '},{name:"utr5",value:' + k[5] + '},{name:"utr3",value:' + k[6] + '},{name:"intron",value:' + k[7] + '},{name:"upstream",value:' + k[8] + '},{name:"downstream",value:' + k[9] + '},{name:"ncRNA",value:' + k[10] + '},{name:"intergenic",value:' + k[11] + '},{name:"unknown",value:' + k[12] + '}],'

                ##circos pic   add by chenyanli
            if ('cnv_call_cnvnator' in ANALYSIS) and ('cnv_call_freec' not in ANALYSIS):
                cir_dir = os.path.join(projdir, 'SV', m2, 'cnvnator', 'Circos')
                cir_images_dir = os.path.join(odir, 'src/pictures/circos')
                #assert not os.system('cp %s %s' % (os.path.join(cir_dir,m2+'.png'),os.path.join(cir_images_dir,m2+'.png')))
                assert not os.system(
                    'convert -resize 800 %s %s' %
                    (os.path.join(cir_dir, m2 + '.png'),
                     os.path.join(cir_images_dir, m2 + '.jpg')))
                assert not os.system(
                    'convert -resize 90 %s %s' %
                    (os.path.join(cir_dir, m2 + '.png'),
                     os.path.join(cir_images_dir, m2 + '.JPEG')))
                context['fig_cnv_circos'].append([
                    '"' + os.path.join('src/pictures/circos', m2 + '.jpg') +
                    '"', '"' +
                    os.path.join('src/pictures/circos', m2 + '.JPEG') + '"', m2
                ])
            if 'cnv_call_freec' in ANALYSIS:
                cir_dir = os.path.join(projdir, 'SV', m2, 'freec', 'Circos')
                cir_images_dir = os.path.join(odir, 'src/pictures/Circos')
                #assert not os.system('cp %s %s' % (os.path.join(cir_dir,m2+'.png'),os.path.join(cir_images_dir,m2+'.png')))
                context['fig_CNV'].append([
                    '"' + os.path.join('src/pictures/Circos',
                                       m2 + '.Chr_CNV.jpg') + '"', m2
                ])
                assert not os.system(
                    'convert -resize 800 %s %s' %
                    (os.path.join(cir_dir, m2 + '.png'),
                     os.path.join(cir_images_dir, m2 + '.jpg')))
                assert not os.system(
                    'convert -resize 90 %s %s' %
                    (os.path.join(cir_dir, m2 + '.png'),
                     os.path.join(cir_images_dir, m2 + '.JPEG')))
                context['fig_cnv_circos'].append([
                    '"' + os.path.join('src/pictures/Circos', m2 + '.jpg') +
                    '"', '"' +
                    os.path.join('src/pictures/Circos', m2 + '.JPEG') + '"', m2
                ])
            ###added by yincan 20171209
            ############################# WES cnv_stat.xls
            if 'cnv_call_conifer' in ANALYSIS:
                cnv_stat_file = os.path.join(projdir, 'SV', m2, 'conifer',
                                             m2 + '.conifer.stat.xls')
                cnv_stat.append(cnv_stat_file)
                cnv_stat_tmp = svtxts2tab(cnv_stat)
                context['table_cnv'] = cnv_stat_tmp
                context['WES_cnv_CoNIFER'] = True
                with open(cnv_stat_file, 'r') as f:
                    for k in islice(f, 1, None):
                        k = k.strip()
                        k = k.split('\t')
                        if k[1] == 'dup':
                            cnv_datagai[
                                m2] = ':[{name:"dup",value:' + k[2] + '},'
                            cnv_gai[
                                m2] = ':[{name:"cds",value:' + k[3] + '},{name:"splicing",value:' + k[4] + '},{name:"utr5",value:' + k[5] + '},{name:"utr3",value:' + k[6] + '},   {name:"intron",value:' + k[7] + '},{name:"upstream",value:' + k[8] + '},{name:"downstream",value:' + k[9] + '},{name:"ncRNA",value:' + k[10] + '},{name:"intergenic",value:' + k[11] + '},{name:"unknown",value:' + k[12] + '}],'
                        elif k[1] == 'del':
                            cnv_datalos[
                                m2] = '{name:"del",value:' + k[2] + '}],'
                            cnv_los[
                                m2] = ':[{name:"cds",value:' + k[3] + '},{name:"splicing",value:' + k[4] + '},{name:"utr5",value:' + k[5] + '},{name:"utr3",value:' + k[6] + '},   {name:"intron",value:' + k[7] + '},{name:"upstream",value:' + k[8] + '},{name:"downstream",value:' + k[9] + '},{name:"ncRNA",value:' + k[10] + '},{name:"intergenic",value:' + k[11] + '},{name:"unknown",value:' + k[12] + '}],'
            ###########################
            if 'cnv_call_conifer' in ANALYSIS:
                Conifer_dir = os.path.join(projdir, 'SV', m2, 'conifer')
                Conifer_pic = os.path.join(Conifer_dir, 'pics')
                cir_images_dir = os.path.join(odir, 'src/pictures/Circos')
                cnv_anno_file = os.path.join(
                    Conifer_dir, m2 + '.conifer.hg19_multianno.xls')
                cnv_anno.append(cnv_anno_file)
                # assert not os.system(
                #     "python %s/Varition/CNV/CoNIFER/conifer_v0.2.2/cnv_chrom_plot.py %s %s %s"
                #     % (moduledir, cnv_anno_file, ref, sam))

                context['WES_cnv_CoNIFER'] = True
                if os.path.exists(
                        os.path.join(Conifer_dir, m2 + '.Chr_CNV.png')):
                    assert not os.system(
                        'convert -resize 800 %s %s' %
                        (os.path.join(Conifer_dir, m2 + '.Chr_CNV.png'),
                         os.path.join(cir_images_dir, m2 + '.Chr_CNV.jpg')))
                    context['fig_CoNIFER'].append([
                        '"' + os.path.join('src/pictures/Circos',
                                           m2 + '.Chr_CNV.jpg') + '"', m2
                    ])
                    if len(txts2tab([cnv_anno[-1]])[:6]) > len(
                            context['table_cnv_anno']):
                        context['table_cnv_anno'] = txts2tab([cnv_anno[-1]])[:6]
                    ###added 20171209
                    else:
                        context['table_cnv_anno'] = txts2tab([cnv_anno[0]])[:6]
                    ###########################################################
                    #context['table_cnv_anno_head']=['染色体','CNV起始位置','CNV终止位置','基因名称','CNV所在的区域','CNV相关的转录本',' 变异情况描述','Gencode注释','CpG岛名称','染色体区段','变异相关的microRNA和snoRNA','预测microRNA的靶点','变异保守性预测','转录因子结合位点保守分值','重复片段','dgv数据库','gwasCatalog数据库','重复序列','基因组中区域的功能注释','基因组中区域的功能注释','基因组中区域的功能注释','基因组中区域的功能注释','基因组中区域的功能注释','基因组中区域的功能注释','发生CNV的区域的大小','CNV的编号','CNV的类型']

samplename = []
samplelane = {}
for i in context['table_qc']:

    samplename.append(i[0][0])
    samplelane[i[0][0]] = []
    for j in i:
        samplelane[i[0][0]].append(j[0] + '_' + j[1] + '_' + j[2])
for i in context['table_qc']:
    qu = len(i)
    if qu != 1:
        context['qukey'] = ""
    else:
        context['qukey'] = "qua"
samplename1 = samplename[0]
samplelane1 = samplelane[samplename[0]][0]

context['FilterSV_anno'] = []
context['FilterCNV_anno'] = []

## advance
## not ready
if rep_ty == 'advance':
    if 'filter_db' in ANALYSIS:
        filter_stat = []
        vcf_dir = os.path.join(advanceDir, 'Merged_vcf')
        vcf_bri_dir = os.path.join(advanceBriDir, 'FilterDB')
        filter_stat_file = os.path.join(
            vcf_dir, 'Filter/SNP/snp.merged.filter.stat.xls')
        filter_stat.append(filter_stat_file)
        context['filterDB'] = True
        context['table_snp_filter'] = txts2tab(filter_stat)
        context['table_merge_anno'] = txts2tab([
            os.path.join(vcf_bri_dir,
                         "snp.merged.freq.func.syn.deleterious.brief.xls")
        ])[:6]

    if 'filter_sv' in ANALYSIS:
        for i in samplename:
            fsv = [
                os.path.join(advanceBriDir, 'FilterSV',
                             i + '.LikelyDeleterious.SV.brief.xls')
            ]
            context['FilterSV'] = True
            if os.path.exists(fsv[0]):
                if len(txts2tab(fsv)[:6]) > len(context['FilterSV_anno']):
                    context['FilterSV_anno'] = txts2tab(fsv)[:6]

    if 'filter_cnv' in ANALYSIS:
        for i in samplename:
            fcnv = [
                os.path.join(advanceBriDir, 'FilterCNV',
                             i + '.LikelyDeleterious.CNV.brief.xls')
            ]
            context['FilterCNV'] = True
            if os.path.exists(fcnv[0]):
                if len(txts2tab(fcnv)[:6]) > len(context['FilterCNV_anno']):
                    context['FilterCNV_anno'] = txts2tab(fcnv)[:6]
                    # print fcnv
                    # print context['FilterCNV_anno'][0]
                    # exit()
    if 'ppi' in ANALYSIS:
        ppi_geneint_file = [
            os.path.join(advanceDir, 'PPI', 'Gene_interactions.xls')
        ]
        context['ppinter'] = True
        if os.path.exists(ppi_geneint_file[0]):
            context['table_ppi_inter'] = txts2tab(ppi_geneint_file)[:6]

    if 'hpa' in ANALYSIS:
        context['HPA'] = True
        if hpa_mode == 'brief':
            context['HPA_brief'] = True
            hpa_file = os.path.join(advanceDir, 'HPA', 'CandidateGene.Expression.Annotation.Brief.xls')
        else:
            context['HPA_complex'] = True
            hpa_file = os.path.join(advanceDir, 'HPA', 'CandidateGene.Expression.Annotation.Brief.xls')
        if os.path.exists(hpa_file):
            context['table_hpa']=txts2tab(hpa_file)[:6]

    if 'filter_drug' in ANALYSIS:
        drug_anno_file = [
            os.path.join(advanceDir, 'Pharmacogenomics',
                         'Pharmacogenetics.xls')
        ]
        context['drug'] = True
        if os.path.exists(drug_anno_file[0]):
            context['drug_anno'] = txts2tab(drug_anno_file)[:6]

    if 'site_association' in ANALYSIS:
        #context['asso']=True
        context['asso_site'] = True
        asso_file = [
            os.path.join(advanceDir, 'SiteAS', 'AssoAna_Allele.final.xls')
        ]
        if os.path.exists(asso_file[0]):
            context['table_asso'] = txts2tab(asso_file)[:6]
        context['fig_asso_QQ'] = []
        context['fig_asso_Manhattan'] = []
        if not os.path.exists(odir + '/src/pictures/SiteAS/'):
            os.mkdir(odir + '/src/pictures/SiteAS/')
        context['fig_asso_QQ'].append([
            "src/pictures/SiteAS/AssoAna_Allele_qqplot_Pvalue_Allele.png", "QQ"
        ])
        context['fig_asso_Manhattan'].append([
            "src/pictures/SiteAS/AssoAna_Allele_ManhattanPlot_Pvalue_Allele.png",
            "Manhattan"
        ])

        assert not os.system('cp ' + advanceDir +
                             '/SiteAS/AssoAna_Allele_qqplot_Pvalue_Allele.png '
                             + odir + '/src/pictures/SiteAS/')
        assert not os.system(
            'cp ' + advanceDir +
            '/SiteAS/AssoAna_Allele_ManhattanPlot_Pvalue_Allele.png ' + odir +
            '/src/pictures/SiteAS/')

    if 'gene_association' in ANALYSIS:
        context['GeneAS'] = True
        GeneAS_file = [
            os.path.join(advanceDir, 'GeneAS', 'snp.burden.result.xls')
        ]
        if os.path.exists(GeneAS_file[0]):
            context['table_GeneAS_anno'] = txts2tab(GeneAS_file)[:6]
        context['fig_GeneAS_QQ'] = []
        os.popen(
            'convert ' + advanceDir + '/GeneAS/QQ_plot.snp.burden.snp.pdf ' +
            advanceDir + '/GeneAS/QQ_plot.snp.burden.snp.png')
        if not os.path.exists(odir + '/src/pictures/GeneAS/'):
            os.mkdir(odir + '/src/pictures/GeneAS/')
        context['fig_GeneAS_QQ'].append(
            ["src/pictures/GeneAS/QQ_plot.snp.burden.snp.png", "QQ"])

        assert not os.system(
            'cp ' + advanceDir + '/GeneAS/QQ_plot.snp.burden.snp.png ' + odir +
            '/src/pictures/GeneAS/')

    if 'pathway' in ANALYSIS:
        context['fig_pathway_go'] = []
        context['fig_pathway_kegg'] = []
        if not os.path.exists(odir + '/src/pictures/Pathway/'):
            os.mkdir(odir + '/src/pictures/Pathway/')
            os.mkdir(odir + '/src/pictures/Pathway/png')
        if os.popen('ls ' + advanceDir + '/Pathway/KEGG_maps/*.svg').readlines(
        ) != []:
            keggsvg = os.popen('ls ' + advanceDir + '/Pathway/KEGG_maps/*.svg'
                               ).readlines()[-1].split('/')[-1].strip()
            keggpng = os.popen(
                'ls ' + advanceDir + '/Pathway/KEGG_maps/png/*.png').readlines(
                )[-1].split('/')[-1].strip()

            context['svgke'] = keggsvg
            context['pngke'] = keggpng
            cmd = '''
                cp {advanceDir}/Pathway/KEGG_maps/{keggsvg} {odir}/src/pictures/Pathway/
                cp {advanceDir}/Pathway/KEGG_maps/png/{keggpng} {odir}/src/pictures/Pathway/png/
                '''.format(**locals())
            # print cmd
            assert not os.system(cmd)
            context['Pathway_kegg'] = True
        else:
            context['Pathway_kegg'] = False
        context['fig_pathway_kegg'].append(
            ["src/pictures/Pathway/Pathway_kegg_dotplot.png", "KEGG"])
        context['fig_pathway_go'].append(
            ["src/pictures/Pathway/Pathway_cc_dotplot.png", "GO_CC"])
        context['fig_pathway_go'].append(
            ["src/pictures/Pathway/Pathway_bp_dotplot.png", "GO_BP"])
        context['fig_pathway_go'].append(
            ["src/pictures/Pathway/Pathway_mf_dotplot.png", "GO_MF"])

        assert not os.system(
            'cp ' + advanceDir + '/Pathway/Pathway/Pathway_kegg_dotplot.png ' +
            odir + '/src/pictures/Pathway/')
        assert not os.system(
            'cp ' + advanceDir + '/Pathway/Pathway/Pathway_cc_dotplot.png ' +
            odir + '/src/pictures/Pathway/')
        assert not os.system(
            'cp ' + advanceDir + '/Pathway/Pathway/Pathway_bp_dotplot.png ' +
            odir + '/src/pictures/Pathway/')
        assert not os.system(
            'cp ' + advanceDir + '/Pathway/Pathway/Pathway_mf_dotplot.png ' +
            odir + '/src/pictures/Pathway/')

with open(odir + "/qcstat.xls", "w") as f:
    f.write(
        'Sample name' + '\t' + 'Novo ID' + '\t' + 'Lane' + '\t' + 'Raw reads' +
        '\t' + 'Raw data(G)' + '\t' + 'Raw depth(x)' + '\t' + 'Effective' +
        '\t' + 'Error' + '\t' + 'Q20' + '\t' + 'Q30' + '\t' + 'GC' + '\n')
    for i in context['table_qc']:
        for j in i:
            f.write(j[0] + '\t' + j[1] + '\t' + j[2] + '\t' +
                    '\t'.join(qc_sta[j[0] + '_' + j[1] + '_' + j[2]]) + '\n')

names = ''
for i in samplename:
    names += '"' + i + '",'
names = '[' + names[:-1] + ']'
# print names
divlen = len('*'.join(samplename))
context['divlen'] = (divlen / 100 + 2) * 25
sampnum = len(samplename)
datanow2 = ''
ermd = ''
erzd = ''
gczd = ''
chrdata = ''
legendata = ''
depthcumu = ''
depthdist = ''
legendallclose = 'e.setOption({legend:{selected:{'
legendallopen = 'e.setOption({legend:{selected:{'
deptha = {}
cova = {}
if rep_ty != 'qc':
    context['table_dep'] = ((np.array(context['table_depth'])).T).tolist()
    with open(odir + "/mappingstat.xls", 'w') as f:
        for i in context['table_dep']:
            for n in i:
                f.write(n + '\t')
            f.write('\n')
    with open(odir + "/src/pictures/depth/covByChr.png.1.txt", 'r') as f:
        for k in islice(f, 1, None):
            k = k.strip()
            name1 = k.split('\t')[0]
            depth1 = '","'.join(k.split('\t')[1:])
            deptha[
                name1] = '{"name":"' + name1 + '","type":"bar","yAxisIndex": "0","data":["' + depth1 + '"]},'
    with open(odir + "/src/pictures/depth/covByChr.png.2.txt", 'r') as f:
        for k in islice(f, 1, None):
            k = k.strip()
            name1 = k.split('\t')[0]
            cov1 = '","'.join(k.split('\t')[1:])
            cova[
                name1] = '{"name":"' + name1 + '","type":"line","yAxisIndex": "1","data":["' + cov1 + '"]},'
    for i in samplename:
        chrdata += cova[i] + deptha[i]
        legendata += '{name:"' + i + '",type:"scatter",data:0},'
    chrdata = '[' + chrdata[:-1] + ']'
    legendata = '[' + legendata[:-1] + ']'
    for i in samplename:
        legendallclose += '"' + i + '":false,'
        legendallopen += '"' + i + '":true,'
        cumudep = ''
        distdep = ''
        with open(projdir + '/Alnstat/' + i + '/cumu.xls', 'r') as f:
            for k in islice(f, 1, None):
                k = k.strip().split('\t')
                if (string.atoi(k[0]) > depthavg):
                    break
                else:
                    cumudep += '[' + k[0] + ',' + k[1] + '],'
        with open(projdir + '/Alnstat/' + i + '/depth_frequency.xls',
                  'r') as f:
            for k in f:
                k = k.strip().split('\t')
                if (string.atoi(k[0]) > depthavg):
                    break
                else:
                    distdep += '[' + k[0] + ',' + k[1] + '],'
        if sampnum > 10:
            depthcumu += '{name:"' + i + '",type:"scatter",symbolSize:7,large:1,largeThreshold:100,data:[' + cumudep[:
                                                                                                                     -1] + ']},'
            depthdist += '{name:"' + i + '",type:"scatter",symbolSize:7,large:1,largeThreshold:100,data:[' + distdep[:
                                                                                                                     -1] + ']},'
        else:
            depthcumu += '{name:"' + i + '",type:"scatter",symbolSize:7,large:1,largeThreshold:1000,data:[' + cumudep[:
                                                                                                                      -1] + ']},'
            depthdist += '{name:"' + i + '",type:"scatter",symbolSize:7,large:1,largeThreshold:1000,data:[' + distdep[:
                                                                                                                      -1] + ']},'
    legendallclose = legendallclose[:-1] + '}}});'
    legendallopen = legendallopen[:-1] + '}}});'
    depthcumu = '[' + depthcumu[:-1] + ']'
    depthdist = '[' + depthdist[:-1] + ']'

snvexoniz1 = ''
indelexoniz1 = ''
snvfuncz1 = ''
indelfuncz1 = ''
indel_newoldz = ''
indel_homhetz = ''
snp_newoldz = ''
snp_homhetz = ''
snp_traitz = ''
snp_tstvz = ''
snp_newtstvz = ''
sv_dataz = ''
sv_data2z = ''
cnv_dataz = ''
sv_insz = ''
sv_invz = ''
sv_delz = ''
sv_traz = ''
sv_ins2z = ''
sv_inv2z = ''
sv_del2z = ''
sv_itx2z = ''
sv_ctx2z = ''
sv_tdp2z = ''
cnv_gaiz = ''
cnv_losz = ''
indel_lenz = ''
indel_lennz = ''
if rep_ty == 'primary' or rep_ty == 'advance':
    if 'snpindel_call' in ANALYSIS:
        if seq_ty != 'TS':
            for i in samplename:
                if ('snpindel_call_samtools' in ANALYSIS
                        or 'snpindel_call_samtools_multi' in ANALYSIS):
                    ge_vcf_file = os.path.join(
                        projdir, "Mutation", i + '.samtools',
                        i + '.samtools.indel.annovar.hg19_multianno.xls.gz')
                elif 'snpindel_call_gatk' in ANALYSIS:
                    ge_vcf_file = os.path.join(
                        projdir, "Mutation", i + '.gatk',
                        i + '.gatk.indel.annovar.hg19_multianno.xls.gz')
                elif 'snpindel_call_sentieon' in ANALYSIS:
                    ge_vcf_file = os.path.join(
                        projdir, "Mutation", i + '.sentieon',
                        i + '.sentieon.indel.annovar.hg19_multianno.xls.gz')

                A = utils.safe_open(ge_vcf_file, 'r')
                cds, ncds = {}, {}
                for n in A:
                    if n.startswith('#'): continue
                    l = n.strip().split('\t')
                    x = len(l[5]) - len(l[4])
                    m = l[10]
                    if m == '.': continue
                    pos = m.split(';')
                    if 'exonic' in pos or 'splicing' in pos:
                        if x not in cds:
                            cds[x] = 1
                        else:
                            cds[x] += 1
                    else:
                        if x not in ncds:
                            ncds[x] = 1
                        else:
                            ncds[x] += 1
                A.close()

                for k in xrange(-30, 31):
                    if k not in cds:
                        cds[k] = 0
                    if k not in ncds:
                        ncds[k] = 0

                nsum, csum = 0, 0
                for k in xrange(-30, 31):
                    if k in ncds:
                        nsum += ncds[k]
                    if k in cds:
                        csum += cds[k]
                x, g_y, e_y = [], [], []
                for k in xrange(-30, 31):
                    if k in ncds:
                        x.append(k)
                        if nsum == 0:
                            g_y_value == 0
                            g_y.append(g_y_value)
                        else:
                            g_y_value = float(ncds[k]) / nsum
                            g_y.append(g_y_value)
                        if csum == 0:
                            e_y_value == 0
                            e_y.append(e_y_value)
                        else:
                            e_y_value = float(cds[k]) / csum
                            e_y.append(e_y_value)
                indel_lenz += '"' + i + '":[' + ','.join(map(str, e_y)) + '],'
                indel_lennz += '"' + i + '":[' + ','.join(map(str, g_y)) + '],'
            indel_lenz = '{' + indel_lenz[:-1] + '}'
            indel_lennz = '{' + indel_lennz[:-1] + '}'
        for i in context['table_exfunc_snp'][1:]:
            snvexoniz1 += "'" + i[0] + "':[{name:'" + context['table_exfunc_snp'][0][1] + "',value:" + i[1] + "},{name:'" + context['table_exfunc_snp'][0][2] + "',value:" + i[2] + "},{name:'" + context['table_exfunc_snp'][0][3] + "',value:" + i[3] + "},{name:'" + context['table_exfunc_snp'][0][4] + "',value:" + i[4] + "},{name:'" + context['table_exfunc_snp'][0][5] + "',value:" + i[5] + "}],"
        snvexoniz1 = '{' + snvexoniz1[:-1] + '}'
        for i in context['table_exfunc_indel'][1:]:
            indelexoniz1 += "'" + i[0] + "':[{name:'" + context['table_exfunc_indel'][0][1] + "',value:" + i[1] + "},{name:'" + context['table_exfunc_indel'][0][2] + "',value:" + i[2] + "},{name:'" + context['table_exfunc_indel'][0][3] + "',value:" + i[3] + "},{name:'" + context['table_exfunc_indel'][0][4] + "',value:" + i[4] + "},{name:'" + context['table_exfunc_indel'][0][5] + "',value:" + i[5] + "},{name:'" + context['table_exfunc_indel'][0][6] + "',value:" + i[6] + "},{name:'" + context['table_exfunc_indel'][0][7] + "',value:" + i[7] + "}],"
        indelexoniz1 = '{' + indelexoniz1[:-1] + '}'
        for i in context['table_func_snp'][1:]:
            snvfuncz1 += "'" + i[0] + "':[{name:'" + context['table_func_snp'][0][1] + "',value:" + i[1] + "},{name:'" + context['table_func_snp'][0][2] + "',value:" + i[2] + "},{name:'" + context['table_func_snp'][0][3] + "',value:" + i[3] + "},{name:'" + context['table_func_snp'][0][4] + "',value:" + i[4] + "},{name:'" + context['table_func_snp'][0][5] + "',value:" + i[5] + "},{name:'" + context['table_func_snp'][0][6] + "',value:" + i[6] + "},{name:'" + context['table_func_snp'][0][7] + "',value:" + i[7] + "},{name:'" + context['table_func_snp'][0][8] + "',value:" + i[8] + "},{name:'" + context['table_func_snp'][0][9] + "',value:" + i[9] + "},{name:'" + context['table_func_snp'][0][10] + "',value:" + i[10] + "},{name:'" + context['table_func_snp'][0][11] + "',value:" + i[11] + "}],"
        snvfuncz1 = '{' + snvfuncz1[:-1] + '}'
        for i in context['table_func_indel'][1:]:
            indelfuncz1 += "'" + i[0] + "':[{name:'" + context['table_func_indel'][0][1] + "',value:" + i[1] + "},{name:'" + context['table_func_indel'][0][2] + "',value:" + i[2] + "},{name:'" + context['table_func_indel'][0][3] + "',value:" + i[3] + "},{name:'" + context['table_func_indel'][0][4] + "',value:" + i[4] + "},{name:'" + context['table_func_indel'][0][5] + "',value:" + i[5] + "},{name:'" + context['table_func_indel'][0][6] + "',value:" + i[6] + "},{name:'" + context['table_func_indel'][0][7] + "',value:" + i[7] + "},{name:'" + context['table_func_indel'][0][8] + "',value:" + i[8] + "},{name:'" + context['table_func_indel'][0][9] + "',value:" + i[9] + "},{name:'" + context['table_func_indel'][0][10] + "',value:" + i[10] + "},{name:'" + context['table_func_indel'][0][11] + "',value:" + i[11] + "}],"
        indelfuncz1 = '{' + indelfuncz1[:-1] + '}'
        for i in context['table_trait_indel'][1:]:
            indel_homhetz += '"' + i[0] + '":[{name:"Hom",value:' + i[3] + '},{name:"Het",value:' + i[2] + '}],'
            indel_newoldz += '"' + i[0] + '":[{name:"novel",value:' + i[4] + '},{name:"in dbSNP",value:' + str(
                string.atoi(i[1]) - string.atoi(i[4])) + '}],'
        indel_homhetz = '{' + indel_homhetz[:-1] + '}'
        indel_newoldz = '{' + indel_newoldz[:-1] + '}'
        for i in context['table_trait_snp'][1:]:
            snp_homhetz += '"' + i[0] + '":[{name:"Hom",value:' + i[3] + '},{name:"Het",value:' + i[2] + '}],'
            snp_newoldz += '"' + i[0] + '":[{name:"novel",value:' + i[4] + '},{name:"in dbSNP",value:' + str(
                string.atoi(i[1]) - string.atoi(i[4])) + '}],'
        snp_homhetz = '{' + snp_homhetz[:-1] + '}'
        snp_newoldz = '{' + snp_newoldz[:-1] + '}'
        for i in context['table_ts_tv'][1:]:
            snp_tstvz += '"' + i[0] + '":[{name:"ts",value:' + i[4] + '}, {name:"tv",value:' + i[6] + '}],'
            snp_newtstvz += '"' + i[0] + '":[{name:"novel ts",value:' + i[1] + '}, {name:"novel tv",value:' + i[3] + '}],'
        snp_tstvz = '{' + snp_tstvz[:-1] + '}'
        snp_newtstvz = '{' + snp_newtstvz[:-1] + '}'

    if 'sv_call' in ANALYSIS:
        if 'sv_call_breakdancer' in ANALYSIS:
            for i in samplename:
                sv_dataz += '"' + i + '"' + sv_datatra[i] + sv_datains[i] + sv_datadel[i] + sv_datainv[i]
                sv_insz += '"' + i + '"' + sv_ins[i]
                sv_invz += '"' + i + '"' + sv_inv[i]
                sv_delz += '"' + i + '"' + sv_del[i]
                sv_traz += '"' + i + '"' + sv_tra[i]
            sv_dataz = '{' + sv_dataz[:-1] + '}'
            sv_traz = '{' + sv_traz[:-1] + '}'
            sv_insz = '{' + sv_insz[:-1] + '}'
            sv_delz = '{' + sv_delz[:-1] + '}'
            sv_invz = '{' + sv_invz[:-1] + '}'
        if 'sv_call_crest' in ANALYSIS:
            for i in samplename:
                sv_data2z += '"' + i + '"' + sv_datains2[i] + sv_datadel2[i] + sv_datainv2[i] + sv_dataitx2[i] + sv_datactx2[i]
                sv_ins2z += '"' + i + '"' + sv_ins2[i]
                sv_del2z += '"' + i + '"' + sv_del2[i]
                sv_inv2z += '"' + i + '"' + sv_inv2[i]
                sv_itx2z += '"' + i + '"' + sv_itx2[i]
                sv_ctx2z += '"' + i + '"' + sv_ctx2[i]
            sv_data2z = '{' + sv_data2z[:-1] + '}'
            sv_itx2z = '{' + sv_itx2z[:-1] + '}'
            sv_ctx2z = '{' + sv_ctx2z[:-1] + '}'
            sv_ins2z = '{' + sv_ins2z[:-1] + '}'
            sv_del2z = '{' + sv_del2z[:-1] + '}'
            sv_inv2z = '{' + sv_inv2z[:-1] + '}'
        if 'sv_call_lumpy' in ANALYSIS:
            for i in samplename:
                sv_data2z += '"' + i + '"' + sv_datains2[i] + sv_datadel2[i] + sv_datainv2[i] + sv_dataitx2[i] + sv_datactx2[i]
                sv_ins2z += '"' + i + '"' + sv_ins2[i]
                sv_del2z += '"' + i + '"' + sv_del2[i]
                sv_inv2z += '"' + i + '"' + sv_inv2[i]
                sv_itx2z += '"' + i + '"' + sv_itx2[i]
                sv_ctx2z += '"' + i + '"' + sv_ctx2[i]
            sv_data2z = '{' + sv_data2z[:-1] + '}'
            sv_itx2z = '{' + sv_itx2z[:-1] + '}'
            sv_ctx2z = '{' + sv_ctx2z[:-1] + '}'
            sv_ins2z = '{' + sv_ins2z[:-1] + '}'
            sv_del2z = '{' + sv_del2z[:-1] + '}'
            sv_inv2z = '{' + sv_inv2z[:-1] + '}'
        if 'sv_call_breakmer' in ANALYSIS:
            for i in samplename:
                sv_data2z += '"' + i + '"' + sv_datains2[i] + sv_datadel2[i] + sv_datainv2[i] + sv_dataitx2[i] + sv_datactx2[i] + sv_datatdp2[i]
                sv_ins2z += '"' + i + '"' + sv_ins2[i]
                sv_del2z += '"' + i + '"' + sv_del2[i]
                sv_inv2z += '"' + i + '"' + sv_inv2[i]
                sv_itx2z += '"' + i + '"' + sv_itx2[i]
                sv_ctx2z += '"' + i + '"' + sv_ctx2[i]
                sv_tdp2z += '"' + i + '"' + sv_tdp2[i]
            sv_data2z = '{' + sv_data2z[:-1] + '}'
            sv_itx2z = '{' + sv_itx2z[:-1] + '}'
            sv_ctx2z = '{' + sv_ctx2z[:-1] + '}'
            sv_ins2z = '{' + sv_ins2z[:-1] + '}'
            sv_del2z = '{' + sv_del2z[:-1] + '}'
            sv_inv2z = '{' + sv_inv2z[:-1] + '}'
            sv_tdp2z = '{' + sv_tdp2z[:-1] + '}'
    if 'cnv_call' in ANALYSIS:
        for i in samplename:
            cnv_dataz += '"' + i + '"' + cnv_datagai[i] + cnv_datalos[i]
            cnv_gaiz += '"' + i + '"' + cnv_gai[i]
            cnv_losz += '"' + i + '"' + cnv_los[i]
        cnv_dataz = '{' + cnv_dataz[:-1] + '}'
        cnv_gaiz = '{' + cnv_gaiz[:-1] + '}'
        cnv_losz = '{' + cnv_losz[:-1] + '}'

for i in samplename:
    for j in samplelane[i]:
        num = ['0']
        halfline = 0

        cmd = string.Template('''
            tail -1 ${projdir}/QC/${i}/${j}.stat |
            awk -F '\\t' -v OFS='\\t' '{print $2, $3, $4, $5, $6}' \\
            > ${odir}/src/qcclean.stat
        ''').safe_substitute(**locals())

        # print cmd
        assert not os.system(cmd)

        with open(odir + "/src/qcclean.stat", 'r') as f:
            for k in f:
                k = k.strip()
                k = k.split('\t')
                n = "'" + j + "':[{name:'Clean reads',value:" + k[1] + "},{name:'Containing N',value:" + k[2] + "},{name:'Low quality',value:" + k[3] + "},{name:'Adapter related',value:" + k[4] + "}],"
                datanow2 += n
        with open(projdir + "/QC/" + i + "/clean_" + j + ".QM") as f:
            for k in f:
                halfline += 1
                num.append(str(halfline))

        assert not os.system(
            "cd " + projdir + "/QC/" + i + "&& awk -F '\t' '{print $3}' clean_"
            + j + ".QM|tr -s '\\n' ',' >" + odir + "/src/qcclean.stat")
        with open(odir + "/src/qcclean.stat", 'r') as f:
            for k in f:
                k = k.strip()
        ermd += '"' + j + '":[' + k[:-1] + '],'

        assert not os.system(
            "cd " + projdir + "/QC/" + i + "&& awk -F '\t' '{print $2}' clean_"
            + j + ".QM|tr -s '\\n' ',' >" + odir + "/src/qcclean.stat")
        with open(odir + "/src/qcclean.stat", 'r') as f:
            for k in f:
                k = k.strip()
        erzd += '"' + j + '":[' + k[:-1] + '],'

        assert not os.system(
            "cd " + projdir + "/QC/" + i + "&& awk -F '\t' '{print $4}' clean_"
            + j + ".GC|tr -s '\\n' ',' >" + odir + "/src/qcclean.stat")
        with open(odir + "/src/qcclean.stat", 'r') as f:
            for k in f:
                k = k.strip()
        gczA = "{name:'A',type:'line',data:[" + k[:-1] + "]},"

        assert not os.system(
            "cd " + projdir + "/QC/" + i + "&& awk -F '\t' '{print $7}' clean_"
            + j + ".GC|tr -s '\\n' ',' >" + odir + "/src/qcclean.stat")
        with open(odir + "/src/qcclean.stat", 'r') as f:
            for k in f:
                k = k.strip()
        gczT = "{name:'T',type:'line',data:[" + k[:-1] + "]},"

        assert not os.system(
            "cd " + projdir + "/QC/" + i + "&& awk -F '\t' '{print $10}' clean_"
            + j + ".GC|tr -s '\\n' ',' >" + odir + "/src/qcclean.stat")
        with open(odir + "/src/qcclean.stat", 'r') as f:
            for k in f:
                k = k.strip()
        gczG = "{name:'G',type:'line',data:[" + k[:-1] + "]},"

        assert not os.system(
            "cd " + projdir + "/QC/" + i + "&& awk -F '\t' '{print $13}' clean_"
            + j + ".GC|tr -s '\\n' ',' >" + odir + "/src/qcclean.stat")
        with open(odir + "/src/qcclean.stat", 'r') as f:
            for k in f:
                k = k.strip()
        gczC = "{name:'C',type:'line',data:[" + k[:-1] + "]},"

        assert not os.system(
            "cd " + projdir + "/QC/" + i + "&& awk -F '\t' '{print $16}' clean_"
            + j + ".GC|tr -s '\\n' ',' >" + odir + "/src/qcclean.stat")
        with open(odir + "/src/qcclean.stat", 'r') as f:
            for k in f:
                k = k.strip()
        gczN = "{name:'N',type:'line',data:[" + k[:-1] + "]},"
        gczd += '"' + j + '":[' + gczA + gczT + gczG + gczC + gczN[:-1] + '],'
        os.system("rm -f " + odir + "/src/qcclean.stat")

gczd = "{" + gczd[:-1] + "}"
ermd = "{" + ermd[:-1] + "}"
erzd = "{" + erzd[:-1] + "}"
halfline = str(halfline)
datanow2 = "{" + datanow2[:-1] + "}"

if ('phenolyzer' in ANALYSIS) and (diseases != ''):
    net_dir = os.path.join(advanceDir, 'Network')
    assert not os.system(
        'cp -f ' + net_dir + '/network.js ' + odir + '/src/js')
    gene_dir = os.path.join(advanceDir, 'IntegrateResults',
                            'Total.candidate.gene.xls')

    DisGeNet_dir = os.path.join(advanceDir, 'Network')
    DisGeNet_file = [os.path.join(DisGeNet_dir, 'DisGeNet_shared_gene.xls')][0]

    disgenetfile = len(utils.safe_open(DisGeNet_file, 'r').readlines())
    if disgenetfile > 1:
        os.popen('awk -F "\\t" \'{print $1,$3,$5,$6,$7,$8,$9,$10}\' OFS="\\t" '
                 '' + DisGeNet_dir + '/DisGeNet_shared_gene.xls > '
                 '' + DisGeNet_dir + '/DisGeNet_final')
        DisGeNet_anno = os.path.join(advanceDir, 'Network', 'DisGeNet_final')
        context['table_DisGeNet_anno'] = txts2tab([DisGeNet_anno])[:6]
        context['disgetnet'] = True
    else:
        context['disgetnet'] = False

    if os.path.exists(gene_dir):
        genelist = ''
        with open(gene_dir, 'r') as f:
            for i in f:
                genelist += i.strip() + ','
        genelist = genelist[:-1]

    if genelist != '':
        context['genelist'] = True
    else:
        context['genelist'] = False
num = ','.join(num)

with open(odir + '/src/js/network.js', 'a') as datajs:
    datajs.write('datanow2=' + datanow2 + '\n')
    datajs.write('num=[' + num + ']\n')
    datajs.write('nbdata2 = datanow2["' + samplelane1 + '"]\n')
    datajs.write('ermd=' + ermd + '\n')
    datajs.write('ermdata = ermd["' + samplelane1 + '"]\n')
    datajs.write('erzd=' + erzd + '\n')
    datajs.write('erzdata = erzd["' + samplelane1 + '"]\n')
    datajs.write('gczd1=' + gczd + '\n')
    datajs.write('gczdata1 = gczd1["' + samplelane1 + '"]\n')
    if rep_ty != 'qc':
        datajs.write('depthcumu=' + depthcumu + '\n')
        datajs.write('depthdist=' + depthdist + '\n')
        datajs.write('chrdata=' + chrdata + '\n')
        datajs.write('legendata=' + legendata + '\n')
        datajs.write('names=' + names + '\n')
        datajs.write('depthavg=' + str(int(depthavg)) + '\n')
        datajs.write(
            'var legendallclose=function(e){' + legendallclose + '};\n')
        datajs.write('var legendallopen=function(e){' + legendallopen + '};\n')
        if rep_ty == 'primary' or rep_ty == 'advance':
            if 'snpindel_call' in ANALYSIS:
                datajs.write('snvexoniz2=' + snvexoniz1 + '\n')
                datajs.write('snvexoni2 = snvexoniz2["' + samplename1 + '"]\n')
                datajs.write('indel_lenz=' + indel_lenz + '\n')
                if seq_ty != "TS":
                    datajs.write(
                        'indel_len=indel_lenz["' + samplename1 + '"]\n')
                    datajs.write('indel_lennz=' + indel_lennz + '\n')
                    datajs.write(
                        'indel_lenn=indel_lennz["' + samplename1 + '"]\n')
                datajs.write('indelexoniz2=' + indelexoniz1 + '\n')
                datajs.write(
                    'indelexoni2 = indelexoniz2["' + samplename1 + '"]\n')
                datajs.write('indelfuncz2=' + indelfuncz1 + '\n')
                datajs.write(
                    'indelfunc2 = indelfuncz2["' + samplename1 + '"]\n')
                datajs.write('snvfuncz2=' + snvfuncz1 + '\n')
                datajs.write('snvfunc2 = snvfuncz2["' + samplename1 + '"]\n')
                datajs.write('snp_tstvz=' + snp_tstvz + '\n')
                datajs.write('snp_tstv=snp_tstvz["' + samplename1 + '"]\n')
                datajs.write('snp_newtstvz=' + snp_newtstvz + '\n')
                datajs.write(
                    'snp_newtstv=snp_newtstvz["' + samplename1 + '"]\n')
                datajs.write('snp_homhetz=' + snp_homhetz + '\n')
                datajs.write('snp_homhet=snp_homhetz["' + samplename1 + '"]\n')
                datajs.write('indel_homhetz=' + indel_homhetz + '\n')
                datajs.write(
                    'indel_homhet=indel_homhetz["' + samplename1 + '"]\n')
                datajs.write('snp_newoldz=' + snp_newoldz + '\n')
                datajs.write('snp_newold=snp_newoldz["' + samplename1 + '"]\n')
                datajs.write('indel_newoldz=' + indel_newoldz + '\n')
                datajs.write(
                    'indel_newold=indel_newoldz["' + samplename1 + '"]\n')
            if 'cnv_call' in ANALYSIS:
                datajs.write('cnv_dataz2=' + cnv_dataz + '\n')
                datajs.write('cnv_data2=cnv_dataz2["' + samplename1 + '"]\n')
                datajs.write('cnv_gaiz=' + cnv_gaiz + '\n')
                datajs.write('cnv_gai=cnv_gaiz["' + samplename1 + '"]\n')
                datajs.write('cnv_losz=' + cnv_losz + '\n')
                datajs.write('cnv_los=cnv_losz["' + samplename1 + '"]\n')
            if ('sv_call_breakdancer' in ANALYSIS) and ('sv_call_crest' not in ANALYSIS):
                datajs.write('sv_dataz2=' + sv_dataz + '\n')
                datajs.write('sv_data2=sv_dataz2["' + samplename1 + '"]\n')
                datajs.write('sv_traz=' + sv_traz + '\n')
                datajs.write('sv_tra=sv_traz["' + samplename1 + '"]\n')
                datajs.write('sv_insz=' + sv_insz + '\n')
                datajs.write('sv_ins=sv_insz["' + samplename1 + '"]\n')
                datajs.write('sv_delz=' + sv_delz + '\n')
                datajs.write('sv_del=sv_delz["' + samplename1 + '"]\n')
                datajs.write('sv_invz=' + sv_invz + '\n')
                datajs.write('sv_inv=sv_invz["' + samplename1 + '"]\n')
            if 'sv_call_crest' in ANALYSIS:
                datajs.write('sv_dataz2=' + sv_data2z + '\n')
                datajs.write('sv_data2=sv_dataz2["' + samplename1 + '"]\n')
                datajs.write('sv_itx2z=' + sv_itx2z + '\n')
                datajs.write('sv_itx2=sv_itx2z["' + samplename1 + '"]\n')
                datajs.write('sv_ins2z=' + sv_ins2z + '\n')
                datajs.write('sv_ins2=sv_ins2z["' + samplename1 + '"]\n')
                datajs.write('sv_del2z=' + sv_del2z + '\n')
                datajs.write('sv_del2=sv_del2z["' + samplename1 + '"]\n')
                datajs.write('sv_inv2z=' + sv_inv2z + '\n')
                datajs.write('sv_inv2=sv_inv2z["' + samplename1 + '"]\n')
                datajs.write('sv_ctx2z=' + sv_ctx2z + '\n')
                datajs.write('sv_ctx2=sv_ctx2z["' + samplename1 + '"]\n')
            if 'sv_call_lumpy' in ANALYSIS:
                datajs.write('sv_dataz2=' + sv_data2z + '\n')
                datajs.write('sv_data2=sv_dataz2["' + samplename1 + '"]\n')
                datajs.write('sv_itx2z=' + sv_itx2z + '\n')
                datajs.write('sv_itx2=sv_itx2z["' + samplename1 + '"]\n')
                datajs.write('sv_ins2z=' + sv_ins2z + '\n')
                datajs.write('sv_ins2=sv_ins2z["' + samplename1 + '"]\n')
                datajs.write('sv_del2z=' + sv_del2z + '\n')
                datajs.write('sv_del2=sv_del2z["' + samplename1 + '"]\n')
                datajs.write('sv_inv2z=' + sv_inv2z + '\n')
                datajs.write('sv_inv2=sv_inv2z["' + samplename1 + '"]\n')
                datajs.write('sv_ctx2z=' + sv_ctx2z + '\n')
                datajs.write('sv_ctx2=sv_ctx2z["' + samplename1 + '"]\n')
            if 'sv_call_breakmer' in ANALYSIS:
                datajs.write('sv_dataz2=' + sv_data2z + '\n')
                datajs.write('sv_data2=sv_dataz2["' + samplename1 + '"]\n')
                datajs.write('sv_itx2z=' + sv_itx2z + '\n')
                datajs.write('sv_itx2=sv_itx2z["' + samplename1 + '"]\n')
                datajs.write('sv_ins2z=' + sv_ins2z + '\n')
                datajs.write('sv_ins2=sv_ins2z["' + samplename1 + '"]\n')
                datajs.write('sv_del2z=' + sv_del2z + '\n')
                datajs.write('sv_del2=sv_del2z["' + samplename1 + '"]\n')
                datajs.write('sv_inv2z=' + sv_inv2z + '\n')
                datajs.write('sv_inv2=sv_inv2z["' + samplename1 + '"]\n')
                datajs.write('sv_ctx2z=' + sv_ctx2z + '\n')
                datajs.write('sv_ctx2=sv_ctx2z["' + samplename1 + '"]\n')
                datajs.write('sv_tdp2z=' + sv_tdp2z + '\n')
                datajs.write('sv_tdp2=sv_tdp2z["' + samplename1 + '"]\n')

#显隐性模式筛选
a = []
b = []
c = []
sam = args['sample']
if sam != 'Null':
    samf = open(sam, 'r')
    for line in samf:
        if line.startswith('#'): continue
        a = line.strip().split()
        if a[0] not in b:
            b.append(a[0])
        c.append(a[0] + '_' + a[1])

for each in b:
    modelD_anno_snp_stat = []
    modelD_anno_indel_stat = []

    if 'filter_acmg' in ANALYSIS:
        damlevel_stat = []
        damlevel_dir = os.path.join(advanceDir, 'ACMG')
        damlevel_bri_dir = os.path.join(advanceBriDir, 'ACMG')
        damlevel_stat_file = os.path.join(damlevel_dir, each, each,
                                          each + '.snp.indel.stat.xls')
        damlevel_anno_file = [
            os.path.join(damlevel_bri_dir, each,
                         each + '.snp.indel.Pathogenic.brief.xls')
        ]
        damlevel_stat.append(damlevel_stat_file)
        context['filterDB'] = True
        context['DamLevel'] = True
        context['table_damlevel_filter'] = txts2tab(damlevel_stat)
        if os.path.exists(damlevel_anno_file[0]):
            
            context['damlevel_anno'] = txts2tab(damlevel_anno_file)[:6]

    if 'model_dominant' in ANALYSIS:
        MF_Dir = os.path.join(advanceDir, 'ModelF')
        MF_bri_Dir = os.path.join(advanceBriDir, 'ModelF')
        modelD_anno_snp_file = os.path.join(MF_bri_Dir, each,
                                            each + '.snp.dominant.brief.xls')
        modelD_anno_indel_file = os.path.join(
            MF_bri_Dir, each, each + '.indel.dominant.brief.xls')
        model_snp_Dir = os.path.join(
            advanceDir, 'ModelF', each,
            each + '.snp.dominant.xls.CandidateGene.xls')
        model_indel_Dir = os.path.join(
            advanceDir, 'ModelF', each,
            each + '.indel.dominant.xls.CandidateGene.xls')
        modelD_anno_snp_stat.append(modelD_anno_snp_file)
        modelD_anno_indel_stat.append(modelD_anno_indel_file)
        snpanno = len(utils.safe_open(modelD_anno_snp_file, 'r').readlines())
        indelanno = len(utils.safe_open(modelD_anno_indel_file, 'r').readlines())
        if snpanno > 1 and indelanno > 1:
            context['ModelD'] = True
            context['table_modelD_anno'] = txts2tab(
                [modelD_anno_snp_stat[0]])[:6]
            modelD_snplist = len(utils.safe_open(model_snp_Dir, 'r').readlines())
            modelD_indellist = len(utils.safe_open(model_indel_Dir, 'r').readlines())
            context['familyid'] = each
            context['modelD_listsnp'] = modelD_snplist
            context['modelD_listindel'] = modelD_indellist
        elif snpanno > 1 and 1 >= indelanno:
            context['ModelD'] = True
            context['table_modelD_anno'] = txts2tab(
                [modelD_anno_snp_stat[0]])[:6]
            modelD_snplist = len(utils.safe_open(model_snp_Dir, 'r').readlines())
            context['familyid'] = each
            context['modelD_listsnp'] = modelD_snplist
            context['modelD_listindel'] = 0
        elif 1 >= snpanno and indelanno > 1:
            context['ModelD'] = True
            context['table_modelD_anno'] = txts2tab(
                [modelD_anno_indel_stat[0]])[:6]
            modelD_indellist = len(utils.safe_open(model_indel_Dir, 'r').readlines())
            context['familyid'] = each
            context['modelD_listsnp'] = 0
            context['modelD_listindel'] = modelD_indellist
        else:
            context['ModelD'] = False

    modelR_anno_snp_stat = []
    modelR_anno_indel_stat = []
    modelC_anno_stat = []
    if 'model_recessive' in ANALYSIS:
        MF_Dir = os.path.join(advanceDir, 'ModelF')
        MF_bri_Dir = os.path.join(advanceBriDir, 'ModelF')
        modelR_anno_snp_file = os.path.join(MF_bri_Dir, each, 
                                            each + '.snp.recessive.brief.xls')
        modelR_anno_indel_file = os.path.join(
            MF_bri_Dir, each, each + '.indel.recessive.brief.xls')
        modelC_anno_file = os.path.join(
            MF_bri_Dir, each, each + '.snp.indel.compound_heterozygous.brief.xls')
        modelR_snp_Dir = os.path.join(
            advanceDir, 'ModelF', each,
            each + '.snp.recessive.xls.CandidateGene.xls')
        modelR_indel_Dir = os.path.join(
            advanceDir, 'ModelF', each,
            each + '.indel.recessive.xls.CandidateGene.xls')
        modelC_snp_Dir = os.path.join(
            advanceDir, 'ModelF', each,
            each + '.snp.indel.compound_heterozygous.xls.CandidateGene.xls')
        modelR_anno_snp_stat.append(modelR_anno_snp_file)
        modelR_anno_indel_stat.append(modelR_anno_indel_file)
        modelC_anno_stat.append(modelC_anno_file)
        snpanno = len(utils.safe_open(modelR_anno_snp_file, 'r').readlines())
        indelanno = len(utils.safe_open(modelR_anno_indel_file, 'r').readlines())
        canno = len(utils.safe_open(modelC_anno_file, 'r').readlines())
        if snpanno > 1 and indelanno > 1 and canno > 1:
            context['familyid'] = each
            context['ModelR'] = True
            context['table_modelR_anno'] = txts2tab(
                [modelR_anno_snp_stat[0]])[:6]
            context['ModelC'] = True
            context['table_modelC_anno'] = txts2tab([modelC_anno_stat[0]])[:6]
            modelR_snplist = len(utils.safe_open(modelR_snp_Dir, 'r').readlines())
            modelR_indellist = len(
                utils.safe_open(modelR_indel_Dir, 'r').readlines())
            context['modelR_listsnp'] = modelR_snplist
            context['modelR_listindel'] = modelR_indellist
            modelC_snplist = len(utils.safe_open(modelC_snp_Dir, 'r').readlines())
            context['modelC_listsnp'] = modelC_snplist
        elif snpanno > 1 and indelanno > 1 and 1 >= canno:
            context['familyid'] = each
            context['ModelR'] = True
            context['table_modelR_anno'] = txts2tab(
                [modelR_anno_snp_stat[0]])[:6]
            context['ModelC'] = False
            modelR_snplist = len(utils.safe_open(modelR_snp_Dir, 'r').readlines())
            modelR_indellist = len(
                utils.safe_open(modelR_indel_Dir, 'r').readlines())
            context['modelR_listsnp'] = modelR_snplist
            context['modelR_listindel'] = modelR_indellist
            context['modelC_listsnp'] = 0
        elif snpanno > 1 and 1 >= indelanno and canno > 1:
            context['familyid'] = each
            context['ModelR'] = True
            context['table_modelR_anno'] = txts2tab(
                [modelR_anno_snp_stat[0]])[:6]
            context['ModelC'] = True
            context['table_modelC_anno'] = txts2tab([modelC_anno_stat[0]])[:6]
            modelR_snplist = len(utils.safe_open(modelR_snp_Dir, 'r').readlines())
            context['modelR_listsnp'] = modelR_snplist
            context['modelR_listindel'] = 0
            modelC_snplist = len(utils.safe_open(modelC_snp_Dir, 'r').readlines())
            context['modelC_listsnp'] = modelC_snplist
        elif snpanno > 1 and 1 >= indelanno and 1 >= canno:
            context['familyid'] = each
            context['ModelR'] = True
            context['table_modelR_anno'] = txts2tab(
                [modelR_anno_snp_stat[0]])[:6]
            context['ModelC'] = False
            modelR_snplist = len(utils.safe_open(modelR_snp_Dir, 'r').readlines())
            context['modelR_listsnp'] = modelR_snplist
            context['modelR_listindel'] = 0
            context['modelC_listsnp'] = 0
        elif 1 >= snpanno and indelanno > 1 and canno > 1:
            context['familyid'] = each
            context['ModelR'] = True
            context['table_modelR_anno'] = txts2tab(
                [modelR_anno_indel_stat[0]])[:6]
            context['ModelC'] = True
            context['table_modelC_anno'] = txts2tab([modelC_anno_stat[0]])[:6]
            modelR_indellist = len(
                utils.safe_open(modelR_indel_Dir, 'r').readlines())
            context['modelR_listsnp'] = 0
            context['modelR_listindel'] = modelR_indellist
            modelC_snplist = len(utils.safe_open(modelC_snp_Dir, 'r').readlines())
            context['modelC_listsnp'] = modelC_snplist
        elif 1 >= snpanno and indelanno > 1 and 1 >= canno:
            context['familyid'] = each
            context['ModelR'] = True
            context['table_modelR_anno'] = txts2tab(
                [modelR_anno_indel_stat[0]])[:6]
            context['ModelC'] = False
            modelR_indellist = len(
                utils.safe_open(modelR_indel_Dir, 'r').readlines())
            context['modelR_listsnp'] = 0
            context['modelR_listindel'] = modelR_indellist
            context['modelC_listsnp'] = 0
        elif 1 >= snpanno and 1 >= indelanno and canno > 1:
            context['familyid'] = each
            context['ModelR'] = True
            context['ModelC'] = True
            context['table_modelC_anno'] = txts2tab([modelC_anno_stat[0]])[:6]
            context['modelR_listsnp'] = 0
            context['modelR_listindel'] = 0
            modelC_snplist = len(utils.safe_open(modelC_snp_Dir, 'r').readlines())
            context['modelC_listsnp'] = modelC_snplist
        else:
            context['ModelR'] = False
            context['ModelC'] = False


    if 'share_compare' in ANALYSIS:
        shareCompare_stat = []
        shareCompare_file = os.path.join(
            advanceBriDir, 'Share',
            'snp.indel.filter.ShareGene.PatientShare.NotInNormal.brief.xls')
        if not os.path.exists(shareCompare_file):
            shareCompare_file = os.path.join(
                advanceBriDir, 'Share',
                'snp.indel.filter.ShareGene.PatientShare.brief.xls')
        shareCompare_stat.append(shareCompare_file)
        context['ShareCompare'] = True
        context['table_shareCompare'] = txts2tab([shareCompare_stat[0]])[:6]

    # ========================= Denovo ============================
    print 'denovo start'
    # only show intersect results
    if 'denovo' in ANALYSIS:

        denovo_rate_dir = os.path.join(advanceDir, 'Denovo', 'DenovoRate')
        denovo_intersect_dir = os.path.join(advanceDir, 'Denovo', 'Intersect')
        denovo_intersect_brief_dir = os.path.join(advanceBriDir, 'Denovo')

        # did means FamilyID_SampleID for denovo sample
        for tid in os.listdir(denovo_intersect_brief_dir):
            if tid in c:
                did = tid
                break     # only show one sample

        # 1 Filter Stat
        denovo_intersect_snp_stat = os.path.join(denovo_intersect_dir, did, 'SNP', did + '.denovo.snp.intersect.filter.stat.xls')
        context['table_denovo_intersect_snp_stat'] = txts2tab([denovo_intersect_snp_stat])[:6]

        denovo_intersect_indel_stat = os.path.join(denovo_intersect_dir, did, 'INDEL', did + '.denovo.indel.intersect.filter.stat.xls')
        context['table_denovo_intersect_indel_stat'] = txts2tab([denovo_intersect_indel_stat])[:6]

        # 2 Denovo Result [brief]
        denovo_intersect_snp_result = os.path.join(denovo_intersect_brief_dir, did, did + '.denovo.snp.intersect.freq.func.syn.deleterious.brief.xls')
        context['table_denovo_intersect_snp_result'] = txts2tab([denovo_intersect_snp_result])[:6]

        denovo_intersect_indel_result = os.path.join(denovo_intersect_brief_dir, did, did + '.denovo.indel.intersect.freq.func.syn.deleterious.brief.xls')
        context['table_denovo_intersect_indel_result'] = txts2tab([denovo_intersect_indel_result])[:6]

        # 3 Denovo Rate
        denovo_rate_snp = os.path.join(denovo_rate_dir, 'denovo_rate_snp.xls')
        context['table_denovo_rate_snp'] = txts2tab([denovo_rate_snp])[:6]

        denovo_rate_indel = os.path.join(denovo_rate_dir, 'denovo_rate_indel.xls')
        context['table_denovo_rate_indel'] = txts2tab([denovo_rate_indel])[:6]


    # if 'denovo_samtools' in ANALYSIS:
    #     denovos_indel_stat = []
    #     denovos_snp_anno = []
    #     denovos_anno_stat = []
    #     denovos_dir = os.path.join(advanceDir, 'Denovo', 'DenovoSam')
    #     for tid in os.listdir(denovos_dir):
    #         if tid in c:
    #             did = tid
    #             break
    #     denovos_snp_anno_file = os.path.join(
    #         denovos_dir, did, did + '.denovo.snp.annovar.hg19_multianno.xls')
    #     denovos_anno_stat_file = os.path.join(
    #         denovos_dir, did, did + '.denovo.snp.filter.stat.xls')
    #     denovos_indel_stat_file = os.path.join(
    #         denovos_dir, did, did + '.denovo.indel.filter.stat.xls')
    #     denovos_snp_anno.append(denovos_snp_anno_file)
    #     denovos_anno_stat.append(denovos_anno_stat_file)
    #     denovos_indel_stat.append(denovos_indel_stat_file)
        
    #     context['Denovoid'] = did
    #     context['denovo_samtools'] = True
    #     context['table_denovos_snp_anno'] = txts2tab([denovos_snp_anno[0]])[:6]
    #     context['table_denovos_anno_stat'] = txts2tab(
    #         [denovos_anno_stat[0]])[:6]
    #     context['table_denovos_indel_stat'] = txts2tab(
    #         [denovos_indel_stat[0]])[:6]

    # denovog_snp_anno = []
    # denovog_anno_stat = []
    # if 'denovo_denovogear' in ANALYSIS:
    #     denovog_dir = os.path.join(advanceDir, 'Denovo', 'DenovoGear')
    #     for tid in os.listdir(denovog_dir):
    #         if tid in c:
    #             did = tid
    #             break
    #     denovog_snp_anno_file = os.path.join(
    #         denovog_dir, did, 'SNP',
    #         did + '.denovo.snp.annovar.hg19_multianno.xls')
    #     denovog_anno_stat_file = os.path.join(
    #         denovog_dir, did, 'SNP', did + '.denovo.snp.filter.stat.xls')
    #     denovog_snp_anno.append(denovog_snp_anno_file)
    #     denovog_anno_stat.append(denovog_anno_stat_file)
    #     context['Denovoid'] = did
    #     context['denovo_denovogear'] = True
    #     context['table_denovog_snp_anno'] = txts2tab([denovog_snp_anno[0]])[:6]
    #     context['table_denovog_anno_stat'] = txts2tab(
    #         [denovog_anno_stat[0]])[:6]

    # denovo_filter_indel_stat = []
    # denovo_filter_snp_anno = []
    # denovo_filter_anno_stat = []
    # if 'denovo_filter' in ANALYSIS:
    #     denovo_filter_dir = os.path.join(advanceDir, 'Denovo', 'denovo_filter')
    #     for tid in os.listdir(denovo_filter_dir):
    #         if tid in c:
    #             did = tid
    #             break
    #     denovo_filter_snp_anno_file = os.path.join(
    #         denovo_filter_dir, did, 'SNP',
    #         did + '.trio.merge.snp.denovo.annovar.hg19_multianno.xls')
    #     denovo_filter_anno_stat_file = os.path.join(
    #         denovo_filter_dir, did, 'SNP', 'Denovo.' + did + '.snp.filter.stat.xls')
    #     denovo_filter_indel_stat_file = os.path.join(
    #         denovo_filter_dir, did, 'INDEL',
    #         'Denovo.' + did + '.indel.filter.stat.xls')
    #     denovo_filter_snp_anno.append(denovo_filter_snp_anno_file)
    #     denovo_filter_anno_stat.append(denovo_filter_anno_stat_file)
    #     denovo_filter_indel_stat.append(denovo_filter_indel_stat_file)
    #     context['Denovoid'] = did
    #     context['denovo_filter'] = True
    #     context['table_denovo_filter_snp_anno'] = txts2tab([denovo_filter_snp_anno[0]])[:6]
    #     context['table_denovo_filter_anno_stat'] = txts2tab(
    #         [denovo_filter_anno_stat[0]])[:6]
    #     context['table_denovo_filter_indel_stat'] = txts2tab(
    #         [denovo_filter_indel_stat[0]])[:6]

    # denovor_snp_anno = []
    # denovor_anno_stat = []
    # if 'denovo_raw' in ANALYSIS:
    #     denovor_dir = os.path.join(advanceDir, 'Denovo', 'DenovoR')
    #     for tid in os.listdir(denovor_dir):
    #         if tid in c:
    #             did = tid
    #             break
    #     denovor_snp_anno_file = os.path.join(
    #         denovor_dir, did, 'SNP',
    #         did + '.trio.merge.snp.denovo.final.annovar.hg19_multianno.xls')
    #     denovor_anno_stat_file = os.path.join(
    #         denovor_dir, did, 'SNP', 'Denovo.' + did + '.snp.filter.stat.xls')
    #     denovor_snp_anno.append(denovor_snp_anno_file)
    #     denovor_anno_stat.append(denovor_anno_stat_file)
    #     context['Denovoid'] = did
    #     context['DenovoR'] = True
    #     context['table_denovor_snp_anno'] = txts2tab([denovor_snp_anno[0]])[:6]
    #     context['table_denovor_anno_stat'] = txts2tab(
    #         [denovor_anno_stat[0]])[:6]

    # denovosaf_snp_anno = []
    # denovorate_snp = []
    # denovosaf_indel_anno = []
    # denovorate_indel = []
    # denovosaf_snp_brief = []
    # if ('denovo_samtools' in ANALYSIS) and ('denovo_filter' in ANALYSIS):
    #     denovos_dir = os.path.join(advanceDir, 'Denovo', 'DenovoSam')
    #     denovosaf_dir = os.path.join(advanceDir, 'Denovo', 'Intersect')
    #     denovosaf_bri_dir = os.path.join(advanceBriDir, 'Denovo')
    #     denovorate_dir = os.path.join(advanceDir, 'Denovo')
    #     for tid in os.listdir(denovos_dir):
    #         if tid in c:
    #             did = tid
    #             break
    #     denovosaf_snp_anno_file = os.path.join(
    #         denovosaf_dir, did, did + '.denovo.snp.intersect.filter.stat.xls')
    #     denovosaf_snp_anno.append(denovosaf_snp_anno_file)
    #     denovorate_snp_file = os.path.join(denovorate_dir, 'DenovoSNPRate.xls')
    #     denovorate_snp.append(denovorate_snp_file)
    #     denovosaf_indel_anno_file = os.path.join(
    #         denovosaf_dir, did,
    #         did + '.denovo.indel.intersect.filter.stat.xls')
    #     denovosaf_indel_anno.append(denovosaf_indel_anno_file)
    #     denovorate_indel_file = os.path.join(denovorate_dir,
    #                                          'DenovoInDelRate.xls')
    #     denovorate_indel.append(denovorate_indel_file)
    #     denovorate_anno_brief = os.path.join(
    #         denovosaf_bri_dir,
    #         did + '.denovo.snp.intersect.freq.func.syn.deleterious.brief.xls')
    #     denovosaf_snp_brief.append(denovorate_anno_brief)
    #     context['denovo_filter'] = True
    #     context['denovo_samtools'] = True
    #     context['table_denovosaf_snp'] = txts2tab([denovosaf_snp_anno[0]])[:6]
    #     context['table_denovorate'] = txts2tab([denovorate_snp[0]])[1:]
    #     context['table_denovosaf_indel'] = txts2tab(
    #         [denovosaf_indel_anno[0]])[:6]
    #     context['table_denovorate_indel'] = txts2tab([denovorate_indel[0]])[1:]
    #     context['table_denovorate_brief'] = txts2tab(
    #         [denovosaf_snp_brief[0]])[:6]

    denovosv_anno = []
    if 'denovo_sv' in ANALYSIS:
        denovosv_dir = os.path.join(advanceDir, 'Denovo', 'DenovoSV')
        denovosv_bri_dir = os.path.join(advanceBriDir, 'Denovo')
        for tid in os.listdir(denovosv_dir):
            if tid in c:
                did = tid
                break
        context['svid'] = did
        denovosv_anno_file = os.path.join(denovosv_dir, did,
                                          did + '.denovoSV.hg19_multianno.xls')
        denovosv_anno.append(denovosv_anno_file)
        context['table_denovosv'] = txts2tab([denovosv_anno[0]])[:6]

    denovocnv_anno = []
    if 'denovo_cnv' in ANALYSIS:
        denovocnv_dir = os.path.join(advanceDir, 'Denovo', 'DenovoCNV')
        denovocnv_bri_dir = os.path.join(advanceBriDir, 'Denovo')
        for tid in os.listdir(denovocnv_dir):
            if tid in c:
                did = tid
                break
        context['cnvid'] = did
        denovocnv_anno_file = os.path.join(
            denovocnv_dir, did, did + '.denovoCNV.hg19_multianno.xls')
        denovocnv_anno.append(denovocnv_anno_file)
        context['table_denovocnv'] = txts2tab([denovocnv_anno[0]])[:6]

    #####非编码区筛选,changed by yincan 20171221
    if 'filter_noncoding' in ANALYSIS:
        noncoding_stat = []
        NC_Dir = os.path.join(advanceDir, 'Noncoding')
        NC_bri_Dir = os.path.join(advanceBriDir, 'Noncoding')
        NC_snp = [
            os.path.join(
                NC_Dir, each, each +
                '.snp.anno.noncoding.GnomAD_EAS_AF.conserv.epigenome.xls.gz')
        ]

        noncoding_stat_file_snp = os.path.join(
            NC_Dir, each, each + '.snp.anno.noncoding.filter.stat.xls.gz')
        noncoding_stat.append(noncoding_stat_file_snp)

        context['table_snp_noncoding'] = txts2tab(noncoding_stat)
        if os.path.exists(NC_snp[0]):
            context['Noncoding_anno'] = txts2tab(NC_snp)[:6]

sam = args['sample']
if sam != 'Null':
    samf = open(sam, 'r')
    for line in samf:
        if line.lower().startswith('#familyid'):
            lines = line.lower().strip().split('\t')
            continue
        if line.startswith('#'): continue
        line = line.strip().split('\t')
        if 'data' in lines:
            datapos = lines.index('data')
            if len(line) >= (datapos + 1):
                if line[datapos] == "1" or line[datapos] == "3":
                    each = line[0]
                    linkage_stat = []
                    context['fig_linkage'] = []
                    if 'linkage' in ANALYSIS:
                        linkage_dir = os.path.join(advanceDir, 'MerLinkage')
                        linkage_stat_file = os.path.join(
                            linkage_dir, 'Linkage_' + each + '_me', 'report',
                            'LinkageAnalysis_' + each + '.xls')
                        linkage_stat.append(linkage_stat_file)
                        linkagefile = len(
                            utils.safe_open(linkage_stat_file, 'r').readlines())
                        if linkagefile > 1:
                            context['Linkage'] = True
                            context['table_linkage_stat'] = txts2tab(
                                [linkage_stat[0]])[:6]
                            linkage_pic_dir = os.path.join(
                                linkage_dir, 'Linkage_' + each + '_me',
                                'report')
                            linkage_images_dir = os.path.join(
                                odir, 'src/pictures/Other')
                            assert not os.system(
                                'cp %s %s' %
                                (os.path.join(linkage_pic_dir, 'merlin_R_' +
                                              each + '-nonparametric.png'),
                                 os.path.join(linkage_images_dir, 'merlin_R_' +
                                              each + '-nonparametric.png')))
                            assert not os.system(
                                'convert -resize 800 %s %s' %
                                (os.path.join(linkage_pic_dir, 'merlin_R_' +
                                              each + '-nonparametric.png'),
                                 os.path.join(linkage_images_dir, 'merlin_R_' +
                                              each + '-nonparametric.png')))
                            assert not os.system(
                                'convert -resize 90 %s %s' %
                                (os.path.join(linkage_pic_dir, 'merlin_R_' +
                                              each + '-nonparametric.png'),
                                 os.path.join(linkage_images_dir, 'merlin_R_' +
                                              each + '-nonparametric.JPEG')))
                            context['fig_linkage'].append([
                                '"' + os.path.join('src/pictures/Other',
                                                   'merlin_R_' + each +
                                                   '-nonparametric.png') + '"',
                                '"' + os.path.join('src/pictures/Other',
                                                   'merlin_R_' + each +
                                                   '-nonparametric.JPEG') + '"'
                            ])

sam = args['sample']
if sam != 'Null':
    samf = open(sam, 'r')
    for line in samf:
        if line.lower().startswith('#familyid'):
            lines = line.lower().strip().split('\t')
            continue
        if line.startswith('#'): continue
        line = line.strip().split('\t')
        if 'data' in lines:
            datapos = lines.index('data')
            if len(line) >= (datapos + 1):
                if line[datapos] != "0":
                    line1 = line[0]
                    line2 = line[1]
                    roh_snp_anno = []
                    if 'roh' in ANALYSIS:
                        roh_dir = os.path.join(advanceDir, 'ROH')
                        roh_bri_dir = os.path.join(advanceBriDir, 'ROH')
                        roh_snp_file = os.path.join(
                            roh_bri_dir, line2 + '_roh_anno_snp.brief.xls')
                        #roh_snp_file=os.path.join(roh_dir,line2,line2+'_roh_anno_snp.xls')
                        roh_snp_anno.append(roh_snp_file)
                        context['ROH'] = True
                        context['table_roh_anno'] = txts2tab(
                            [roh_snp_anno[0]])[:6]
            else:
                line1 = line[0]
                line2 = line[1]
                roh_snp_anno = []
                if 'roh' in ANALYSIS:
                    roh_dir = os.path.join(advanceDir, 'ROH')
                    roh_bri_dir = os.path.join(advanceBriDir, 'ROH')
                    roh_snp_file = os.path.join(
                        roh_bri_dir, line2 + '_roh_anno_snp.brief.xls')
                    #roh_snp_file=os.path.join(roh_dir,line2,line2+'_roh_anno_snp.xls')
                    roh_snp_anno.append(roh_snp_file)
                    context['ROH'] = True
                    context['table_roh_anno'] = txts2tab([roh_snp_anno[0]])[:6]
        else:
            line1 = line[0]
            line2 = line[1]
            roh_snp_anno = []
            if 'roh' in ANALYSIS:
                roh_dir = os.path.join(advanceDir, 'ROH')
                roh_bri_dir = os.path.join(advanceBriDir, 'ROH')
                roh_snp_file = os.path.join(roh_bri_dir,
                                            line2 + '_roh_anno_snp.brief.xls')
                roh_snp_anno.append(roh_snp_file)
                context['ROH'] = True
                context['table_roh_anno'] = txts2tab([roh_snp_anno[0]])[:6]

    #func_enrichement=[]
    #if Funcenrichement:
    #    Func_dir=os.path.join(advanceDir,'Functional_enrichment')
    #    print Func_dir
    #    Func_enrichment_file=os.path.join(Func_dir,'Functional_enrichment.xls')
    #    func_enrichement.append(Func_enrichment_file)
    #    context['Funcenrichement']=True
    #    context['table_func_enrichement']=txts2tab([func_enrichement[0]])[:6]

odir = odir.rstrip('/')

nonnum = 2
fsvnum = 2
fcnvnum = 2

if 'filter_acmg' in ANALYSIS:
    nonnum += 1
    fsvnum += 1
    fcnvnum += 1
if 'filter_noncoding' in ANALYSIS:
    fsvnum += 1
    fcnvnum += 1
if 'filter_sv' in ANALYSIS:
    fcnvnum += 1

sharenum = 1
denovonum = 1
linkagenum = 1
rohnum = 1
ppinum = 1
hpanum = 1
if ('model_dominant' in ANALYSIS) or ('model_recessive' in ANALYSIS):
    sharenum += 1
    denovonum += 1
    linkagenum += 1
    rohnum += 1
    ppinum += 1
    hpanum += 1
if 'share_compare' in ANALYSIS:
    denovonum += 1
    linkagenum += 1
    rohnum += 1
    ppinum += 1
    hpanum += 1
if 'denovo' in ANALYSIS:
    linkagenum += 1
    rohnum += 1
    ppinum += 1
    hpanum += 1
if 'linkage' in ANALYSIS:
    rohnum += 1
    ppinum += 1
    hpanum += 1
#if
context['nonnum'] = nonnum
context['fsvnum'] = fsvnum

context['sharenum'] = sharenum
context['denovonum'] = denovonum
context['linkagenum'] = linkagenum
context['rohnum'] = rohnum
context['hpanum'] = hpanum

if WES_xten == 'Y':
    context['WES_xten'] = True

# ====================
# =  Render to HTML  =
# ====================
gc_tpl = env.get_template('GCcheck.html')
gc_html = gc_tpl.render(context)

with utils.safe_open(os.path.join(odir, 'src', 'pictures', 'GC', 'GCcheck.html'), 'w') as out:
    out.write(gc_html)

if english == 'Y':
    context['ER'] = True
    template = env.get_template('Report_template_disease_English.html')
else:
    template = env.get_template('index.html')


reportDir = os.path.split(odir)[0]
reportF = reportDir + '/' + hetongnum + '-' + bg + '-5.tar'

if rep_ty == 'qc':
    ReportType = 'QC'
    outfile = os.path.join(odir, 'QC_Report.html')
elif rep_ty == 'mapping':
    ReportType = 'Mapping'
    outfile = os.path.join(odir, 'Mapping_Report.html')
elif rep_ty == 'primary':
    ReportType = 'Primary'
    outfile = os.path.join(odir, 'Primary_Report.html')
elif rep_ty == 'advance':
    ReportType = 'Advance'
    outfile = os.path.join(odir, 'Advance_Report.html')
else:
    print "wrong rep_ty\n"

#  PDF Report
pdf_report_href = False
if pdf_report == 'Y':
    temp_html = outfile.replace('html', 'temp.html')
    with open(temp_html, 'w') as temp:
        context['pdf_report'] = True
        pdf_html = template.render(context)
        temp.write(pdf_html)

    context['pdf_report_href'] = 'src/{ReportType}_Report.pdf'.format(**locals())

    cmd = '''
    wkhtmltopdf \\
        -n \\
        --print-media-type \\
        --images \\
        --page-size A4 \\
        --debug-javascript \\
        --enable-javascript \\
        --no-stop-slow-scripts \\
        --javascript-delay 15000 \\
        {temp_html} \\
        {ReportType}/src/{ReportType}_Report.pdf

    rm -f {temp_html}
    '''.format(**locals())

    print cmd
    os.system(cmd)

# Disease Report
find_doid_bg_html = False
if diseases and args['disease'] == 'Y' and rep_ty == 'advance':
    
    database = '{}/project/DisGeNet.json'.format(config.CONFIG.get('software', 'soft_dir'))

    doids = utils.get_disease_id(diseases, database)

    for doid in doids.split(';'):
        doid_bg_html = '{disease_db}/Disease_BackGround/{doid}__*html'.format(**locals())
        print 'looking for file:', doid_bg_html
        if glob.glob(doid_bg_html):
            find_doid_bg_html = True
            context['doid_bg_html'] = os.path.basename(glob.glob(doid_bg_html)[0])
            break

    if find_doid_bg_html:

        disease_report_gene = []
        disease_report_table = [] 

        disease_result = '{projdir}/Advance/{newjob}/Disease/Integrate.disease.xls'.format(**locals())
        with utils.safe_open(disease_result) as f:
            for line in f:
                linelist = line.strip().split('\t')
                disease_report_table.append(linelist)
                if linelist[0] != "GeneCandidate":
                    disease_report_gene += linelist[0].split(',')
        disease_report_gene = set(disease_report_gene)

        context['disease_report_table'] = disease_report_table
        context['disease_report_gene'] = disease_report_gene
    else:
        print '[warn] this disease is not supported to do disease report yet:', diseases

# Generate HTML
context['pdf_report'] = False
html = template.render(context)

with utils.safe_open(outfile, 'w') as out:
    out.write(html)

# Package report result
pubreportdir = config.CONFIG.get('software', 'report_store')
if 'WES' in seq_ty:
    projectreportdir = os.path.join(
        pubreportdir, 'WES' + '.' + hetongnum + '.' + 'Reports', newjob)
    tmpname = 'WES' + '.' + hetongnum + '.' + 'Reports'
else:
    projectreportdir = os.path.join(
        pubreportdir + '/' + seq_ty + '.' + hetongnum + '.' + 'Reports', newjob)
    tmpname = seq_ty + '.' + hetongnum + '.' + 'Reports'

remote_report_dir = os.path.join(
    remote_dir, tmpname, newjob, os.path.split(odir)[1])

projectdir = os.path.join(pubreportdir, tmpname)
print 'The report was stored to: ', projectdir

if not os.path.exists(projectdir):
    os.mkdir(projectdir)
    os.system('chmod -R 775 %s' % projectdir)
if not os.path.exists(projectreportdir):
    os.mkdir(projectreportdir)
    os.system('chmod -R 775 %s' % projectreportdir)
os.system('ln -sf %s %s' % (odir, projectreportdir))


odir_base = os.path.basename(odir)
analy_array = args['analy_array']
report_date = utils.get_now_time('%Y%m%d')
hetong = '{}-{}'.format(
    pn.split()[0], bg)

cmd = '''
    source  ~/.bash_profile

    set -e
    cd {odir}/..

    tar -hcf {hetong}.{report_date}.{ReportType}.tar {ReportType}
    ln -sf {hetong}.{report_date}.{ReportType}.tar {hetong}-5.tar

    python2 {REPORT_DIR}/report_upload.py {odir} {remote_report_dir}

    sendEmail -f humaninfo@novogene.com -t {mail} \\
        -u "[{projdir}] {ReportType} report" \\
        -o message-content-type=text \\
        -o message-charset=utf8 \\
        -xu humaninfo@novogene.com \\
        -xp DhumanB0206 \\
        -s 183.57.48.39 \\
        -o tls=no -m "
    Congratulations! Your {ReportType} report of {projdir} has done!\n
    Please visit http://{report_host}/{tmpname}/{newjob}/{odir_base}/{ReportType}_Report.html for the full report.\n
    and visit http://{report_host}/{tmpname}/{newjob}/{odir_base}/src/pictures/GC/GCcheck.html for GC check information."
'''

if datastat == 'Y':
    cmd += '''
    python {moduledir}/ProjStat/mongoinsert/mongoupdate.py \\
        --ref {ref} \\
        --project {proj_n} \\
        --projdir {projdir} \\
        --sample {sam} \\
        --suffix {suffix} \\
        --analy_array {analy_array} \\
        --mail "{mail}" \\
        --seqsty {seq_ty} \\
        --repsty {rep_ty}
    '''

cmd = cmd.format(**locals())
print textwrap.dedent('\033[33;40m' + cmd + '\033[0m')

if 'NJ' in ROOT_DIR:
    host = 'njlogin04'
else:
    host = 'login04'

# os.system("ssh {} '{}'".format(host, cmd))

print cumudep
print '---' * 10
print distdep

# import json
# print json.dumps(context, indent=2, ensure_ascii=False)
