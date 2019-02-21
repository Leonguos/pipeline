#!/usr/bin/env python
# -*- coding=utf-8 -*-
'''
Project associate tools
'''
import os
from collections import defaultdict, OrderedDict

from .common import get_index, get_value, safe_open, get_flowcell_lib


class Project(object):

    def __init__(self, analydir, sample_info, sample_info_done, sample_list,
                 qc_list, qc_status, mapping_status, is_advance=True):
        self.analydir = analydir
        self.sample_info = sample_info
        self.sample_info_done = sample_info_done
        self.sample_list = sample_list
        self.qc_list = qc_list
        self.qc_status = qc_status
        self.mapping_status = mapping_status
        self.is_advance = is_advance

    @staticmethod
    def get_sample_infos(sample_list, sample_info, sample_info_done=None, is_advance=True):

        fenqi = 'B1'
        tissue = None
        disease_name = None
        sample_infos = OrderedDict()       # 只包含有测序数据的样本
        sample_infos_all = OrderedDict()   # 包含全部样本（包括为例连锁分析添加的无测序数据样本）

        samples_in_list = set()
        with open(sample_list) as f:
            for line in f:
                linelist = line.strip().split('\t')
                if not line.strip():
                    continue
                elif line.startswith('#') and 'sampleid' in line.lower():
                    headerlist = map(str.lower, linelist)
                    continue
                elif line.startswith('#'):
                    continue

                sampleid = linelist[headerlist.index('sampleid')]
                samples_in_list.add(sampleid)

        with open(sample_info) as f:
            for line in f:
                if line.lower().startswith('#b'):
                    fenqi = line.strip('#').strip()
                elif line.lower().startswith('#disease'):
                    disease_name = line.strip().split(':')[-1].strip().strip('"').strip("'")
                elif line.lower().startswith('#tissue'):
                    tissue = line.strip().split(':')[-1].strip().strip('"').strip("'")
                elif line.lower().startswith('#familyid'):
                    title = [t.lower() for t in line.lstrip('#').strip().split('\t')]
                elif line.startswith('#') or not line.strip():
                    continue
                else:
                    linelist = line.strip().split('\t')
                    # print linelist
                    sampleID = get_value(linelist, get_index('sampleid', title))

                    familyid = get_value(linelist, get_index('familyid', title))
                    novoid = get_value(linelist, get_index('novoid', title))
                    sex = get_value(linelist, get_index('sex', title))
                    phenotype = get_value(linelist, get_index('normal/patient', title))
                    pn = get_value(linelist, get_index('pn', title))
                    diseaseid = get_value(linelist, get_index('disease', title))
                    pa = get_value(linelist, get_index('pa', title))
                    ma = get_value(linelist, get_index('ma', title))

                    # 只有data列标0时表示该样本无测序数据
                    data = get_value(linelist, get_index('data', title))

                    if is_advance and familyid in ('.', None):
                        print '[error] advance analysis needs a FamilyID'

                    sample_infos_all[sampleID] = {
                        'familyid': familyid,
                        'sampleid': sampleID,
                        'novoid': novoid,
                        'sex': sex,
                        'phenotype': phenotype,
                        'diseaseid': diseaseid,
                        'data': data,
                        'pa': pa,
                        'ma': ma,
                    }

                    # 跳过sample_info中没有测序数据的样本
                    if (sampleID not in samples_in_list) or (data == '0'):
                        continue

                    sample_infos[sampleID] = {
                        'familyid': familyid,
                        'sampleid': sampleID,
                        'novoid': novoid,
                        'sex': sex,
                        'phenotype': phenotype,
                        'diseaseid': diseaseid,
                        'data': data,
                        'pa': pa,
                        'ma': ma,
                    }

        # 不重新分析sample_info_done中的样本
        sample_done = []
        if sample_info_done and os.path.isfile(sample_info_done):
            with open(sample_info_done) as f:
                for line in f:
                    linelist = line.strip().split('\t')
                    if line.startswith('#') or not line.strip():
                        if 'sampleid' in line.lower():
                            headerlist = map(str.lower, linelist)
                            # print headerlist
                        continue
                    sample_done.append(linelist[headerlist.index('sampleid')])

                # print sample_done

        return fenqi, tissue, disease_name, sample_infos, sample_infos_all, sample_done

    @property
    def get_sample_lists(self):

        sample_lists = defaultdict(dict)
        samples_in_info = self.get_sample_infos(
            self.sample_list, self.sample_info, self.sample_info_done, self.is_advance)[-3].keys()

        with open(self.sample_list) as f:
            for line in f:
                if not line.strip():
                    continue
                elif line.startswith('#'):
                    title = [t.lower() for t in line.lstrip('#').strip().split('\t')]
                    # print title
                else:
                    linelist = line.strip().split('\t')
                    # print linelist
                    # 跳过sample_list中有，sample_info中没有的样本
                    sampleID = get_value(linelist, get_index('sampleid', title))
                    if sampleID not in samples_in_info:
                        print '[warn] sample {} not in sample_info, skipped'.format(sampleID)
                        continue

                    lane = get_value(linelist, get_index('lane', title, 'ori_lane'))[-1]

                    patientID = get_value(linelist, get_index('patientid', title))
                    libID = get_value(linelist, get_index('libid', title))
                    novoid = get_value(linelist, get_index('novoid', title))
                    index = get_value(linelist, get_index('index', title))
                    path = get_value(linelist, get_index('path', title)).rstrip('/')

                    # print path

                    flowcell, libID = get_flowcell_lib(path, libID, index, lane)

                    if not flowcell:
                        flowcell = 'NoFC'

                    # Require rawdata but file not exists
                    if self.qc_status == 'waiting' and flowcell == 'NoFC':
                        print '  [error] rawdata not exists for sample {} lane {}: {}'.format(sampleID, lane, os.path.join(path, libID))
                        exit(1)

                    if 'lanes' not in sample_lists[sampleID]:
                        sample_lists[sampleID]['lanes'] = []

                    items = {
                        'patientID': patientID,
                        'path': path,
                        'libID': libID,
                        'novoid': novoid,
                        'index': index,
                        'lane': lane,
                        'flowcell': flowcell,
                        'flowcell_lane': '{}_L{}'.format(flowcell, lane)
                    }
                    if items not in sample_lists[sampleID]['lanes']:
                        sample_lists[sampleID]['lanes'].append(items)
        # print dict(sample_lists)
        return dict(sample_lists)

    @staticmethod
    def get_qc_lists(qc_list):

        qc_lists = defaultdict(dict)

        with open(qc_list) as f:
            for line in f:
                if not line.strip():
                    continue
                elif line.startswith('#'):
                    title = [t.lower() for t in line[1:].strip().split('\t')]
                    # print title
                else:
                    linelist = line.strip().split('\t')
                    # print linelist
                    title = title or []
                    flowcell_lane = get_value(linelist, get_index('flowcell_lane', title) or 0)
                    patientID = get_value(linelist, get_index('patientid', title) or 1)
                    sampleID = get_value(linelist, get_index('sampleid', title) or 2)
                    libID = get_value(linelist, get_index('libid', title) or 3)
                    novoid = get_value(linelist, get_index('novoid', title) or 4)
                    index = get_value(linelist, get_index('index', title) or 5)
                    path = get_value(linelist, get_index('path', title) or 6)

                    if 'lanes' not in qc_lists[sampleID]:
                        qc_lists[sampleID]['lanes'] = []

                    items = {
                        'patientID': patientID,
                        'path': path,
                        'libID': libID,
                        'novoid': novoid,
                        'index': index,
                        'flowcell_lane': flowcell_lane
                    }
                    if items not in qc_lists[sampleID]['lanes']:
                        qc_lists[sampleID]['lanes'].append(items)
        return dict(qc_lists)

    def update_qc_list(self):

        sample_lists = self.get_sample_lists

        if not os.path.exists(self.qc_list):
            print 'generate qc_list ...'
            with open(self.qc_list, 'w') as out:
                title = '#Flowcell_Lane PatientID SampleID LibID NovoID Index Path'.split()
                out.write('\t'.join(title) + '\n')
                for sampleID, items in sample_lists.items():
                    for item in items['lanes']:
                        flowcell = item['flowcell']
                        lane = item['lane']
                        flowcell_lane = '{}_L{}'.format(flowcell, lane)
                        patientID = item['patientID']
                        libID = item['libID']
                        novoid = item['novoid']
                        index = item['index']
                        path = item['path']
                        linelist = [flowcell_lane, patientID, sampleID, libID, novoid, index, path]
                        if flowcell:
                            out.write('\t'.join(linelist) + '\n')
        else:
            print 'update qc_list ...'
            qc_lists = self.get_qc_lists(self.qc_list)
            # print qc_lists
            jiace_samples = []
            with open(self.qc_list, 'a') as out:
                for sampleid, items in sample_lists.items():
                    for item in items['lanes']:
                        flowcell = item['flowcell']

                        if flowcell == 'NoFC':
                            continue

                        lane = item['lane']
                        flowcell_lane = '{}_L{}'.format(flowcell, lane)
                        patientID = item['patientID']
                        libID = item['libID']
                        novoid = item['novoid']
                        index = item['index']
                        path = item['path']
                        linelist = [flowcell_lane, patientID, sampleid, libID, novoid, index, path]
                        if sampleid in qc_lists:
                            old_flowcell_lanes = [each.get('flowcell_lane') for each in qc_lists[sampleid]['lanes']]
                            # print sampleid, flowcell_lane, old_flowcell_lanes
                            # exit()
                            # print sampleid, flowcell_lane
                            if flowcell_lane not in old_flowcell_lanes:
                                self.check_sort_bam(patientID, sampleid)
                                jiace_samples.append(sampleid)
                                sample_lists[sampleid].update({'jiace': True})
                        out.write('\t'.join(linelist) + '\n')
            if jiace_samples:
                print '  jiace samples: \033[31m{}\033[0m'.format(jiace_samples)
            cmd = 'sort -u {qc_list} > temp && mv temp {qc_list}'.format(qc_list=self.qc_list)
            os.system(cmd)
            print '  updated qc_list ...'
            
        return sample_lists

    def check_sort_bam(self, patientID, sampleID):

        sort_bam = '{analydir}/Mapping/{patientID}.{sampleID}/{sampleID}.sort.bam'.format(
            analydir=self.analydir, patientID=patientID, sampleID=sampleID
        )

        if not os.path.exists(sort_bam):
            print '[error] sort.bam not exists for jiace sample: {}'.format(sampleID)
            exit(1)

# the end
